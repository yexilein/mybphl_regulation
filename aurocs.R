
library(tidyverse)

#' predictors is a matrix where each column is a predictor and each row is a sample.
#' label_matrix is a binary matrix where columns are labels and each row is a sample.
#' 1 indicates that the sample on this row belongs to the label on this column.
compute_aurocs = function(predictors, label_matrix, return_tie_correction = FALSE) {
    label_matrix = as.matrix(label_matrix)
    n_samples = nrow(label_matrix)
    n_positives = colSums(label_matrix)
    n_negatives = n_samples - n_positives
    if (is(predictors, "dgCMatrix")) {
        # we shift all ranks after the matrix multiplication to keep
        # the predictor matrix sparse
        ranks = rank_sparse(predictors)
        sum_of_positive_ranks = as.matrix(Matrix::crossprod(label_matrix, ranks)) +
            outer(n_positives, rank_zero(predictors))
    } else {
        predictors = as.matrix(predictors)
        ranks = matrixStats::colRanks(predictors, ties.method = "average", preserveShape=TRUE)
        sum_of_positive_ranks = crossprod(label_matrix, ranks)
        colnames(sum_of_positive_ranks) = colnames(predictors)
    }
    if (return_tie_correction) {
      tie_correction = compute_tie_correction(ranks)
    }
    result = (sum_of_positive_ranks / n_positives - (n_positives+1)/2) / n_negatives
    if (return_tie_correction) {
      return(list(aurocs = result, tie_corrections = tie_correction))
    } else {
      return(result)
    }
}

design_matrix = function(cell_type, scale_columns=FALSE) {
  factors = levels(as.factor(cell_type))
  if (length(factors) > 1) {
    result = model.matrix(~as.factor(cell_type)-1)
  } else {
    result = matrix(1, nrow = length(cell_type), ncol = 1)
  }
  colnames(result) = factors
  if (scale_columns) {
    result = scale(result, center = FALSE, scale = colSums(result))
  }
  return(result)
}

# these 2 functions only rank non-zeros, implicitly shifting the matrix of ranks
# to keep the matrix sparse according to the formula:
#   true_ranks(M) = rank_zero(M) + rank_sparse(M), where:
#     + rank_zero(M) = #negatives + (#zeros + 1)/2 = #neg + #zeros/2 + 0.5
#     + rank_sparse(M) =
#        if positive: (rank(nnz) + #zeros) - rank_zero = rank(nnz) - #neg + #zeros/2 - 0.5
#        if negative: rank(nnz) - rank_zero = rank(nnz) - #neg - #zeros/2 - 0.5
#     + note that #zeros = nrow(M) - diff(M@p)
# faster than solution using base::tapply, probably data.table would be faster
rank_sparse = function(M) {
    nnz = diff(M@p)
    ranks = tibble(x = M@x, j = rep.int(1:ncol(M), nnz)) %>%
        group_by(j) %>%
        mutate(rank_ = as.vector(colRanks(as.matrix(x), ties.method = "average")))
    R = M
    offset = rep.int((nrow(M)-nnz)/2, nnz) # #zeros/2
    n_neg = rep.int(Matrix::colSums(M<0), nnz)
    is_positive = 2*(M@x>0) - 1 # +1 if positive, -1 if negative
    R@x = ranks$rank_ - n_neg + is_positive*offset - 0.5
    return(R)
}

# #negatives + (#zeros+1)/2
#   #zeros = nrow(M) - diff(M@p)
#   #negatives = 
rank_zero = function(M) {
    return(Matrix::colSums(M<0) + ((nrow(M) - diff(M@p)) + 1) / 2)
}

# For the following two functions, see
#
#   https://en.wikipedia.org/wiki/Mann%E2%80%93Whitney_U_test#Normal_approximation_and_tie_correction
#
# The tie correction effectively computes lost variance because of ties (compared to discrete uniform).
# Computing the Wikipedia formula naively is slow, this method is equivalent and fast.
compute_tie_correction = function(ranks) {
    if (is(ranks, "dgCMatrix")) {
        observed_var = colVars_sparse(ranks)
    } else {
        observed_var = matrixStats::colVars(as.matrix(ranks))
    }
    max_var = var(seq_len(nrow(ranks)))
    return((max_var-observed_var) * 12 / nrow(ranks))
}

auroc_p_value = function(aurocs, label_matrix, two_tailed = TRUE, tie_correction = 0, log.p = FALSE) {
    p = colSums(label_matrix)
    n = nrow(label_matrix) - p
  
    # Careful: NAs arise from tie_correction (predictor with 0 variance)
    if (length(tie_correction) > 1) {
        Z = (aurocs - 0.5) * sqrt(12*n*p)
        Z = t(t(Z) / sqrt(nrow(label_matrix)+1-tie_correction))
    } else {
        Z = (aurocs - 0.5) / sqrt((nrow(label_matrix)+1-tie_correction)/(12*n*p))
    }
  
    result = Z
    if (two_tailed) {
        is_not_na = !is.na(Z)
        result[Z<=0 & is_not_na] = pnorm(Z[Z<=0 & is_not_na], log.p = log.p) * 2
        result[Z>0 & is_not_na] = pnorm(Z[Z>0 & is_not_na], lower.tail=FALSE, log.p = log.p) * 2
    } else {
        result = pnorm(Z, lower.tail=FALSE, log.p = log.p)
    }
    return(result)
}

colVars_sparse = function(M) {
    result = (Matrix::colMeans(M**2) - Matrix::colMeans(M)**2)*nrow(M)/(nrow(M)-1)
    return(result)
}

partial_auc = function(predictors, label_matrix, top_f=0.05, seed=17) {
    set.seed(seed)
    label_matrix = as.matrix(label_matrix)
    n_samples = nrow(label_matrix)
    n_positives = colSums(label_matrix)
    n_negatives = n_samples - n_positives
    min_rank = round((1-top_f) * n_samples)
    if (is(predictors, "dgCMatrix")) {
        # we only consider the ranks of non-zero values
        ranks = rank_sparse_no_offset(predictors)
        ranks@x = ranks@x - min_rank
        ranks@x[ranks@x < 0] = 0
        partial_positives = as.matrix(Matrix::crossprod(label_matrix, ranks>0))
        sum_of_positive_ranks = as.matrix(Matrix::crossprod(label_matrix, ranks))
        # check if the fraction of zeros > top_f
        n_zeros = n_samples - diff(predictors@p)
        n_missing = n_zeros - min_rank
        missing_values = which(n_missing > 0)
        if (length(missing_values) > 0) {
            warning("Partial AUC includes 0 values in ", length(missing_values), "/", ncol(predictors), " predictors.")
            n_missing = ifelse(n_missing<0, 0, n_missing)
            # correction: we ignored all 0-values that should contribute a little bit to the AUC (top part)
            #  - on average we add by chance #new = #missing * #remaining positives / # zeros
            new_positives = scale(n_positives-partial_positives, center=FALSE, scale = n_zeros / n_missing)
            #  - the missing AUC is #missing / #negatives * (#new/2) / #positives
            # (we divide by #negatives*#positives later, so we omit it here)
            correction_zeros = scale(new_positives/2, center=FALSE, scale = 1/n_missing)
        }
    } else {
        predictors = as.matrix(predictors)
        ranks = matrixStats::colRanks(predictors, ties.method = "random", preserveShape=TRUE) - min_rank
        ranks[ranks < 0] = 0
        sum_of_positive_ranks = crossprod(label_matrix, ranks)
        colnames(sum_of_positive_ranks) = colnames(predictors)
    }
    result = (sum_of_positive_ranks - partial_positives*(partial_positives-1)/2 + correction_zeros) / (n_negatives*n_positives)
    return(as.matrix(result))
}

rank_sparse_no_offset = function(M) {
    nnz = diff(M@p)
    ranks = tibble(x = M@x, j = rep.int(1:ncol(M), nnz)) %>%
        group_by(j) %>%
        mutate(rank_ = as.numeric(colRanks(as.matrix(x), ties.method = "random")) + nrow(M) - n())
    R = M
    R@x = ranks$rank_
    return(R)
}


