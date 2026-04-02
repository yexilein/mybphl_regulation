
library(tidyverse)
library(Seurat)
library(Signac)


####################################
#              RNA                 #
####################################

find_all_markers_rna = function(srt, group_by, condition=NULL, cells = colnames(srt), min_cells=20, ...) {
    keep_cell = colnames(srt) %in% cells
    if (!is.null(condition)) {
        cond = as.factor(srt[[]][[condition]][keep_cell])
        message("Adjusting for conditional differences according to ", condition,
                " with ", levels(cond)[1], " as the control condition (reorder factor levels to change the control condition). ",
                "For groups with fewer than ", min_cells, " control cells, the next available level will be treated as the control.")
    }
    group = as.factor(srt[[]][[group_by]][keep_cell])
    result = lapply(set_names(levels(group)), function(g) {
        message("Looking for differential expression in ", g, "...")
        find_markers_rna(srt=srt, group_by=group_by, goi=g, condition=condition, cells=cells, min_cells=min_cells, ...)
    })
    result = result[!sapply(result, is.null)]
    if (is.null(condition)) {
        result = bind_rows(result, .id = "group")
    } else {
        result = purrr::transpose(result)
        result$group = bind_rows(result$group, .id = "group")
        result$condition = bind_rows(result$condition, .id = "group")
    }
    return(result)
}

find_markers_rna = function(srt, group_by, goi, assay_name="RNA", ...) {
    # create layer with CP10K info
    srt[[assay_name]] = NormalizeData(srt[[assay_name]], normalization.method="RC")
    result = find_markers(srt=srt, group_by, goi, assay_name=assay_name, ...)
    return(result)
}

find_markers = function(srt, group_by, goi, condition=NULL, cells = colnames(srt), features = NULL, min_cells=20, assay_name = "RNA", family="gaussian", ranked=family=="gaussian", fc_pseudocount=10, max_batch_size=20000, .valid_cond=NULL) {
    if (is.null(features)) {
        features = rownames(srt[[assay_name]])
    }
    
    # check that there are enough cells
    keep_cell = colnames(srt) %in% cells
    is_group = srt[[]][[group_by]][keep_cell] == goi
    if (sum(is_group) < min_cells) {
        message("Not enough cells in group, skipping...")
        return(NULL)
    }
    
    # change reference condition level if not enough cells in the control condition
    if (!is.null(condition)) {
        cond = as.factor(srt[[]][[condition]][keep_cell])
        if (is.null(.valid_cond)) {
            valid_cond = levels(cond)[tabulate(cond[is_group], nlevels(cond)) >= min_cells]
            invalid_cond = levels(cond)[!(levels(cond) %in% valid_cond)]
            if (length(valid_cond) == 0) {
                message("No condition matching required number of cells, skipping...")
                return(NULL)
            }
            message("Control condition: ", valid_cond[1], " (conditions with insufficient cells: ",
                    ifelse(length(invalid_cond)==0, "none", paste(invalid_cond, collapse=", ")), ").")
        } else {
            valid_cond = .valid_cond
            invalid_cond = levels(cond)[!(levels(cond) %in% valid_cond)]
        }
        cond = factor(cond, c(valid_cond, invalid_cond))
    }
    
    # cut into smaller batches if necessary
    if (length(features) > max_batch_size) {
        n_batches = ceiling(length(features) / max_batch_size)
        batch_index = round(seq(0, length(features), length.out = n_batches+1))
        result = lapply(1:n_batches, function(i) {
            find_markers(
                srt=srt, group_by=group_by, goi=goi, condition=condition, cells=cells, features=features[(batch_index[i]+1):batch_index[i+1]],
                min_cells=min_cells, assay_name=assay_name, family=family, ranked=ranked, fc_pseudocount=fc_pseudocount, max_batch_size=max_batch_size, .valid_cond=valid_cond
            )
        })
        if (is.null(condition)) {
            return(bind_rows(result))
        } else {
            result = purrr::transpose(result)
            result$group = bind_rows(result$group)
            result$condition = bind_rows(result$condition)
            return(result)
        }
    }
    
    # dependent variable
    full_expr = LayerData(srt[[assay_name]], "data")[features,which(keep_cell),drop=FALSE]
    full_counts = LayerData(srt[[assay_name]], "counts")[features,which(keep_cell),drop=FALSE]
    full_detection = full_expr > 0
    keep_feature = rowSums(full_detection, na.rm=TRUE)>=min_cells
    if (family=="gaussian") {
        full_counts = full_expr
        if (ranked) {
            # break ties by adding random number? -> penalty against weakly expressed genes
            full_counts = matrixStats::rowRanks(as.matrix(full_counts), ties.method = "average", preserveShape=TRUE) / ncol(full_counts)
        }
    } else if (family=="binomial") {
        # logistic regression -> binarize data
        full_counts = full_detection
    } else if (family=="poisson") {
        full_counts = full_counts
    } else {
        stop("Invalid model family:", family)
    }
    full_expr = full_expr[keep_feature,,drop=FALSE]
    full_counts = full_counts[keep_feature,,drop=FALSE]
    full_detection = full_detection[keep_feature,,drop=FALSE]

    # independent variables 
    is_group = srt[[]][[group_by]][keep_cell] == goi
    if (is.null(condition)) {
        model_matrix = model.matrix(~ is_group)  
    } else {
        model_matrix = model.matrix(~ is_group * cond)
    }
    # gaussian uses normalized data -> no need to add # counts as a covariable
    if (family != "gaussian") {
        log_ncounts = scale(log(srt[[paste0("nCount_", assay_name)]]))
        model_matrix = cbind(model_matrix, model.matrix(~ 0 + log_ncounts))
    }
    
    # run model
    result = apply(full_counts, 1, function(counts) {
        is_not_na = !is.na(1*counts)
        my_glm = fastglm::fastglm(model_matrix[is_not_na,], 1*counts[is_not_na], family = get(family)())
        # stat is stored in 4th column (t value or z value), p-value is stored in 5th column (Pr(>|t|) or Pr(>|z|))
        subresult = as_tibble(summary(my_glm)$coefficients, rownames = "variable") %>%
            mutate(group = sub("^is_group([^:]+).*", "\\1", variable),
                   group = ifelse(group==variable, NA, group),
                   condition = sub("^is_group[^:]+:cond(\\w)", "\\1", variable),
                   condition = ifelse(condition==variable, NA, condition)) %>%
            filter(!is.na(group) | !is.na(condition)) %>%
            select(condition, avg_diff = Estimate, 4) 
        return(subresult)
    })
    result = bind_rows(result, .id = "gene")
    if ("t value" %in% names(result)) {
        result = result %>%
            mutate(log_pval = log(2)+pt(abs(`t value`), ncol(full_counts)-2, lower.tail = FALSE, log.p = TRUE)) %>% # if more accurate p-value is needed
            select(-`t value`)
    } else {
        result = result %>%
            mutate(log_pval = log(2)+pnorm(abs(`z value`), lower.tail = FALSE, log.p = TRUE)) %>% # if more accurate p-value is needed
            select(-`z value`) 
    }
    
    # format results
    group_result = result %>%
        filter(is.na(condition)) %>%
        select(-condition) %>%
        mutate(log_fdr = my_fdr(log_pval, log.p=TRUE),
               fold_change = (MatrixGenerics::rowMeans2(full_expr, cols=which(is_group))+fc_pseudocount/sum(is_group)) /
                               (MatrixGenerics::rowMeans2(full_expr, cols=which(!is_group))+fc_pseudocount/sum(!is_group)),
               n1 = MatrixGenerics::rowSums2(full_detection, cols=which(is_group)),
               n2 = MatrixGenerics::rowSums2(full_detection, cols=which(!is_group)),
               pct.1 = n1/sum(is_group),
               pct.2 = n2/sum(!is_group)) %>%
        arrange(log_pval)
    if (is.null(condition)) {
        return(group_result)
    } else {
        # /!\ first condition is used as reference level -> not included in the results table
        group_result = mutate(group_result, condition = levels(cond)[1], .after = 1)
        pos_table = design_matrix(cond)[,-1,drop=FALSE]
        neg_table = 1 - pos_table
        # only keep cells in the group of interest
        pos_table[!is_group,] = 0
        neg_table[!is_group,] = 0
        cond_result = result %>%
            filter(!is.na(condition)) %>%
            mutate(log_fdr = my_fdr(log_pval, log.p=TRUE),
                   fold_change = as.vector(t(mean_expression(full_expr, pos_table, pseudocount=fc_pseudocount) / mean_expression(full_expr, neg_table, pseudocount=fc_pseudocount))),
                   n1 = as.vector(t(sum_expression(full_detection, pos_table))),
                   n2 = as.vector(t(sum_expression(full_detection, neg_table))),
                   pct.1 = as.vector(t(mean_expression(full_detection, pos_table))),
                   pct.2 = as.vector(t(mean_expression(full_detection, neg_table)))) %>%
            arrange(log_pval) %>%
            filter(condition %in% valid_cond)
        return(list(group = group_result, condition = cond_result))
    }
}

my_fdr = function (p_values, log.p = FALSE) {
    i = length(p_values):1L
    o = order(p_values, decreasing = TRUE)
    n = length(p_values)
    result = rep(0, n)
    if (log.p) {
        result[o] = pmin(0, cummin(log(n) - log(i) + p_values[o]))
    }
    else {
        result[o] = pmin(1, cummin(n/i * p_values[o]))
    }
    return(result)
}

# what’s a reasonable pseudo-count?
# pseudo-count=1 -> detected in 2 vs 0 -> fold_change of 3 (too much)
# pseudo-count=10 -> detected in 20 vs 0 -> fold_change of 3 (OK)
mean_expression = function(expression, design_matrix, pseudocount=0) {
    return(scale(sum_expression(expression, design_matrix, pseudocount), center = FALSE, scale = colSums(design_matrix)))
}

sum_expression = function(expression, design_matrix, pseudocount=0) {
    return(as.matrix(expression %*% design_matrix)+pseudocount)
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

make_volcano = function(markers, fc_label="fold_change", p_label="log_fdr", group_label="group", log2.fc=FALSE, log.p=TRUE, max_log10fdr=Inf, max_log2fc=Inf, topn=10, fdr_threshold=0.05, log2fc_threshold=log2(1.5)) {
    markers$log_fc = if(log2.fc) markers[[fc_label]] else log2(markers[[fc_label]])
    markers$log_fdr = if(log.p) markers[[p_label]]/log(10) else log10(markers[[p_label]])
    lapply(sort(unique(markers[[group_label]])), function(ct) {
        markers %>%
            filter(.data[[group_label]] == ct) %>%
            mutate(log_fdr = ifelse(-log_fdr>max_log10fdr, -max_log10fdr, log_fdr)) %>%
            mutate(log_fc = sign(log_fc)*ifelse(abs(log_fc)>max_log2fc, max_log2fc, abs(log_fc))) %>%
            mutate(is_significant = log_fdr<log10(fdr_threshold) & abs(log_fc) > log2fc_threshold) %>%
            mutate(r_fc = rank(ifelse(is_significant, log_fc, 0)), r_p = rank(ifelse(is_significant, sign(log_fc)*log_fdr, 0))) %>%
            mutate(highlight = (r_fc<=topn | r_fc>(n()-topn) | r_p<=topn | r_p > (n()-topn)) & is_significant) %>%
            ggplot(aes(x=log_fc,y=-log_fdr)) +
            geom_point(aes(col=is_significant), size=0.1, show.legend=FALSE) +
            ggrepel::geom_text_repel(data=. %>% filter(highlight), aes(label=gene), max.overlaps = Inf) +
#            geom_hline(yintercept = -log10(fdr_threshold), linetype="dashed") +
#            geom_vline(xintercept = c(log2fc_threshold, -log2fc_threshold), linetype = "dashed") +
            labs(title=ct, x="log2(FC)", y="-log10(FDR)") +
            scale_color_manual(values=c("gray60","black"))
    })
}


####################################
#              ATAC                #
####################################

find_all_markers_atac = function(srt, group_by, condition=NULL, cells = colnames(srt), min_cells=20, ...) {
    keep_cell = colnames(srt) %in% cells
    if (!is.null(condition)) {
        cond = as.factor(srt[[]][[condition]][keep_cell])
        message("Adjusting for conditional differences according to ", condition,
                " with ", levels(cond)[1], " as the control condition (reorder factor levels to change the control condition). ",
                "For groups with fewer than ", min_cells, " control cells, the next available level will be treated as the control.")
    }
    group = as.factor(srt[[]][[group_by]][keep_cell])
    result = lapply(set_names(levels(group)), function(g) {
        message("Looking for differential accessibility in ", g, "...")
        find_markers_atac(srt=srt, group_by=group_by, goi=g, condition=condition, cells=cells, min_cells=min_cells, ...)
    })
    result = result[!sapply(result, is.null)]
    if (is.null(condition)) {
        result = bind_rows(result, .id = "group")
    } else {
        result = purrr::transpose(result)
        result$group = bind_rows(result$group, .id = "group")
        result$condition = bind_rows(result$condition, .id = "group")
    }
    return(result)
}

find_markers_atac = function(srt, group_by, goi, assay_name="ATAC", family="binomial", ...) {
    result = find_markers(srt, group_by, goi, assay_name=assay_name, family=family, ...)
    return(result)
}


####################################
#            ChromVAR              #
####################################
run_chromvar = function(srt, peak_width=250, assay_name = "ATAC", genome=BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10) {
    # restrict peak width and re-count ~ 10 minutes on VPN
    if (!is.null(peak_width) & peak_width > 0) {
        new_peaks = IRanges::resize(granges(srt[[assay_name]]), width = peak_width, fix = "center")    
        new_matrix = FeatureMatrix(Fragments(srt[[assay_name]]), features = new_peaks, process_n = 10000)
        new_ca = Signac::CreateChromatinAssay(new_matrix, annotation = Annotation(srt[[assay_name]]), fragments = Fragments(srt[[assay_name]]))
    } else {
        new_ca = srt[[assay_name]]
    }    

    # find motifs ~ 2 minutes on VPN
    if (is.null(Motifs(new_ca))) {
        new_ca = AddMotifs(new_ca, pfm = jaspar_pfm(), genome = genome)
    }

    # only use 4 cores to limit memory usage ~ 7 minutes on VPN
    default_param = BiocParallel::registered()[[1]]
    BiocParallel::register(BiocParallel::MulticoreParam(4))
    result = RunChromVAR(new_ca, genome = genome)
    BiocParallel::register(default_param)
    gc()
    return(result)
}

compute_all_dtfa = function(srt, group_by, cells = colnames(srt), ...) {
    keep_cell = colnames(srt) %in% cells
    all_groups = levels(as.factor(srt[[]][[group_by]][keep_cell]))
    result = lapply(set_names(all_groups), function(g) {
        message(paste0("Looking for differential TF activity in ", g, "..."))
        compute_dtfa(srt=srt, group_by=group_by, goi=g, cells=cells, ...)
    })
    result = bind_rows(result, .id = "group")
    return(result)
}

compute_dtfa = function(srt, group_by, goi, cells=colnames(srt), assay_name="chromvar", correct_fragments = TRUE) {
    keep_cell = colnames(srt) %in% cells
    is_group = srt[[]][[group_by]][keep_cell] == goi
    n_fragments = srt$nCount_ATAC[keep_cell]
    tfa = srt[[assay_name]]@data[, cells, drop=FALSE]
    motif_to_tf = unlist(Motifs(srt[["ATAC"]])@motif.names)
    
    my_data = tibble(as_tibble(t(tfa)), n_fragments=n_fragments, group=ifelse(is_group, "goi", "background")) |>
        pivot_longer(c(-n_fragments, -group), names_to = "motif", values_to = "tfa") |>
        drop_na() |>
        group_by(motif) |>
        filter(length(unique(group)) > 1) |>
        ungroup()
    if (correct_fragments) {
        lm_result = my_data %>%
            nest_by(motif) %>%
            reframe(as_tibble(summary(lm(tfa ~ scale(log(n_fragments)) * group, data = data))$coefficients, rownames="variable"))
    } else {
        lm_result = my_data %>%
            nest_by(motif) %>%
            reframe(as_tibble(summary(lm(tfa ~ group, data = data))$coefficients, rownames="variable"))
    }
    result = lm_result %>%
        filter(variable == "groupgoi") %>%
        select(motif, avg_diff=`Estimate`, p_val = `Pr(>|t|)`) %>% # se=`Std. Error`, t=`t value`
        mutate(fdr = p.adjust(p_val, "fdr")) %>%
        arrange(p_val)
    motif_to_intercept = filter(lm_result, variable == "(Intercept)") %>%
        select(motif, Estimate) %>%
        deframe()
    result$mean_background = motif_to_intercept[result$motif]
    result$motif.name = motif_to_tf[result$motif]
    return(result)
}

compute_condition_dtfa = function(srt, coi="mut", condition_name = "condition", assay_name="chromvar") {
    cond = setNames(as.factor(srt[[]][[condition_name]]), colnames(srt))
    condition_markers = lapply(unique(Idents(srt)), function(ct) {
        message("Looking for conditional differences in ", ct)
        ct_cells = WhichCells(srt, idents = ct)
        # only proceed in at least 20 cells per condition
        if (min(tabulate(cond[ct_cells], nbins = nlevels(cond))) < 20) { return(NULL) }
        subresult = compute_dtfa(srt, condition_name, coi, cells = ct_cells, assay_name=assay_name) %>%
            mutate(group = ct)
    })
    condition_markers = bind_rows(condition_markers)
    return(condition_markers)
}
