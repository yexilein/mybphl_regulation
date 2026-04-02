
library(Seurat)
library(matrixStats)
library(GenomicRanges)
source("aurocs.R")


jaspar_pfm = function(filename="data/jaspar2024.rds", force_download=FALSE) {
    if (!file.exists(filename) | force_download) {
        result = TFBSTools::getMatrixSet(
            JASPAR2024::JASPAR2024()@db,
            opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
        )
        saveRDS(result, filename)
    } else {
        result=readRDS(filename)
    }
    return(result)
}

all_tfs = function(cache_file = "data/tf_list.txt") {
    if (file.exists(cache_file)) {
        result = scan(cache_file, what = "string")
    } else {
        library(org.Mm.eg.db)
        library(AnnotationDbi)
        #org.Mm.eg.db
        #columns(org.Mm.eg.db)

        my_goid = "GO:0003700"
        df = AnnotationDbi::select(org.Mm.eg.db, keys = my_goid, columns = c("SYMBOL", "MGI"), keytype = "GOALL")
        result = unique(df$SYMBOL)
        write(result, cache_file)
    } 
    return(result)
}

human_to_mouse = function(gene_list, data_directory="data") {
    hgnc_to_mgi = mouse_to_human_mapping(data_directory) %>% dplyr::select(symbol_human, symbol_mouse) %>% deframe
    result = setNames(hgnc_to_mgi[gene_list], gene_list)
    return(result)
}
5
mouse_to_human_mapping = function(data_directory="data", force_recompute=FALSE) {
    mapping_file = file.path(data_directory, "mouse_to_human.csv")
    if (!file.exists(mapping_file) | force_recompute) {
        date = gsub("-","_",Sys.Date(),fixed = TRUE)
        mgi_file = file.path(data_directory, paste0("mgi_homology_", date,".rpt"))
        download.file("https://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt", mgi_file)
        mgi = read_delim(mgi_file, delim="\t", show_col_types = FALSE)
        one_to_one = mgi %>%
            dplyr::select(db_class_key=`DB Class Key`, organism=`Common Organism Name`, symbol=`Symbol`, entrez_id=`EntrezGene ID`, mgi_id=`Mouse MGI ID`, hgnc_id=`HGNC ID`) %>%
            mutate(organism = ifelse(organism=="mouse, laboratory", "mouse", organism)) %>%
            group_by(db_class_key) %>%
            dplyr::filter(n() == 2 & length(unique(organism)) == 2) %>%
            ungroup() %>%
            pivot_wider(db_class_key, names_from = organism, values_from = c(symbol, entrez_id, mgi_id, hgnc_id)) %>%
            dplyr::select(-mgi_id_human, -hgnc_id_mouse)
        write_csv(one_to_one, mapping_file)
    }
    return(read_csv(mapping_file, show_col_types = FALSE))
}

smooth_data = function(srt, assay_name="RNA", layer="data", graph_name=paste0(assay_name, "_nn"), features=NULL, normalize_by_row=TRUE) {
    expr = SeuratObject::LayerData(srt, assay=assay_name, layer=layer)
    if (!is.null(features)) {
        expr = expr[features,,drop=FALSE]
    }
    knn = Graphs(srt, graph_name)[colnames(expr),colnames(expr)]
    # /!\ not symmetric -> by default, row 1 = 20 neighbors of cell 1, etc.
    # cell = average of its neighbors -> normalize by row (overrepresentation of well-connected cells)
    # cell = propagated equally to its neighbors -> normalize by column (conservation of overall content)
    # each is considered to be its own neighbor (good!)
    knn = as(knn, "dgCMatrix")
    if (normalize_by_row) {
        knn = knn / rowSums(knn)
    } else {
        knn = knn / colSums(knn)
    }
    # S = expr.t(knn) = tcrossprod(expr,knn)
    result = Matrix::tcrossprod(expr,knn)
    return(result)
}

smooth_feature = function(srt, features, assay_name="RNA", layer="data", graph_name=paste0(assay_name, "_nn"), normalize_by_row=TRUE) {
    expr = t(as.matrix(FetchData(srt, features)))
    knn = Graphs(srt, graph_name)[colnames(expr),colnames(expr)]
    # /!\ not symmetric -> by default, row 1 = 20 neighbors of cell 1, etc.
    # cell = average of its neighbors -> normalize by row (overrepresentation of well-connected cells)
    # cell = propagated equally to its neighbors -> normalize by column (conservation of overall content)
    # each is considered to be its own neighbor (good!)
    knn = as(knn, "dgCMatrix")
    if (normalize_by_row) {
        knn = knn / rowSums(knn)
    } else {
        knn = knn / colSums(knn)
    }
    # S = expr.t(knn) = tcrossprod(expr,knn)
    result = Matrix::tcrossprod(expr,knn)
    return(result)
}


#################
# GRN initiation
#################

compute_gene_stats = function(srt, nbins=100) {
    # basic stats: detection + average expression
    pct = rowMeans(SeuratObject::LayerData(srt[["RNA"]])>0)
    avg = rowMeans(SeuratObject::LayerData(srt[["RNA"]]))

    # advanced stats: %variance of detection patterns explained by cell types (ANOVA style)
    n_genes = nrow(srt[["RNA"]])
    all_ct = as.factor(srt$cell_type)
    n_ct = tabulate(all_ct)
    mean_global = pct
    ss_global = MatrixGenerics::rowVars(SeuratObject::LayerData(srt[["RNA"]])>0) * (n_genes-1)
    mean_ct = t(as.matrix(Pando::aggregate_matrix(t(SeuratObject::LayerData(srt[["RNA"]])>0), all_ct, fun="mean")))
    ss_ct = colSums(n_ct*t((mean_ct-mean_global)^2))
    f_stat = (ss_ct/(length(n_ct)-1)) / ((ss_global - ss_ct) / (n_genes - length(n_ct)))
    pval = 1-pf(f_stat, df1 = length(n_ct)-1, df2 = n_genes-length(n_ct))
    
    # compute standardized variance (~quantile-normalization according to overall detection rate)
    my_df = tibble(gene = rownames(srt[["RNA"]]), pct, avg, ss_ct, f_stat, pval) |>
        mutate(pct_bin = as.numeric(cut(pct, breaks = unique(quantile(pct, probs = seq(0,1,l=nbins+1))), include.lowest = TRUE))) |>
        group_by(pct_bin) |>
        mutate(normalized_var = (log(ss_ct) - mean(log(ss_ct)))/sd(log(ss_ct))) |>
        ungroup()
    result = dplyr::select(my_df, gene, avg, pct, normalized_var, f_stat, pval)

    return(result)
}

select_targets = function(gene_stats, min_pct=0.05, min_normalized_var = 1) {
    result = gene_stats |> filter(pct >= min_pct, normalized_var >= min_normalized_var) |> pull(gene) |> unique()
    return(result)
}

select_regulators = function(gene_stats, min_pct=0.05, min_tf_pct = 0.01, tf_list = all_tfs()) {
    result = gene_stats |> filter(pct>min_pct | (gene %in% tf_list & pct > min_tf_pct)) |> pull(gene) |> unique()
    return(result)
}

# methods 
#  - Signac: use linkPeaks (adjust for accessibility, peak width, %GC). /!\ can’t specify upstream and downstream values properly, issue for long genes
#  - custom: adjust for accessibility /!\ seems to favor lowly accessible peaks a bit too much...
#  - Pando: pick all peaks near genes
find_ccres = function(srt, my_targets, method = c("custom","signac","pando"), upstream=1e5, downstream=1e5, rank_=FALSE, region_layer="counts", ...) {
    # TODO: XGB/SCENIC+-type scoring
    message("Identifying cCREs...")
    method = match.arg(method)
    
    # advanced stats
    if (method == "custom") { # ~5min for 2k genes
        peak_score = my_link_peaks(srt, my_targets, upstream=upstream, downstream=downstream, rank_=rank_, region_layer=region_layer, ...)
    } else if (method == "signac") { # ~30min-1h30 for 2k genes
        peak_score = LinkPeaks(srt, peak.assay = "ATAC", expression.assay = "RNA", peak.slot = region_layer,
                               method = ifelse(rank_, "spearman", "pearson"),
                               genes.use = my_targets, pvalue_cutoff=1, score_cutoff=0.01)[["ATAC"]] |>
            Links() |>
            as.data.frame() |>
            dplyr::select(target = gene, region = peak, peak_score = score, peak_zscore = zscore, pval = pvalue)
    } else if (method == "pando") { # ~seconds
        gene_annot = Signac:::CollapseToLongestTranscript(Annotation(srt[["ATAC"]]))
        gene_annot = gene_annot[gene_annot$gene_name %in% my_targets]
        peaks = Pando:::find_peaks_near_genes(granges(srt[["ATAC"]]), genes = gene_annot, upstream=upstream, downstream=downstream)
        peak_score = tibble(target = rep.int(colnames(peaks), diff(peaks@p)), region = rownames(peaks)[peaks@i]) |>
            mutate(peak_score = cor_gene_region_matched(srt, target, region, rank_=rank_, region_layer=region_layer)$rho,
                   peak_zscore = NA, pval = NA)
    } else {
        stop("Unknown cCRE scoring method: ", method)
    }

    # basic stats
    peak_score$pct_peak = rowMeans(SeuratObject::LayerData(srt[["ATAC"]])>0)[peak_score$region]
    #peak_score$avg = rowMeans(SeuratObject::LayerData(srt[["ATAC"]]))[peak_score$region]
    
    # add TSS info
    tss_pos = GetTSSPositions(Annotation(srt[["ATAC"]]), biotypes = NULL)
    tss_pos = tss_pos[match(peak_score$target,tss_pos$gene_name)]
    peak_score$tss_distance = IRanges::distance(StringToGRanges(peak_score$region), tss_pos)

    # add gene info
    gene_pos = Signac:::CollapseToLongestTranscript(Annotation(srt[["ATAC"]]))
    gene_pos = gene_pos[match(peak_score$target,gene_pos$gene_name)]
    peak_score$gene_distance = IRanges::distance(StringToGRanges(peak_score$region), gene_pos)
    return(peak_score)
}

my_link_peaks = function(srt, my_targets, upstream=1e5, downstream=1e5, rank_=rank_, region_layer="counts", ...) {
    # Step 0: make background, 90s for ~2k genes
    bg = make_gene_region_background(srt, my_targets, rank_=rank_, region_layer=region_layer, ...)
    all_peaks = StringToGRanges(names(bg$peak_bin))
    all_targets = unique(bg$bg_matrix$target)
    gene_pos = Signac:::CollapseToLongestTranscript(Annotation(srt[["ATAC"]]))
    gene_pos = gene_pos[match(all_targets, gene_pos$gene_name)]
    
    # Step 1: peak to gene association, 40s for ~2k genes
    peaks = Pando:::find_peaks_near_genes(all_peaks, genes = gene_pos, upstream = upstream, downstream = downstream)
    result = tibble(target = rep.int(colnames(peaks), diff(peaks@p)), region = rownames(peaks)[peaks@i+1]) |>
        mutate(pbin = bg$peak_bin[region])

    # Step 2: compute correlation, 4min for ~2k genes
    result = result |>
        mutate(peak_score = cor_gene_region_matched(srt, target, region, rank_=rank_, region_layer=region_layer)$rho) |>
        left_join(bg$bg_matrix, by=join_by(pbin,target)) |>
        mutate(peak_zscore = (peak_score-mu)/sd, pval = pnorm(-abs(peak_zscore)), fdr = p.adjust(pval, "fdr")) |>
        dplyr::select(target, region, peak_score, peak_zscore, pval, fdr)

    return(result)
}

make_gene_region_background = function(srt, my_targets, min_cells=10, n_background=1000, bg_seed=17, bg_pbin=20, rank_=FALSE, region_layer="counts") {
    peak_counts = rowSums(SeuratObject::LayerData(srt[["ATAC"]])>0)
    all_peak_names = names(peak_counts)[peak_counts>min_cells]
    peak_counts = peak_counts[all_peak_names]
    all_peaks = StringToGRanges(all_peak_names)
    gene_counts = rowSums(SeuratObject::LayerData(srt[["RNA"]], features=my_targets)>0)
    all_targets = names(gene_counts)[gene_counts>min_cells]
    
    gene_pos = Signac:::CollapseToLongestTranscript(Annotation(srt[["ATAC"]]))
    gene_pos = gene_pos[match(all_targets, gene_pos$gene_name)]
    
    # Background: ~90s with 2k genes and default parameters (1k peaks / bin)
    n_cells = ncol(srt)
    all_peaks$peak_bin = as.numeric(cut(peak_counts, unique(quantile(peak_counts, probs = seq(0,1,l=bg_pbin+1))), include.lowest = TRUE))
    peak_bin = setNames(all_peaks$peak_bin, GRangesToString(all_peaks))
    set.seed(bg_seed)
    bg_matrix = bind_rows(lapply(seq_along(unique(all_peaks$peak_bin)), function(my_peak_bin) {
        bg_peaks = sample(all_peaks[all_peaks$peak_bin == my_peak_bin], n_background)
        bg_rho = cor_gene_region(srt, all_targets, GRangesToString(bg_peaks), rank_=rank_, region_layer=region_layer)
        is_same_chr = sapply(as.character(seqnames(gene_pos)), function(chr) { as.vector(seqnames(bg_peaks) == chr) })
        bg_rho[is_same_chr] = NA
        return(tibble(
            target = all_targets,
            pbin = my_peak_bin,
            mu = colMeans(bg_rho, na.rm=TRUE),
            sd = colSds(bg_rho, na.rm=TRUE)
        ))
    }))

    return(list(bg_matrix=bg_matrix, peak_bin = peak_bin))
}

make_gene_region_background_old = function(srt, min_cells=10, n_background=1000, bg_seed=17, bg_pbin=20, bg_gbin=10) {
    peak_counts = rowSums(SeuratObject::LayerData(srt[["ATAC"]])>0)
    all_peak_names = names(peak_counts)[peak_counts>min_cells]
    peak_counts = peak_counts[all_peak_names]
    all_peaks = StringToGRanges(all_peak_names)
    gene_counts = rowSums(SeuratObject::LayerData(srt[["RNA"]])>0)
    all_genes = names(gene_counts)[gene_counts>min_cells]
    gene_counts = gene_counts[all_genes]
    gene_pos = Signac:::CollapseToLongestTranscript(Annotation(srt[["ATAC"]]))
    gene_pos = gene_pos[match(all_genes, gene_pos$gene_name)]
    
    # Background: ~90s with 2k genes and default parameters (1k peaks / bin)
    n_cells = ncol(srt)
    all_peaks$peak_bin = as.numeric(cut(peak_counts, unique(quantile(peak_counts, probs = seq(0,1,l=bg_pbin+1))), include.lowest = TRUE))
    peak_bin = setNames(all_peaks$peak_bin, GRangesToString(all_peaks))
    gene_bin = as.numeric(cut(gene_counts, unique(quantile(gene_counts, probs = seq(0,1,l=bg_gbin+1))), include.lowest = TRUE))
    gene_bin = setNames(gene_bin, all_genes)
    set.seed(bg_seed)
    bg_matrix = bind_rows(lapply(seq_along(unique(all_peaks$peak_bin)), function(my_peak_bin) {
        bg_peaks = sample(all_peaks[all_peaks$peak_bin == my_peak_bin], n_background)
        bg_rho = cor_gene_region(srt, all_genes, GRangesToString(bg_peaks))
        is_same_chr = sapply(as.character(seqnames(gene_pos)), function(chr) { as.vector(seqnames(bg_peaks) == chr) })
        bg_rho[is_same_chr] = NA
        return(tibble(
            gbin = unique(gene_bin),
            pbin = my_peak_bin,
            mu = sapply(gbin, function(i) { mean(bg_rho[,gene_bin == i], na.rm=TRUE) }),
            sd = sapply(gbin, function(i) { sd(bg_rho[,gene_bin == i], na.rm=TRUE) })
        ))
    }))

    return(list(gene_bin = gene_bin, peak_bin = peak_bin, bg_matrix=bg_matrix))
}

# /!\ this function take as input a table of *matched* genes and regions
# e.g., row 1 = gene_1 region_1 -> output = cor(gene_1, region_1)
cor_gene_region_matched = function(srt, genes, regions, rank_=TRUE, region_layer="data") {
    g_to_r = split(regions, factor(genes, unique(genes)))
    g = SeuratObject::LayerData(srt[["RNA"]])[names(g_to_r),,drop=FALSE]
    r = SeuratObject::LayerData(srt[["ATAC"]], layer=region_layer)[unique(regions),,drop=FALSE]
    # performance: SeuratObject::LayerData >> [,,drop=FALSE] >> my_tscale > crossprod
    # performance: accessing the data is very costly... not clear why. Even SeuratObject::LayerData is very expensive -> lots of overhead
    my_cor = bind_rows(lapply(seq_along(g_to_r), function(i) {
        scaled_rna = my_tscale(g[names(g_to_r)[i],,drop=FALSE], rank_=rank_)
        scaled_atac = my_tscale(r[g_to_r[[i]],,drop=FALSE], rank_ = rank_)
        tibble(gene=names(g_to_r)[i], region = g_to_r[[i]], rho = as.vector(crossprod(scaled_rna, scaled_atac)))
    }))
    # reorder rows to match input order
    row_order = match(paste0(genes,regions), paste0(my_cor$gene, my_cor$region))
    my_cor = my_cor[row_order,]
    return(my_cor)
}

# /!\ this function take as input a table of *unmatched* genes and regions
# e.g., (g1, ..., gN) (r1, ..., rM) -> output: all possible combinations g1xr1, g1xr2, ...
cor_gene_region = function(srt, genes, regions, rank_=TRUE, gene_assay="RNA", region_assay="ATAC", region_layer="data") {
    scaled_rna = my_tscale(SeuratObject::LayerData(srt[[gene_assay]])[genes,,drop=FALSE], rank_=rank_)
    scaled_atac = my_tscale(SeuratObject::LayerData(srt[[region_assay]], layer=region_layer)[regions,,drop=FALSE], rank_ = rank_)
    my_cor = crossprod(scaled_atac, scaled_rna)
    return(my_cor)
}

# cor(x,y) = crossprod(scale(x)/sqrt(nrow(x)-1), scale(y)/sqrt(nrow(y)-1))
my_tscale = function(x, rank_=FALSE) {
    if (rank_) {
        result = scale(t(rowRanks(as.matrix(x), useNames = TRUE, ties.method = "average")))
    } else {
        result = scale(t(as.matrix(x)))
    }
    result = result / sqrt(nrow(result)-1)
    return(result)
}

filter_ccres = function(ccres, topn=20, min_zscore=1.5, min_pct=0.01, rank_by="score") {
    average_score = ccres |>
        group_by(target, peak_id) |>
        filter(abs(peak_zscore) > min_zscore, pct_peak>min_pct) |>
        summarize(n=n(),
                  zscore = mean(peak_zscore, na.rm = TRUE),
                  score = mean(peak_score),
                  tss_distance = mean(tss_distance),
                  gene_distance = mean(gene_distance),
                  pct_peak = exp(mean(log(pct_peak))),
                  tissue = paste(tissue, collapse = ",")) |>
        ungroup()
    pos_ccre = average_score |> filter(zscore>min_zscore) |> slice_max(.data[[rank_by]], by=target, n = topn)
    neg_ccre = average_score |> filter(zscore< -min_zscore) |> slice_max(.data[[rank_by]], by=target, n = topn)
    return(bind_rows(pos_ccre, neg_ccre))
}

find_motifs = function(ccres, pfm=jaspar_pfm(), p.cutoff=5e-4, genome = BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10, cache_dir = NULL, filename="motifs.csv.gz") {
    if (is.null(cache_dir)) {
        cache_file = NULL
    } else {
        cache_file = file.path(cache_dir, filename)
    }
    if (is.null(cache_file) || !file.exists(cache_file)) {
        message("Finding motifs...")
        all_ccres = StringToGRanges(unique(ccres$peak_id))
        result = motifmatchr::matchMotifs(pwms = pfm, subject = all_ccres, out="scores", genome = genome, bg="genome", p.cutoff=p.cutoff)
        result = motifmatchr::motifScores(result)
        result = tibble(region = GRangesToString(all_ccres)[result@i+1],
                        motif = rep.int(colnames(result), diff(result@p)),
                        motif_score = result@x)
        ## alternative find all positions to have the precise position of the motif
        #result = motifmatchr::matchMotifs(pwms = pfm, subject = all_ccres, out="positions", genome = genome, bg="genome", p.cutoff=p.cutoff) 
        #result = find_best_motifs(result, all_ccres)
        result = add_motif_pval(result, cache_dir = cache_dir)
        result = add_tf_name(result, pfm)
        if (!is.null(cache_file)) {
            write_csv(result, cache_file)
        }
    } else {
        result= read_csv(cache_file, show_col_types = FALSE)
    }
    return(result)
}

find_best_motifs = function(motif_score, my_query) {
    region_name = GRangesToString(my_query)
    result = bind_rows(lapply(motif_score, function(ms) {
        motif_to_peak = IRanges::findOverlaps(ms, my_query)    
        result = tibble(motif_position = GRangesToString(ms),
                        motif_score = ms$score,
                        region = region_name[subjectHits(motif_to_peak)]) |>
            group_by(region) |>
            slice_max(motif_score, n=1, with_ties = FALSE) |>
            ungroup()
        return(result)
    }), .id = "motif")
    return(result)
}

add_motif_pval = function(motif_score, cache_dir = NULL) {
    motif_pval_threshold = find_motif_pval(cache_dir = cache_dir)
    result = bind_rows(lapply(unique(motif_score$motif), function(m) {
        score_to_pval = filter(motif_pval_threshold, motif == m) |>
            arrange(threshold) |>
            group_by(threshold) |> 
            slice_min(pval) |>
            ungroup()
        subresult = filter(motif_score, motif == m) |>
            mutate(motif_pval = cut(motif_score, breaks = c(score_to_pval$threshold, Inf), include.lowest = TRUE, labels = score_to_pval$pval)) |>
            mutate(motif_pval = as.numeric(as.character(motif_pval)))
        return(subresult)
    }))
    return(result)
}

find_motif_pval = function(all_pfm = jaspar_pfm(), cache_dir=NULL) {
    if (is.null(cache_dir)) {
        cache_file = NULL
    } else {
        cache_file = file.path(cache_dir, "motif_pval_threshold.csv")
    }
    if (is.null(cache_file) || !file.exists(cache_file)) {
        # see matchMotifs, matchMotifs_helper code
        pwms <- do.call(TFBSTools::PWMatrixList, lapply(all_pfm, TFBSTools::toPWM))
        bg <- motifmatchr:::get_bg(bg_method="genome", subject=NULL, genome=BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10)
        motif_mats <- motifmatchr:::convert_pwms(pwms, bg)
        
        # default p-value threshold 5e-5
        my_pval = sort(c(c(1:9)*10**(-3), c(1:9)*10**(-4), c(1:9)*10**(-5)))
        system.time({
        pval_thresholds = lapply(my_pval, motifmatchr:::get_thresholds, mats = motif_mats, nuc_freqs = bg)
        })
        names(pval_thresholds) = as.character(my_pval)
        thr_forward = lapply(pval_thresholds, function(p) { p[1:length(motif_mats)] })
        thr_reverse = lapply(pval_thresholds, function(p) { p[(length(motif_mats)+1):(2*length(motif_mats))] })
        
        # in our case forward and reverse have the same values
        result = bind_rows(lapply(thr_forward, function(t) { tibble(motif = names(all_pfm), threshold = t) }), .id = "pval")
        result$pval = as.numeric(result$pval)
        if (!is.null(cache_file)) {
            write_csv(result, cache_file)
        }
    } else {
        result = read_csv(cache_file, show_col_types = FALSE)
    }
    return(result)
}

add_tf_name = function(result, pfm, data_directory="data") {
    motif2tf = tibble(motif=names(pfm), tf=sapply(pfm, function(m) { m@name })) %>%
            mutate(tf_mouse = human_to_mouse(tf, data_directory)) %>%
            mutate(tf = ifelse(is.na(tf_mouse), tf, tf_mouse)) %>%
            dplyr::select(motif, tf)
    return(left_join(result, motif2tf))
}

run_chromvar = function(srt, peak_width=250, assay_name = "ATAC", genome=BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10, ncores=4, ...) {
    # restrict peak width and re-count ~ 10 minutes on VPN
    if (!is.null(peak_width) && peak_width > 0) {
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
    BiocParallel::register(BiocParallel::MulticoreParam(ncores))
    result = RunChromVAR(new_ca, genome = genome, ...)
    BiocParallel::register(default_param)
    gc()
    return(result)
}

#####################
# Regulon stats
####################
# output: cell-level
cell_target_auc = function(srt, target_list, gene_universe=rownames(srt), partial_auc=1) {
    gene_universe = gene_universe[gene_universe %in% rownames(srt)]
    regulon_matrix = list_to_matrix(target_list, gene_universe)
    if (partial_auc < 1) {
        result = partial_auc(GetAssayData(srt)[gene_universe, ], regulon_matrix, top_f = partial_auc)
    } else {
        result = compute_aurocs(GetAssayData(srt)[gene_universe, ], regulon_matrix)   
    }
    return(result)
}

# output: cell-level
cell_enhancer_auc = function(sgc, enhancer_list, assay_name = "eregulon") {
    #if (is.null(GetAssay(sgc)@bias)) {
    #    sgc = InsertionBias(sgc, genome=BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10)
    #}
    #background_peaks = chromVAR::getBackgroundPeaks(GetAssayData(sgc))
    #result = chromVAR::computeDeviations()
    # convert e-regulons to matrix format (peak x regulon matrix)
    regulon_matrix = list_to_matrix(enhancer_list, rownames(sgc))
    default_param = BiocParallel::registered()[[1]]
    BiocParallel::register(BiocParallel::MulticoreParam(4))
    sgc = RunChromVAR(sgc, genome=BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10, motif.matrix=regulon_matrix, new.assay.name=assay_name, verbose=FALSE)
    BiocParallel::register(default_param)
    gc()
    return(sgc)
}

list_to_matrix = function(feature_list, all_features=NULL) {
    if (!is.list(feature_list) | is.null(names(feature_list)) | any(duplicated(names(feature_list)))) {
        stop("Feature list must be a named list with unique names!")
    }
    features = factor(unlist(feature_list), all_features)
    col_names = factor(rep(names(feature_list), sapply(feature_list, length)), names(feature_list))
    result = Matrix::sparseMatrix(i=as.numeric(features),j=as.numeric(col_names),
                                  dims = c(length(levels(features)), length(levels(col_names))),
                                  dimnames = list(levels(features), levels(col_names)))
    return(result)
}

## expand GRN
find_regulators = function(srt, target, extend_by=5e5, smooth = FALSE) {
    my_range = Signac:::FindRegion(srt, assay="ATAC", target)
    start(my_range) = start(my_range) - extend_by
    end(my_range) = end(my_range) + extend_by
    
    putative_cre = rownames(srt[["ATAC"]])[overlapsAny(StringToGRanges(rownames(srt[["ATAC"]])), my_range)]
    cre_to_tf = region_to_tf(srt, putative_cre)
    regulations = tibble(
        target = target, region=rownames(cre_to_tf)[cre_to_tf@i+1], tf=rep.int(colnames(cre_to_tf), diff(cre_to_tf@p))
    )
    if (smooth) {
        rna = t(smooth_data(srt, features = unique(c(regulations$target, regulations$tf))))
        atac = t(smooth_data(srt, features = unique(regulations$region), assay_name = "ATAC", graph_name = "RNA_nn"))
    } else {
        rna = t(SeuratObject::LayerData(srt[["RNA"]]))
        atac = t(SeuratObject::LayerData(srt[["ATAC"]]))
    }
    result = score_regulations(regulations, rna, atac)

    # add TSS info
    tss_pos = GetTSSPositions(Annotation(srt[["ATAC"]]))
    tss_pos = tss_pos[tss_pos$gene_name == target]
    result$tss_distance = IRanges::distance(StringToGRanges(result$region), tss_pos)

    # add Signac score for peaks
    peak_score = Links(LinkPeaks(srt, peak.assay = "ATAC", expression.assay = "RNA", genes.use = target, distance = extend_by, pvalue_cutoff=1, score_cutoff=-1)[["ATAC"]]) |>
        as.data.frame() |>
        select(region = peak, peak_score = score, peak_zscore = zscore)
    result = left_join(result, peak_score)
        
    return(result)
}

score_regulations = function(regulations, rna, atac) {
    # gather initial candidate regulatory elements and TFs
    all_regions = unique(regulations$region)
    all_tfs = unique(regulations$tf)
    all_targets = unique(regulations$target)
    scaled_tf = scale(as.matrix(rna[,all_tfs,drop=FALSE]))
    scaled_region = scale(as.matrix(atac[,all_regions,drop=FALSE]))
    scaled_target = scale(as.matrix(rna[,all_targets,drop=FALSE]))
    
    # compute correlation with target and filter
    tf_target = (crossprod(scaled_tf, scaled_target) / (nrow(rna)-1)) %>%
        as_tibble(rownames = "tf") %>%
        pivot_longer(-tf, names_to = "target", values_to = "tf_target")
    region_target = (crossprod(scaled_region, scaled_target) / (nrow(rna)-1)) %>%
        as_tibble(rownames = "region") %>%
        pivot_longer(-region, names_to = "target", values_to = "region_target")
    tf_region = (crossprod(scaled_tf, scaled_region) / (nrow(rna)-1)) %>%
        as_tibble(rownames = "tf") %>%
        pivot_longer(-tf, names_to = "region", values_to = "tf_region")
    
    # create result table
    result = select(regulations, tf, region, target) |>
        left_join(tf_target) |>
        left_join(region_target) |>
        left_join(tf_region)

    # compute 3-way "correlation"
    scaled_tf_region = scale(as.matrix(rna[,result$tf,drop=FALSE])*as.matrix(atac[,result$region,drop=FALSE]))
    result$tf_region_target = as.vector(crossprod(scaled_tf_region, scaled_target) / (nrow(rna)-1))
    return(result)
}

region_to_tf = function(srt, regions, keep_known_tfs =TRUE) {
    result = GetMotifData(srt[["ATAC"]])[unique(regions),, drop=FALSE]
    colnames(result) = tf_names(srt)
    if (keep_known_tfs) {
        result = result[, colnames(result) %in% rownames(srt[["RNA"]])]
    }
    return(result)
}

tf_names = function(srt) {
    tf_names = ConvertMotifID(srt[["ATAC"]], id = colnames(Motifs(srt[["ATAC"]])))
    tf_names = human_to_mouse(tf_names)
    tf_names = ifelse(is.na(tf_names), names(tf_names), tf_names)
    return(tf_names)
}

