library(tidyverse)
library(GenomicRanges)
source("dataset.R")
source("grn.R")


main = function() {
    args = commandArgs(trailingOnly=TRUE)
    all_stages = !is.na(args[1]) & args[1] == "TRUE"
    message("Using all stages ", all_stages)
    grn_dir = ifelse(all_stages, "results/grn_all_stages", "results/grn")
    dir.create(grn_dir, showWarnings = FALSE, recursive = TRUE)
    
    # pre-select target genes and regulators of interest
    gene_stats = all_gene_stats(cache_dir = grn_dir, all_stages=all_stages)
    my_targets = select_targets(gene_stats)
    additional_targets = gene_stats$gene[startsWith(gene_stats$gene,"Myh")]
    my_targets = unique(c(my_targets, additional_targets))
    my_regulators = select_regulators(gene_stats)
    my_regulators = intersect(my_regulators, all_tfs())
    
    # pre-select cCREs, ~15 minutes
    ccres = preselect_ccres(my_targets, cache_dir=grn_dir, all_stages=all_stages)
    
    # find motifs in cCREs, ~10 minutes
    motifs = find_motifs(ccres, p.cutoff=1e-3, cache_dir=grn_dir)
    
    # compute TF-cCRE association, ~5 minutes
    cache_file = file.path(grn_dir, "tf_ccre.rds")
    if (!file.exists(cache_file)) {
        message("Identifying regulators (TF -> cCRE)...")
        pb = load_pseudobulk(all_stages = all_stages)
        mc = load_metacells(all_stages = all_stages)
        
        pb_cor = cor_peak_gene(pb, my_regulators, unique(ccres$peak_id), suffix = "_pb")
        mc_cor = cor_peak_gene(mc, my_regulators, unique(ccres$peak_id), suffix = "_mc")
        
        assertthat::are_equal(pb_cor$gene, mc_cor$gene)
        assertthat::are_equal(pb_cor$peak, mc_cor$peak)
        result = pb_cor
        result = cbind(result, mc_cor[,c(-1,-2)])
        result = dplyr::rename(result, tf = gene)
        
        # keep highly correlated TFs or top 100 TFs
        min_cor = 0.05
        nbins=10
        tf_counts = rowSums(SeuratObject::LayerData(pb[["RNA"]], "counts", features = my_regulators))
        tf_bin = as.numeric(cut(tf_counts, unique(quantile(tf_counts, probs = seq(0,1,l=nbins+1))), include.lowest = TRUE))
        tf_bin = setNames(tf_bin, my_regulators)
        result = result |>
            mutate(bg_bin = tf_bin[tf]) |>
            mutate(tf_rank = data.table::frank(-pooled_mc), tf_quantile = tf_rank / n(), .by = "peak") |>
            #drop_na() |>
            mutate(tf_z = (pooled_mc - mean(pooled_mc))/sd(pooled_mc), .by = c("peak","bg_bin")) |>
            filter(tf_rank <= 100 | pooled_mc > min_cor, tf_z > 0) |>
            select(-bg_bin)
        saveRDS(result, cache_file)
    }

    # compate chromVAR-cCRE association, ~ 20 minutes
    cache_file = file.path(grn_dir, "cv_ccre.rds")
    if (!file.exists(cache_file)) {
        message("Identifying regulators from ChromVAR...")
        all_motifs = readRDS("srt/motifs.rds")
        
        pb = load_pseudobulk()
        Motifs(pb[["ATAC"]]) = all_motifs
        pb[["chromvar"]] = run_chromvar(pb, peak_width = NULL)

        mc = load_metacells()
        # rule of thumb ~ 1s per motif (50 iterations, 1 core), ~1.5s per motif (100 iterations)
        # loops over motifs with bplapply and uses matrix multiplication on background peaks
        # -> increased efficiency for high number of background peaks, lower efficiency with more motifs
        # -> it would probably be much more efficient to loop over the background and do the matrix multiplication over motifs
        Motifs(mc[["ATAC"]]) = all_motifs
        # 4 cores, 900 motifs -> 7 min
        mc[["chromvar"]] = run_chromvar(mc, peak_width = NULL, ncores=4)

        
        pb_cor = cor_peak_gene(pb, rownames(pb[["chromvar"]]), unique(ccres$peak_id), suffix = "_chromvar_pb", gene_assay = "chromvar") |>
            dplyr::rename(motif=gene) |>
            add_tf_name(jaspar_pfm()) |>
            dplyr::select(-motif) |>
            dplyr::select(peak, gene=tf, everything())
        
        mc_cor = cor_peak_gene(mc, rownames(mc[["chromvar"]]), unique(ccres$peak_id), suffix = "_chromvar_mc", gene_assay = "chromvar") |>
            dplyr::rename(motif=gene) |>
            add_tf_name(jaspar_pfm()) |>
            dplyr::select(-motif) |>
            dplyr::select(peak, gene=tf, everything())
        
        assertthat::are_equal(pb_cor$gene, mc_cor$gene)
        assertthat::are_equal(pb_cor$peak, mc_cor$peak)
        result = pb_cor
        result = cbind(result, mc_cor[,c(-1,-2)])
        result = dplyr::rename(result, tf = gene)
        
        # only keep positively correlated TFs and top TFs
        min_cor = 0.05
        result = result |>
            group_by(peak) |>
            mutate(chromvar_rank = data.table::frank(-pooled_chromvar_mc),
                   chromvar_quantile = chromvar_rank / n(),
                   chromvar_z = (pooled_chromvar_mc - mean(pooled_chromvar_mc)) / sd(pooled_chromvar_mc)) |>
            ungroup() |>
            filter(chromvar_rank <= 100 | pooled_chromvar_mc > min_cor, chromvar_z > 0)

        saveRDS(result, cache_file)
    }

    # compute TF-target association, ~ 5 minutes
    cache_file = file.path(grn_dir, "tf_target.rds")
    if (!file.exists(cache_file)) {
        message("Identifying regulators (TF -> target)...")
        pb = load_pseudobulk()
        mc = load_metacells()

        pb_cor = cor_peak_gene(pb, my_regulators, my_targets, suffix = "_pb", region_assay="RNA")
        mc_cor = cor_peak_gene(mc, my_regulators, my_targets, suffix = "_mc", region_assay="RNA")
        
        assertthat::are_equal(pb_cor$gene, mc_cor$gene)
        assertthat::are_equal(pb_cor$peak, mc_cor$peak)
        result = pb_cor
        result = cbind(result, mc_cor[,c(-1,-2)])
        result = dplyr::rename(result, tf = gene)
        
        # keep highly correlated TFs or top 100 TFs
        min_cor = 0.05
        nbins=10
        tf_counts = rowSums(LayerData(pb[["RNA"]], "counts", features = my_regulators))
        tf_bin = as.numeric(cut(tf_counts, unique(quantile(tf_counts, probs = seq(0,1,l=nbins+1))), include.lowest = TRUE))
        tf_bin = setNames(tf_bin, my_regulators)
        result = result |>
            mutate(bg_bin = tf_bin[tf]) |>
            mutate(tf_rank = data.table::frank(-pooled_mc), tf_quantile = tf_rank / n(), .by = "peak") |>
            mutate(tf_z = (pooled_mc - mean(pooled_mc))/sd(pooled_mc), .by = c("peak","bg_bin")) |>
            filter(tf_rank <= 100, pooled_mc > min_cor, tf_z > 0) |>
            select(-bg_bin)
        result = dplyr::rename(result, target = peak)
        saveRDS(result, cache_file)
    }
}

all_gene_stats = function(cache_dir="results/grn", all_stages=FALSE) {
    cache_file = file.path(cache_dir,"gene_stats.csv.gz")
    if (!file.exists(cache_file)) {
        result = bind_rows(lapply(set_names(all_grn_tissues(all_stages=all_stages)), function(my_tissue) {
            compute_gene_stats(load_tissue(my_tissue))
        }), .id="tissue")
        result = filter(result, !startsWith(gene, "Gm"), !endsWith(gene, "Rik"))
        write_csv(result, cache_file)
    } else {
        result = read_csv(cache_file, show_col_types = FALSE)
    }
    return(result)
}

preselect_ccres = function(my_targets, cache_dir = "results/grn", all_stages = FALSE) {
    cache_file = file.path(cache_dir, "ccres.csv.gz")
    if (!file.exists(cache_file)) {
        ccres = lapply(set_names(all_grn_tissues(all_stages = all_stages)), function(my_tissue) {
            subcache_file = file.path(cache_dir, glue::glue("ccres_{my_tissue}.csv.gz"))
            if (!file.exists(subcache_file)) {
                message(my_tissue)
                srt = load_tissue(my_tissue)
                t = system.time({ c = find_ccres(srt, my_targets) })
                message(t["elapsed"])
                write_csv(c, subcache_file)
            }
            return(read_csv(subcache_file, show_col_types=FALSE))
        })
        result = bind_rows(ccres, .id="tissue")
        result = add_peak_id(result)
        result = filter_ccres(result)
        write_csv(result, cache_file)
    }
    return(read_csv(cache_file))
}

add_peak_id = function(result) {
    all_peaks = unique(result$region)
    common_peaks = rtracklayer::import("srt/combined_peaks.bed")
    hits = GenomicRanges::findOverlaps(StringToGRanges(all_peaks), common_peaks)
    peak_to_id = tibble(region = all_peaks[queryHits(hits)],
                        peak_id = GRangesToString(common_peaks[subjectHits(hits)]))
    result = left_join(result, peak_to_id)
}

cor_peak_gene = function(srt, genes, regions, suffix="", ...) {
    result = tibble(peak=rep(regions, length(genes)), gene=rep(genes, each=length(regions)))
    result[[glue::glue("pooled{suffix}")]] = c(cor_gene_region(srt, genes, regions, ...))
    gc()
    for (my_tissue in unique(srt$tissue)) {
        result[[glue::glue("{my_tissue}{suffix}")]] = c(cor_gene_region(srt[, srt$tissue == my_tissue], genes, regions, ...))
        gc()
    }
    return(result)
}

if (sys.nframe() == 0) {
    main()
}