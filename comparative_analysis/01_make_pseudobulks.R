
library(Seurat)
library(Signac)
library(tidyverse)
source("dataset.R")
source("grn.R")


main = function() {
    if (!file.exists("srt/combined_peaks.bed")) {
        rtracklayer::export(make_peak_set(), "srt/combined_peaks.bed")
    }
    peak_set = rtracklayer::import("srt/combined_peaks.bed")

    # find standard motifs
    cache_file = "srt/motifs.rds"
    if (file.exists(cache_file)) {
        all_motifs = readRDS(cache_file)
    } else {
        peak_set = rtracklayer::import("srt/combined_peaks.bed")
        all_motifs = CreateMotifMatrix(peak_set, pwm = jaspar_pfm(), genome=BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10)
        all_motifs = CreateMotifObject(all_motifs)
        saveRDS(all_motifs, cache_file)
    }

    # 10-25min per dataset
    make_pseudobulks("eom", peak_set, all_motifs, mc_resolution=10) # 54 cells / metacell
    make_pseudobulks("limb", peak_set, all_motifs, mc_resolution=5) # 62 cells / metacell
    make_pseudobulks("limb_all_stages", peak_set, all_motifs, mc_resolution=5) # 60 cells / metacell
    make_pseudobulks("larynx", peak_set, all_motifs, mc_resolution=10) # 72 cells / metacell
    make_pseudobulks("heart", peak_set, all_motifs, mc_resolution=5) # 61 cells / metacell
}

make_peak_set = function() {
    peak_sets = lapply(all_tissues(), function(t) {
        peak_set_name = if(t == "eom") { "adult_peaks.bed" } else { "combined_peaks.bed" }
        message(t, peak_set_name)
        load_peak_set(glue::glue("{t}_multiome"), peak_set_name)
    })
    result = Reduce(combine_peaks, peak_sets)
    result = filter_standard_chr(result)
    return(result)
}

combine_peaks = function(peaks1, peaks2, min_width=20, max_width=10000) {
    # Create a unified set of peaks to quantify in each dataset
    combined_peaks = IRanges::reduce(x = c(peaks1, peaks2))
    # Filter out bad peaks based on length
    peakwidths = IRanges::width(combined_peaks)
    combined_peaks = combined_peaks[peakwidths < max_width & peakwidths > min_width]
    return(combined_peaks)
}

filter_standard_chr = function(peak_set) {
    keep_peak = as.vector(seqnames(peak_set)) %in% GenomeInfoDb::standardChromosomes(peak_set)
    peak_set = peak_set[keep_peak, ]
    return(peak_set)
}

make_pseudobulks = function(tissue_name, peak_set, all_motifs, mc_resolution = 10,
                            pb_cache = glue::glue("srt/{tissue_name}_pseudobulk.rds"),
                            mc_cache = glue::glue("srt/{tissue_name}_metacell.rds")) {
    if (file.exists(pb_cache) & file.exists(mc_cache)) {
        return()
    }
    
    srt = load_tissue(tissue_name)
    # recompute ATAC counts (~5-10 minutes)
    atac = FeatureMatrix(Fragments(srt[["ATAC"]]), features = peak_set, cells = colnames(srt),  process_n = 10000)
    srt[["ATAC"]] = CreateChromatinAssay(atac, annotation=Annotation(srt[["ATAC"]]), motifs = all_motifs)

    ## compute pseudobulk
    if (!file.exists(pb_cache)) {
        pb = AggregateExpression(srt, assays = c("RNA","ATAC"), return.seurat = TRUE)
        pb[["ATAC"]] = CreateChromatinAssay(counts = LayerData(pb[["ATAC"]], "counts"), annotation = Annotation(srt[["ATAC"]]))
        pb = RunTFIDF(pb, assay = "ATAC")
        saveRDS(pb, pb_cache)
    }

    ## compute metacells -> target ~50-70 cells/metacell
    if (!file.exists(mc_cache)) {
        mc = lapply(unique(srt$orig.ident), function(oi) {
            make_metacells(srt[, srt$orig.ident == oi], resolution = mc_resolution)
        })
        if (length(mc) > 1) {
            mc = merge(mc[[1]], mc[-1])
        } else {
            mc = mc[[1]]
        }
        message("Average # cells / metacell: ", round(mean(mc$n_cells)))
        saveRDS(mc, mc_cache)
    }
}

make_metacells = function(srt, resolution) {
    DefaultAssay(srt) = "RNA"
    srt = FindVariableFeatures(srt)
    srt = ScaleData(srt)
    srt = RunPCA(srt, npcs = 30)
    srt = FindNeighbors(srt, dims = 1:30)
    srt = FindClusters(srt, resolution = resolution)
    mc = AggregateExpression(srt, assays = c("RNA","ATAC"), return.seurat = TRUE)
    mc$label = vote_label(srt$cell_type, srt$seurat_clusters)[colnames(mc)]
    mc$n_cells = as.vector(table(Idents(srt))[gsub("^g","",colnames(mc))])
    mc[["ATAC"]] = CreateChromatinAssay(counts = LayerData(mc[["ATAC"]], "counts"), annotation = Annotation(srt[["ATAC"]]))
    mc = RunTFIDF(mc, assay = "ATAC")
    return(mc)
}

vote_label = function(ct_label, mc_id) {
    annotation = tibble(label = ct_label, metacell = mc_id) |>
            dplyr::count(metacell, label) |>
            slice_max(n, n=1,with_ties = FALSE, by = metacell) |>
            transmute(metacell = paste0("g",metacell), label) |>
            deframe()
    return(annotation)
}

if (sys.nframe() == 0) {
    main()
}