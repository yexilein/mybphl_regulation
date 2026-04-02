
library(Seurat)
library(Signac)
library(tidyverse)
source("dataset.R")

main = function() {
    args = commandArgs(trailingOnly=TRUE)
    all_stages = !is.na(args[1]) & args[1] == "TRUE"
    message("Using all stages ", all_stages)
    
    peak_set = rtracklayer::import("srt/combined_peaks.bed")
    all_motifs = readRDS("srt/motifs.rds")
    counts = lapply(all_grn_tissues(all_stages = all_stages), \(tissue_name) {
        srt = load_tissue(tissue_name)
        result = FeatureMatrix(Fragments(srt[["ATAC"]]), features = peak_set, cells = colnames(srt),  process_n = 10000)
        return(result)
    })
    tissue = rep(all_grn_tissues(all_stages = all_stages), sapply(counts, ncol))
    counts = do.call(cbind, counts)
    srt = CreateSeuratObject(CreateChromatinAssay(counts, motifs = all_motifs))
    srt$tissue = tissue

    BiocParallel::register(BiocParallel::MulticoreParam(4))
    # ~ 20min
    system.time({ srt = RunChromVAR(srt, genome = BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10, niterations=50) })
    file_suffix = ifelse(all_stages, "_all_stages", "")
    saveRDS(srt, glue::glue("srt/all_chromvar{file_suffix}.rds"))
}

if (sys.nframe() == 0) {
    main()
}