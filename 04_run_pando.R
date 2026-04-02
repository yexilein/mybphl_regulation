# basic pipeline: https://quadbiolab.github.io/Pando/articles/getting_started.html
# initiating GRN using Seurat/Signac: https://quadbiolab.github.io/Pando/articles/regions.html
# inferring GRN using predefined region-target associations: https://quadbiolab.github.io/Pando/articles/association.html
# presentation of models that can be used to infer the GRN: https://quadbiolab.github.io/Pando/articles/models.html

library(tidyverse)
library(Pando)
library(Signac)
library(Matrix)
source("grn.R")
source("dataset.R")


main = function() {
    ## Parameters
    # command-line: run_pando.R <RESOLUTION> <TISSUE="all"> <ALL_STAGES=FALSE>
    args = commandArgs(trailingOnly=TRUE)
    ct_resolution = args[1] # "single_cell", "metacell", "pseudobulk"
    my_tissue = ifelse(is.na(args[2]), "all", args[2])
    all_stages = !is.na(args[3]) & args[3] == "TRUE"
    # default parameters
    tf_cor = ifelse(ct_resolution=="single_cell", 0.01, 0.05) # minimum TF-target correlation
    stage_suffix = ifelse(all_stages, "_all_stages", "")
    grn_dir = glue::glue("results/grn{stage_suffix}")
    pando_dir = glue::glue("results/pando{stage_suffix}")
    dir.create(pando_dir, FALSE, TRUE)
    message(glue::glue("Running Pando, tissue: {my_tissue}, resolution: {ct_resolution}, min cor: {tf_cor}, all stages {all_stages}."))
    
    ## GRN initiation
    message("Loading and preparing data...")
    # target <-> cCRE
    ccre_table = read_csv(glue::glue("{grn_dir}/ccres.csv.gz"), show_col_types = FALSE)
    #  /!\ in Pando, gene names are treated as variables -> fails for genes with weird names
    ccre_table = filter(ccre_table, target == make.names(target))
    targets = unique(ccre_table$target)
    ccres = StringToGRanges(unique(ccre_table$peak_id))
    # /!\ Pando requires the target name in a metadata column called "gene_name"
    ccre_to_target = StringToGRanges(ccre_table$peak_id)
    ccre_to_target$gene_name = ccre_table$target 
    # TF <-> cCRE, keeping only TFs with a motif_score > 0
    motifs = read_csv(glue::glue("{grn_dir}/motifs.csv.gz"), show_col_types = FALSE) |>
        dplyr::filter(motif_score > 0) |>
        dplyr::select(region, tf)
    tf_ccre = readRDS(glue::glue("{grn_dir}/tf_ccre.rds")) |>
        dplyr::select(region=peak, tf) |>
        inner_join(motifs)
    # multi-ome fdata
    srt = load_pando_data(my_tissue, ct_resolution, ccres, all_stages = all_stages)
    
    ## Pando pipeline
    message(glue::glue("Running Pando on {ncol(srt)} cells, {length(unique(tf_ccre$tf))} TF motifs and {length(targets)} genes."))
    message("Initiating GRN...")
    grn = initiate_pando(srt, ccres, tf_ccre)
    message("Inferring GRN...")
    grn = infer_grn(grn, method = "xgb", genes = targets, peak_to_gene_domains = ccre_to_target, parallel = FALSE, verbose=2, tf_cor = tf_cor)
    message("Writing results to file...")
    saveRDS(GetGRN(grn), glue::glue("{pando_dir}/{ct_resolution}_{my_tissue}.rds"))
    message("Done!")
}

load_pando_data = function(my_tissue, ct_resolution, ccres, all_stages=all_stages) {
    if (my_tissue == "all") {
        if (ct_resolution == "single_cell") {
            srt = lapply(all_grn_tissues(all_stages), function(t) {
                my_srt = load_tissue(t)
                atac = FeatureMatrix(Fragments(my_srt[["ATAC"]]), features = ccres, cells = colnames(my_srt))
                my_srt[["ATAC"]] = CreateChromatinAssay(atac, annotation = Annotation(my_srt[["ATAC"]]))
                return(my_srt)
            })
        } else {
            srt = lapply(all_grn_tissues(all_stages), function(t) { load_data(t, ct_resolution) })
        }
        srt = merge(srt[[1]], srt[2:length(srt)], add.cell.ids = setdiff(all_tissues(), "heart"))
        srt = JoinLayers(srt)
    } else {
        if (ct_resolution == "single_cell") {
            srt = load_tissue(my_tissue)
            atac = FeatureMatrix(Fragments(srt[["ATAC"]]), features = ccres, cells = colnames(srt))
            srt[["ATAC"]] = CreateChromatinAssay(atac, annotation = Annotation(srt[["ATAC"]]))
        } else {
            srt = load_data(my_tissue, ct_resolution)
        }
    }
    srt = RunTFIDF(srt, assay="ATAC")
    return(srt)
}

initiate_pando = function(srt, ccres, tf_ccre) {
    # mostly nothing happens, *except* that regions overlapping exons are trimmed!
    grn = initiate_grn(srt, peak_assay = "ATAC", regions = ccres, exclude_exons = FALSE) # ~2s
    
    ## option 1: Pando recommended method
    # JASPAR: ~ 5 min for 746 motifs on laptop, ~10 min on cluster
    # JASPAR-CIS: ~ 1h30 for 1,590 motifs on cluster
    # this is a simple wrapper around AddMotifs in Signac -> better do it ourselves
    #grn = find_motifs(grn, pfm = motifs, genome = BSgenome.Mmusculus.UCSC.mm10, motif_tfs = my_motif2tf)
    
    ## option 2: manual motif search
    #cand_ranges = grn@grn@regions@ranges
    #motifs = find_motifs(cand_ranges, p.cutoff=5e-4, cache_file="results/pando/motifs.csv.gz")
    #motifs = dplyr::filter(motifs, motif_score > 10)
    #ccre_to_tf = vec_to_mat(motifs$region, motifs$tf)
    #grn@grn@regions@motifs = CreateMotifObject(ccre_to_tf)
    
    ## option 3: arbitrary link between TFs and cCREs from GRN initiation
    ccre_to_tf = vec_to_mat(tf_ccre$region, tf_ccre$tf)
    # reorder cCREs to match ordering in Pando
    # /!\ Pando may subset some of the peaks -> map from orignal cCREs to Pando cCREs
    ccre_to_pando = findOverlaps(grn@grn@regions@ranges, ccres)
    ccre_to_tf = ccre_to_tf[GRangesToString(ccres)[subjectHits(ccre_to_pando)],]
    grn@grn@regions@motifs = CreateMotifObject(ccre_to_tf)
    tfs = colnames(ccre_to_tf)
    grn@grn@regions@tfs = sparseMatrix(i = seq_along(tfs), j = seq_along(tfs), x=1, dimnames = list(tfs, tfs))
    
    return(grn)
}

vec_to_mat = function(x,y) {
    x = as.factor(x)
    y = as.factor(y)
    return(sparseMatrix(i = as.numeric(x), j = as.numeric(y), x=1, dims = c(nlevels(x), nlevels(y)), dimnames = list(levels(x),levels(y))))
}

if (sys.nframe() == 0) {
    main()
}
