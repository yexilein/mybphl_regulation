
library(tidyverse)
source("../..//processed_data//processed_data.R", chdir = TRUE)


all_tissues = function() {
    c("eom","limb","larynx","heart")
}

all_grn_tissues = function(all_stages = FALSE) {
    if (all_stages) {
        c("eom","limb_all_stages","larynx")   
    } else {
        c("eom","limb","larynx")
    }
}

load_tissue = function(tissue_name, object_name=NULL) {
    dataset_name = c("eom"="eom_multiome","limb"="limb_multiome","limb_all_stages"="limb_multiome","larynx"="larynx_multiome","heart"="heart_multiome")[tissue_name]
    if (is.null(object_name)) {
        object_name = c("eom"="myo_adult","limb"="myo_adult","limb_all_stages"="myo","larynx"="all","heart"="all")[tissue_name]
    }
    result = load_srt(dataset_name = dataset_name, object_name=object_name)
    return(result)
}

load_pseudobulk = function(prefix=".", all_stages=FALSE) {
    suffix = ifelse(all_stages, "_all_stages", "")
    cache_file = glue::glue("{prefix}/srt/all_pseudobulk{suffix}.rds")
    if (!file.exists(cache_file)) {
        result = load_data(all_grn_tissues(all_stages = all_stages),"pseudobulk", prefix=prefix)
        saveRDS(result, cache_file)
    } else {
        result = readRDS(cache_file)
    }
    return(result)
}

load_data = function(tissue, suffix, prefix=".") {
    srt = lapply(set_names(tissue), function(t) { readRDS(glue::glue("{prefix}/srt/{t}_{suffix}.rds"))})
    for (i in seq_along(srt)) {
        srt[[i]]$tissue = names(srt)[i]
    }
    if (length(srt)>1) {
        srt = merge(srt[[1]], srt[-1])
    } else {
        srt = srt[[1]]
    }
    srt = JoinLayers(srt)
    srt = RunTFIDF(srt, assay="ATAC")
    return(srt)
}

load_metacells = function(prefix=".", all_stages=FALSE) {
    suffix = ifelse(all_stages, "_all_stages", "")
    cache_file = glue::glue("{prefix}/srt/all_metacell{suffix}.rds")
    if (!file.exists(cache_file)) {
        result = load_data(all_grn_tissues(all_stages=all_stages),"metacell", prefix=prefix)
        saveRDS(result, cache_file)
    } else {
        result = readRDS(cache_file)
    }
    return(result)
}

adult_ct = function() {
    c("GL Myh1 (Adult)", "GL Myh4 (Adult)", "OL (Adult)", "MTJ (Adult)", "NMJ (Adult)",
      "MTJ", "Myonuclei Myh1/Myh2 (Adult)", "Myonuclei Myh4 (Adult)", "NMJ")
}
