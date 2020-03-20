library(tidyverse)
library(bulletxtrctr)
library(x3ptools)
library(randomForest)
library(foreach)
library(doParallel)

#args <- commandArgs(TRUE)
#fold_id <- args[1] # should be a fold number between 1 and 5

#data_filename <- paste0("/media/Raven/Variability/processed_data_storage/blue/blue_pairwise_fold", fold_id, ".rda")
data_filename <- "/media/Raven/Variability/processed_data_storage/blue/blue_pairwise_all.rda"
make_comparison_grid <- function(scan_ids, self_comparisons = T){
  combos <- combn(scan_ids, 2) %>% t()
  combos_df <- data.frame(land1 = combos[,1], land2 = combos[,2]) 
  if(self_comparisons == T){
    self_combos <- data.frame(land1 = scan_ids, land2 = scan_ids)
    combos_df <- rbind(combos_df, self_combos) %>%
      mutate(pairing_id = paste0(land1, "_", land2))
  }
  else(combos_df <- combos_df %>% 
         mutate(pairing_id = paste0(land1, "_", land2)))
  return(combos_df)
}


blue_sigs <- readRDS("data/variability_scans/Blue_all_sigs_updated.rda")

blue_sigs <- blue_sigs %>% filter(barrel == "Barrel Blue")


blue_lands <- unique(blue_sigs$scan_id)

#blue_comparisons <- data.frame(
#expand.grid(land1 = blue_lands, land2 = blue_lands), stringsAsFactors = FALSE)

blue_comparisons <- make_comparison_grid(scan_ids = blue_lands, self_comparisons = T)

#blue_comparisons$fold <- sort(rep(1:5, length.out = nrow(blue_comparisons)))

#blue_comparisons <- blue_comparisons %>% filter(fold == fold_id)

## here need to split this into chunks
numCores <- 12 # on the server
registerDoParallel(cores = numCores)



list_out = foreach(m = 1:nrow(blue_comparisons)) %dopar% {
  
  
  blue_comps <- blue_comparisons[m,]
  blue_comps <- blue_comps %>% mutate(
    aligned = purrr::map2(.x = land1, .y = land2, .f = function(xx, yy) {
      land1 <- blue_sigs$sigs[blue_sigs$scan_id == xx][[1]]
      land2 <- blue_sigs$sigs[blue_sigs$scan_id == yy][[1]]
      land1$bullet <- "first-land"
      land2$bullet <- "second-land"
      
      sig_align(land1$sig, land2$sig)
    })
  )
  
  #blue_comps <- blue_comps %>% mutate(
  #  ccf0 = aligned %>% 
  #    purrr::map_dbl(.f = function(x) extract_feature_ccf(x$lands)),
  #  lag0 = aligned %>% 
  #    purrr::map_dbl(.f = function(x) extract_feature_lag(x$lands)),
  #  D0 = aligned %>% 
  #    purrr::map_dbl(.f = function(x) extract_feature_D(x$lands)),
  #  length0 = aligned %>% 
  #    purrr::map_dbl(.f = function(x) extract_feature_length(x$lands)),
  #  overlap0 = aligned %>% 
  #    purrr::map_dbl(.f = function(x) extract_feature_overlap(x$lands))
  #)
  
  blue_comps <- blue_comps %>% mutate(
    striae = aligned %>% purrr::map(.f = sig_cms_max, span = 75) 
  )
  #blue_comps <- blue_comps %>% mutate(
  #  cms_per_mm = purrr::map2(striae, aligned, .f = function(s, a) {
  #    extract_feature_cms_per_mm(s$lines, a$lands, resolution=0.645)
  #  }),
  #  matches0 = striae %>% purrr::map_dbl(.f = function(s) {
  #    bulletxtrctr:::extract_helper_feature_n_striae(s$lines, type = "peak", match = TRUE)
  #  }),
  #  mismatches0 = striae %>% purrr::map_dbl(.f = function(s) {
  #    bulletxtrctr:::extract_helper_feature_n_striae(s$lines, type = "peak", match = FALSE)
  #  })
  #)
  
  blue_comps <- blue_comps %>% mutate(
    features = purrr::map2(.x = aligned, .y = striae, .f = extract_features_all, resolution = 0.645)
  )
  
  #blue_comps <- blue_comps %>% mutate(
  #  legacy_features = purrr::map(striae, extract_features_all_legacy, resolution = 0.645)
  #)
  
  #blue_comps <- blue_comps %>% tidyr::unnest(legacy_features) 
  
  ## THIS STEP IS IMPORTANT
  blue_comps <- blue_comps %>% select(aligned, features) %>% 
    tidyr::unnest(features) %>%
    select(ccf, rough_cor, D, sd_D, matches_per_mm, mismatches_per_mm, cms_per_mm, non_cms_per_mm, sum_peaks)
  
  names(blue_comps) <- c("ccf", "rough_cor", "D", "sd_D", "matches", "mismatches", "cms", "non_cms", "sum_peaks")
  # scale features before using them in the random forest, legacy features can be used out of the box
  blue_comps$rfscore <- predict(bulletxtrctr::rtrees, newdata = blue_comps, type = "prob")[,2]
  
  #rfscores_blue <- blue_comps %>% select(land1, land2, rfscore)
  
  #rfscore <- blue_comps %>% pull(rfscore)
  return(blue_comps)
  ## this is the matrix that is returned as list element m, in the list list_out.
}

blue_comparisons$rf_feats_and_score <- list_out
blue_comparisons <- blue_comparisons %>% unnest(rf_feats_and_score)
#blue_comparisons <- blue_comparisons %>% mutate(rfscore = purrr::map_dbl(rfscore, .f = function(x) x[[1]]))
#rfscores_fold <- unlist(list_out)

saveRDS(blue_comparisons, file = data_filename)
rm(list_out)
gc()