library(tidyverse)
library(bulletxtrctr)
library(x3ptools)
library(randomForest)


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


pink_sigs <- readRDS("data/variability_scans/Pink_all_sigs.rda")

pink_lands <- unique(pink_sigs$scan_id) 
#pink_comparisons <- data.frame(
#expand.grid(land1 = pink_lands, land2 = pink_lands), stringsAsFactors = FALSE)
pink_comparisons <- make_comparison_grid(scan_ids = pink_lands, self_comparisons = T)
pink_comparisons <- pink_comparisons %>% mutate(
  aligned = purrr::map2(.x = land1, .y = land2, .f = function(xx, yy) {
    land1 <- pink_sigs$sigs[pink_sigs$scan_id == xx][[1]]
    land2 <- pink_sigs$sigs[pink_sigs$scan_id == yy][[1]]
    land1$bullet <- "first-land"
    land2$bullet <- "second-land"
    
    sig_align(land1$sig, land2$sig)
  })
)

pink_comparisons <- pink_comparisons %>% mutate(
  ccf0 = aligned %>% 
    purrr::map_dbl(.f = function(x) extract_feature_ccf(x$lands)),
  lag0 = aligned %>% 
    purrr::map_dbl(.f = function(x) extract_feature_lag(x$lands)),
  D0 = aligned %>% 
    purrr::map_dbl(.f = function(x) extract_feature_D(x$lands)),
  length0 = aligned %>% 
    purrr::map_dbl(.f = function(x) extract_feature_length(x$lands)),
  overlap0 = aligned %>% 
    purrr::map_dbl(.f = function(x) extract_feature_overlap(x$lands))
)

pink_comparisons <- pink_comparisons %>% mutate(
  striae = aligned %>% purrr::map(.f = sig_cms_max, span = 75) 
)
pink_comparisons <- pink_comparisons %>% mutate(
  cms_per_mm = purrr::map2(striae, aligned, .f = function(s, a) {
    extract_feature_cms_per_mm(s$lines, a$lands, resolution=0.645)
  }),
  matches0 = striae %>% purrr::map_dbl(.f = function(s) {
    bulletxtrctr:::extract_helper_feature_n_striae(s$lines, type = "peak", match = TRUE)
  }),
  mismatches0 = striae %>% purrr::map_dbl(.f = function(s) {
    bulletxtrctr:::extract_helper_feature_n_striae(s$lines, type = "peak", match = FALSE)
  })
)

pink_comparisons <- pink_comparisons %>% mutate(
  features = purrr::map2(.x = aligned, .y = striae, .f = extract_features_all, resolution = 0.645)
)

pink_comparisons <- pink_comparisons %>% mutate(
  legacy_features = purrr::map(striae, extract_features_all_legacy, resolution = 0.645)
)

pink_comparisons <- pink_comparisons %>% tidyr::unnest(legacy_features) 
# scale features before using them in the random forest, legacy features can be used out of the box
pink_comparisons$rfscore <- predict(bulletxtrctr::rtrees, newdata = pink_comparisons, type = "prob")[,2]

rfscores_pink <- pink_comparisons %>% select(land1, land2, rfscore)
saveRDS(rfscores_pink, "data/variability_scans/pink_rfscores.rda")

