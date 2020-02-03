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


green_sigs <- readRDS("data/variability_scans/Green_all_sigs.rda")

green_lands <- unique(green_sigs$scan_id) 
#green_comparisons <- data.frame(
#expand.grid(land1 = green_lands, land2 = green_lands), stringsAsFactors = FALSE)
green_comparisons <- make_comparison_grid(scan_ids = green_lands, self_comparisons = T)
green_comparisons <- green_comparisons %>% mutate(
  aligned = purrr::map2(.x = land1, .y = land2, .f = function(xx, yy) {
    land1 <- green_sigs$sigs[green_sigs$scan_id == xx][[1]]
    land2 <- green_sigs$sigs[green_sigs$scan_id == yy][[1]]
    land1$bullet <- "first-land"
    land2$bullet <- "second-land"
    
    sig_align(land1$sig, land2$sig)
  })
)

green_comparisons <- green_comparisons %>% mutate(
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

green_comparisons <- green_comparisons %>% mutate(
  striae = aligned %>% purrr::map(.f = sig_cms_max, span = 75) 
)
green_comparisons <- green_comparisons %>% mutate(
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

green_comparisons <- green_comparisons %>% mutate(
  features = purrr::map2(.x = aligned, .y = striae, .f = extract_features_all, resolution = 0.645)
)

green_comparisons <- green_comparisons %>% mutate(
  legacy_features = purrr::map(striae, extract_features_all_legacy, resolution = 0.645)
)

green_comparisons <- green_comparisons %>% tidyr::unnest(legacy_features) 
# scale features before using them in the random forest, legacy features can be used out of the box
green_comparisons$rfscore <- predict(bulletxtrctr::rtrees, newdata = green_comparisons, type = "prob")[,2]

rfscores_green <- green_comparisons %>% select(land1, land2, rfscore)
saveRDS(rfscores_green, "data/variability_scans/green_rfscores.rda")

