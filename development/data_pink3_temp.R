library(tidyverse)
library(bulletxtrctr)
library(x3ptools)
library(gridExtra)
library(randomForest)

pink_sigs <- readRDS("data/variability_scans/pink_123_sigs.rda")


pink_bullet3 <- pink_sigs %>% filter(bullet == "Bullet 3") 
pink_b3_lands <- unique(pink_bullet3$scan_id) 
pink_b3_comparisons <- data.frame(
  expand.grid(land1 = pink_b3_lands, land2 = pink_b3_lands), stringsAsFactors = FALSE)
pink_b3_comparisons <- pink_b3_comparisons %>% mutate(
  aligned = purrr::map2(.x = land1, .y = land2, .f = function(xx, yy) {
    land1 <- pink_bullet3$sigs[pink_bullet3$scan_id == xx][[1]]
    land2 <- pink_bullet3$sigs[pink_bullet3$scan_id == yy][[1]]
    land1$bullet <- "first-land"
    land2$bullet <- "second-land"
    
    sig_align(land1$sig, land2$sig)
  })
)

pink_b3_comparisons <- pink_b3_comparisons %>% mutate(
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

pink_b3_comparisons <- pink_b3_comparisons %>% mutate(
  striae = aligned %>% purrr::map(.f = sig_cms_max, span = 75) 
)
pink_b3_comparisons <- pink_b3_comparisons %>% mutate(
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

pink_b3_comparisons <- pink_b3_comparisons %>% mutate(
  features = purrr::map2(.x = aligned, .y = striae, .f = extract_features_all, resolution = 0.645)
)

pink_b3_comparisons <- pink_b3_comparisons %>% mutate(
  legacy_features = purrr::map(striae, extract_features_all_legacy, resolution = 0.645)
)

pink_b3_comparisons <- pink_b3_comparisons %>% tidyr::unnest(legacy_features) 
# scale features before using them in the random forest, legacy features can be used out of the box
pink_b3_comparisons$rfscore <- predict(bulletxtrctr::rtrees, newdata = pink_b3_comparisons, type = "prob")[,2]


rfscores_pink_b3 <- pink_b3_comparisons %>% select(land1, land2, rfscore)
saveRDS(rfscores_pink_b3, "data/variability_scans/pink_b3_rfscores.rda")

