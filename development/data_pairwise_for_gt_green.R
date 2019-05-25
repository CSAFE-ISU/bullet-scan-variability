library(tidyverse)
library(bulletxtrctr)
library(x3ptools)
library(gridExtra)
library(randomForest)

green_sigs <- readRDS("data/variability_scans/green_123_sigs.rda")

# COMPARE BULLETS 1, 2, 3 FIRST to get "reference"
green_bullets_123 <- green_sigs %>% filter(operator == "Connor", round == "Round 1", machine == "Sneox1") 
green_b123_lands <- unique(green_bullets_123$scan_id) 
green_b123_comparisons <- data.frame(
  expand.grid(land1 = green_b123_lands, land2 = green_b123_lands), stringsAsFactors = FALSE)
green_b123_comparisons <- green_b123_comparisons %>% mutate(
  aligned = purrr::map2(.x = land1, .y = land2, .f = function(xx, yy) {
    land1 <- green_bullets_123$sigs[green_bullets_123$scan_id == xx][[1]]
    land2 <- green_bullets_123$sigs[green_bullets_123$scan_id == yy][[1]]
    land1$bullet <- "first-land"
    land2$bullet <- "second-land"
    
    sig_align(land1$sig, land2$sig)
  })
)
green_b123_comparisons <- green_b123_comparisons %>% mutate(
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

green_b123_comparisons <- green_b123_comparisons %>% mutate(
  striae = aligned %>% purrr::map(.f = sig_cms_max, span = 75) 
)
green_b123_comparisons <- green_b123_comparisons %>% mutate(
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


green_b123_comparisons <- green_b123_comparisons %>% mutate(
  features = purrr::map2(.x = aligned, .y = striae, .f = extract_features_all, resolution = 0.645)
)

green_b123_comparisons <- green_b123_comparisons %>% mutate(
  legacy_features = purrr::map(striae, extract_features_all_legacy, resolution = 0.645)
)

green_b123_comparisons <- green_b123_comparisons %>% tidyr::unnest(legacy_features) 
# scale features before using them in the random forest, legacy features can be used out of the box
green_b123_comparisons$rfscore <- predict(bulletxtrctr::rtrees, newdata = green_b123_comparisons, type = "prob")[,2]

rfscores_green_b123 <- green_b123_comparisons %>% select(land1, land2, rfscore)
saveRDS(rfscores_green_b123, "data/variability_scans/green_b123_rfscores.rda")

# NOW DO ALL OF BULLET 1  

green_bullet1 <- green_sigs %>% filter(bullet == "Bullet 1") 
green_b1_lands <- unique(green_bullet1$scan_id) 
green_b1_comparisons <- data.frame(
  expand.grid(land1 = green_b1_lands, land2 = green_b1_lands), stringsAsFactors = FALSE)
green_b1_comparisons <- green_b1_comparisons %>% mutate(
  aligned = purrr::map2(.x = land1, .y = land2, .f = function(xx, yy) {
    land1 <- green_bullet1$sigs[green_bullet1$scan_id == xx][[1]]
    land2 <- green_bullet1$sigs[green_bullet1$scan_id == yy][[1]]
    land1$bullet <- "first-land"
    land2$bullet <- "second-land"
    
    sig_align(land1$sig, land2$sig)
  })
)

green_b1_comparisons <- green_b1_comparisons %>% mutate(
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

green_b1_comparisons <- green_b1_comparisons %>% mutate(
  striae = aligned %>% purrr::map(.f = sig_cms_max, span = 75) 
)
green_b1_comparisons <- green_b1_comparisons %>% mutate(
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

green_b1_comparisons <- green_b1_comparisons %>% mutate(
  features = purrr::map2(.x = aligned, .y = striae, .f = extract_features_all, resolution = 0.645)
)

green_b1_comparisons <- green_b1_comparisons %>% mutate(
  legacy_features = purrr::map(striae, extract_features_all_legacy, resolution = 0.645)
)

green_b1_comparisons <- green_b1_comparisons %>% tidyr::unnest(legacy_features) 
# scale features before using them in the random forest, legacy features can be used out of the box
green_b1_comparisons$rfscore <- predict(bulletxtrctr::rtrees, newdata = green_b1_comparisons, type = "prob")[,2]

rfscores_green_b1 <- green_b1_comparisons %>% select(land1, land2, rfscore)
saveRDS(rfscores_green_b1, "data/variability_scans/green_b1_rfscores.rda")

# NOW DO ALL OF BULLET 2  

green_bullet2 <- green_sigs %>% filter(bullet == "Bullet 2") 
green_b2_lands <- unique(green_bullet2$scan_id) 
green_b2_comparisons <- data.frame(
  expand.grid(land1 = green_b2_lands, land2 = green_b2_lands), stringsAsFactors = FALSE)
green_b2_comparisons <- green_b2_comparisons %>% mutate(
  aligned = purrr::map2(.x = land1, .y = land2, .f = function(xx, yy) {
    land1 <- green_bullet2$sigs[green_bullet2$scan_id == xx][[1]]
    land2 <- green_bullet2$sigs[green_bullet2$scan_id == yy][[1]]
    land1$bullet <- "first-land"
    land2$bullet <- "second-land"
    
    sig_align(land1$sig, land2$sig)
  })
)

green_b2_comparisons <- green_b2_comparisons %>% mutate(
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

green_b2_comparisons <- green_b2_comparisons %>% mutate(
  striae = aligned %>% purrr::map(.f = sig_cms_max, span = 75) 
)
green_b2_comparisons <- green_b2_comparisons %>% mutate(
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

green_b2_comparisons <- green_b2_comparisons %>% mutate(
  features = purrr::map2(.x = aligned, .y = striae, .f = extract_features_all, resolution = 0.645)
)

green_b2_comparisons <- green_b2_comparisons %>% mutate(
  legacy_features = purrr::map(striae, extract_features_all_legacy, resolution = 0.645)
)

green_b2_comparisons <- green_b2_comparisons %>% tidyr::unnest(legacy_features) 
# scale features before using them in the random forest, legacy features can be used out of the box
green_b2_comparisons$rfscore <- predict(bulletxtrctr::rtrees, newdata = green_b2_comparisons, type = "prob")[,2]

#head(comparisons_user)

#saveRDS(comparisons_user, "../data/user_meas_comparisons.rda")

rfscores_green_b2 <- green_b2_comparisons %>% select(land1, land2, rfscore)
saveRDS(rfscores_green_b2, "data/variability_scans/green_b2_rfscores.rda")

# NOW DO ALL OF BULLET 3  

green_bullet3 <- green_sigs %>% filter(bullet == "Bullet 3") 
green_b3_lands <- unique(green_bullet3$scan_id) 
green_b3_comparisons <- data.frame(
  expand.grid(land1 = green_b3_lands, land2 = green_b3_lands), stringsAsFactors = FALSE)
green_b3_comparisons <- green_b3_comparisons %>% mutate(
  aligned = purrr::map2(.x = land1, .y = land2, .f = function(xx, yy) {
    land1 <- green_bullet3$sigs[green_bullet3$scan_id == xx][[1]]
    land2 <- green_bullet3$sigs[green_bullet3$scan_id == yy][[1]]
    land1$bullet <- "first-land"
    land2$bullet <- "second-land"
    
    sig_align(land1$sig, land2$sig)
  })
)

green_b3_comparisons <- green_b3_comparisons %>% mutate(
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

green_b3_comparisons <- green_b3_comparisons %>% mutate(
  striae = aligned %>% purrr::map(.f = sig_cms_max, span = 75) 
)
green_b3_comparisons <- green_b3_comparisons %>% mutate(
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

green_b3_comparisons <- green_b3_comparisons %>% mutate(
  features = purrr::map2(.x = aligned, .y = striae, .f = extract_features_all, resolution = 0.645)
)

green_b3_comparisons <- green_b3_comparisons %>% mutate(
  legacy_features = purrr::map(striae, extract_features_all_legacy, resolution = 0.645)
)

green_b3_comparisons <- green_b3_comparisons %>% tidyr::unnest(legacy_features) 
# scale features before using them in the random forest, legacy features can be used out of the box
green_b3_comparisons$rfscore <- predict(bulletxtrctr::rtrees, newdata = green_b3_comparisons, type = "prob")[,2]


rfscores_green_b3 <- green_b3_comparisons %>% select(land1, land2, rfscore)
saveRDS(rfscores_green_b3, "data/variability_scans/green_b3_rfscores.rda")

