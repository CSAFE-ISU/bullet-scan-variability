library(tidyverse)
library(bulletxtrctr)
library(x3ptools)
library(randomForest)

green_sigs <- readRDS("data/variability_scans/green_123_sigs.rda")
green_sigs <- green_sigs %>% 
  select(scan_id, sigs)

green_gt <- read_csv("data/variability_scans/green_gt.csv")
green_sigs <- full_join(green_sigs, green_gt)

green_ops <- read_csv("data/variability_scans/operators_blind.csv")
green_sigs <- left_join(green_sigs, green_ops)

# Barrel-Land G-1
green_l1 <- green_sigs %>% filter(unique_id == "BL G-1") 
green_lands1 <- unique(green_l1$scan_id) 
green_l1_comparisons <- data.frame(
  expand.grid(land1 = green_lands1, land2 = green_lands1), stringsAsFactors = FALSE)
green_l1_comparisons <- green_l1_comparisons %>% mutate(
  aligned = purrr::map2(.x = land1, .y = land2, .f = function(xx, yy) {
    land1 <- green_l1$sigs[green_l1$scan_id == xx][[1]]
    land2 <- green_l1$sigs[green_l1$scan_id == yy][[1]]
    land1$bullet <- "first-land"
    land2$bullet <- "second-land"
    
    sig_align(land1$sig, land2$sig)
  })
)
green_l1_comparisons <- green_l1_comparisons %>% mutate(
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

green_l1_comparisons <- green_l1_comparisons %>% mutate(
  striae = aligned %>% purrr::map(.f = sig_cms_max, span = 75) 
)
green_l1_comparisons <- green_l1_comparisons %>% mutate(
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


green_l1_comparisons <- green_l1_comparisons %>% mutate(
  features = purrr::map2(.x = aligned, .y = striae, .f = extract_features_all, resolution = 0.645)
)

green_l1_comparisons <- green_l1_comparisons %>% mutate(
  legacy_features = purrr::map(striae, extract_features_all_legacy, resolution = 0.645)
)

green_l1_comparisons <- green_l1_comparisons %>% tidyr::unnest(legacy_features) 
# scale features before using them in the random forest, legacy features can be used out of the box
green_l1_comparisons$rfscore <- predict(bulletxtrctr::rtrees, newdata = green_l1_comparisons, type = "prob")[,2]

rfscores_green_l1 <- green_l1_comparisons %>% select(land1, land2, rfscore)
saveRDS(rfscores_green_l1, "data/variability_scans/green_l1_rfscores.rda")


# Barrel-Land G-2
green_l2 <- green_sigs %>% filter(unique_id == "BL G-2") 
green_lands2 <- unique(green_l2$scan_id) 
green_l2_comparisons <- data.frame(
  expand.grid(land1 = green_lands2, land2 = green_lands2), stringsAsFactors = FALSE)
green_l2_comparisons <- green_l2_comparisons %>% mutate(
  aligned = purrr::map2(.x = land1, .y = land2, .f = function(xx, yy) {
    land1 <- green_l2$sigs[green_l2$scan_id == xx][[1]]
    land2 <- green_l2$sigs[green_l2$scan_id == yy][[1]]
    land1$bullet <- "first-land"
    land2$bullet <- "second-land"
    
    sig_align(land1$sig, land2$sig)
  })
)
green_l2_comparisons <- green_l2_comparisons %>% mutate(
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

green_l2_comparisons <- green_l2_comparisons %>% mutate(
  striae = aligned %>% purrr::map(.f = sig_cms_max, span = 75) 
)
green_l2_comparisons <- green_l2_comparisons %>% mutate(
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


green_l2_comparisons <- green_l2_comparisons %>% mutate(
  features = purrr::map2(.x = aligned, .y = striae, .f = extract_features_all, resolution = 0.645)
)

green_l2_comparisons <- green_l2_comparisons %>% mutate(
  legacy_features = purrr::map(striae, extract_features_all_legacy, resolution = 0.645)
)

green_l2_comparisons <- green_l2_comparisons %>% tidyr::unnest(legacy_features) 
# scale features before using them in the random forest, legacy features can be used out of the box
green_l2_comparisons$rfscore <- predict(bulletxtrctr::rtrees, newdata = green_l2_comparisons, type = "prob")[,2]

rfscores_green_l2 <- green_l2_comparisons %>% select(land1, land2, rfscore)
saveRDS(rfscores_green_l2, "data/variability_scans/green_l2_rfscores.rda")

# Barrel-Land G-3
green_l3 <- green_sigs %>% filter(unique_id == "BL G-3") 
green_lands3 <- unique(green_l3$scan_id) 
green_l3_comparisons <- data.frame(
  expand.grid(land1 = green_lands3, land2 = green_lands3), stringsAsFactors = FALSE)
green_l3_comparisons <- green_l3_comparisons %>% mutate(
  aligned = purrr::map2(.x = land1, .y = land2, .f = function(xx, yy) {
    land1 <- green_l3$sigs[green_l3$scan_id == xx][[1]]
    land2 <- green_l3$sigs[green_l3$scan_id == yy][[1]]
    land1$bullet <- "first-land"
    land2$bullet <- "second-land"
    
    sig_align(land1$sig, land2$sig)
  })
)
green_l3_comparisons <- green_l3_comparisons %>% mutate(
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

green_l3_comparisons <- green_l3_comparisons %>% mutate(
  striae = aligned %>% purrr::map(.f = sig_cms_max, span = 75) 
)
green_l3_comparisons <- green_l3_comparisons %>% mutate(
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


green_l3_comparisons <- green_l3_comparisons %>% mutate(
  features = purrr::map2(.x = aligned, .y = striae, .f = extract_features_all, resolution = 0.645)
)

green_l3_comparisons <- green_l3_comparisons %>% mutate(
  legacy_features = purrr::map(striae, extract_features_all_legacy, resolution = 0.645)
)

green_l3_comparisons <- green_l3_comparisons %>% tidyr::unnest(legacy_features) 
# scale features before using them in the random forest, legacy features can be used out of the box
green_l3_comparisons$rfscore <- predict(bulletxtrctr::rtrees, newdata = green_l3_comparisons, type = "prob")[,2]

rfscores_green_l3 <- green_l3_comparisons %>% select(land1, land2, rfscore)
saveRDS(rfscores_green_l3, "data/variability_scans/green_l3_rfscores.rda")


# Barrel-Land G-4
green_l4 <- green_sigs %>% filter(unique_id == "BL G-4") 
green_lands4 <- unique(green_l4$scan_id) 
green_l4_comparisons <- data.frame(
  expand.grid(land1 = green_lands4, land2 = green_lands4), stringsAsFactors = FALSE)
green_l4_comparisons <- green_l4_comparisons %>% mutate(
  aligned = purrr::map2(.x = land1, .y = land2, .f = function(xx, yy) {
    land1 <- green_l4$sigs[green_l4$scan_id == xx][[1]]
    land2 <- green_l4$sigs[green_l4$scan_id == yy][[1]]
    land1$bullet <- "first-land"
    land2$bullet <- "second-land"
    
    sig_align(land1$sig, land2$sig)
  })
)
green_l4_comparisons <- green_l4_comparisons %>% mutate(
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

green_l4_comparisons <- green_l4_comparisons %>% mutate(
  striae = aligned %>% purrr::map(.f = sig_cms_max, span = 75) 
)
green_l4_comparisons <- green_l4_comparisons %>% mutate(
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


green_l4_comparisons <- green_l4_comparisons %>% mutate(
  features = purrr::map2(.x = aligned, .y = striae, .f = extract_features_all, resolution = 0.645)
)

green_l4_comparisons <- green_l4_comparisons %>% mutate(
  legacy_features = purrr::map(striae, extract_features_all_legacy, resolution = 0.645)
)

green_l4_comparisons <- green_l4_comparisons %>% tidyr::unnest(legacy_features) 
# scale features before using them in the random forest, legacy features can be used out of the box
green_l4_comparisons$rfscore <- predict(bulletxtrctr::rtrees, newdata = green_l4_comparisons, type = "prob")[,2]

rfscores_green_l4 <- green_l4_comparisons %>% select(land1, land2, rfscore)
saveRDS(rfscores_green_l4, "data/variability_scans/green_l4_rfscores.rda")


# Barrel-Land G-5
green_l5 <- green_sigs %>% filter(unique_id == "BL G-5") 
green_lands5 <- unique(green_l5$scan_id) 
green_l5_comparisons <- data.frame(
  expand.grid(land1 = green_lands5, land2 = green_lands5), stringsAsFactors = FALSE)
green_l5_comparisons <- green_l5_comparisons %>% mutate(
  aligned = purrr::map2(.x = land1, .y = land2, .f = function(xx, yy) {
    land1 <- green_l5$sigs[green_l5$scan_id == xx][[1]]
    land2 <- green_l5$sigs[green_l5$scan_id == yy][[1]]
    land1$bullet <- "first-land"
    land2$bullet <- "second-land"
    
    sig_align(land1$sig, land2$sig)
  })
)
green_l5_comparisons <- green_l5_comparisons %>% mutate(
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

green_l5_comparisons <- green_l5_comparisons %>% mutate(
  striae = aligned %>% purrr::map(.f = sig_cms_max, span = 75) 
)
green_l5_comparisons <- green_l5_comparisons %>% mutate(
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


green_l5_comparisons <- green_l5_comparisons %>% mutate(
  features = purrr::map2(.x = aligned, .y = striae, .f = extract_features_all, resolution = 0.645)
)

green_l5_comparisons <- green_l5_comparisons %>% mutate(
  legacy_features = purrr::map(striae, extract_features_all_legacy, resolution = 0.645)
)

green_l5_comparisons <- green_l5_comparisons %>% tidyr::unnest(legacy_features) 
# scale features before using them in the random forest, legacy features can be used out of the box
green_l5_comparisons$rfscore <- predict(bulletxtrctr::rtrees, newdata = green_l5_comparisons, type = "prob")[,2]

rfscores_green_l5 <- green_l5_comparisons %>% select(land1, land2, rfscore)
saveRDS(rfscores_green_l5, "data/variability_scans/green_l5_rfscores.rda")


# Barrel-Land G-6
green_l6 <- green_sigs %>% filter(unique_id == "BL G-6") 
green_lands6 <- unique(green_l6$scan_id) 
green_l6_comparisons <- data.frame(
  expand.grid(land1 = green_lands6, land2 = green_lands6), stringsAsFactors = FALSE)
green_l6_comparisons <- green_l6_comparisons %>% mutate(
  aligned = purrr::map2(.x = land1, .y = land2, .f = function(xx, yy) {
    land1 <- green_l6$sigs[green_l6$scan_id == xx][[1]]
    land2 <- green_l6$sigs[green_l6$scan_id == yy][[1]]
    land1$bullet <- "first-land"
    land2$bullet <- "second-land"
    
    sig_align(land1$sig, land2$sig)
  })
)
green_l6_comparisons <- green_l6_comparisons %>% mutate(
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

green_l6_comparisons <- green_l6_comparisons %>% mutate(
  striae = aligned %>% purrr::map(.f = sig_cms_max, span = 75) 
)
green_l6_comparisons <- green_l6_comparisons %>% mutate(
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


green_l6_comparisons <- green_l6_comparisons %>% mutate(
  features = purrr::map2(.x = aligned, .y = striae, .f = extract_features_all, resolution = 0.645)
)

green_l6_comparisons <- green_l6_comparisons %>% mutate(
  legacy_features = purrr::map(striae, extract_features_all_legacy, resolution = 0.645)
)

green_l6_comparisons <- green_l6_comparisons %>% tidyr::unnest(legacy_features) 
# scale features before using them in the random forest, legacy features can be used out of the box
green_l6_comparisons$rfscore <- predict(bulletxtrctr::rtrees, newdata = green_l6_comparisons, type = "prob")[,2]

rfscores_green_l6 <- green_l6_comparisons %>% select(land1, land2, rfscore)
saveRDS(rfscores_green_l6, "data/variability_scans/green_l6_rfscores.rda")





