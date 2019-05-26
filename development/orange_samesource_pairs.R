library(tidyverse)
library(bulletxtrctr)
library(x3ptools)
library(randomForest)

orange_sigs <- readRDS("data/variability_scans/orange_123_sigs.rda")
orange_sigs <- orange_sigs %>% 
  select(scan_id, sigs)

orange_gt <- read_csv("data/variability_scans/orange_gt.csv")
orange_sigs <- full_join(orange_sigs, orange_gt)

orange_ops <- read_csv("data/variability_scans/operators_blind.csv")
orange_sigs <- left_join(orange_sigs, orange_ops)

# Barrel-Land O-1
orange_l1 <- orange_sigs %>% filter(unique_id == "BL O-1") 
orange_lands1 <- unique(orange_l1$scan_id) 
orange_l1_comparisons <- data.frame(
  expand.grid(land1 = orange_lands1, land2 = orange_lands1), stringsAsFactors = FALSE)
orange_l1_comparisons <- orange_l1_comparisons %>% mutate(
  aligned = purrr::map2(.x = land1, .y = land2, .f = function(xx, yy) {
    land1 <- orange_l1$sigs[orange_l1$scan_id == xx][[1]]
    land2 <- orange_l1$sigs[orange_l1$scan_id == yy][[1]]
    land1$bullet <- "first-land"
    land2$bullet <- "second-land"
    
    sig_align(land1$sig, land2$sig)
  })
)
orange_l1_comparisons <- orange_l1_comparisons %>% mutate(
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

orange_l1_comparisons <- orange_l1_comparisons %>% mutate(
  striae = aligned %>% purrr::map(.f = sig_cms_max, span = 75) 
)
orange_l1_comparisons <- orange_l1_comparisons %>% mutate(
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


orange_l1_comparisons <- orange_l1_comparisons %>% mutate(
  features = purrr::map2(.x = aligned, .y = striae, .f = extract_features_all, resolution = 0.645)
)

orange_l1_comparisons <- orange_l1_comparisons %>% mutate(
  legacy_features = purrr::map(striae, extract_features_all_legacy, resolution = 0.645)
)

orange_l1_comparisons <- orange_l1_comparisons %>% tidyr::unnest(legacy_features) 
# scale features before using them in the random forest, legacy features can be used out of the box
orange_l1_comparisons$rfscore <- predict(bulletxtrctr::rtrees, newdata = orange_l1_comparisons, type = "prob")[,2]

rfscores_orange_l1 <- orange_l1_comparisons %>% select(land1, land2, rfscore)
saveRDS(rfscores_orange_l1, "data/variability_scans/orange_l1_rfscores.rda")


# Barrel-Land O-2
orange_l2 <- orange_sigs %>% filter(unique_id == "BL O-2") 
orange_lands2 <- unique(orange_l2$scan_id) 
orange_l2_comparisons <- data.frame(
  expand.grid(land1 = orange_lands2, land2 = orange_lands2), stringsAsFactors = FALSE)
orange_l2_comparisons <- orange_l2_comparisons %>% mutate(
  aligned = purrr::map2(.x = land1, .y = land2, .f = function(xx, yy) {
    land1 <- orange_l2$sigs[orange_l2$scan_id == xx][[1]]
    land2 <- orange_l2$sigs[orange_l2$scan_id == yy][[1]]
    land1$bullet <- "first-land"
    land2$bullet <- "second-land"
    
    sig_align(land1$sig, land2$sig)
  })
)
orange_l2_comparisons <- orange_l2_comparisons %>% mutate(
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

orange_l2_comparisons <- orange_l2_comparisons %>% mutate(
  striae = aligned %>% purrr::map(.f = sig_cms_max, span = 75) 
)
orange_l2_comparisons <- orange_l2_comparisons %>% mutate(
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


orange_l2_comparisons <- orange_l2_comparisons %>% mutate(
  features = purrr::map2(.x = aligned, .y = striae, .f = extract_features_all, resolution = 0.645)
)

orange_l2_comparisons <- orange_l2_comparisons %>% mutate(
  legacy_features = purrr::map(striae, extract_features_all_legacy, resolution = 0.645)
)

orange_l2_comparisons <- orange_l2_comparisons %>% tidyr::unnest(legacy_features) 
# scale features before using them in the random forest, legacy features can be used out of the box
orange_l2_comparisons$rfscore <- predict(bulletxtrctr::rtrees, newdata = orange_l2_comparisons, type = "prob")[,2]

rfscores_orange_l2 <- orange_l2_comparisons %>% select(land1, land2, rfscore)
saveRDS(rfscores_orange_l2, "data/variability_scans/orange_l2_rfscores.rda")

# Barrel-Land O-3
orange_l3 <- orange_sigs %>% filter(unique_id == "BL O-3") 
orange_lands3 <- unique(orange_l3$scan_id) 
orange_l3_comparisons <- data.frame(
  expand.grid(land1 = orange_lands3, land2 = orange_lands3), stringsAsFactors = FALSE)
orange_l3_comparisons <- orange_l3_comparisons %>% mutate(
  aligned = purrr::map2(.x = land1, .y = land2, .f = function(xx, yy) {
    land1 <- orange_l3$sigs[orange_l3$scan_id == xx][[1]]
    land2 <- orange_l3$sigs[orange_l3$scan_id == yy][[1]]
    land1$bullet <- "first-land"
    land2$bullet <- "second-land"
    
    sig_align(land1$sig, land2$sig)
  })
)
orange_l3_comparisons <- orange_l3_comparisons %>% mutate(
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

orange_l3_comparisons <- orange_l3_comparisons %>% mutate(
  striae = aligned %>% purrr::map(.f = sig_cms_max, span = 75) 
)
orange_l3_comparisons <- orange_l3_comparisons %>% mutate(
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


orange_l3_comparisons <- orange_l3_comparisons %>% mutate(
  features = purrr::map2(.x = aligned, .y = striae, .f = extract_features_all, resolution = 0.645)
)

orange_l3_comparisons <- orange_l3_comparisons %>% mutate(
  legacy_features = purrr::map(striae, extract_features_all_legacy, resolution = 0.645)
)

orange_l3_comparisons <- orange_l3_comparisons %>% tidyr::unnest(legacy_features) 
# scale features before using them in the random forest, legacy features can be used out of the box
orange_l3_comparisons$rfscore <- predict(bulletxtrctr::rtrees, newdata = orange_l3_comparisons, type = "prob")[,2]

rfscores_orange_l3 <- orange_l3_comparisons %>% select(land1, land2, rfscore)
saveRDS(rfscores_orange_l3, "data/variability_scans/orange_l3_rfscores.rda")


# Barrel-Land O-4
orange_l4 <- orange_sigs %>% filter(unique_id == "BL O-4") 
orange_lands4 <- unique(orange_l4$scan_id) 
orange_l4_comparisons <- data.frame(
  expand.grid(land1 = orange_lands4, land2 = orange_lands4), stringsAsFactors = FALSE)
orange_l4_comparisons <- orange_l4_comparisons %>% mutate(
  aligned = purrr::map2(.x = land1, .y = land2, .f = function(xx, yy) {
    land1 <- orange_l4$sigs[orange_l4$scan_id == xx][[1]]
    land2 <- orange_l4$sigs[orange_l4$scan_id == yy][[1]]
    land1$bullet <- "first-land"
    land2$bullet <- "second-land"
    
    sig_align(land1$sig, land2$sig)
  })
)
orange_l4_comparisons <- orange_l4_comparisons %>% mutate(
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

orange_l4_comparisons <- orange_l4_comparisons %>% mutate(
  striae = aligned %>% purrr::map(.f = sig_cms_max, span = 75) 
)
orange_l4_comparisons <- orange_l4_comparisons %>% mutate(
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


orange_l4_comparisons <- orange_l4_comparisons %>% mutate(
  features = purrr::map2(.x = aligned, .y = striae, .f = extract_features_all, resolution = 0.645)
)

orange_l4_comparisons <- orange_l4_comparisons %>% mutate(
  legacy_features = purrr::map(striae, extract_features_all_legacy, resolution = 0.645)
)

orange_l4_comparisons <- orange_l4_comparisons %>% tidyr::unnest(legacy_features) 
# scale features before using them in the random forest, legacy features can be used out of the box
orange_l4_comparisons$rfscore <- predict(bulletxtrctr::rtrees, newdata = orange_l4_comparisons, type = "prob")[,2]

rfscores_orange_l4 <- orange_l4_comparisons %>% select(land1, land2, rfscore)
saveRDS(rfscores_orange_l4, "data/variability_scans/orange_l4_rfscores.rda")


# Barrel-Land O-5
orange_l5 <- orange_sigs %>% filter(unique_id == "BL O-5") 
orange_lands5 <- unique(orange_l5$scan_id) 
orange_l5_comparisons <- data.frame(
  expand.grid(land1 = orange_lands5, land2 = orange_lands5), stringsAsFactors = FALSE)
orange_l5_comparisons <- orange_l5_comparisons %>% mutate(
  aligned = purrr::map2(.x = land1, .y = land2, .f = function(xx, yy) {
    land1 <- orange_l5$sigs[orange_l5$scan_id == xx][[1]]
    land2 <- orange_l5$sigs[orange_l5$scan_id == yy][[1]]
    land1$bullet <- "first-land"
    land2$bullet <- "second-land"
    
    sig_align(land1$sig, land2$sig)
  })
)
orange_l5_comparisons <- orange_l5_comparisons %>% mutate(
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

orange_l5_comparisons <- orange_l5_comparisons %>% mutate(
  striae = aligned %>% purrr::map(.f = sig_cms_max, span = 75) 
)
orange_l5_comparisons <- orange_l5_comparisons %>% mutate(
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


orange_l5_comparisons <- orange_l5_comparisons %>% mutate(
  features = purrr::map2(.x = aligned, .y = striae, .f = extract_features_all, resolution = 0.645)
)

orange_l5_comparisons <- orange_l5_comparisons %>% mutate(
  legacy_features = purrr::map(striae, extract_features_all_legacy, resolution = 0.645)
)

orange_l5_comparisons <- orange_l5_comparisons %>% tidyr::unnest(legacy_features) 
# scale features before using them in the random forest, legacy features can be used out of the box
orange_l5_comparisons$rfscore <- predict(bulletxtrctr::rtrees, newdata = orange_l5_comparisons, type = "prob")[,2]

rfscores_orange_l5 <- orange_l5_comparisons %>% select(land1, land2, rfscore)
saveRDS(rfscores_orange_l5, "data/variability_scans/orange_l5_rfscores.rda")


# Barrel-Land O-6
orange_l6 <- orange_sigs %>% filter(unique_id == "BL O-6") 
orange_lands6 <- unique(orange_l6$scan_id) 
orange_l6_comparisons <- data.frame(
  expand.grid(land1 = orange_lands6, land2 = orange_lands6), stringsAsFactors = FALSE)
orange_l6_comparisons <- orange_l6_comparisons %>% mutate(
  aligned = purrr::map2(.x = land1, .y = land2, .f = function(xx, yy) {
    land1 <- orange_l6$sigs[orange_l6$scan_id == xx][[1]]
    land2 <- orange_l6$sigs[orange_l6$scan_id == yy][[1]]
    land1$bullet <- "first-land"
    land2$bullet <- "second-land"
    
    sig_align(land1$sig, land2$sig)
  })
)
orange_l6_comparisons <- orange_l6_comparisons %>% mutate(
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

orange_l6_comparisons <- orange_l6_comparisons %>% mutate(
  striae = aligned %>% purrr::map(.f = sig_cms_max, span = 75) 
)
orange_l6_comparisons <- orange_l6_comparisons %>% mutate(
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


orange_l6_comparisons <- orange_l6_comparisons %>% mutate(
  features = purrr::map2(.x = aligned, .y = striae, .f = extract_features_all, resolution = 0.645)
)

orange_l6_comparisons <- orange_l6_comparisons %>% mutate(
  legacy_features = purrr::map(striae, extract_features_all_legacy, resolution = 0.645)
)

orange_l6_comparisons <- orange_l6_comparisons %>% tidyr::unnest(legacy_features) 
# scale features before using them in the random forest, legacy features can be used out of the box
orange_l6_comparisons$rfscore <- predict(bulletxtrctr::rtrees, newdata = orange_l6_comparisons, type = "prob")[,2]

rfscores_orange_l6 <- orange_l6_comparisons %>% select(land1, land2, rfscore)
saveRDS(rfscores_orange_l6, "data/variability_scans/orange_l6_rfscores.rda")





