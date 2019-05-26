library(tidyverse)
library(bulletxtrctr)
library(x3ptools)
library(randomForest)

pink_sigs <- readRDS("data/variability_scans/pink_123_sigs.rda")
pink_sigs <- pink_sigs %>% 
  select(scan_id, sigs)

pink_gt <- read_csv("data/variability_scans/pink_gt.csv")
pink_sigs <- full_join(pink_sigs, pink_gt)

pink_ops <- read_csv("data/variability_scans/operators_blind.csv")
pink_sigs <- left_join(pink_sigs, pink_ops)

# Barrel-Land P-1
pink_l1 <- pink_sigs %>% filter(unique_id == "BL P-1") 
pink_lands1 <- unique(pink_l1$scan_id) 
pink_l1_comparisons <- data.frame(
  expand.grid(land1 = pink_lands1, land2 = pink_lands1), stringsAsFactors = FALSE)
pink_l1_comparisons <- pink_l1_comparisons %>% mutate(
  aligned = purrr::map2(.x = land1, .y = land2, .f = function(xx, yy) {
    land1 <- pink_l1$sigs[pink_l1$scan_id == xx][[1]]
    land2 <- pink_l1$sigs[pink_l1$scan_id == yy][[1]]
    land1$bullet <- "first-land"
    land2$bullet <- "second-land"
    
    sig_align(land1$sig, land2$sig)
  })
)
pink_l1_comparisons <- pink_l1_comparisons %>% mutate(
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

pink_l1_comparisons <- pink_l1_comparisons %>% mutate(
  striae = aligned %>% purrr::map(.f = sig_cms_max, span = 75) 
)
pink_l1_comparisons <- pink_l1_comparisons %>% mutate(
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


pink_l1_comparisons <- pink_l1_comparisons %>% mutate(
  features = purrr::map2(.x = aligned, .y = striae, .f = extract_features_all, resolution = 0.645)
)

pink_l1_comparisons <- pink_l1_comparisons %>% mutate(
  legacy_features = purrr::map(striae, extract_features_all_legacy, resolution = 0.645)
)

pink_l1_comparisons <- pink_l1_comparisons %>% tidyr::unnest(legacy_features) 
# scale features before using them in the random forest, legacy features can be used out of the box
pink_l1_comparisons$rfscore <- predict(bulletxtrctr::rtrees, newdata = pink_l1_comparisons, type = "prob")[,2]

rfscores_pink_l1 <- pink_l1_comparisons %>% select(land1, land2, rfscore)
saveRDS(rfscores_pink_l1, "data/variability_scans/pink_l1_rfscores.rda")


# Barrel-Land P-2
pink_l2 <- pink_sigs %>% filter(unique_id == "BL P-2") 
pink_lands2 <- unique(pink_l2$scan_id) 
pink_l2_comparisons <- data.frame(
  expand.grid(land1 = pink_lands2, land2 = pink_lands2), stringsAsFactors = FALSE)
pink_l2_comparisons <- pink_l2_comparisons %>% mutate(
  aligned = purrr::map2(.x = land1, .y = land2, .f = function(xx, yy) {
    land1 <- pink_l2$sigs[pink_l2$scan_id == xx][[1]]
    land2 <- pink_l2$sigs[pink_l2$scan_id == yy][[1]]
    land1$bullet <- "first-land"
    land2$bullet <- "second-land"
    
    sig_align(land1$sig, land2$sig)
  })
)
pink_l2_comparisons <- pink_l2_comparisons %>% mutate(
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

pink_l2_comparisons <- pink_l2_comparisons %>% mutate(
  striae = aligned %>% purrr::map(.f = sig_cms_max, span = 75) 
)
pink_l2_comparisons <- pink_l2_comparisons %>% mutate(
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


pink_l2_comparisons <- pink_l2_comparisons %>% mutate(
  features = purrr::map2(.x = aligned, .y = striae, .f = extract_features_all, resolution = 0.645)
)

pink_l2_comparisons <- pink_l2_comparisons %>% mutate(
  legacy_features = purrr::map(striae, extract_features_all_legacy, resolution = 0.645)
)

pink_l2_comparisons <- pink_l2_comparisons %>% tidyr::unnest(legacy_features) 
# scale features before using them in the random forest, legacy features can be used out of the box
pink_l2_comparisons$rfscore <- predict(bulletxtrctr::rtrees, newdata = pink_l2_comparisons, type = "prob")[,2]

rfscores_pink_l2 <- pink_l2_comparisons %>% select(land1, land2, rfscore)
saveRDS(rfscores_pink_l2, "data/variability_scans/pink_l2_rfscores.rda")

# Barrel-Land P-3
pink_l3 <- pink_sigs %>% filter(unique_id == "BL P-3") 
pink_lands3 <- unique(pink_l3$scan_id) 
pink_l3_comparisons <- data.frame(
  expand.grid(land1 = pink_lands3, land2 = pink_lands3), stringsAsFactors = FALSE)
pink_l3_comparisons <- pink_l3_comparisons %>% mutate(
  aligned = purrr::map2(.x = land1, .y = land2, .f = function(xx, yy) {
    land1 <- pink_l3$sigs[pink_l3$scan_id == xx][[1]]
    land2 <- pink_l3$sigs[pink_l3$scan_id == yy][[1]]
    land1$bullet <- "first-land"
    land2$bullet <- "second-land"
    
    sig_align(land1$sig, land2$sig)
  })
)
pink_l3_comparisons <- pink_l3_comparisons %>% mutate(
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

pink_l3_comparisons <- pink_l3_comparisons %>% mutate(
  striae = aligned %>% purrr::map(.f = sig_cms_max, span = 75) 
)
pink_l3_comparisons <- pink_l3_comparisons %>% mutate(
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


pink_l3_comparisons <- pink_l3_comparisons %>% mutate(
  features = purrr::map2(.x = aligned, .y = striae, .f = extract_features_all, resolution = 0.645)
)

pink_l3_comparisons <- pink_l3_comparisons %>% mutate(
  legacy_features = purrr::map(striae, extract_features_all_legacy, resolution = 0.645)
)

pink_l3_comparisons <- pink_l3_comparisons %>% tidyr::unnest(legacy_features) 
# scale features before using them in the random forest, legacy features can be used out of the box
pink_l3_comparisons$rfscore <- predict(bulletxtrctr::rtrees, newdata = pink_l3_comparisons, type = "prob")[,2]

rfscores_pink_l3 <- pink_l3_comparisons %>% select(land1, land2, rfscore)
saveRDS(rfscores_pink_l3, "data/variability_scans/pink_l3_rfscores.rda")


# Barrel-Land P-4
pink_l4 <- pink_sigs %>% filter(unique_id == "BL P-4") 
pink_lands4 <- unique(pink_l4$scan_id) 
pink_l4_comparisons <- data.frame(
  expand.grid(land1 = pink_lands4, land2 = pink_lands4), stringsAsFactors = FALSE)
pink_l4_comparisons <- pink_l4_comparisons %>% mutate(
  aligned = purrr::map2(.x = land1, .y = land2, .f = function(xx, yy) {
    land1 <- pink_l4$sigs[pink_l4$scan_id == xx][[1]]
    land2 <- pink_l4$sigs[pink_l4$scan_id == yy][[1]]
    land1$bullet <- "first-land"
    land2$bullet <- "second-land"
    
    sig_align(land1$sig, land2$sig)
  })
)
pink_l4_comparisons <- pink_l4_comparisons %>% mutate(
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

pink_l4_comparisons <- pink_l4_comparisons %>% mutate(
  striae = aligned %>% purrr::map(.f = sig_cms_max, span = 75) 
)
pink_l4_comparisons <- pink_l4_comparisons %>% mutate(
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


pink_l4_comparisons <- pink_l4_comparisons %>% mutate(
  features = purrr::map2(.x = aligned, .y = striae, .f = extract_features_all, resolution = 0.645)
)

pink_l4_comparisons <- pink_l4_comparisons %>% mutate(
  legacy_features = purrr::map(striae, extract_features_all_legacy, resolution = 0.645)
)

pink_l4_comparisons <- pink_l4_comparisons %>% tidyr::unnest(legacy_features) 
# scale features before using them in the random forest, legacy features can be used out of the box
pink_l4_comparisons$rfscore <- predict(bulletxtrctr::rtrees, newdata = pink_l4_comparisons, type = "prob")[,2]

rfscores_pink_l4 <- pink_l4_comparisons %>% select(land1, land2, rfscore)
saveRDS(rfscores_pink_l4, "data/variability_scans/pink_l4_rfscores.rda")


# Barrel-Land P-5
pink_l5 <- pink_sigs %>% filter(unique_id == "BL P-5") 
pink_lands5 <- unique(pink_l5$scan_id) 
pink_l5_comparisons <- data.frame(
  expand.grid(land1 = pink_lands5, land2 = pink_lands5), stringsAsFactors = FALSE)
pink_l5_comparisons <- pink_l5_comparisons %>% mutate(
  aligned = purrr::map2(.x = land1, .y = land2, .f = function(xx, yy) {
    land1 <- pink_l5$sigs[pink_l5$scan_id == xx][[1]]
    land2 <- pink_l5$sigs[pink_l5$scan_id == yy][[1]]
    land1$bullet <- "first-land"
    land2$bullet <- "second-land"
    
    sig_align(land1$sig, land2$sig)
  })
)
pink_l5_comparisons <- pink_l5_comparisons %>% mutate(
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

pink_l5_comparisons <- pink_l5_comparisons %>% mutate(
  striae = aligned %>% purrr::map(.f = sig_cms_max, span = 75) 
)
pink_l5_comparisons <- pink_l5_comparisons %>% mutate(
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


pink_l5_comparisons <- pink_l5_comparisons %>% mutate(
  features = purrr::map2(.x = aligned, .y = striae, .f = extract_features_all, resolution = 0.645)
)

pink_l5_comparisons <- pink_l5_comparisons %>% mutate(
  legacy_features = purrr::map(striae, extract_features_all_legacy, resolution = 0.645)
)

pink_l5_comparisons <- pink_l5_comparisons %>% tidyr::unnest(legacy_features) 
# scale features before using them in the random forest, legacy features can be used out of the box
pink_l5_comparisons$rfscore <- predict(bulletxtrctr::rtrees, newdata = pink_l5_comparisons, type = "prob")[,2]

rfscores_pink_l5 <- pink_l5_comparisons %>% select(land1, land2, rfscore)
saveRDS(rfscores_pink_l5, "data/variability_scans/pink_l5_rfscores.rda")


# Barrel-Land P-6
pink_l6 <- pink_sigs %>% filter(unique_id == "BL P-6") 
pink_lands6 <- unique(pink_l6$scan_id) 
pink_l6_comparisons <- data.frame(
  expand.grid(land1 = pink_lands6, land2 = pink_lands6), stringsAsFactors = FALSE)
pink_l6_comparisons <- pink_l6_comparisons %>% mutate(
  aligned = purrr::map2(.x = land1, .y = land2, .f = function(xx, yy) {
    land1 <- pink_l6$sigs[pink_l6$scan_id == xx][[1]]
    land2 <- pink_l6$sigs[pink_l6$scan_id == yy][[1]]
    land1$bullet <- "first-land"
    land2$bullet <- "second-land"
    
    sig_align(land1$sig, land2$sig)
  })
)
pink_l6_comparisons <- pink_l6_comparisons %>% mutate(
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

pink_l6_comparisons <- pink_l6_comparisons %>% mutate(
  striae = aligned %>% purrr::map(.f = sig_cms_max, span = 75) 
)
pink_l6_comparisons <- pink_l6_comparisons %>% mutate(
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


pink_l6_comparisons <- pink_l6_comparisons %>% mutate(
  features = purrr::map2(.x = aligned, .y = striae, .f = extract_features_all, resolution = 0.645)
)

pink_l6_comparisons <- pink_l6_comparisons %>% mutate(
  legacy_features = purrr::map(striae, extract_features_all_legacy, resolution = 0.645)
)

pink_l6_comparisons <- pink_l6_comparisons %>% tidyr::unnest(legacy_features) 
# scale features before using them in the random forest, legacy features can be used out of the box
pink_l6_comparisons$rfscore <- predict(bulletxtrctr::rtrees, newdata = pink_l6_comparisons, type = "prob")[,2]

rfscores_pink_l6 <- pink_l6_comparisons %>% select(land1, land2, rfscore)
saveRDS(rfscores_pink_l6, "data/variability_scans/pink_l6_rfscores.rda")
