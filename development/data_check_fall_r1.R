library(tidyverse)
library(bulletxtrctr)
library(x3ptools)
library(gridExtra)
library(randomForest)
library(zoo)

mark_crosscuts <- readRDS("data/variability_scans/Mark_round1_cc.rda")
connor_crosscuts <- readRDS("data/variability_scans/Connor_round1_cc.rda")

crosscuts <- rbind(connor_crosscuts, mark_crosscuts)

crosscuts <- crosscuts %>% 
  mutate(scan_id = paste0(barrel, "-", land, "-", bullet, "-", 
                          operator, "-", machine, "-", round))

crosscuts <- crosscuts %>% mutate(
  grooves = ccdata %>% 
    purrr::map(.f = cc_locate_grooves, method = "lassofull", 
               adjust = 10, return_plot = TRUE)
)



ccdata <- crosscuts$ccdata
grooves <- crosscuts$grooves



library(shiny)

#if (interactive()) {
shinyApp(
  
  ui = fluidPage(
    selectInput("k","Investigate kth plot:", selected = 1,
                choices=1:length(grooves)),
    textOutput("groovelocations"),
    actionButton("confirm", "Confirm"),
    actionButton("save", "Save"),
    plotOutput("groovePlot", click = "plot_click"),
    verbatimTextOutput("info")
  ),
  
  server = function(input, output, session) {
    output$groovePlot <- renderPlot({
      k <- as.numeric(input$k)
      p <- grooves[[k]]$plot 
      
      p
    })
    output$groovelocations <- renderText({
      paste("Left Groove: ",grooves[[as.numeric(input$k)]]$groove[1], 
            " Right Groove: ",grooves[[as.numeric(input$k)]]$groove[2])
    })
    observeEvent(input$confirm,{
      cat(str(input$k))
      updateSelectInput(session, "k","Investigate kth plot:", 
                        selected = as.numeric(input$k)+1,
                        choices=1:length(grooves))
    })
    observeEvent(input$save,{
      saveRDS(grooves, file="data/variability_scans/Mark_check_grooves.rda")
      cat("groove data saved\n")
    })
    
    observeEvent(input$plot_click,{
      k <- as.numeric(input$k)
      xloc <- input$plot_click$x
      
      gr <- grooves[[k]]$groove
      if (abs(gr[1]-xloc) < abs(gr[2]-xloc)) {
        grooves[[k]]$groove[1] <<- xloc
      } else {
        grooves[[k]]$groove[2] <<- xloc
      }
      output$groovePlot <- renderPlot({ 
        k <- as.numeric(input$k)
        p <- grooves[[k]]$plot + 
          geom_vline(xintercept = grooves[[k]]$groove[1], colour="green") +
          geom_vline(xintercept = grooves[[k]]$groove[2], colour="green")
        
        p
      })
      
    })
    output$info <- renderText({
      paste0("x=", input$plot_click$x, "\ny=", input$plot_click$y)
    })
    
  },
  
  options = list(height = 500)
)
#}
#saveRDS(grooves, file="data/variability_scans/Mark_check_grooves.rda")

grooves <- readRDS("data/variability_scans/Mark_check_grooves.rda")
crosscuts$grooves <- grooves


crosscuts <- crosscuts %>% mutate(
  sigs = purrr::map2(
    .x = ccdata, .y = grooves, 
    .f = function(x, y) {
      cc_get_signature(
        ccdata = x, grooves = y, span1 = 0.75, span2 = 0.03)
    })
)

sigs <- crosscuts %>% select(-ccdata)
saveRDS(sigs, "data/variability_scans/Mark_check_sigs.rda")


#### PINK ####
pink_sigs <- sigs %>% filter(barrel == "Barrel Pink")

pink_lands <- unique(pink_sigs$scan_id) 
pink_comparisons <- data.frame(
  expand.grid(land1 = pink_lands, land2 = pink_lands), stringsAsFactors = FALSE)
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



#### ORANGE #### 
orange_sigs <- sigs %>% filter(barrel == "Barrel Orange")

orange_lands <- unique(orange_sigs$scan_id) 
orange_comparisons <- data.frame(
  expand.grid(land1 = orange_lands, land2 = orange_lands), stringsAsFactors = FALSE)
orange_comparisons <- orange_comparisons %>% mutate(
  aligned = purrr::map2(.x = land1, .y = land2, .f = function(xx, yy) {
    land1 <- orange_sigs$sigs[orange_sigs$scan_id == xx][[1]]
    land2 <- orange_sigs$sigs[orange_sigs$scan_id == yy][[1]]
    land1$bullet <- "first-land"
    land2$bullet <- "second-land"
    
    sig_align(land1$sig, land2$sig)
  })
)

orange_comparisons <- orange_comparisons %>% mutate(
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

orange_comparisons <- orange_comparisons %>% mutate(
  striae = aligned %>% purrr::map(.f = sig_cms_max, span = 75) 
)
orange_comparisons <- orange_comparisons %>% mutate(
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

orange_comparisons <- orange_comparisons %>% mutate(
  features = purrr::map2(.x = aligned, .y = striae, .f = extract_features_all, resolution = 0.645)
)

orange_comparisons <- orange_comparisons %>% mutate(
  legacy_features = purrr::map(striae, extract_features_all_legacy, resolution = 0.645)
)

orange_comparisons <- orange_comparisons %>% tidyr::unnest(legacy_features) 
# scale features before using them in the random forest, legacy features can be used out of the box
orange_comparisons$rfscore <- predict(bulletxtrctr::rtrees, newdata = orange_comparisons, type = "prob")[,2]

rfscores_orange <- orange_comparisons %>% select(land1, land2, rfscore)




rfscores <- rbind(rfscores_pink, rfscores_orange)
saveRDS(rfscores, "data/variability_scans/mark_check_rfscores.rda")


rfscores <- rfscores %>% 
  mutate(id_S1 = land1, id_S2 = land2) %>% 
  separate(col = id_S1, 
           into = c("barrel_S1", "land_S1", "bullet_S1", 
                    "operator_S1", "machine_S1", "round_S1"), sep = "-") %>% 
  separate(col = id_S2, 
           into = c("barrel_S2", "land_S2", "bullet_S2", 
                    "operator_S2", "machine_S2", "round_S2"), sep = "-")

rfscores %>%
  filter(barrel_S1 == "Barrel Pink", bullet_S1 == "Bullet 1") %>%
  filter(operator_S1 == "Connor" & operator_S2 == "Mark") %>% 
  filter(machine_S1 == "Sneox1") %>% 
  ggplot(aes(x = land_S2, y = land1, fill = rfscore)) +
  geom_tile() +
  scale_fill_gradient2(low = "grey80", high = "darkorange", 
                       midpoint = .5) +
  facet_wrap(~machine_S2+bullet_S2) +
  #xlab("Land A") +
  #ylab("Land B") +
  theme(aspect.ratio = 1, axis.text.x = element_text(angle =20))

rfscores %>%
  filter(barrel_S1 == "Barrel Orange", bullet_S1 == "Bullet 1") %>%
  filter(operator_S1 == "Connor" & operator_S2 == "Mark") %>% 
  filter(machine_S1 == "Sneox1") %>% 
  ggplot(aes(x = land_S2, y = land1, fill = rfscore)) +
  geom_tile() +
  scale_fill_gradient2(low = "grey80", high = "darkorange", 
                       midpoint = .5) +
  facet_wrap(~machine_S2+bullet_S2) +
  #xlab("Land A") +
  #ylab("Land B") +
  theme(aspect.ratio = 1, axis.text.x = element_text(angle =20))

