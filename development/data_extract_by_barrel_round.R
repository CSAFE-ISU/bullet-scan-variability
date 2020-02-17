library(tidyverse)
library(bulletxtrctr)
library(x3ptools)
library(zoo)


args <- commandArgs(TRUE)
barrel_id <- args[1] # should be "Orange", "Green", or "Pink"
round_id <- args[2] # should be 1, 2, 3, 4 or 5

data_filename <- paste0(barrel_id, "_round", round_id, "_cc_update.rda")

ccdata_ten_crosscut <- function(ccdata, y = NULL, range = 1e-05){
  #x3pdat <- bulletxtrctr::check_x3p(x3p)
  x3p_df <- na.trim(ccdata)
  ys <- unique(x3p_df$y)
  if (is.null(y))
    y <- median(ys)
  ind <- which.min(abs(y - ys))
  y_min <- ys[ind - 5]
  y_max <- ys[ind + 4]
  lower <- min(y_min, y_max)
  upper <- max(y_min, y_max)
  x3p_df_fix <- x3p_df[x3p_df$y >= lower & x3p_df$y <= upper,]
  return(na.omit(x3p_df_fix))
}


round_id_full <- paste0("Round ", as.character(round_id))
barrel_id_full <- paste0("Barrel ", barrel_id)

filenames <- list.files("/media/Raven/Variability/", pattern = "*.x3p", full.names = T, recursive = T)
scans <- data.frame(source = filenames)

scans <- scans %>%
  mutate(operator = purrr::map_chr(as.character(source), .f = function(filename){
    strsplit(strsplit(filename, "//")[[1]][2], "/")[[1]][1]
  }),
  round = purrr::map_chr(as.character(source), .f = function(filename){
    str_trim(strsplit(filename, "-")[[1]][2])
  }), 
  machine = purrr::map_chr(as.character(source), .f = function(filename){
    str_trim(strsplit(filename, "-")[[1]][6])
  }), 
  barrel = purrr::map_chr(as.character(source), .f = function(filename){
    str_trim(strsplit(filename, "-")[[1]][3])
  }), 
  bullet = purrr::map_chr(as.character(source), .f = function(filename){
    str_trim(strsplit(filename, "-")[[1]][4])
  }), 
  land = purrr::map_chr(as.character(source), .f = function(filename){
    str_trim(strsplit(filename, "-")[[1]][5])
  }))

scans <- scans %>% filter(barrel == barrel_id_full,round == round_id_full)

scans <- scans %>% mutate(x3p = purrr::map(as.character(source), .f = function(filename){
  x3ptools::read_x3p(filename)
}))


scans <- scans %>% mutate(sizeX = purrr::map_dbl(x3p, .f = function(x3p){
  as.numeric(x3p$matrix.info$MatrixDimension$SizeX[[1]][1])
}), 
sizeY = purrr::map_dbl(x3p, .f = function(x3p){
  as.numeric(x3p$matrix.info$MatrixDimension$SizeY[[1]][1])
}), 
date = purrr::map_chr(x3p, .f = function(x3p){
  strsplit(x3p$general.info$Date[[1]][1], "T")[[1]][1]
}),
time = purrr::map_chr(x3p, .f = function(x3p){
  strsplit(x3p$general.info$Date[[1]][1], "T")[[1]][2]
}))


scans <- scans %>% mutate(
  x3p_increment_small = purrr::map_lgl(x3p, .f = function(x) x$header.info$incrementX == 0.000000645)
)

scans <- scans %>% mutate(
  x3p = purrr::map2(.x = x3p, .y = x3p_increment_small, .f = function(x,y){
    x3p <- if(y == T) {x3p_m_to_mum(x)} else{x}
    return(x3p)
  }))

scans <- scans %>% mutate(ccdata = purrr::map(x3p, .f = function(x3p){
  x3ptools::x3p_to_df(x3p)
}))


crosscuts <- scans %>%
  mutate(
    crosscut_safe = purrr::map(x3p, .f = purrr::safely(x3p_crosscut_optimize)),
    crosscut = purrr::map_dbl(crosscut_safe, "result"))#,
#crosscut = purrr::map_dbl(crosscut, function(x) ifelse(is.null(x), 200, x))
#)

crosscut_by_op <- crosscuts %>% group_by(operator) %>% summarise(median_cc = median(crosscut, na.rm = T))

crosscuts <- crosscuts %>% mutate(crosscut = purrr::map2_dbl(.x = crosscut, .y = operator, .f = function(x,y){
  crosscut_new <- if(is.na(x)){crosscut_by_op[crosscut_by_op$operator == y,]$median_cc} else{x}
  return(crosscut_new)
}))

crosscuts <- crosscuts %>% select(-x3p)

crosscuts <- crosscuts %>% mutate(
  ccdata = purrr::map2(.x = ccdata, .y = crosscut, 
                       .f = ccdata_ten_crosscut)
)


crosscuts <- crosscuts %>% mutate(
  ccdata = purrr::map(ccdata, .f = function(ccdata){
    ccdata %>% group_by(x) %>% summarise(y = mean(y, na.rm = T), value = mean(value, na.rm = T))
  })
)

saveRDS(crosscuts, file = paste0("data/variability_scans/", data_filename)) 
