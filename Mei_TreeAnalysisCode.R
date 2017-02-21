library(dplyr)
library(tidyr)
library(magrittr)
library(zoo)
rm(list = ls()); gc(reset = TRUE)
source("iFormFunctionalMapping_OLS.R")

mei <- read.csv("data.meitree.raw.longformat.csv", stringsAsFactors = FALSE)
markers <- read.csv("MeiTrees_imputed_merged.csv")
markers <- markers[,c(2,9:5562)]


mei <- mei[grepl("^L",mei$progeny), ]
mei$Tree <- substr(mei$progeny, 2, 5)
mei$Date <- factor(mei$Date)

mei <- mei %>%
  group_by(Tree, Date) %>%
  summarise(HT = mean(as.numeric(HT), na.rm = TRUE), DIA = mean(as.numeric(DIA), na.rm = TRUE)) %>%
  mutate(time = as.numeric(Date))
mei$HT <- na.locf(mei$HT)
mei$DIA <- na.locf(mei$DIA)

dat <- merge(mei, markers, by = "Tree")
dat <- dat[, -2]
dat <- dat %>% arrange(Tree, time)


markers_want <- c("AATTC_nn_np_2517_a", "AATTC_nn_np_2815_a", "CATG_nn_np_3479_a", "CATG_nn_np_1284_a",
                  "CATG_lm_ll_3153_a", "AATTC_lm_ll_3034_a", "AATTC_hk_hk_278_a", "AATTC_hk_hk_479_d",
                  "AATTC_nn_np_554_a", "AATTC_nn_np_1615_a", "AATTC_nn_np_929_a", "CATG_hk_hk_648_a")

sort(match(markers_want, names(dat)))

samp <- sample(names(dat), 10)
test_marks <- union(markers_want, samp)

ht_run <- c("Tree", "time", "HT", test_marks)
dat_ht <- dat[,ht_run]

system.time({
  
fit_ht <- iForm_FunctionalMap(formula = HT ~ .,
                      data = dat_ht,
                      id_col = "Tree",
                      time_col = "time",
                      heredity = "strong",
                      higher_order = FALSE,
                      poly_num = 6)
  
})

lapply(markers_want, function(x) grepl(x, names(coef(fit_ht$fit))))


dia_run <- c("Tree", "time", "DIA", markers_want)
dat_dia <- dat[,dia_run]

system.time({
  
  fit_dia <- iForm_FunctionalMap(formula = DIA ~ .,
                                data = dat_dia,
                                id = "Tree",
                                time_col = "time",
                                heredity = "weak",
                                higher_order = FALSE,
                                poly_num = 5)
  
})