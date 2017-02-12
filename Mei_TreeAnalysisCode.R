library(dplyr)
library(tidyr)
library(magrittr)
library(zoo)
rm(list = ls()); gc(reset = TRUE)

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

system.time({
  
fit_ht <- iForm_FunctionalMap(formula = HT ~ .,
                      data = dat[,-3],
                      id = "Tree",
                      time_col = "time",
                      heredity = "weak",
                      higher_order = FALSE)
  
})


system.time({
  
  fit_dia <- iForm_FunctionalMap(formula = DIA ~ .,
                                data = dat[,-2],
                                id = "Tree",
                                time_col = "time",
                                heredity = "weak",
                                higher_order = FALSE)
  
})