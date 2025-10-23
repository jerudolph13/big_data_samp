
################################################################################
#
# Project: Big Data Sampling
#
# Purpose: Summarize across the repeat D&R analyses
#
# Last Update: 05 Aug 2024
#
################################################################################

library("tidyverse")

res1 <- read_csv("~/sampling/results/res_cox_divide10p_repeat.csv") %>%
  summarize(hr.mean = mean(HR))

res2 <- read_csv("~/sampling/results/res_cox_divide50p_repeat.csv") %>%
  summarize(hr.mean = mean(HR))

res3 <- read_csv("~/sampling/results/res_cif_divide10p_repeat.csv") %>%
  summarize(rr1.mean = mean(RR1),
            rr2.mean = mean(RR2),
            rr3.mean = mean(RR3),
            rr4.mean = mean(RR4),
            rr5.mean = mean(RR5))

res4 <- read_csv("~/sampling/results/res_cif_divide50p_repeat.csv") %>%
  summarize(rr1.med = median(RR1),
            rr2.med = median(RR2),
            rr3.med = median(RR3),
            rr4.med = median(RR4),
            rr5.med = median(RR5))

