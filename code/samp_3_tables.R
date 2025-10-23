
################################################################################
#
# Project: Big Data Sampling
#
# Purpose: Combine results into tables
#
# Last Update: 25 Jul 2025
#
################################################################################

library("tidyverse")
library("tidyselect")

analysis <- "cif" # "cif", "cox", "poisson


# Read in RR results ------------------------------------------------------

read_res <- function (samp) {
  
  if (samp %in% c("full", "sc25", "sc10", "cc25", "cc10")) {
    
    res.pe <- read_csv(paste0("../results/res_cif_", samp, "_pointest.csv")) %>% 
      rename(point_est_time = comp_time) %>% 
      select(-rep)
    
    res.b <- read_csv(paste0("../results/res_cif_", samp, "_500.csv")) %>% 
      summarize(boot_time = mean(comp_time),
                SE1 = sd(RR1), SE2 = sd(RR2), SE3 = sd(RR3), SE4 = sd(RR4), SE5 = sd(RR5))
    
    res <- bind_cols(res.pe, res.b)
    
  }  else if (samp %in% c("divide10", "divide20", "divide50")) {

    res.pe <- read_csv(paste0("../results/res_cif_", samp, "_pointest.csv")) %>% 
      select(analysis, RR1, RR2, RR3, RR4, RR5, comp_time) %>%
      rename(point_est_time = comp_time) 
    
    res.b <- read_csv(paste0("../results/res_cif_", samp, ".csv")) %>% 
      select(SE1, SE2, SE3, SE4, SE5, comp_time) %>%
      rename(boot_time = comp_time)
    
    res <- bind_cols(res.pe, res.b)
    
  } else if (samp=="divide10p") {
    
    res <- read_csv(paste0("../results/res_cif_", samp, "_pointest.csv")) %>% 
      select(-starts_with("n_invalid")) %>% 
      rename(point_est_time = comp_time) %>% 
      mutate(boot_time = NA)
  
  } else if (samp %in% c("divide20p", "divide50p")) {
    
    res.pe1 <- read_csv(paste0("../results/res_cif_", samp, "_pointest.csv")) %>% 
      select(-starts_with("n_invalid")) %>% 
      rename(point_est_time = comp_time) %>% 
      mutate(boot_time = NA,
             analysis = paste0(samp, "_10cores"))
    
    res.pe2 <- read_csv(paste0("../results/res_cif_", samp, "_pointest2.csv")) %>% 
      select(-starts_with("n_invalid")) %>% 
      rename(point_est_time = comp_time) %>% 
      mutate(boot_time = NA,
             analysis = paste0(samp, "_mcores"))
    
    res <- bind_rows(res.pe1, res.pe2)
    
  }
  
  return(res)
}


res <- lapply(c("full", "divide10", "divide10p", "divide20", "divide20p", "divide50", "divide50p",
                "sc25", "sc10", "cc25", "cc10"), function(x){read_res(x)})
res <- bind_rows(res)

write_csv(res, paste0("../results/res_cif.csv"))


# Read in HR results ------------------------------------------------------

read_res <- function (samp) {
  
  res <- read_csv(paste0("../results/res_cox_", samp, ".csv")) 
    
  return(res)
}


res <- lapply(c("full", "divide10", "divide10p", "divide20", "divide20p", "divide50", "divide50p",
                "sc25", "sc10", "cc25", "cc10"), function(x){read_res(x)})
res <- bind_rows(res)

write_csv(res, paste0("../results/res_cox.csv"))


# Read in IRR results -----------------------------------------------------

read_res <- function (samp) {
  
  if (samp %in% c("full", "sc25", "sc10", "cc25", "cc10")) {
    
    res.pe <- read_csv(paste0("../results/res_poisson_", samp, "_pointest.csv")) %>% 
      rename(point_est_time = comp_time) %>% 
      select(-rep)
    
    res.b <- read_csv(paste0("../results/res_poisson_", samp, "_500.csv")) %>% 
      summarize(boot_time = mean(comp_time),
                SE = sd(IRR))
    
    res <- bind_cols(res.pe, res.b)
    
  }  else if (samp %in% c("divide10", "divide20", "divide50")) {
    
    res.pe <- read_csv(paste0("../results/res_poisson_", samp, "_pointest.csv")) %>% 
      select(analysis, IRR, comp_time) %>%
      rename(point_est_time = comp_time) 
    
    res.b <- read_csv(paste0("../results/res_poisson_", samp, ".csv")) %>% 
      select(SE, comp_time) %>%
      rename(boot_time = comp_time)
    
    res <- bind_cols(res.pe, res.b)
    
  } else if (samp=="divide10p") {
    
    res <- read_csv(paste0("../results/res_poisson_", samp, "_pointest.csv")) %>% 
      rename(point_est_time = comp_time) %>% 
      mutate(boot_time = NA)
    
  } else if (samp %in% c("divide20p", "divide50p")) {
    
    res.pe1 <- read_csv(paste0("../results/res_poisson_", samp, "_pointest.csv")) %>% 
      rename(point_est_time = comp_time) %>% 
      mutate(boot_time = NA,
             analysis = paste0(samp, "_10cores"))
    
    res.pe2 <- read_csv(paste0("../results/res_poisson_", samp, "_pointest2.csv")) %>% 
      rename(point_est_time = comp_time) %>% 
      mutate(boot_time = NA,
             analysis = paste0(samp, "_mcores"))
    
    res <- bind_rows(res.pe1, res.pe2)
    
  }
  
  return(res)
}


res <- lapply(c("full", "divide10", "divide10p", "divide20", "divide20p", "divide50", "divide50p",
                "sc25", "sc10", "cc25", "cc10"), function(x){read_res(x)})
res <- bind_rows(res)

write_csv(res, paste0("../results/res_poisson.csv"))

