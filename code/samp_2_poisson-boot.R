
################################################################################
#
# Project: Big Data Sampling
#
# Purpose: Run risk analysis
#
# Last Update: 03 Jun 2025
#
################################################################################

lib <- "~/R/4.3"
packages <- c("tidyverse", "splines", "survival", "broom", "meta", "parallel")
for (package in packages){
  if (!(package %in% installed.packages(lib=lib))) {install.packages(package, lib=lib)}
  suppressPackageStartupMessages(library(package, character.only=T, lib.loc=lib, quietly=T))
}

args <- commandArgs(trailingOnly=TRUE)
analysis <- as.character(args[1])
boot_start <- 0
boot_end <- 500
cores <- 10


# Read in data ------------------------------------------------------------

dat <- read_csv("/cms01/data/dua/57285/users/c-jrudolp9-57285/samp_data.csv") %>%
  mutate(across(c(hiv, sex, race, state, enroll_period, n_cond), ~ as.factor(.x)),
         delta = as.numeric(first_event=="lung"))

start_time <- Sys.time()


# Bootstrap resample ------------------------------------------------------

boot_rep <- function(r) {
  
  set.seed(r)
  
  samp <- table(dat[sample(1:nrow(dat), nrow(dat), replace=T), (names(dat)=="bene_id")])
  
  boot <- NULL
  if (r==0) {
    boot <- dat %>%
      rename(bid = bene_id)
  } else {
    for(zzz in 1:max(samp)){
      cc <- dat[dat$bene_id %in% names(samp[samp %in% c(zzz:max(samp))]),]
      cc$bid <- paste0(cc$bene_id, zzz)
      boot <- rbind(boot, cc)
    }
    boot <- select(boot, -bene_id)
  }
  
  
  # Sample and analyze ------------------------------------------------------
  
  
  # Function for modeling
  models <- function(data) {
    
    est <- data.frame(rep = NA,
                      IRR = NA)
    
    # Weighted Poisson model
    if (analysis %in% c("cc10", "cc25")) {
      
      unique <- filter(data, !duplicated(bid))
      
      mod_hiv <- glm(hiv ~ bs(age, df=3) + race + sex + state + enroll_period,
                     family=binomial(link="logit"), weights=wt, data=unique)$fitted.values
      den_hiv <- ifelse(unique$hiv==1, mod_hiv, 1 - mod_hiv)
      
      mod_hiv <- glm(hiv ~ 1,
                     family=binomial(link="logit"), weights=wt, data=unique)$fitted.values
      num_hiv <- ifelse(unique$hiv==1, mod_hiv, 1 - mod_hiv)
      
      unique$wt2 <- num_hiv/den_hiv
      data <- left_join(data, select(unique, bid, wt2), by="bid") %>%
        mutate(wt3 = wt*wt2)
      
      mod <- glm(I(first_event=="lung") ~ hiv + offset(log(person_yrs)), weights=wt3, family=poisson(link="log"), data=data)
      res <- tidy(mod)
      
    } else {
      
      mod_hiv <- glm(hiv ~ bs(age, df=3) + race + sex + state + enroll_period,
                     family=binomial(link="logit"), weights=wt, data=data)$fitted.values
      den_hiv <- ifelse(data$hiv==1, mod_hiv, 1 - mod_hiv)
      
      mod_hiv <- glm(hiv ~ 1,
                     family=binomial(link="logit"), weights=wt, data=data)$fitted.values
      num_hiv <- ifelse(data$hiv==1, mod_hiv, 1 - mod_hiv)
      
      data$wt2 <- (num_hiv/den_hiv)*data$wt    
      
      mod <- glm(I(first_event=="lung") ~ hiv + offset(log(person_yrs)), weights=wt2, family=poisson(link="log"), data=data)
      res <- tidy(mod)
      
    }
    
    est$rep <- r
    est$IRR <- res$estimate[2]
    
    return(est)
    
  }
  
  
  # Run analysis
  if (analysis=="full") {
    
    samp <- boot %>%
      mutate(wt = 1)
    
    res <- models(samp)
    
  } else if (analysis=="sc25") {
    
    p1 <- 1
    p0 <- 0.25
    
    dat1 <- filter(boot, hiv==1) %>%
      mutate(wt = 1/p1)
    dat0 <- filter(boot, hiv==0)[rbinom(nrow(filter(boot, hiv==0)), 1, p0)==1, ] %>%
      mutate(wt = 1/p0)
    
    samp <- bind_rows(dat1, dat0)
    
    res <- models(samp)
    
  } else if (analysis=="sc10") {
    
    p1 <- 1
    p0 <- 0.10
    
    dat1 <- filter(boot, hiv==1) %>%
      mutate(wt = 1/p1)
    dat0 <- filter(boot, hiv==0)[rbinom(nrow(filter(boot, hiv==0)), 1, p0)==1, ] %>%
      mutate(wt = 1/p0)
    
    samp <- bind_rows(dat1, dat0)
    
    res <- models(samp)
    
  } else if (analysis=="cc25") {
    
    p1 <- 1
    p0 <- 0.25
    
    dat1 <- filter(boot, hiv==1) %>%
      mutate(wt = 1/p1,
             end = ifelse(first_event=="lung", person_yrs - 0.001, person_yrs),
             start = 0)
    dat0 <- filter(boot, hiv==0)[rbinom(nrow(filter(boot, hiv==0)), 1, p0)==1, ] %>%
      mutate(wt = 1/p0,
             end = ifelse(first_event=="lung", person_yrs - 0.001, person_yrs),
             start = 0)
    cases <- filter(boot, first_event=="lung") %>%
      mutate(wt = 1,
             end = person_yrs,
             start = person_yrs - 0.001)
    
    samp <- bind_rows(dat1, dat0, cases)
    
    res <- models(samp)
    
  } else if (analysis=="cc10") {
    
    p1 <- 1
    p0 <- 0.10
    
    dat1 <- filter(boot, hiv==1) %>%
      mutate(wt = 1/p1,
             end = ifelse(first_event=="lung", person_yrs - 0.001, person_yrs),
             start = 0)
    dat0 <- filter(boot, hiv==0)[rbinom(nrow(filter(boot, hiv==0)), 1, p0)==1, ] %>%
      mutate(wt = 1/p0,
             end = ifelse(first_event=="lung", person_yrs - 0.001, person_yrs),
             start = 0)
    cases <- filter(boot, first_event=="lung") %>%
      mutate(wt = 1,
             end = person_yrs,
             start = person_yrs - 0.001)
    
    samp <- bind_rows(dat1, dat0, cases)
    
    res <- models(samp)
    
  }
  
  return(res)
  
}

#all.res <- lapply(0:nboot, function(x) {boot_rep(x)})
all.res <- mclapply(boot_start:boot_end, function(x) {boot_rep(x)}, mc.cores=cores, mc.set.seed=F)
all.res <- bind_rows(all.res)

all.res$analysis <- analysis
end_time <- Sys.time()
all.res$comp_time <- difftime(end_time, start_time, units="hours")


# Output results ----------------------------------------------------------

if (boot_end==0) {
  write_csv(all.res, paste0("./results/res_poisson_", analysis, "_pointest.csv"))
} else {
  write_csv(all.res, paste0("./results/res_poisson_", analysis, "_", boot_end, ".csv"))
}


