
################################################################################
#
# Project: Big Data Sampling
#
# Purpose: Run risk analysis (full, sub-cohort, case-cohort)
#
# Last Update: 26 Jun 2025
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
boot_end <- 0
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
    
    est <- data.frame(rep = r,
                      RR1 = NA, RR2 = NA, RR3 = NA, RR4 = NA, RR5 = NA)
    
    # Weighted cumulative incidence function
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
      
      mod.cif <- survfit(Surv(start, end, delta) ~  hiv,
                         weights=wt3, id=bid, data=data, conf.type="none", se.fit=F, timefix=F)
      
    } else {
      
      mod_hiv <- glm(hiv ~ bs(age, df=3) + race + sex + state + enroll_period,
                     family=binomial(link="logit"), weights=wt, data=data)$fitted.values
      den_hiv <- ifelse(data$hiv==1, mod_hiv, 1 - mod_hiv)
      
      mod_hiv <- glm(hiv ~ 1,
                     family=binomial(link="logit"), weights=wt, data=data)$fitted.values
      num_hiv <- ifelse(data$hiv==1, mod_hiv, 1 - mod_hiv)
      
      data$wt2 <- (num_hiv/den_hiv)*data$wt    
      
      mod.cif <- survfit(Surv(person_yrs, delta) ~  hiv,
                         weights=wt2, data=data, conf.type="none", se.fit=F, timefix=F)
    }
    
    fit <- summary(mod.cif)
    res <- data.frame(time=fit$time,
                      estimate=fit$surv,
                      strata=fit$strata)
    
    est.cif1 <- res %>% 
      filter(time<=1) %>%
      filter(!duplicated(strata, fromLast=T)) %>%
      select(estimate, strata)
    
    est$RR1 <- log((1 - est.cif1$estimate[est.cif1$strata=="hiv=1"]) / 
                     (1 - est.cif1$estimate[est.cif1$strata=="hiv=0"]))
    
    est.cif2 <- res %>% 
      filter(time<=2) %>%
      filter(!duplicated(strata, fromLast=T)) %>%
      select(estimate, strata)
    
    est$RR2 <- log((1 - est.cif2$estimate[est.cif2$strata=="hiv=1"]) / 
                     (1 - est.cif2$estimate[est.cif2$strata=="hiv=0"]))
    
    est.cif3 <- res %>% 
      filter(time<=3) %>%
      filter(!duplicated(strata, fromLast=T)) %>%
      select(estimate, strata)
    
    est$RR3 <- log((1 - est.cif3$estimate[est.cif3$strata=="hiv=1"]) / 
                     (1 - est.cif3$estimate[est.cif3$strata=="hiv=0"]))
    
    est.cif4 <- res %>% 
      filter(time<=4) %>%
      filter(!duplicated(strata, fromLast=T)) %>%
      select(estimate, strata)
    
    est$RR4 <- log((1 - est.cif4$estimate[est.cif4$strata=="hiv=1"]) / 
                     (1 - est.cif4$estimate[est.cif4$strata=="hiv=0"]))
    
    est.cif5 <- res %>% 
      filter(time<=5) %>%
      filter(!duplicated(strata, fromLast=T)) %>%
      select(estimate, strata)
    
    est$RR5 <- log((1 - est.cif5$estimate[est.cif5$strata=="hiv=1"]) / 
                     (1 - est.cif5$estimate[est.cif5$strata=="hiv=0"]))
    
    return(est)
    
  }
  
  # Function for divide & recombine sampling
  samp <- function(s) {
    
    samp <- filter(boot, s==split) %>%
      mutate(wt = 1)
    
    res.k <- models(samp)
    
    return(res.k)
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

#all.res <- lapply(boot_start:boot_end, function(x) {boot_rep(x)})
all.res <- mclapply(boot_start:boot_end, function(x) {boot_rep(x)}, mc.cores=cores, mc.set.seed=F)
all.res <- bind_rows(all.res)

all.res$analysis <- analysis
end_time <- Sys.time()
all.res$comp_time <- difftime(end_time, start_time, units="hours")


# Output results ----------------------------------------------------------

if (boot_end==0) {
  write_csv(all.res, paste0("./results/res_cif_", analysis, "_pointest.csv"))
} else {
  write_csv(all.res, paste0("./results/res_cif_", analysis, "_", boot_end, ".csv"))
}

