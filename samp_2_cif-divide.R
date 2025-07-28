
################################################################################
#
# Project: Big Data Sampling
#
# Purpose: Run risk analysis (divide and recombine)
#
# Last Update: 15 Jul 2025
#
################################################################################

lib <- "~/R/4.3"
packages <- c("tidyverse", "splines", "survival", "broom", "meta", "parallel")
for (package in packages){
  if (!(package %in% installed.packages(lib=lib))) {install.packages(package, lib=lib)}
  suppressPackageStartupMessages(library(package, character.only=T, lib.loc=lib, quietly=T))
}

set.seed(123)

args <- commandArgs(trailingOnly=TRUE)
analysis <- as.character(args[1])
boot_start <- 0
boot_end <- 0
cores <- 50


# Read in data ------------------------------------------------------------

dat <- read_csv("/cms01/data/dua/57285/users/c-jrudolp9-57285/samp_data.csv") %>%
  mutate(across(c(hiv, sex, race, state, enroll_period, n_cond), ~ as.factor(.x)),
         delta = as.numeric(first_event=="lung"))

start_time <- Sys.time()


# Sample and analyze ------------------------------------------------------

# Function for modeling
models <- function(data) {
  
  # For the divide & recombine, we estimated bootstrap SE within each subsample
  # in order to use the meta-regression to recombine
  boot_rep <- function(r) {
    
    est <- data.frame(rep = r,
                      RR1 = NA, RR2 = NA, RR3 = NA, RR4 = NA, RR5 = NA)
    
    set.seed(r)
    b.samp <- table(data[sample(1:nrow(data), nrow(data), replace=T), (names(data)=="bene_id")])
    
    boot <- NULL
    if (r==0) {
      boot <- data %>%
        rename(bid = bene_id)
    } else {
      for(zzz in 1:max(b.samp)){
        cc <- data[data$bene_id %in% names(b.samp[b.samp %in% c(zzz:max(b.samp))]),]
        cc$bid <- paste0(cc$bene_id, zzz)
        boot <- rbind(boot, cc)
      }
      boot <- select(boot, -bene_id)
    }
    
    # Weighted cumulative incidence function
    mod_hiv <- glm(hiv ~ bs(age, df=3) + race + sex + state + enroll_period,
                   family=binomial(link="logit"), weights=wt, data=boot)$fitted.values
    den_hiv <- ifelse(boot$hiv==1, mod_hiv, 1 - mod_hiv)
    
    mod_hiv <- glm(hiv ~ 1,
                   family=binomial(link="logit"), weights=wt, data=boot)$fitted.values
    num_hiv <- ifelse(boot$hiv==1, mod_hiv, 1 - mod_hiv)
    
    boot$wt2 <- (num_hiv/den_hiv)*boot$wt    
    
    mod.cif <- survfit(Surv(person_yrs, delta) ~  hiv,
                       weights=wt2, data=boot, conf.type="none", se.fit=F, timefix=F)
    
    fit <- summary(mod.cif)
    res <- data.frame(time=fit$time,
                      estimate=fit$surv,
                      strata=fit$strata)
    
    est.cif1 <- res %>% 
      filter(time<=1) %>%
      filter(!duplicated(strata, fromLast=T)) %>%
      select(estimate, strata)
    
    if (nrow(est.cif1)<2) { # In case there are no events by t 
      est$RR1 <- NA
    } else {
      est$RR1 <- log((1 - est.cif1$estimate[est.cif1$strata=="hiv=1"]) / 
                       (1 - est.cif1$estimate[est.cif1$strata=="hiv=0"]))
    }
    
    est.cif2 <- res %>% 
      filter(time<=2) %>%
      filter(!duplicated(strata, fromLast=T)) %>%
      select(estimate, strata)
    
    if (nrow(est.cif2)<2) { # In case there are no events by t
      est$RR2 <- NA
    } else {
      est$RR2 <- log((1 - est.cif2$estimate[est.cif2$strata=="hiv=1"]) / 
                       (1 - est.cif2$estimate[est.cif2$strata=="hiv=0"]))
    }
    
    est.cif3 <- res %>% 
      filter(time<=3) %>%
      filter(!duplicated(strata, fromLast=T)) %>%
      select(estimate, strata)
    
    if (nrow(est.cif3)<2) { # In case there are no events by t
      est$RR3 <- NA
    } else {
      est$RR3 <- log((1 - est.cif3$estimate[est.cif3$strata=="hiv=1"]) / 
                       (1 - est.cif3$estimate[est.cif3$strata=="hiv=0"]))
    }
    
    est.cif4 <- res %>% 
      filter(time<=4) %>%
      filter(!duplicated(strata, fromLast=T)) %>%
      select(estimate, strata)
    
    if (nrow(est.cif4)<2) { # In case there are no events by t
      est$RR4 <- NA
    } else {
      est$RR4 <- log((1 - est.cif4$estimate[est.cif4$strata=="hiv=1"]) / 
                       (1 - est.cif4$estimate[est.cif4$strata=="hiv=0"]))
    }
    
    est.cif5 <- res %>% 
      filter(time<=5) %>%
      filter(!duplicated(strata, fromLast=T)) %>%
      select(estimate, strata)
    
    if (nrow(est.cif5)<2) { # In case there are no events by t
      est$RR5 <- NA
    } else {
      est$RR5 <- log((1 - est.cif5$estimate[est.cif5$strata=="hiv=1"]) / 
                       (1 - est.cif5$estimate[est.cif5$strata=="hiv=0"]))
    }
    
    return(est)
  }
  
  #all.res <- lapply(boot_start:boot_end, function(x) {boot_rep(x)})
  all.res <- mclapply(boot_start:boot_end, function(x) {boot_rep(x)}, mc.cores=cores, mc.set.seed=F)
  all.res <- bind_rows(all.res)
  
  # Summarize across bootsraps, discarding replicates in which there were no events in PWH
  res <- data.frame(RR1 = all.res$RR1[all.res$rep==0],
                    SE1 = sd(all.res$RR1, na.rm=T),
                    RR2 = all.res$RR2[all.res$rep==0],
                    SE2 = sd(all.res$RR2, na.rm=T),
                    RR3 = all.res$RR3[all.res$rep==0],
                    SE3 = sd(all.res$RR3, na.rm=T),
                    RR4 = all.res$RR4[all.res$rep==0],
                    SE4 = sd(all.res$RR4, na.rm=T),
                    RR5 = all.res$RR5[all.res$rep==0],
                    SE5 = sd(all.res$RR5, na.rm=T))
  
  return(res)
}


# Function for divide & recombine sampling
samp <- function(s) {
  
  samp <- filter(dat, s==split) %>%
    mutate(wt = 1)
  
  res.k <- models(samp)
  
  return(res.k)
}

# Run analysis
if (analysis=="divide10") {
  
  k <- 10
  dat$split <- sample(rep(1:k, ceiling(nrow(dat)/k))[1:nrow(dat)])
  
  res.samp <- lapply(1:k, function(x){samp(x)})
  res.samp <- bind_rows(res.samp)
  
  if (boot_end==0) {
    res <- data.frame(RR1 = mean(res.samp$RR1, na.rm=T),
                      SE1 = NA,
                      n_invalid1 = sum(is.na(res.samp$RR1)),
                      RR2 = mean(res.samp$RR2, na.rm=T),
                      SE2 = NA,
                      n_invalid2 = sum(is.na(res.samp$RR2)),
                      RR3 = mean(res.samp$RR3, na.rm=T),
                      SE3 = NA,
                      n_invalid3 = sum(is.na(res.samp$RR3)),
                      RR4 = mean(res.samp$RR4, na.rm=T),
                      SE4 = NA,
                      n_invalid4 = sum(is.na(res.samp$RR4)),
                      RR5 = mean(res.samp$RR5, na.rm=T),
                      SE5 = NA,
                      n_invalid5 = sum(is.na(res.samp$RR5)))
    
  } else {   
    # If a split has invalid RR (due to no events by t), drop it from metaregression 
    res.pool1 <- metagen(RR1, SE1, data=filter(res.samp, !is.na(RR1)))
    res.pool2 <- metagen(RR2, SE2, data=filter(res.samp, !is.na(RR2)))
    res.pool3 <- metagen(RR3, SE3, data=filter(res.samp, !is.na(RR3)))
    res.pool4 <- metagen(RR4, SE4, data=filter(res.samp, !is.na(RR4)))
    res.pool5 <- metagen(RR5, SE5, data=filter(res.samp, !is.na(RR5)))
    
    res <- data.frame(RR1 = res.pool1$TE.fixed,
                      SE1 = res.pool1$seTE.fixed,
                      RR2 = res.pool2$TE.fixed,
                      SE2 = res.pool2$seTE.fixed,
                      RR3 = res.pool3$TE.fixed,
                      SE3 = res.pool3$seTE.fixed,
                      RR4 = res.pool4$TE.fixed,
                      SE4 = res.pool4$seTE.fixed,
                      RR5 = res.pool5$TE.fixed,
                      SE5 = res.pool5$seTE.fixed)
  }
  
} else if (analysis=="divide10p") {
  
  k <- 10
  dat$split <- sample(rep(1:k, ceiling(nrow(dat)/k))[1:nrow(dat)])
  
  res.samp <- mclapply(1:k, function(x){samp(x)}, mc.cores=cores, mc.set.seed=F)
  res.samp <- bind_rows(res.samp)
  
  if (boot_end==0) {
    res <- data.frame(RR1 = mean(res.samp$RR1, na.rm=T),
                      SE1 = NA,
                      n_invalid1 = sum(is.na(res.samp$RR1)),
                      RR2 = mean(res.samp$RR2, na.rm=T),
                      SE2 = NA,
                      n_invalid2 = sum(is.na(res.samp$RR2)),
                      RR3 = mean(res.samp$RR3, na.rm=T),
                      SE3 = NA,
                      n_invalid3 = sum(is.na(res.samp$RR3)),
                      RR4 = mean(res.samp$RR4, na.rm=T),
                      SE4 = NA,
                      n_invalid4 = sum(is.na(res.samp$RR4)),
                      RR5 = mean(res.samp$RR5, na.rm=T),
                      SE5 = NA,
                      n_invalid5 = sum(is.na(res.samp$RR5)))
    
  } else {
    # If a split has invalid RR (due to no events by t), drop it from metaregression
    res.pool1 <- metagen(RR1, SE1, data=filter(res.samp, !is.na(RR1)))
    res.pool2 <- metagen(RR2, SE2, data=filter(res.samp, !is.na(RR2)))
    res.pool3 <- metagen(RR3, SE3, data=filter(res.samp, !is.na(RR3)))
    res.pool4 <- metagen(RR4, SE4, data=filter(res.samp, !is.na(RR4)))
    res.pool5 <- metagen(RR5, SE5, data=filter(res.samp, !is.na(RR5)))
    
    res <- data.frame(RR1 = res.pool1$TE.fixed,
                      SE1 = res.pool1$seTE.fixed,
                      RR2 = res.pool2$TE.fixed,
                      SE2 = res.pool2$seTE.fixed,
                      RR3 = res.pool3$TE.fixed,
                      SE3 = res.pool3$seTE.fixed,
                      RR4 = res.pool4$TE.fixed,
                      SE4 = res.pool4$seTE.fixed,
                      RR5 = res.pool5$TE.fixed,
                      SE5 = res.pool5$seTE.fixed)
  }
  
} else if (analysis=="divide20") {
  
  k <- 20
  dat$split <- sample(rep(1:k, ceiling(nrow(dat)/k))[1:nrow(dat)])
  
  res.samp <- lapply(1:k, function(x){samp(x)})
  res.samp <- bind_rows(res.samp)
  
  if (boot_end==0) {
    res <- data.frame(RR1 = mean(res.samp$RR1, na.rm=T),
                      SE1 = NA,
                      n_invalid1 = sum(is.na(res.samp$RR1)),
                      RR2 = mean(res.samp$RR2, na.rm=T),
                      SE2 = NA,
                      n_invalid2 = sum(is.na(res.samp$RR2)),
                      RR3 = mean(res.samp$RR3, na.rm=T),
                      SE3 = NA,
                      n_invalid3 = sum(is.na(res.samp$RR3)),
                      RR4 = mean(res.samp$RR4, na.rm=T),
                      SE4 = NA,
                      n_invalid4 = sum(is.na(res.samp$RR4)),
                      RR5 = mean(res.samp$RR5, na.rm=T),
                      SE5 = NA,
                      n_invalid5 = sum(is.na(res.samp$RR5)))
    
  } else {
    # If a split has invalid RR (due to no events by t), drop it from metaregression
    res.pool1 <- metagen(RR1, SE1, data=filter(res.samp, !is.na(RR1)))
    res.pool2 <- metagen(RR2, SE2, data=filter(res.samp, !is.na(RR2)))
    res.pool3 <- metagen(RR3, SE3, data=filter(res.samp, !is.na(RR3)))
    res.pool4 <- metagen(RR4, SE4, data=filter(res.samp, !is.na(RR4)))
    res.pool5 <- metagen(RR5, SE5, data=filter(res.samp, !is.na(RR5)))
    
    res <- data.frame(RR1 = res.pool1$TE.fixed,
                      SE1 = res.pool1$seTE.fixed,
                      RR2 = res.pool2$TE.fixed,
                      SE2 = res.pool2$seTE.fixed,
                      RR3 = res.pool3$TE.fixed,
                      SE3 = res.pool3$seTE.fixed,
                      RR4 = res.pool4$TE.fixed,
                      SE4 = res.pool4$seTE.fixed,
                      RR5 = res.pool5$TE.fixed,
                      SE5 = res.pool5$seTE.fixed)
  }
  
} else if (analysis=="divide20p") {
  
  k <- 20
  dat$split <- sample(rep(1:k, ceiling(nrow(dat)/k))[1:nrow(dat)])
  
  res.samp <- mclapply(1:k, function(x){samp(x)}, mc.cores=cores, mc.set.seed=F)
  res.samp <- bind_rows(res.samp)
  
  if (boot_end==0) {
    res <- data.frame(RR1 = mean(res.samp$RR1, na.rm=T),
                      SE1 = NA,
                      n_invalid1 = sum(is.na(res.samp$RR1)),
                      RR2 = mean(res.samp$RR2, na.rm=T),
                      SE2 = NA,
                      n_invalid2 = sum(is.na(res.samp$RR2)),
                      RR3 = mean(res.samp$RR3, na.rm=T),
                      SE3 = NA,
                      n_invalid3 = sum(is.na(res.samp$RR3)),
                      RR4 = mean(res.samp$RR4, na.rm=T),
                      SE4 = NA,
                      n_invalid4 = sum(is.na(res.samp$RR4)),
                      RR5 = mean(res.samp$RR5, na.rm=T),
                      SE5 = NA,
                      n_invalid5 = sum(is.na(res.samp$RR5)))
    
  } else {
    # If a split has invalid RR (due to no events by t), drop it from metaregression
    res.pool1 <- metagen(RR1, SE1, data=filter(res.samp, !is.na(RR1)))
    res.pool2 <- metagen(RR2, SE2, data=filter(res.samp, !is.na(RR2)))
    res.pool3 <- metagen(RR3, SE3, data=filter(res.samp, !is.na(RR3)))
    res.pool4 <- metagen(RR4, SE4, data=filter(res.samp, !is.na(RR4)))
    res.pool5 <- metagen(RR5, SE5, data=filter(res.samp, !is.na(RR5)))
    
    res <- data.frame(RR1 = res.pool1$TE.fixed,
                      SE1 = res.pool1$seTE.fixed,
                      RR2 = res.pool2$TE.fixed,
                      SE2 = res.pool2$seTE.fixed,
                      RR3 = res.pool3$TE.fixed,
                      SE3 = res.pool3$seTE.fixed,
                      RR4 = res.pool4$TE.fixed,
                      SE4 = res.pool4$seTE.fixed,
                      RR5 = res.pool5$TE.fixed,
                      SE5 = res.pool5$seTE.fixed)
  }
  
} else if (analysis=="divide50") {
  
  k <- 50
  dat$split <- sample(rep(1:k, ceiling(nrow(dat)/k))[1:nrow(dat)])
  
  res.samp <- lapply(1:k, function(x){samp(x)})
  res.samp <- bind_rows(res.samp)
  
  if (boot_end==0) {
    res <- data.frame(RR1 = mean(res.samp$RR1, na.rm=T),
                      SE1 = NA,
                      n_invalid1 = sum(is.na(res.samp$RR1)),
                      RR2 = mean(res.samp$RR2, na.rm=T),
                      SE2 = NA,
                      n_invalid2 = sum(is.na(res.samp$RR2)),
                      RR3 = mean(res.samp$RR3, na.rm=T),
                      SE3 = NA,
                      n_invalid3 = sum(is.na(res.samp$RR3)),
                      RR4 = mean(res.samp$RR4, na.rm=T),
                      SE4 = NA,
                      n_invalid4 = sum(is.na(res.samp$RR4)),
                      RR5 = mean(res.samp$RR5, na.rm=T),
                      SE5 = NA,
                      n_invalid5 = sum(is.na(res.samp$RR5)))
    
  } else {
    # If a split has invalid RR (due to no events by t), drop it from metaregression
    res.pool1 <- metagen(RR1, SE1, data=filter(res.samp, !is.na(RR1)))
    res.pool2 <- metagen(RR2, SE2, data=filter(res.samp, !is.na(RR2)))
    res.pool3 <- metagen(RR3, SE3, data=filter(res.samp, !is.na(RR3)))
    res.pool4 <- metagen(RR4, SE4, data=filter(res.samp, !is.na(RR4)))
    res.pool5 <- metagen(RR5, SE5, data=filter(res.samp, !is.na(RR5)))
    
    res <- data.frame(RR1 = res.pool1$TE.fixed,
                      SE1 = res.pool1$seTE.fixed,
                      RR2 = res.pool2$TE.fixed,
                      SE2 = res.pool2$seTE.fixed,
                      RR3 = res.pool3$TE.fixed,
                      SE3 = res.pool3$seTE.fixed,
                      RR4 = res.pool4$TE.fixed,
                      SE4 = res.pool4$seTE.fixed,
                      RR5 = res.pool5$TE.fixed,
                      SE5 = res.pool5$seTE.fixed)
  }
  
} else if (analysis=="divide50p") {
  
  k <- 50
  dat$split <- sample(rep(1:k, ceiling(nrow(dat)/k))[1:nrow(dat)])
  
  res.samp <- mclapply(1:k, function(x){samp(x)}, mc.cores=cores, mc.set.seed=F)
  res.samp <- bind_rows(res.samp)
  
  if (boot_end==0) {
    res <- data.frame(RR1 = mean(res.samp$RR1, na.rm=T),
                      SE1 = NA,
                      n_invalid1 = sum(is.na(res.samp$RR1)),
                      RR2 = mean(res.samp$RR2, na.rm=T),
                      SE2 = NA,
                      n_invalid2 = sum(is.na(res.samp$RR2)),
                      RR3 = mean(res.samp$RR3, na.rm=T),
                      SE3 = NA,
                      n_invalid3 = sum(is.na(res.samp$RR3)),
                      RR4 = mean(res.samp$RR4, na.rm=T),
                      SE4 = NA,
                      n_invalid4 = sum(is.na(res.samp$RR4)),
                      RR5 = mean(res.samp$RR5, na.rm=T),
                      SE5 = NA,
                      n_invalid5 = sum(is.na(res.samp$RR5)))
    
  } else {
    # If a split has invalid RR (due to no events by t), drop it from metaregression
    res.pool1 <- metagen(RR1, SE1, data=filter(res.samp, !is.na(RR1)))
    res.pool2 <- metagen(RR2, SE2, data=filter(res.samp, !is.na(RR2)))
    res.pool3 <- metagen(RR3, SE3, data=filter(res.samp, !is.na(RR3)))
    res.pool4 <- metagen(RR4, SE4, data=filter(res.samp, !is.na(RR4)))
    res.pool5 <- metagen(RR5, SE5, data=filter(res.samp, !is.na(RR5)))
    
    res <- data.frame(RR1 = res.pool1$TE.fixed,
                      SE1 = res.pool1$seTE.fixed,
                      RR2 = res.pool2$TE.fixed,
                      SE2 = res.pool2$seTE.fixed,
                      RR3 = res.pool3$TE.fixed,
                      SE3 = res.pool3$seTE.fixed,
                      RR4 = res.pool4$TE.fixed,
                      SE4 = res.pool4$seTE.fixed,
                      RR5 = res.pool5$TE.fixed,
                      SE5 = res.pool5$seTE.fixed)
  }
} 


res$analysis <- analysis
end_time <- Sys.time()
res$comp_time <- difftime(end_time, start_time, units="hours")


# Output results ----------------------------------------------------------

if (boot_end==0) {
  write_csv(res, paste0("./results/res_cif_", analysis, "_pointest2.csv"))
} else {
  write_csv(res, paste0("./results/res_cif_", analysis, ".csv"))
}
