
################################################################################
#
# Project: Big Data Sampling
#
# Purpose: Run rate analysis (divide and recombine)
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
boot_end <- 500
cores <- 10


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
                      IRR = NA)
    
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
    
    # Weighted Poisson model
    mod_hiv <- glm(hiv ~ bs(age, df=3) + race + sex + state + enroll_period,
                   family=binomial(link="logit"), weights=wt, data=boot)$fitted.values
    den_hiv <- ifelse(boot$hiv==1, mod_hiv, 1 - mod_hiv)
    
    mod_hiv <- glm(hiv ~ 1,
                   family=binomial(link="logit"), weights=wt, data=boot)$fitted.values
    num_hiv <- ifelse(boot$hiv==1, mod_hiv, 1 - mod_hiv)
    
    boot$wt2 <- (num_hiv/den_hiv)*boot$wt
    
    mod <- glm(I(first_event=="lung") ~ hiv + offset(log(person_yrs)), 
               weights=wt2, family=poisson(link="log"), data=boot)
    res <- tidy(mod)
    
    est$rep <- r
    est$IRR <- res$estimate[2]
    
    return(est)
  }
  
  #all.res <- lapply(boot_start:boot_end, function(x) {boot_rep(x)})
  all.res <- mclapply(boot_start:boot_end, function(x) {boot_rep(x)}, mc.cores=cores, mc.set.seed=F)
  all.res <- bind_rows(all.res)
  
  # Summarize across bootsraps
  res <- data.frame(IRR = all.res$IRR[all.res$rep==0],
                    SE = sd(all.res$IRR))
  
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
    res <- data.frame(IRR = mean(res.samp$IRR),
                      SE = NA)
    
  } else {   
    res.pool <- metagen(IRR, SE, data=res.samp)
    
    res <- data.frame(IRR = res.pool$TE.fixed,
                      SE = res.pool$seTE.fixed)
    
  }
  
} else if (analysis=="divide10p") {
  
  k <- 10
  dat$split <- sample(rep(1:k, ceiling(nrow(dat)/k))[1:nrow(dat)])
  
  res.samp <- mclapply(1:k, function(x){samp(x)}, mc.cores=cores, mc.set.seed=F)
  res.samp <- bind_rows(res.samp)
  
  if (boot_end==0) {
    res <- data.frame(IRR = mean(res.samp$IRR),
                      SE = NA)
    
  } else {
    res.pool <- metagen(IRR, SE, data=res.samp)
    
    res <- data.frame(IRR = res.pool$TE.fixed,
                      SE = res.pool$seTE.fixed)
    
  }
  
} else if (analysis=="divide20") {
  
  k <- 20
  dat$split <- sample(rep(1:k, ceiling(nrow(dat)/k))[1:nrow(dat)])
  
  res.samp <- lapply(1:k, function(x){samp(x)})
  res.samp <- bind_rows(res.samp)
  
  if (boot_end==0) {
    res <- data.frame(IRR = mean(res.samp$IRR), 
                      SE = NA)
    
  } else {
    res.pool <- metagen(IRR, SE, data=res.samp)
    
    res <- data.frame(IRR = res.pool$TE.fixed, 
                      SE = res.pool$seTE.fixed)
    
  }      
  
} else if (analysis=="divide20p") {
  
  k <- 20
  dat$split <- sample(rep(1:k, ceiling(nrow(dat)/k))[1:nrow(dat)])
  
  res.samp <- mclapply(1:k, function(x){samp(x)}, mc.cores=cores, mc.set.seed=F)
  res.samp <- bind_rows(res.samp)
  
  if (boot_end==0) {
    res <- data.frame(IRR = mean(res.samp$IRR),
                      SE = NA)
    
  } else {
    res.pool <- metagen(IRR, SE, data=res.samp)
    
    res <- data.frame(IRR = res.pool$TE.fixed,
                      SE = res.pool$seTE.fixed)
    
  }      
  
} else if (analysis=="divide50") {
  
  k <- 50
  dat$split <- sample(rep(1:k, ceiling(nrow(dat)/k))[1:nrow(dat)])
  
  res.samp <- lapply(1:k, function(x){samp(x)})
  res.samp <- bind_rows(res.samp)
  
  if (boot_end==0) {
    res <- data.frame(IRR = mean(res.samp$IRR),
                      SE = NA)
    
  } else {
    res.pool <- metagen(IRR, SE, data=res.samp)
    
    res <- data.frame(IRR = res.pool$TE.fixed,
                      SE = res.pool$seTE.fixed)
    
  }
  
} else if (analysis=="divide50p") {
  
  k <- 50
  dat$split <- sample(rep(1:k, ceiling(nrow(dat)/k))[1:nrow(dat)])
  
  res.samp <- mclapply(1:k, function(x){samp(x)}, mc.cores=cores, mc.set.seed=F)
  res.samp <- bind_rows(res.samp)
  
  if (boot_end==0) {
    res <- data.frame(IRR = mean(res.samp$IRR),
                      SE = NA)
    
  } else {
    res.samp <- res.samp %>%
      filter(!is.na(SE)) %>%
      filter(SE!=0)
    
    res.pool <- metagen(IRR, SE, data=res.samp)
    
    res <- data.frame(IRR = res.pool$TE.fixed,
                      SE = res.pool$seTE.fixed)
    
  }      
  
} 


res$analysis <- analysis
end_time <- Sys.time()
res$comp_time <- difftime(end_time, start_time, units="hours")


# Output results ----------------------------------------------------------

if (boot_end==0) {
  write_csv(res, paste0("./results/res_poisson_", analysis, "_pointest.csv"))
} else {
  write_csv(res, paste0("./results/res_poisson_", analysis, ".csv"))
}
