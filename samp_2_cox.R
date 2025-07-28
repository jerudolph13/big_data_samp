
################################################################################
#
# Project: Big Data Sampling
#
# Purpose: Run Cox analysis
#
# Last Update: 30 Sep 2024
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
cores <- 50

# Read in data ------------------------------------------------------------

dat <- read_csv("/cms01/data/dua/57285/users/c-jrudolp9-57285/samp_data.csv") %>%
  mutate(across(c(hiv, sex, race, state, enroll_period, n_cond), ~ as.factor(.x)))


# Sample and analyze ------------------------------------------------------

start_time <- Sys.time()

# Function for modeling
models <- function(data) {
  
  # Weighted Cox model
  if (analysis %in% c("cc10", "cc25")) {
    
    unique <- filter(data, !duplicated(bene_id))	  
    
    mod_hiv <- glm(hiv ~ bs(age, df=3) + race + sex + state + enroll_period,
                   family=binomial(link="logit"), weights=wt,  data=unique)$fitted.values
    den_hiv <- ifelse(unique$hiv==1, mod_hiv, 1 - mod_hiv)
    
    mod_hiv <- glm(hiv ~ 1,
                   family=binomial(link="logit"), weights=wt, data=unique)$fitted.values
    num_hiv <- ifelse(unique$hiv==1, mod_hiv, 1 - mod_hiv)
    
    unique$wt2 <- num_hiv/den_hiv
    
    data <- left_join(data, select(unique, bene_id, wt2), by="bene_id") %>%
      mutate(wt3 = wt*wt2)
    
    mod.cox <- coxph(Surv(start, end, as.numeric(first_event=="lung")) ~  hiv,
                     weights=wt3, id=bene_id, data=data)
    est.cox <- tidy(mod.cox)
    print(est.cox)
    
  } else {
    
    mod_hiv <- glm(hiv ~ bs(age, df=3) + race + sex + state + enroll_period,
                   family=binomial(link="logit"), weights=wt, data=data)$fitted.values
    den_hiv <- ifelse(data$hiv==1, mod_hiv, 1 - mod_hiv)
    
    mod_hiv <- glm(hiv ~ 1,
                   family=binomial(link="logit"), weights=wt, data=data)$fitted.values
    num_hiv <- ifelse(data$hiv==1, mod_hiv, 1 - mod_hiv)
    
    data$wt2 <- (num_hiv/den_hiv)*data$wt
    print(summary(data$wt))
    
    mod.cox <- coxph(Surv(person_yrs, as.numeric(first_event=="lung")) ~  hiv,
                     weights=wt2, id=bene_id, robust=T, data=data)
    est.cox <- tidy(mod.cox)
    print(est.cox)
    
  }
  
  est <- data.frame(HR = est.cox$estimate,
                    SE = est.cox$robust.se)
  
  return(est)
  
}

# Function for sampling
samp <- function(s) {
  
  samp <- filter(dat, s==split) %>%
    mutate(wt = 1)
  
  res.k <- models(samp)
  
  return(res.k)
}

# Run analysis
if (analysis=="full") {
  samp <- dat %>%
    mutate(wt = 1)
  
  res <- models(samp)
  
} else if (analysis=="divide10") {
  
  k <- 10
  dat$split <- sample(rep(1:k, ceiling(nrow(dat)/k))[1:nrow(dat)])
  
  res.samp <- lapply(1:k, function(x){samp(x)})
  res.samp <- bind_rows(res.samp)
  
  res.pool <- metagen(HR, SE, data=res.samp)
  res <- data.frame(HR = res.pool$TE.fixed,
                    SE = res.pool$seTE.fixed,
                    HR2 = res.pool$TE.random,
                    SE2 = res.pool$seTE.random)
  
} else if (analysis=="divide10p") {
  
  k <- 10
  dat$split <- sample(rep(1:k, ceiling(nrow(dat)/k))[1:nrow(dat)])
  
  res.samp <- mclapply(1:k, function(x){samp(x)}, mc.cores=cores)
  res.samp <- bind_rows(res.samp)
  
  res.pool <- metagen(HR, SE, data=res.samp)
  res <- data.frame(HR = res.pool$TE.fixed,
                    SE = res.pool$seTE.fixed,
                    HR2 = res.pool$TE.random,
                    SE2 = res.pool$seTE.random)
  
} else if (analysis=="divide20") {
  
  k <- 20
  dat$split <- sample(rep(1:k, ceiling(nrow(dat)/k))[1:nrow(dat)])
  
  res.samp <- lapply(1:k, function(x){samp(x)})
  res.samp <- bind_rows(res.samp)
  
  res.pool <- metagen(HR, SE, data=res.samp)
  res <- data.frame(HR = res.pool$TE.fixed,
                    SE = res.pool$seTE.fixed,
                    HR2 = res.pool$TE.random,
                    SE2 = res.pool$seTE.random)
  
} else if (analysis=="divide20p") {
  
  k <- 20
  dat$split <- sample(rep(1:k, ceiling(nrow(dat)/k))[1:nrow(dat)])
  
  res.samp <- mclapply(1:k, function(x){samp(x)}, mc.cores=cores)
  res.samp <- bind_rows(res.samp)
  
  res.pool <- metagen(HR, SE, data=res.samp)
  res <- data.frame(HR = res.pool$TE.fixed,
                    SE = res.pool$seTE.fixed,
                    HR2 = res.pool$TE.random,
                    SE2 = res.pool$seTE.random)
  
} else if (analysis=="divide50") {
  
  k <- 50
  dat$split <- sample(rep(1:k, ceiling(nrow(dat)/k))[1:nrow(dat)])
  
  res.samp <- lapply(1:k, function(x){samp(x)})
  res.samp <- bind_rows(res.samp)
  
  res.pool <- metagen(HR, SE, data=res.samp)
  res <- data.frame(HR = res.pool$TE.fixed,
                    SE = res.pool$seTE.fixed,
                    HR2 = res.pool$TE.random,
                    SE2 = res.pool$seTE.random)
  
} else if (analysis=="divide50p") {
  
  k <- 50
  dat$split <- sample(rep(1:k, ceiling(nrow(dat)/k))[1:nrow(dat)])
  
  res.samp <- mclapply(1:k, function(x){samp(x)}, mc.cores=cores)
  res.samp <- bind_rows(res.samp)
  
  res.pool <- metagen(HR, SE, data=res.samp)
  res <- data.frame(HR = res.pool$TE.fixed,
                    SE = res.pool$seTE.fixed,
                    HR2 = res.pool$TE.random,
                    SE2 = res.pool$seTE.random)
  
} else if (analysis=="sc25") {
  
  p1 <- 1
  p0 <- 0.25
  
  dat1 <- filter(dat, hiv==1) %>%
    mutate(wt = 1/p1)
  dat0 <- filter(dat, hiv==0)[rbinom(nrow(filter(dat, hiv==0)), 1, p0)==1, ] %>%
    mutate(wt = 1/p0)
  
  samp <- bind_rows(dat1, dat0)
  
  res <- models(samp)
  
} else if (analysis=="sc10") {
  
  p1 <- 1
  p0 <- 0.10
  
  dat1 <- filter(dat, hiv==1) %>%
    mutate(wt = 1/p1)
  dat0 <- filter(dat, hiv==0)[rbinom(nrow(filter(dat, hiv==0)), 1, p0)==1, ] %>%
    mutate(wt = 1/p0)
  
  samp <- bind_rows(dat1, dat0)
  
  res <- models(samp)
  
} else if (analysis=="cc25") {
  
  p1 <- 1
  p0 <- 0.25
  
  dat1 <- filter(dat, hiv==1) %>%
    mutate(wt = 1/p1,
           end = ifelse(first_event=="lung", person_yrs - 0.001, person_yrs),
           start = 0)
  dat0 <- filter(dat, hiv==0)[rbinom(nrow(filter(dat, hiv==0)), 1, p0)==1, ] %>%
    mutate(wt = 1/p0,
           end = ifelse(first_event=="lung", person_yrs - 0.001, person_yrs),
           start = 0)
  cases <- filter(dat, first_event=="lung") %>%
    mutate(wt = 1,
           end = person_yrs,
           start = person_yrs - 0.001)
  
  samp <- bind_rows(dat1, dat0, cases)
  
  res <- models(samp)
  
} else if (analysis=="cc10") {
  
  p1 <- 1
  p0 <- 0.10
  
  dat1 <- filter(dat, hiv==1) %>%
    mutate(wt = 1/p1,
           end = ifelse(first_event=="lung", person_yrs - 0.001, person_yrs),
           start = 0)
  dat0 <- filter(dat, hiv==0)[rbinom(nrow(filter(dat, hiv==0)), 1, p0)==1, ] %>%
    mutate(wt = 1/p0,
           end = ifelse(first_event=="lung", person_yrs - 0.001, person_yrs),
           start = 0)
  cases <- filter(dat, first_event=="lung") %>%
    mutate(wt = 1,
           end = person_yrs,
           start = person_yrs - 0.001)
  
  samp <- bind_rows(dat1, dat0, cases)
  
  res <- models(samp)
  
}

res$analysis <- analysis
end_time <- Sys.time()
res$comp_time <- difftime(end_time, start_time, units="hours")


# Output results ----------------------------------------------------------

write_csv(res, paste0("./results/res_cox_", analysis, ".csv"))

