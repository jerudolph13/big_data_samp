
################################################################################
#
# Project: Big Data Sampling
#
# Purpose: Run non-clustered analysis
#
# Last Update: 05 Aug 2024
#
################################################################################

lib <- "~/R/4.3"
packages <- c("tidyverse", "splines", "sandwich", "broom", "meta", "parallel")
for (package in packages){
  if (!(package %in% installed.packages(lib=lib))) {install.packages(package, lib=lib)}
  suppressPackageStartupMessages(library(package, character.only=T, lib.loc=lib, quietly=T))
}

set.seed(123)

args <- commandArgs(trailingOnly=TRUE)
analysis <- as.character(args[1])
cores <- 10

# Read in data ------------------------------------------------------------

dat <- read_csv("/cms01/data/dua/57285/users/c-jrudolp9-57285/samp_data.csv") %>%
  mutate(across(c(hiv, sex, race, state, enroll_period, n_cond), ~ as.factor(.x)))


# Sample and analyze ------------------------------------------------------

start_time <- Sys.time()

# Function for modeling
models <- function(data) {
  
  if (analysis %in% c("cc10", "cc25")) {
    
    unique <- filter(data, !duplicated(bene_id))
    
    mod_hiv <- glm(hiv ~ bs(age, df=3) + race + sex + state + enroll_period,
                   family=binomial(link="logit"), weights=wt, data=unique)$fitted.values
    den_hiv <- ifelse(unique$hiv==1, mod_hiv, 1 - mod_hiv)
    
    mod_hiv <- glm(hiv ~ 1,
                   family=binomial(link="logit"), weights=wt, data=unique)$fitted.values
    num_hiv <- ifelse(unique$hiv==1, mod_hiv, 1 - mod_hiv)
    
    unique$wt2 <- num_hiv/den_hiv
    data <- left_join(data, select(unique, bene_id, wt2), by="bene_id") %>%
      mutate(wt3 = wt*wt2)
    
    mod <- glm(I(first_event=="lung") ~ hiv + offset(log(person_yrs)), weights=wt3, family=poisson(link="log"), data=data)
    est <- tidy(mod)
    print(est)
    
  } else {
    
    mod_hiv <- glm(hiv ~ bs(age, df=3) + race + sex + state + enroll_period,
                   family=binomial(link="logit"), weights=wt, data=data)$fitted.values
    den_hiv <- ifelse(data$hiv==1, mod_hiv, 1 - mod_hiv)
    
    mod_hiv <- glm(hiv ~ 1,
                   family=binomial(link="logit"), weights=wt, data=data)$fitted.values
    num_hiv <- ifelse(data$hiv==1, mod_hiv, 1 - mod_hiv)
    
    data$wt2 <- (num_hiv/den_hiv)*data$wt
    
    mod <- glm(I(first_event=="lung") ~ hiv + offset(log(person_yrs)), weights=wt2, family=poisson(link="log"), data=data)
    est <- tidy(mod)
    print(est)
    
  }
  
  est2 <- data.frame(IRR = est$estimate[2],
                     SE = sqrt(vcovHC(mod, type="HC1")[2,2]))
  print(est2)
  
  return(est2)
  
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
  
  res.pool <- metagen(IRR, SE, data=res.samp)
  res <- data.frame(IRR = res.pool$TE.fixed,
                    SE = res.pool$seTE.fixed,
                    IRR2 = res.pool$TE.random,
                    SE2 = res.pool$seTE.random)
  
} else if (analysis=="divide10p") {
  
  k <- 10
  dat$split <- sample(rep(1:k, ceiling(nrow(dat)/k))[1:nrow(dat)])
  
  res.samp <- mclapply(1:k, function(x){samp(x)}, mc.cores=cores)
  res.samp <- bind_rows(res.samp)
  
  res.pool <- metagen(IRR, SE, data=res.samp)
  res <- data.frame(IRR = res.pool$TE.fixed,
                    SE = res.pool$seTE.fixed,
                    IRR2 = res.pool$TE.random,
                    SE2 = res.pool$seTE.random)
  
} else if (analysis=="divide20") {
  
  k <- 20
  dat$split <- sample(rep(1:k, ceiling(nrow(dat)/k))[1:nrow(dat)])
  
  res.samp <- lapply(1:k, function(x){samp(x)})
  res.samp <- bind_rows(res.samp)
  
  res.pool <- metagen(IRR, SE, data=res.samp)
  res <- data.frame(IRR = res.pool$TE.fixed,
                    SE = res.pool$seTE.fixed,
                    IRR2 = res.pool$TE.random,
                    SE2 = res.pool$seTE.random)
  
} else if (analysis=="divide20p") {
  
  k <- 20
  dat$split <- sample(rep(1:k, ceiling(nrow(dat)/k))[1:nrow(dat)])
  
  res.samp <- mclapply(1:k, function(x){samp(x)}, mc.cores=cores)
  res.samp <- bind_rows(res.samp)
  
  res.pool <- metagen(IRR, SE, data=res.samp)
  res <- data.frame(IRR = res.pool$TE.fixed,
                    SE = res.pool$seTE.fixed,
                    IRR2 = res.pool$TE.random,
                    SE2 = res.pool$seTE.random)
  
} else if (analysis=="divide50") {
  
  k <- 50
  dat$split <- sample(rep(1:k, ceiling(nrow(dat)/k))[1:nrow(dat)])
  
  res.samp <- lapply(1:k, function(x){samp(x)})
  res.samp <- bind_rows(res.samp)
  
  res.pool <- metagen(IRR, SE, data=res.samp)
  res <- data.frame(IRR = res.pool$TE.fixed,
                    SE = res.pool$seTE.fixed,
                    IRR2 = res.pool$TE.random,
                    SE2 = res.pool$seTE.random)
  
} else if (analysis=="divide50p") {
  
  k <- 50
  dat$split <- sample(rep(1:k, ceiling(nrow(dat)/k))[1:nrow(dat)])
  
  res.samp <- mclapply(1:k, function(x){samp(x)}, mc.cores=cores)
  res.samp <- bind_rows(res.samp)
  
  res.pool <- metagen(IRR, SE, data=res.samp)
  res <- data.frame(IRR = res.pool$TE.fixed,
                    SE = res.pool$seTE.fixed,
                    IRR2 = res.pool$TE.random,
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
           person_yrs = ifelse(first_event=="lung", person_yrs - 0.001, person_yrs))
  dat0 <- filter(dat, hiv==0)[rbinom(nrow(filter(dat, hiv==0)), 1, p0)==1, ] %>%
    mutate(wt = 1/p0,
           person_yrs = ifelse(first_event=="lung", person_yrs - 0.001, person_yrs))
  cases <- filter(dat, first_event=="lung") %>%
    mutate(wt = 1,
           person_yrs = 0.001)
  
  samp <- bind_rows(dat1, dat0, cases)
  
  res <- models(samp)
  
} else if (analysis=="cc10") {
  
  p1 <- 1
  p0 <- 0.10
  
  dat1 <- filter(dat, hiv==1) %>%
    mutate(wt = 1/p1,
           person_yrs = ifelse(first_event=="lung", person_yrs - 0.001, person_yrs))
  dat0 <- filter(dat, hiv==0)[rbinom(nrow(filter(dat, hiv==0)), 1, p0)==1, ] %>%
    mutate(wt = 1/p0,
           person_yrs = ifelse(first_event=="lung", person_yrs - 0.001, person_yrs))
  cases <- filter(dat, first_event=="lung") %>%
    mutate(wt = 1,
           person_yrs = 0.001)
  
  samp <- bind_rows(dat1, dat0, cases)
  
  res <- models(samp)
  
}

res$analysis <- analysis
end_time <- Sys.time()
res$comp_time <- difftime(end_time, start_time, units="hours")


# Output results ----------------------------------------------------------

write_csv(res, paste0("./results/res_poisson_", analysis, ".csv"))

