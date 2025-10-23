
################################################################################
#
# Project: Big Data Sampling
#
# Purpose: Repeat Poisson divide-and-recombine analysis
#
# Last Update: 05 Sep 2025
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
cores <- 10


# Read in data ------------------------------------------------------------

dat <- read_csv("/cms01/data/dua/57285/users/c-jrudolp9-57285/samp_data.csv") %>%
  mutate(across(c(hiv, sex, race, state, enroll_period, n_cond), ~ as.factor(.x)))


# Sample and analyze ------------------------------------------------------

rep_analysis <- function(r) {

  set.seed(123+r)

start_time <- Sys.time()

# Function for modeling
models <- function(data) {

  # Weighted Poisson model

  mod_hiv <- glm(hiv ~ bs(age, df=3) + race + sex + state + enroll_period,
                 family=binomial(link="logit"), weights=wt, data=data)$fitted.values
  den_hiv <- ifelse(data$hiv==1, mod_hiv, 1 - mod_hiv)

  mod_hiv <- glm(hiv ~ 1,
                 family=binomial(link="logit"), weights=wt, data=data)$fitted.values
  num_hiv <- ifelse(data$hiv==1, mod_hiv, 1 - mod_hiv)

  data$wt2 <- (num_hiv/den_hiv)*data$wt

  mod <- glm(I(first_event=="lung") ~ hiv + offset(log(person_yrs)),
             weights=wt2, family=poisson(link="log"), data=data)
  res <- tidy(mod)

  est <- data.frame(IRR = res$estimate[2])

  return(est)

}

# Function for sampling
samp <- function(s) {

  samp <- filter(dat, s==split) %>%
          mutate(wt = 1)

  res.k <- models(samp)

  return(res.k)
}


# Run analysis -----------------------------------------------------------

if (analysis=="divide10p") {

  k <- 10
  dat$split <- sample(rep(1:k, ceiling(nrow(dat)/k))[1:nrow(dat)])

  res.samp <- mclapply(1:k, function(x){samp(x)}, mc.cores=cores, mc.set.seed=F)
  res.samp <- bind_rows(res.samp)

  res <- data.frame(IRR = mean(res.samp$IRR))

} else if (analysis=="divide50p") {

  k <- 50
  dat$split <- sample(rep(1:k, ceiling(nrow(dat)/k))[1:nrow(dat)])

  res.samp <- mclapply(1:k, function(x){samp(x)}, mc.cores=cores)
  res.samp <- bind_rows(res.samp)

  res <- data.frame(IRR = mean(res.samp$IRR))

}

res$analysis <- analysis
end_time <- Sys.time()
res$comp_time <- difftime(end_time, start_time, units="hours")

return(res)

}

all.res <- lapply(0:499, function(x){rep_analysis(x)})
all.res <- bind_rows(all.res)


# Output results ----------------------------------------------------------

write_csv(all.res, paste0("./results/res_poisson_", analysis, "_repeat.csv"))
