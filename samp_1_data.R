
################################################################################
#
# Project: Big Data Sampling
#
# Purpose: Manage raw data
#
# Last Update: 29 Jul 2024
#
################################################################################


lib <- "~/R/4.3"
packages <- c("tidyverse", "lubridate", "readxl")
for (package in packages){
  if (!(package %in% installed.packages(lib=lib))) {install.packages(package, lib=lib)}
  suppressPackageStartupMessages(library(package, character.only=T, lib.loc=lib, quietly=T))
}


# Read in data ------------------------------------------------------------

dat <- read_csv(file="/cms01/data/dua/57285/medican/data_requests/jackie/20240710_jacike_incidences_v7/incid_v7.csv", 
                col_types=paste0("cfffffccDDDDfff", paste(rep("D", 76), collapse=""))) %>%
  rename(start_date = first_elig_period_start_date_65,
         end_date = first_elig_period_end_date_65,
         state = first_elig_state,
         zip = first_elig_zip) %>%
  filter(first_elig_period_rstrct_bnft_65==1) %>% 
  select(bene_id, race, sex, state, zip, dob, dod, start_date, end_date, 
         hiv_cms_typ2_dt, any_cx_typ1_dt, any_cx_typ2_dt, lung_cx_typ2_dt,
         myocard_infarct_typ1_dt, cong_heart_fail_typ1_dt,
         periph_vas_dis_typ1_dt, cereb_vas_dis_typ1_dt, dementia_typ1_dt,
         chron_pulm_dis_typ1_dt, rheumatic_dis_typ1_dt, peptic_ulcer_typ1_dt,
         mild_liv_dis_typ1_dt, mod_sevr_liver_dis_typ1_dt, diabt_wo_chron_complict_typ1_dt,
         diabt_chron_complict_da_typ1_dt, hemiplegia_paraplegia_d_typ1_dt,
         renal_dis_typ1_dt, metastatic_typ1_dt, any_malignancy_typ1_dt)

dat.ruca <- read_excel("../../shared/RUCA2010zipcode.xlsx", sheet = "Data") %>%
  rename(zip = ZIP_CODE) %>%
  select(zip, RUCA1)


# Manage data -------------------------------------------------------------

dat2 <- dat %>% 
  # Define age, HIV status, and comorbidities at baseline
  mutate(baseline = start_date + months(6),
         age = (dob %--% baseline)/years(1),
         hiv = ifelse(is.na(hiv_cms_typ2_dt), 0, 
                      as.numeric(hiv_cms_typ2_dt <= baseline)),
        
         ami = ifelse(is.na(myocard_infarct_typ1_dt), 0, 
                      as.numeric(myocard_infarct_typ1_dt <= baseline)),
         chf = ifelse(is.na(cong_heart_fail_typ1_dt), 0, 
                      as.numeric(cong_heart_fail_typ1_dt <= baseline)),
         pvd = ifelse(is.na(periph_vas_dis_typ1_dt), 0, 
                      as.numeric(periph_vas_dis_typ1_dt <= baseline)),
         cevd = ifelse(is.na(cereb_vas_dis_typ1_dt), 0, 
                       as.numeric(cereb_vas_dis_typ1_dt <= baseline)),
         dementia = ifelse(is.na(dementia_typ1_dt), 0, 
                           as.numeric(dementia_typ1_dt <= baseline)),
         copd = ifelse(is.na(chron_pulm_dis_typ1_dt), 0, 
                       as.numeric(chron_pulm_dis_typ1_dt <= baseline)),
         rheumd = ifelse(is.na(rheumatic_dis_typ1_dt), 0, 
                         as.numeric(rheumatic_dis_typ1_dt <= baseline)),
         pud = ifelse(is.na(peptic_ulcer_typ1_dt), 0, 
                      as.numeric(peptic_ulcer_typ1_dt <= baseline)),
         mld = ifelse(is.na(mild_liv_dis_typ1_dt), 0, 
                      as.numeric(mild_liv_dis_typ1_dt <= baseline)),
         diab = ifelse(is.na(diabt_wo_chron_complict_typ1_dt), 0, 
                       as.numeric(diabt_wo_chron_complict_typ1_dt <= baseline)),
         diabwc = ifelse(is.na(diabt_chron_complict_da_typ1_dt), 0, 
                         as.numeric(diabt_chron_complict_da_typ1_dt <= baseline)),
         hp = ifelse(is.na(hemiplegia_paraplegia_d_typ1_dt), 0, 
                     as.numeric(hemiplegia_paraplegia_d_typ1_dt <= baseline)),
         rend = ifelse(is.na(renal_dis_typ1_dt), 0, 
                       as.numeric(renal_dis_typ1_dt <= baseline)),
         msld = ifelse(is.na(mod_sevr_liver_dis_typ1_dt), 0, 
                       as.numeric(mod_sevr_liver_dis_typ1_dt <= baseline)),
         n_cond = ami + chf + pvd + cevd + dementia + copd + rheumd + pud + mld +
                  diab+ diabwc + hp + rend + msld) %>% 
  # Remove ineligible beneficiaries:
    # Died, had cancer, or dropped out before baseline
  filter((any_cx_typ1_dt > baseline) | is.na(any_cx_typ1_dt), 
         (dod > baseline) | is.na(dod),
         end_date > baseline)

# Determine first event
first <- as_date(apply(dat2[, c("end_date", "lung_cx_typ2_dt", "any_cx_typ2_dt", "dod")],
                       1, min, na.rm=T))
dat2 <- bind_cols(dat2, "event_date"=first)
event_type <- rep(NA, nrow(dat2))
event_type[dat2$event_date == dat2$end_date] <- "censor"
event_type[dat2$event_date == dat2$any_cx_typ2_dt] <- "censor" # Censor at other cancers
event_type[dat2$event_date == dat2$dod] <- "dead"
event_type[dat2$event_date == dat2$lung_cx_typ2_dt] <- "lung"
dat2 <- dat2 %>% 
  bind_cols(., "first_event"=event_type) %>% 
  mutate(person_yrs = (baseline %--% event_date)/years(1))

# Set-up for analysis
dat3 <- dat2 %>% 
  filter(sex!="") %>% 
  mutate(yr_base = year(baseline),
         enroll_period = factor(case_when(
           yr_base %in% 2001:2005 ~ 1,
           yr_base %in% 2006:2010 ~ 2,
           T ~ 3)),
         n_cond = factor(ifelse(n_cond==0, 0,
                                ifelse(n_cond==1, 1, 2)))) %>% 
  select(bene_id, hiv, first_event, person_yrs, age, sex, race, state, zip, enroll_period, n_cond)

# Merge zip code data
dat3$zip <- ifelse(nchar(dat3$zip) == 4 | nchar(dat3$zip) == 8, paste0("0", dat3$zip), dat3$zip)
dat3$zip <- ifelse(nchar(dat3$zip) == 9, substr(dat3$zip, 1, 5), dat3$zip)
dat4 <- left_join(dat3, dat.ruca, by="zip")


# Output data -------------------------------------------------------------

#write_csv(dat3, "/cms01/data/dua/57285/users/c-jrudolp9-57285/samp_data.csv")

# Testing sample
samp <- dat4[rbinom(nrow(dat4), 1, 0.10)==1, ]
write_csv(samp, "/cms01/data/dua/57285/users/c-jrudolp9-57285/samp_data.csv")

