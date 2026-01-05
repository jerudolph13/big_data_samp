This repository contains code for the paper "Sampling for computational efficiency when conducting analyses in big data," In Press at the _American Journal of Epidemiology_ We also include the bash scripts used to run the R programs.
 

CODE DESCRIPTION

Step 1. Set up data for analysis using samp_1_data.R

Step 2. Run the analyses

  RR Estimation
  - samp_2_cif.R : Estimate RR (and bootstrap 95% CIs) in full data, sub-cohort, and case-cohort analyses
  - samp_2_cif-divide.R : Estimate RR (and bootstrap 95% CIs) in divide-and-recombine analyses
  - samp_2_cif-divide-repeat.R : Repeat divide-and-recombine analyses
  
  IRR Estimation
  - samp_2_poisson.R : Estimate IRR with robust SE (note that the SE is invalid for the case-cohort analysis)
  - samp_2_poisson-boot.R : Estimate IRR (and bootstrap 95% CIs) in full data, sub-cohort, and case-cohort analyses
  - samp_2_poisson-divide.R : Estimate IRR (and bootstrap 95% CIs) in divide-and-recombine analyses
  - samp_2_poisson-divide-repeat.R : Repeat divide-and-recombine analyses

  HR Estimation
  - samp_2_cox.R : Estimate HR with robust SE across all sampling approaches
  - samp_2_cox-repeat.R : Repeat divide-and-recombine analyses

Step 3. Summarize results
  - samp_3_repeat-summ.R : Smmarize across repeated divide-and-recombine analyses
  - samp_3_tables.R : Build manuscript tables
  - samp_3_fig.R : Build manuscript figures
