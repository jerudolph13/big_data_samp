The programs are used to run the following analyses:


RR Estimation
- samp_2_cif.R : Estimate RR (and bootstrap 95% CIs) in full data, sub-cohort, and case-cohort analyses
- samp_2_cif-divide.R : Estimate RR (and bootstrap 95% CIs) in divide and recombine analyses

IRR Estimation
- samp_2_poisson.R : Estimate IRR with robust SE (note that the SE is invalid for the case-cohort analysis)
- samp_2_poisson.R : Estimate IRR (and bootstrap 95% CIs) in full data, sub-cohort, and case-cohort analyses
- samp_2_poisson-divide.R : Estimate IRR (and bootstrap 95% CIs) in divide and recombine analyses

HR Estimation
- samp_2_cox.R : Estimate HR with robust SE across all sampling approaches
