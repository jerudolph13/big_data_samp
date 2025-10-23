
################################################################################
#
# Project: Big Data Sampling
#
# Purpose: Visualize results
#
# Last Update: 26 Feb 2025
#
################################################################################

library("tidyverse")
library("readxl")
library("extrafont")
library("patchwork")

# font_import()
# loadfonts(device="win")

#palette <- c(rep(c("#002D72", "#68ACE5"), 6))
palette <- c(rep("black", 12))

thm <- theme_classic() +
  theme(axis.text = element_text(color="black", family="Tahoma"),
        axis.title.x = element_text(color="black", family="Tahoma"),
        axis.title.y = element_blank())


# IRR ---------------------------------------------------------------------

irr <- read_excel("../tables/irr.xlsx") %>% 
  filter(!(Approach %in% c("Divide (10; P; 10)", "Divide (20; P; 10)", "Divide (20; P; 20)", 
                           "Divide (50; P; 10)", "Divide (50; P; 50)"))) %>% 
  mutate(Approach = factor(Approach, 
                           levels=c("Full sample",
                                    "Divide (10; NP)",
                                    "Divide (20; NP)",
                                    "Divide (50; NP)",
                                    "Sub-cohort (25%)",
                                    "Sub-cohort (10%)",
                                    "Case-cohort (25%)",
                                    "Case-cohort (10%)"),
                           labels=c("Full sample",
                                    "Divide (10)",
                                    "Divide (20)",
                                    "Divide (50)",
                                    "Sub-cohort (25%)",
                                    "Sub-cohort (10%)",
                                    "Case-cohort (25%)",
                                    "Case-cohort (10%)")),
         SE = as.numeric(SE))

p1 <- ggplot(data=irr, aes(y=Approach, color=Approach)) +
  labs(tag="D)") +
  geom_vline(aes(xintercept=`Log(IRR)`[Approach=="Full sample"]),
             linewidth=0.5, linetype="dashed") +
  geom_errorbar(aes(xmin = `Log(IRR)` - 1.96*SE, 
                    xmax = `Log(IRR)` + 1.96*SE), linewidth=0.5, width=0.25) +
  geom_point(aes(x=`Log(IRR)`), size=2) +
  scale_x_continuous(limits=c(0.13, 0.91)) +
  scale_y_discrete(limits=rev) +
  scale_color_manual(values=palette) +
  guides(color="none") +
  thm

jpeg("../figures/irr_est.jpeg", height=3.5, width=3, units="in", res=300)
print(p1)
dev.off()


# HR ----------------------------------------------------------------------

hr <- read_excel("../tables/hr.xlsx") %>% 
  filter(!(Approach %in% c("Divide (10; P; 10)", "Divide (20; P; 10)", "Divide (20; P; 20)", 
                           "Divide (50; P; 10)", "Divide (50; P; 50)"))) %>% 
  mutate(Approach = factor(Approach, 
                           levels=c("Full sample",
                                    "Divide (10; NP)",
                                    "Divide (20; NP)",
                                    "Divide (50; NP)",
                                    "Sub-cohort (25%)",
                                    "Sub-cohort (10%)",
                                    "Case-cohort (25%)",
                                    "Case-cohort (10%)"),
                           labels=c("Full sample",
                                    "Divide (10)",
                                    "Divide (20)",
                                    "Divide (50)",
                                    "Sub-cohort (25%)",
                                    "Sub-cohort (10%)",
                                    "Case-cohort (25%)",
                                    "Case-cohort (10%)")),
         Error = abs(Error))

p2 <- ggplot(data=hr, aes(y=Approach, color=Approach)) +
  labs(tag="C)") +
  geom_vline(aes(xintercept=`Log(HR)`[Approach=="Full sample"]),
             linewidth=0.5, linetype="dashed") +
  geom_errorbar(aes(xmin = `Log(HR)` - 1.96*SE, 
                    xmax = `Log(HR)` + 1.96*SE), linewidth=0.5, width=0.25) +
  geom_point(aes(x=`Log(HR)`), size=2) +
  scale_x_continuous(limits=c(0.13, 0.91)) +
  scale_y_discrete(limits=rev) +
  scale_color_manual(values=palette) +
  guides(color="none") +
  thm

jpeg("../figures/hr_est.jpeg", height=3.5, width=3, units="in", res=300)
print(p2)
dev.off()


# RR ----------------------------------------------------------------------

rr <- read_excel("../tables/rr.xlsx") %>% 
  filter(!(Approach %in% c("Divide (10; P; 10)", "Divide (20; P; 10)", "Divide (20; P; 20)", 
                           "Divide (50; P; 10)", "Divide (50; P; 50)"))) %>% 
  mutate(Approach = factor(Approach, 
                           levels=c("Full sample",
                                    "Divide (10; NP)",
                                    "Divide (20; NP)",
                                    "Divide (50; NP)",
                                    "Sub-cohort (25%)",
                                    "Sub-cohort (10%)",
                                    "Case-cohort (25%)",
                                    "Case-cohort (10%)"),
                           labels=c("Full sample",
                                    "Divide (10)",
                                    "Divide (20)",
                                    "Divide (50)",
                                    "Sub-cohort (25%)",
                                    "Sub-cohort (10%)",
                                    "Case-cohort (25%)",
                                    "Case-cohort (10%)")),
         Error = abs(Error))

# t=1
p3 <- ggplot(data=filter(rr, Time==1), aes(y=Approach, color=Approach)) +
  labs(tag="A)") +
  labs(x=expression(paste("log(R", R[1], ")", sep=""))) +
  geom_vline(aes(xintercept=`Log(RR)`[Approach=="Full sample"]),
             linewidth=0.5, linetype="dashed") +
  geom_errorbar(aes(xmin = `Log(RR)` - 1.96*SE, 
                    xmax = `Log(RR)` + 1.96*SE), linewidth=0.5, width=0.25) +
  geom_point(aes(x=`Log(RR)`), size=2) +
  scale_x_continuous(limits=c(0.00, 0.91)) +
  scale_y_discrete(limits=rev) +
  scale_color_manual(values=palette) +
  guides(color="none") +
  thm

jpeg("../figures/rr_est1.jpeg", height=3.5, width=3, units="in", res=300)
print(p3)
dev.off()

# t=3
jpeg("../figures/rr_est3.jpeg", height=3.5, width=3, units="in", res=300)
ggplot(data=filter(rr, Time==3), aes(y=Approach, color=Approach)) +
  labs(x=expression(paste("log(R", R[3], ")", sep=""))) +
  geom_vline(aes(xintercept=`Log(RR)`[Approach=="Full sample"]),
             linewidth=0.5, linetype="dashed") +
  geom_errorbar(aes(xmin = `Log(RR)` - 1.96*SE, 
                    xmax = `Log(RR)` + 1.96*SE), linewidth=0.5, width=0.25) +
  geom_point(aes(x=`Log(RR)`), size=2) +
  scale_x_continuous(limits=c(0.00, 0.91)) +
  scale_y_discrete(limits=rev) +
  scale_color_manual(values=palette) +
  guides(color="none") +
  thm
dev.off()

# t=5
p4 <- ggplot(data=filter(rr, Time==5), aes(y=Approach, color=Approach)) +
  labs(tag="B)",
       x=expression(paste("log(R", R[5], ")", sep=""))) +
  geom_vline(aes(xintercept=`Log(RR)`[Approach=="Full sample"]),
             linewidth=0.5, linetype="dashed") +
  geom_errorbar(aes(xmin = `Log(RR)` - 1.96*SE, 
                    xmax = `Log(RR)` + 1.96*SE), linewidth=0.5, width=0.25) +
  geom_point(aes(x=`Log(RR)`), size=2) +
  scale_x_continuous(limits=c(0.0, 0.91)) +
  scale_y_discrete(limits=rev) +
  scale_color_manual(values=palette) +
  guides(color="none") +
  thm

jpeg("../figures/rr_est5.jpeg", height=3.5, width=3, units="in", res=300)
print(p4)
dev.off()


# Panel plot --------------------------------------------------------------

jpeg("../figures/panel.jpeg", height=5, width=8, units="in", res=300)
p3 + p4 + p2 + p1 + plot_layout(ncol=2, nrow=2, guides="collect") & thm
dev.off()
