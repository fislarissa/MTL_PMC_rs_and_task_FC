##########################################################################################################################

#Longitudinal functional connectivity during rest and task is differentially 
#related to Alzheimer’s pathology and episodic memory in older adults

#Analysis of PREVENT-AD data - longitudinal resting-state and task fMRI with cognitive and PET markers
#Larissa Fischer - Multimodal Neuroimaging Lab, DZNE Magdeburg, 2024/25

#########################################################################################################################

#Import packages:
library(readxl)
library(writexl)
library(dplyr)
library(tidyr)
library(psych)
library(car)
library(lme4)
library(lmtest)
library(sandwich)
library(MASS)
library(ggplot2)
library(sjPlot)
library(sjmisc)
library(sjlabelled)
library(boot.pval)
library(Hmisc)
library(lmerTest)
library(bptest)
library(MuMIn)
library(sjstats)


################################################################################
#set working directory:
setwd("/Users/your/path")

#load data
load("/.../data_analysis_6_FC_all.RData")

################################################################################
#subsets and theme

data_analysis_6_FC_all_non_APOE <- subset(data_analysis_6_FC_all, APOE4_carrier == 0)
data_analysis_6_FC_all_APOE <- subset(data_analysis_6_FC_all, APOE4_carrier == 1)
long_data_all_non_APOE <- subset(long_data_all, APOE4_carrier == 0)
long_data_all_APOE <- subset(long_data_all, APOE4_carrier == 1)

plot_theme <- function() {
  theme(
    panel.border = element_blank(),
    axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
    axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", colour = NA),
    text = element_text(size = 16),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 16)
  )
}


################################################################################
#ANALYSES
################################################################################

################################################################################
# Demographics
################################################################################
#demographics by APOE
summary(data_analysis_6_FC_all$APOE4_carrier)
describeBy(data_analysis_6_FC_all$age_years, data_analysis_6_FC_all$APOE4_carrier)
t.test(age_years ~ APOE4_carrier, data = data_analysis_6_FC_all)
describeBy(data_analysis_6_FC_all$education_years, data_analysis_6_FC_all$APOE4_carrier)
t.test(education_years ~ APOE4_carrier, data = data_analysis_6_FC_all)
summary(data_analysis_6_FC_all_non_APOE$sex) #0 = male
chisq.test(data_analysis_6_FC_all$APOE4_carrier, data_analysis_6_FC_all$sex) 
describeBy(data_analysis_6_FC_all$amyloid_WB, data_analysis_6_FC_all$APOE4_carrier)
t.test(amyloid_WB_boxcox ~ APOE4_carrier, data=data_analysis_6_FC_all)
describeBy(data_analysis_6_FC_all$EC_tau_PET_LR_mean, data_analysis_6_FC_all$APOE4_carrier)
t.test(EC_tau_PET_LR_mean ~ APOE4_carrier, data=data_analysis_6_FC_all) 
#overview amyloid positive
data_analysis_6_FC_all$Amyloid_Group <- factor(ifelse(data_analysis_6_FC_all$amyloid_WB > 1.39, 1, 0)) 
table(data_analysis_6_FC_all$Amyloid_Group, data_analysis_6_FC_all$APOE4_carrier)
chi_test_Abeta_APOE <- chisq.test(data_analysis_6_FC_all$Amyloid_Group, data_analysis_6_FC_all$APOE4_carrier)
print(chi_test_Abeta_APOE)
#overview tau positive
data_analysis_6_FC_all$Tau_Group <- factor(ifelse(data_analysis_6_FC_all$EC_tau_PET_LR_mean > 1.30, 1, 0)) 
table_tau_apoe <- table(data_analysis_6_FC_all$Tau_Group, data_analysis_6_FC_all$APOE4_carrier)
table_tau_apoe
fisher_test_tau_apoe <- fisher.test(table_tau_apoe)
print(fisher_test_tau_apoe)


################################################################################
################################################################################
    ### effects of covariates on FC ###
################################################################################
################################################################################
# repeat for within MTL, within PMC, and between MTL and PMC for rsFC, encoding-FC,
#and retrieval-FC
within_MTL_time_covariates <- lmer(within_MTL~ time_days_scaled*APOE4_carrier + age_years + sex + education_years  + (1 |subj), data=long_data_all, REML = F)
summary(within_MTL_time_covariates)
plot(resid(within_MTL_time_covariates) ~ fitted(within_MTL_time_covariates)) 
qqnorm(resid(within_MTL_time_covariates))
qqline(resid(within_MTL_time_covariates)) 
vif(within_MTL_time_covariates)
tab_model(within_MTL_time_covariates, df.method = "satterthwaite", show.stat = TRUE, show.se = TRUE, show.std = TRUE, file = "within_MTL_rs_time_covariates_time_continuous_int.doc")

#does rsFC slope correlate with enc FC slope?
cor.test(data_analysis_6_FC_all$slope_within_MTL, data_analysis_6_FC_all$slope_within_MTL,  use = "complete.obs")
#repeat for all correlations


################################################################################
### effects of covariates on COGNITION ###
################################################################################
#Change in delayed memory score of the RBANS with covariates:
RBANS_time <- lmer(RBANS_memory_long ~ time_days_scaled+APOE4_carrier + age_years + sex + education_years + (1+time_days_scaled|subj), data=long_data_all, REML = F)
summary(RBANS_time) 
plot(resid(RBANS_time) ~ fitted(RBANS_time)) 
qqnorm(resid(RBANS_time))
qqline(resid(RBANS_time)) 
vif(RBANS_time)
tab_model(RBANS_time,df.method = "satterthwaite" , show.stat = TRUE, show.se = TRUE, show.std = TRUE, file = "RBANS_time.doc") 

plot_RBANS_time <- ggplot(long_data_all, aes(time_days_scaled, RBANS_memory_long, colour = APOE4_carrier)) + 
  geom_point(alpha = 0.6) +
  geom_smooth(aes(time_days_scaled, RBANS_memory_long), method = lm) +
  xlab ("Time") + ylab ("RBANS delayed memory index score") +
  scale_color_manual(labels = c("Non-Carrier", "Carrier"), values = c("#1A85FF", "#D41159")) + labs(color = "APOE4 carrier")+
  scale_x_continuous(breaks = c(-1.3814, -0.6552, 0.05057, 1.4625), labels = c("BL", "12 Months", "24 Months","48 Months")) 
plot_RBANS_time

################################################################################
#Change in corrected hit rate with covariates:
corrhr_time <- lmer(HR_corr ~ time_days_scaled+APOE4_carrier+ age_years + sex + education_years + (1|subj), data=long_data_all, REML = F)
summary(corrhr_time) 
plot(resid(within_MTL_time_covariates) ~ fitted(within_MTL_time_covariates)) 
qqnorm(resid(within_MTL_time_covariates))
qqline(resid(within_MTL_time_covariates)) 
vif(within_MTL_time_covariates)
tab_model(corrhr_time,df.method = "satterthwaite", show.stat = TRUE, show.se = TRUE, show.std = TRUE, file = "corrhr_time.doc") 

plot_corrhr_time <- ggplot(long_data_all, aes(time_days_scaled, HR_corr, colour = APOE4_carrier)) + 
  geom_point(alpha = 0.6) +
  geom_smooth(aes(time_days_scaled, HR_corr), method = lm) +
  xlab ("Time") + ylab ("Corrected hit rate fMRI retrieval task") +
  scale_color_manual(labels = c("Non-Carrier", "Carrier"), values = c("#1A85FF", "#D41159")) + labs(color = "APOE4 carrier")+
  scale_x_continuous(breaks = c(-1.3814, -0.6552, 0.05057, 1.4625), labels = c("BL", "12 Months", "24 Months","48 Months")) 
plot_corrhr_time


################################################################################
################################################################################
### Resting-state ###
################################################################################
################################################################################

################################################################################
#effect of FC  slope on later AD pathology via PET
################################################################################
#Slope (change in) FC on amyloid: 
SLOPEvsAMY <- lm(amyloid_WB_boxcox ~ slope_within_MTL*APOE4_carrier+slope_within_PMC*APOE4_carrier+slope_between_MTL.PMC*APOE4_carrier+age_years +sex+ education_years + time_BL_to_PET, data=data_analysis_6_FC_all)
summary(SLOPEvsAMY) 
qqnorm(resid(SLOPEvsAMY))
qqline(resid(SLOPEvsAMY)) 
bptest(SLOPEvsAMY) 
vif(SLOPEvsAMY) 
tab_model(SLOPEvsAMY, df.method = "satterthwaite", show.stat = TRUE, show.se = TRUE, show.std = TRUE, file = "RS_SLOPEvsAMY_int.doc")
#robust SE
options(scipen = 999)
coeftest(SLOPEvsAMY, vcov. = vcovHC(SLOPEvsAMY, type = "HC1"))
robust_vcov <- vcovHC(SLOPEvsAMY, type = "HC1")
coefci(SLOPEvsAMY, vcov. = robust_vcov, level = 0.95)
#standardized
SLOPEvsAMY_scaled <- lm(scale(amyloid_WB_boxcox) ~ scale(slope_within_MTL)*APOE4_carrier+scale(slope_within_PMC)*APOE4_carrier+scale(slope_between_MTL.PMC)*APOE4_carrier+scale(age_years)+sex+scale(education_years)+ scale(time_BL_to_PET), data=data_analysis_6_FC_all)
options(scipen = 999)
coeftest(SLOPEvsAMY_scaled, vcov. = vcovHC(SLOPEvsAMY_scaled, type = "HC1"))
robust_vcov <- vcovHC(SLOPEvsAMY_scaled, type = "HC1")
coefci(SLOPEvsAMY_scaled, vcov. = robust_vcov, level = 0.95)

#APOE groups
SLOPEvsAMY_APOE <- lm(amyloid_WB_boxcox ~ slope_within_MTL+slope_within_PMC+slope_between_MTL.PMC +age_years+sex+ education_years +time_BL_to_PET, data=data_analysis_6_FC_all_APOE)
summary(SLOPEvsAMY_APOE) 
qqnorm(resid(SLOPEvsAMY_APOE))
qqline(resid(SLOPEvsAMY_APOE)) 
bptest(SLOPEvsAMY_APOE) 
vif(SLOPEvsAMY_APOE)
tab_model(SLOPEvsAMY_APOE, df.method = "satterthwaite", show.stat = TRUE, show.se = TRUE, show.std = TRUE, file = "RS_SLOPEvsAMY_APOE.doc") 

SLOPEvsAMY_non_APOE <- lm(amyloid_WB_boxcox ~ slope_within_MTL+slope_within_PMC+slope_between_MTL.PMC +age_years+ sex+education_years +time_BL_to_PET, data=data_analysis_6_FC_all_non_APOE)
summary(SLOPEvsAMY_non_APOE) 
qqnorm(resid(SLOPEvsAMY_non_APOE))
qqline(resid(SLOPEvsAMY_non_APOE)) 
bptest(SLOPEvsAMY_non_APOE) 
vif(SLOPEvsAMY_non_APOE)
tab_model(SLOPEvsAMY_non_APOE, df.method = "satterthwaite", show.stat = TRUE, show.se = TRUE, show.std = TRUE, file = "RS_SLOPEvsAMY_non_APOE.doc") 

plot_SLOPEvsAMY_line <- ggplot(data_analysis_6_FC_all, aes(slope_within_PMC, amyloid_WB_boxcox, color = APOE4_carrier)) + 
  xlab("Slope of FC strength within PMC during rest") + ylab("Global neocortical amyloid (PET)") +
  geom_point()+geom_smooth(method="lm")+ labs(color = "APOE4 group")+
  scale_color_manual(labels = c("Non-Carrier", "Carrier"), values = c("#1A85FF", "#D41159"))+plot_theme()
plot_SLOPEvsAMY_line 


#Slope (change in) FC on tau: 
SLOPEvsTAU <- lm(EC_tau_PET_LR_mean_boxcox ~ slope_within_MTL*APOE4_carrier+slope_within_PMC*APOE4_carrier+slope_between_MTL.PMC*APOE4_carrier + age_years +sex+ education_years + time_BL_to_PET, data=data_analysis_6_FC_all)
summary(SLOPEvsTAU)
qqnorm(resid(SLOPEvsTAU))
qqline(resid(SLOPEvsTAU)) 
bptest(SLOPEvsTAU)
vif(SLOPEvsTAU)
tab_model(SLOPEvsTAU, df.method = "satterthwaite", show.stat = TRUE, show.se = TRUE, show.std = TRUE, file = "RS_SLOPEvsTAU_int.doc") 


################################################################################
#COGNITION
################################################################################
################################################################################
#Change in delayed memory score on the RBANS with covariates AND FC: 

#FC over time (slope)
SLOPEvsCog <- lm(slope_RBANS ~ slope_within_MTL*APOE4_carrier+slope_within_PMC*APOE4_carrier+slope_between_MTL.PMC*APOE4_carrier+age_years +sex+ education_years, data=data_analysis_6_FC_all)
summary(SLOPEvsCog) 
qqnorm(resid(SLOPEvsCog))
qqline(resid(SLOPEvsCog)) 
bptest(SLOPEvsCog)
vif(SLOPEvsCog)
tab_model(SLOPEvsCog,df.method = "satterthwaite", show.stat = TRUE, show.se = TRUE, show.std = TRUE, file = "RBANS_RS_FC_slope_int.doc") 
# plot
plot_RBANS_FC_slope_line <- ggplot(data_analysis_6_FC_all, aes(slope_between_MTL.PMC, slope_RBANS, color = APOE4_carrier)) + 
  xlab("Slope of FC strength between MTL and PMC during rest") + ylab("Slope of RBANS delayed memory index score") +
  geom_point()+geom_smooth(method="lm")+ labs(color = "APOE4")+
  scale_color_manual(labels = c("Non-Carrier", "Carrier"), values = c("#1A85FF", "#D41159")) +plot_theme()
plot_RBANS_FC_slope_line

################################################################################
#Change in corrected hit rate with covariates AND FC: 

#FC over time (slope)
SLOPEvsCog <- lm(slope_corrhr ~ slope_within_MTL*APOE4_carrier+slope_within_PMC*APOE4_carrier+slope_between_MTL.PMC*APOE4_carrier+age_years +sex+ education_years, data=data_analysis_6_FC_all)
summary(SLOPEvsCog) 
qqnorm(resid(SLOPEvsCog))
qqline(resid(SLOPEvsCog)) 
bptest(SLOPEvsCog)
vif(SLOPEvsCog)
tab_model(SLOPEvsCog,df.method = "satterthwaite", show.stat = TRUE, show.se = TRUE, show.std = TRUE, file = "corrhr_RS_FC_slope_int.doc") 


################################################################################
################################################################################
### ENCODING ###
################################################################################
################################################################################

################################################################################
#effect of FC BL and slope on later AD pathology via PET
################################################################################

#Slope (change in) FC on amyloid: 
SLOPEvsAMY <- lm(amyloid_WB_boxcox ~ slope_within_MTL_enc*APOE4_carrier+slope_within_PMC_enc*APOE4_carrier+slope_between_MTL.PMC_enc*APOE4_carrier+age_years +sex+ education_years + time_BL_to_PET, data=data_analysis_6_FC_all)
summary(SLOPEvsAMY) 
qqnorm(resid(SLOPEvsAMY))
qqline(resid(SLOPEvsAMY)) 
bptest(SLOPEvsAMY) 
vif(SLOPEvsAMY) 
confint(SLOPEvsAMY)
tab_model(SLOPEvsAMY, df.method = "satterthwaite", show.stat = TRUE, show.se = TRUE, show.std = TRUE, file = "ENC_SLOPEvsAMY_int.doc")


#Slope (change in) FC on tau: 
SLOPEvsTAU <- lm(EC_tau_PET_LR_mean_boxcox ~ slope_within_MTL_enc*APOE4_carrier+slope_within_PMC_enc*APOE4_carrier+slope_between_MTL.PMC_enc*APOE4_carrier + age_years +sex+ education_years + time_BL_to_PET, data=data_analysis_6_FC_all)
summary(SLOPEvsTAU) 
qqnorm(resid(SLOPEvsTAU))
qqline(resid(SLOPEvsTAU)) 
bptest(SLOPEvsTAU) 
vif(SLOPEvsTAU)
tab_model(SLOPEvsTAU, df.method = "satterthwaite", show.stat = TRUE, show.se = TRUE, show.std = TRUE, file = "Enc_SLOPEvsTAU.doc_int.doc")
#APOE groups:
SLOPEvsTAU_APOE <- lm(EC_tau_PET_LR_mean_boxcox ~ slope_within_MTL_enc+slope_within_PMC_enc+slope_between_MTL.PMC_enc +age_years+sex+ education_years +time_BL_to_PET, data=data_analysis_6_FC_all_APOE)
summary(SLOPEvsTAU_APOE) 
qqnorm(resid(SLOPEvsTAU_APOE))
qqline(resid(SLOPEvsTAU_APOE)) 
bptest(SLOPEvsTAU_APOE) 
vif(SLOPEvsTAU_APOE)
tab_model(SLOPEvsTAU_APOE, df.method = "satterthwaite", show.stat = TRUE, show.se = TRUE, show.std = TRUE, file = "ENC_SLOPEvsTAU_APOE.doc") 

SLOPEvsTAU_non_APOE <- lm(EC_tau_PET_LR_mean_boxcox ~ slope_within_MTL_enc+slope_within_PMC_enc+slope_between_MTL.PMC_enc +age_years+ sex+education_years +time_BL_to_PET, data=data_analysis_6_FC_all_non_APOE)
summary(SLOPEvsTAU_non_APOE) 
qqnorm(resid(SLOPEvsTAU_non_APOE))
qqline(resid(SLOPEvsTAU_non_APOE)) 
bptest(SLOPEvsTAU_non_APOE) 
vif(SLOPEvsTAU_non_APOE)
tab_model(SLOPEvsTAU_non_APOE, df.method = "satterthwaite", show.stat = TRUE, show.se = TRUE, show.std = TRUE, file = "ENC_SLOPEvsTAU_non_APOE.doc") 

plot_SLOPEvsTAU_line <- ggplot(data_analysis_6_FC_all, aes(slope_within_MTL_enc, EC_tau_PET_LR_mean_boxcox, color = APOE4_carrier)) + 
  xlab("Slope of FC strength within MTL during encoding") + 
  ylab("Entorhinal tau (PET)") + geom_point() + geom_smooth(method="lm") + labs(color = "APOE4 group") +
  scale_color_manual(labels = c("Non-Carrier", "Carrier"), values = c("#1A85FF", "#D41159")) +
  scale_y_continuous(breaks = seq(-0.2, 0.4, by = 0.1)) +
  plot_theme()
plot_SLOPEvsTAU_line


################################################################################
#COGNITION
################################################################################
#Change in delayed memory score on the RBANS with covariates AND FC: 

#FC over time (slope)
SLOPEvsCog <- lm(slope_RBANS ~ slope_within_MTL_enc*APOE4_carrier+slope_within_PMC_enc*APOE4_carrier+slope_between_MTL.PMC_enc*APOE4_carrier+age_years +sex+ education_years, data=data_analysis_6_FC_all)
summary(SLOPEvsCog) 
qqnorm(resid(SLOPEvsCog))
qqline(resid(SLOPEvsCog)) 
bptest(SLOPEvsCog)
vif(SLOPEvsCog)
tab_model(SLOPEvsCog,df.method = "satterthwaite" , show.stat = TRUE, show.se = TRUE, show.std = TRUE, file = "RBANS_ENC_FC_slope_int.doc") 

# plot
plot_RBANS_FC_slope_line <- ggplot(data_analysis_6_FC_all, aes(slope_within_PMC_enc, slope_RBANS, color = APOE4_carrier)) + 
  xlab("FC strength slope within PMC during encoding") + ylab("RBANS delayed memory index score slope") +
  geom_point()+geom_smooth(method="lm")+ labs(color = "APOE4")+
  scale_color_manual(labels = c("Non-Carrier", "Carrier"), values = c("#1A85FF", "#D41159")) 
plot_RBANS_FC_slope_line 

################################################################################
#Change in corrected hit rate with covariates AND FC: 

#FC over time (slope)
SLOPEvsCog <- lm(slope_corrhr ~ slope_within_MTL_enc*APOE4_carrier+slope_within_PMC_enc*APOE4_carrier+slope_between_MTL.PMC_enc*APOE4_carrier+age_years +sex+ education_years, data=data_analysis_6_FC_all)
summary(SLOPEvsCog) 
qqnorm(resid(SLOPEvsCog))
qqline(resid(SLOPEvsCog)) 
bptest(SLOPEvsCog)
vif(SLOPEvsCog)
tab_model(SLOPEvsCog,df.method = "satterthwaite" , show.stat = TRUE, show.se = TRUE, show.std = TRUE, file = "corrhr_ENC_FC_slope_int.doc") 

# plot
plot_corrhr_FC_slope_line <- ggplot(data_analysis_6_FC_all, aes(slope_between_MTL.PMC_enc, slope_corrhr, color = APOE4_carrier)) + 
  xlab("Slope of FC strength between MTL and PMC during encoding") + ylab("Slope of corrected hit rate (fMRI retrieval task)") +
  geom_point()+geom_smooth(method="lm")+ labs(color = "APOE4")+
  scale_color_manual(labels = c("Non-carrier", "Carrier"), values = c("#1A85FF", "#D41159")) +plot_theme()
plot_corrhr_FC_slope_line 


################################################################################
################################################################################
### RETRIEVAL ###
################################################################################
################################################################################

################################################################################
#effect of FC BL and slope on later AD pathology via PET
################################################################################

#Slope (change in) FC on amyloid: 
SLOPEvsAMY <- lm(amyloid_WB_boxcox ~ slope_within_MTL_ret*APOE4_carrier+slope_within_PMC_ret*APOE4_carrier+slope_between_MTL.PMC_ret*APOE4_carrier+age_years +sex+ education_years + time_BL_to_PET, data=data_analysis_6_FC_all)
summary(SLOPEvsAMY) 
qqnorm(resid(SLOPEvsAMY))
qqline(resid(SLOPEvsAMY)) 
bptest(SLOPEvsAMY) 
vif(SLOPEvsAMY) 
confint(SLOPEvsAMY)
tab_model(SLOPEvsAMY, df.method = "satterthwaite", show.stat = TRUE, show.se = TRUE, show.std = TRUE, file = "ret_SLOPEvsAMY_int.doc")


#Slope (change in) FC on tau: 
SLOPEvsTAU <- lm(EC_tau_PET_LR_mean_boxcox ~ slope_within_MTL_ret*APOE4_carrier+slope_within_PMC_ret*APOE4_carrier+slope_between_MTL.PMC_ret*APOE4_carrier + age_years +sex+ education_years + time_BL_to_PET, data=data_analysis_6_FC_all)
summary(SLOPEvsTAU) 
qqnorm(resid(SLOPEvsTAU))
qqline(resid(SLOPEvsTAU)) 
bptest(SLOPEvsTAU)
vif(SLOPEvsTAU)
tab_model(SLOPEvsTAU, df.method = "satterthwaite", show.stat = TRUE, show.se = TRUE, show.std = TRUE, file = "ret_SLOPEvsTAU.doc_int.doc")
#APOE groups:
SLOPEvsTAU_APOE <- lm(EC_tau_PET_LR_mean_boxcox ~ slope_within_MTL_ret+slope_within_PMC_ret+slope_between_MTL.PMC_ret +age_years+sex+ education_years +time_BL_to_PET, data=data_analysis_6_FC_all_APOE)
summary(SLOPEvsTAU_APOE) 
qqnorm(resid(SLOPEvsTAU_APOE))
qqline(resid(SLOPEvsTAU_APOE)) 
bptest(SLOPEvsTAU_APOE)
vif(SLOPEvsTAU_APOE)
tab_model(SLOPEvsTAU_APOE, df.method = "satterthwaite", show.stat = TRUE, show.se = TRUE, show.std = TRUE, file = "ret_SLOPEvsTAU_APOE.doc") 
SLOPEvsTAU_non_APOE <- lm(EC_tau_PET_LR_mean_boxcox ~ slope_within_MTL_ret+slope_within_PMC_ret+slope_between_MTL.PMC_ret +age_years+ sex+education_years +time_BL_to_PET, data=data_analysis_6_FC_all_non_APOE)
summary(SLOPEvsTAU_non_APOE) 
qqnorm(resid(SLOPEvsTAU_non_APOE))
qqline(resid(SLOPEvsTAU_non_APOE)) 
bptest(SLOPEvsTAU_non_APOE)
vif(SLOPEvsTAU_non_APOE)
tab_model(SLOPEvsTAU_non_APOE, df.method = "satterthwaite", show.stat = TRUE, show.se = TRUE, show.std = TRUE, file = "ret_SLOPEvsTAU_non_APOE.doc") 

plot_retSLOPE_PMCvsTAU_line <- ggplot(data_analysis_6_FC_all, aes(slope_within_PMC_ret, EC_tau_PET_LR_mean_boxcox, color = APOE4_carrier)) + 
  xlab("Slope of FC strength within PMC during retrieval") + ylab("Entorhinal tau (PET)") +
  geom_point()+geom_smooth(method="lm")+ labs(color = "APOE4 group")+
  scale_y_continuous(breaks = seq(-0.2, 0.4, by = 0.1)) +
  scale_color_manual(labels = c("Non-Carrier", "Carrier"), values = c("#1A85FF", "#D41159")) +plot_theme()
plot_retSLOPE_PMCvsTAU_line 

plot_retSLOPE_betweenvsTAU_line <- ggplot(data_analysis_6_FC_all, aes(slope_between_MTL.PMC_ret, EC_tau_PET_LR_mean_boxcox, color = APOE4_carrier)) + 
  xlab("Slope of FC strength between MTL and PMC during retrieval") + ylab("Entorhinal tau (PET)") +
  geom_point()+geom_smooth(method="lm")+ labs(color = "APOE4 group")+
  scale_y_continuous(breaks = seq(-0.2, 0.4, by = 0.1)) +
  scale_color_manual(labels = c("Non-Carrier", "Carrier"), values = c("#1A85FF", "#D41159")) +plot_theme()
plot_retSLOPE_betweenvsTAU_line 


################################################################################
#COGNITION
################################################################################
#Change in delayed memory score on the RBANS with covariates AND FC: 

#FC over time (slope)
SLOPEvsCog <- lm(slope_RBANS ~ slope_within_MTL_ret*APOE4_carrier+slope_within_PMC_ret*APOE4_carrier+slope_between_MTL.PMC_ret*APOE4_carrier+age_years +sex+ education_years, data=data_analysis_6_FC_all)
summary(SLOPEvsCog) 
qqnorm(resid(SLOPEvsCog))
qqline(resid(SLOPEvsCog)) 
bptest(SLOPEvsCog)
vif(SLOPEvsCog)
tab_model(SLOPEvsCog,df.method = "satterthwaite" , show.stat = TRUE, show.se = TRUE, show.std = TRUE, file = "RBANS_ret_FC_slope_int.doc") 

# plot
plot_RBANS_FC_slope_line <- ggplot(data_analysis_6_FC_all, aes(slope_within_PMC_ret, slope_RBANS, color = APOE4_carrier)) + 
  xlab("FC strength slope within PMC during retoding") + ylab("RBANS delayed memory index score slope") +
  geom_point()+geom_smooth(method="lm")+ labs(color = "APOE4")+
  scale_color_manual(labels = c("Non-Carrier", "Carrier"), values = c("#1A85FF", "#D41159")) 
plot_RBANS_FC_slope_line 

################################################################################
#Change in corrected hit rate with covariates AND FC: 

#FC over time (slope)
SLOPEvsCog <- lm(slope_corrhr ~ slope_within_MTL_ret*APOE4_carrier+slope_within_PMC_ret*APOE4_carrier+slope_between_MTL.PMC_ret*APOE4_carrier+age_years +sex+ education_years  + time_BL_to_PET, data=data_analysis_6_FC_all)
summary(SLOPEvsCog) 
qqnorm(resid(SLOPEvsCog))
qqline(resid(SLOPEvsCog)) 
bptest(SLOPEvsCog)
vif(SLOPEvsCog)
tab_model(SLOPEvsCog,df.method = "satterthwaite" , show.stat = TRUE, show.se = TRUE, show.std = TRUE, file = "corrhr_ret_FC_slope_int.doc") 

# plot
plot_corrhr_FC_slope_line <- ggplot(data_analysis_6_FC_all, aes(slope_between_MTL.PMC_ret, slope_corrhr, color = APOE4_carrier)) + 
  xlab("FC strength slope between MTL and PMC during retoding") + ylab("Corrected hit rate fMRI retrieval task slope") +
  geom_point()+geom_smooth(method="lm")+ labs(color = "APOE4")+
  scale_color_manual(labels = c("Non-Carrier", "Carrier"), values = c("#1A85FF", "#D41159")) 
plot_corrhr_FC_slope_line


# Additional analyses:
#check if  effect with RBANS attention by using RBANS_attention_long and slope_RBANS_attention in the models above
#run retrieval models for hits only in subgroup of N = 136 with data_analysis_6_FC_all_hits.RData


