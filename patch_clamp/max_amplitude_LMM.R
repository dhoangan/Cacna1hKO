library(lme4)
library(ggplot2)
library(lmerTest)
library(dplyr)
library(MASS)
library(car)
library(MCMCglmm)
library(DHARMa)
library(glmmTMB)



data_table = read.csv("data_patch_clamp_-35mV_15mV.csv")

# response variable: I_pA_IV_diff (T-type current) and I_pA_IV40 (Non-T-type current)
# fixed factor: genotype
# randomfactor: animal


boxplot(I_pA_IV_diff ~ genotype:U_mV, data_table)

df_LMM_1 = data.frame()
df_LMM_2 = data.frame()

df_LMM_TMB = data.frame()

data_subset_1 = subset(data_table,U_mV==-35) # T-type
data_subset_2 = subset(data_table,U_mV==15) # Non-T-type


boxplot(I_pA_IV_diff ~ genotype, data_subset_1, main = paste("I_pA_IV_diff (T-type current) at -35 mV "))
boxplot(I_pA_IV40 ~ genotype, data_subset_2, main = paste("I_pA_IV40 (Non-T-type current) at 15 mV "))


####data_subset_1#######

d1_total.full = lmer(I_pA_IV_diff ~ 1 + (1|animal), data=data_subset_1, REML=FALSE)
d1_total.model_inter = lmer(I_pA_IV_diff ~ genotype + (1|animal) , data=data_subset_1, REML=FALSE)

d1_model_inter = ((anova(d1_total.full, d1_total.model_inter)))
print(d1_model_inter)
d1_p_value = d1_model_inter["d1_total.model_inter","Pr(>Chisq)"]


df_LMM_1 = data.frame(-35,"T-type",d1_p_value)

plot(fitted(d1_total.model_inter),
     resid(d1_total.model_inter, type = "pearson"))
abline(0, 0, col="red")


####data_subset_2#######


d2_total.full = lmer(I_pA_IV40 ~ 1 + (1|animal), data=data_subset_2, REML=FALSE)
d2_total.model_inter = lmer(I_pA_IV40 ~ genotype + (1|animal) , data=data_subset_2, REML=FALSE)

d2_model_inter = ((anova(d2_total.full, d2_total.model_inter)))
print(d2_model_inter)
d2_p_value = d2_model_inter["d2_total.model_inter","Pr(>Chisq)"]


df_LMM_2 = data.frame(15,"Non-T-type",d2_p_value)

plot(fitted(d2_total.model_inter),
     resid(d2_total.model_inter, type = "pearson"))
abline(0, 0, col="red")


sessionInfo()
