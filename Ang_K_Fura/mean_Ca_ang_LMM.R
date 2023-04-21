library(lme4)
library(ggplot2)
library(lmerTest)
library(dplyr)
library(MASS)
library(car)
library(MCMCglmm)
library(DHARMa)
library(glmmTMB)



data_table = read.csv("data_total_ang_dropped_mean_ca.csv")

# response variable: mean_ca_corr
# fixed factor: genotype
# randomfactor: animal, recording


boxplot(mean_ca_corr ~ genotype:potassium:angiotensin, data_table)

pot_list = unique(data_table["potassium"])[[1]]
pot_list_len = length(pot_list)

ang_list = unique(data_table["angiotensin"])[[1]]
ang_list_len = length(ang_list)

df_mean_ang_LMM = data.frame()
df_mean_ang_LMM_TMB = data.frame()

data_table_wt = subset(data_table, genotype=="wt")

mean_per_ang = aggregate(data_table_wt$mean_ca_corr, list(data_table_wt$angiotensin), FUN=mean)
count_per_ang = aggregate(data_table_wt$mean_ca_corr, list(data_table_wt$angiotensin), FUN=length)


# Find the correct underlying probability distribution
p = 4
a = 200
data_subset = subset(data_table, potassium == p& angiotensin==a)
ggplot(data_subset, aes(mean_ca_corr, fill = genotype)) + geom_histogram()
m1.norm <- glmmTMB(mean_ca_corr ~ genotype + (1|recording), 
                   data = data_subset)
#m1.lnorm <- glmmTMB(mean_ca_corr ~ genotype + (1|animal) + (1|animal:recording), 
#                    family = gaussian(link="log"), data = data_table)
#m1.gamma <- glmmTMB(mean_ca_corr ~ genotype + (1|animal) + (1|animal:recording), 
#                    family = Gamma(link="log"), data = data_table)
#m1.tweedie <- glmmTMB(mean_ca_corr ~ genotype + (1|animal) + (1|animal:recording), 
#                      family = tweedie(link="log"), data = data_table)

# Selcet one of the distributions
res <- simulateResiduals(m1.norm, plot = T)
#res <- simulateResiduals(m1.lnorm, plot = T)
#res <- simulateResiduals(m1.gamma, plot = T)
#res <- simulateResiduals(m1.tweedie, plot = T)

#select normal

for (a in ang_list[2:8]){
  print(a)
  
  for (p in pot_list){
    print(p)
    data_subset = subset(data_table, potassium == p& angiotensin==a)
    
    # reduced base model
    model.base <- glmmTMB(mean_ca_corr ~ 1 + (1|recording),
                          data = data_subset)
    
    # full model
    model <- glmmTMB(mean_ca_corr ~ genotype + (1|recording),
                     data = data_subset)
    
    anova_result = anova(model.base, model)
    print(anova_result)
    
    df_con = data.frame(p,a,anova_result$`Pr(>Chisq)`[2])
    df_mean_ang_LMM_TMB = rbind(df_mean_ang_LMM_TMB, df_con)
    
    
  }
}
# date
today = Sys.Date()
print(today)
today = format(today, format="%y%m%d")

#export anova results
write.csv(df_mean_ang_LMM_TMB,paste("anova_result_meanCa_lmmTMB_ang_",today,".csv",sep=""))
write.csv(df_mean_ang_LMM,paste("anova_result_meanCa_LMM_ang_",today,".csv",sep=""))

sessionInfo()
