library(lme4)
library(ggplot2)
library(lmerTest)
library(dplyr)
library(MASS)
library(car)
library(MCMCglmm)
library(DHARMa)
library(glmmTMB)



data_table = read.csv("data_total_pot_dropped_mean_ca.csv")

# response variable: mean_ca_corr
# fixed factor: genotype
# randomfactor: animal, recording


boxplot(mean_ca_corr ~ genotype:potassium:angiotensin, data_table)

pot_list = unique(data_table["potassium"])[[1]]
pot_list_len = length(pot_list)

ang_list = unique(data_table["angiotensin"])[[1]]
ang_list_len = length(ang_list)

df_mean_pot_LMM = data.frame()
df_mean_pot_LMM_TMB = data.frame()

data_table_wt = subset(data_table, genotype=="wt")


mean_per_pot = aggregate(data_table_wt$mean_ca_corr, list(data_table_wt$potassium), FUN=mean)
count_per_pot = aggregate(data_table_wt$mean_ca_corr, list(data_table_wt$potassium), FUN=length)


for (a in ang_list){
  print(a)
  
  for (p in pot_list){
    print(p)
    
    data_subset = subset(data_table, potassium == p& angiotensin==a)
    #data_subset = na.omit(data_subset) # Remove all rows with NAs
    
    
    boxplot(mean_ca_corr ~ genotype, data_subset, main = paste(toString(p),  "mM K+,", toString(a), "pM Ang II", sep = " "))
    
    
    total.full = lmer(mean_ca_corr ~ 1 + (1|animal) + (1|animal:recording), data=data_subset, REML=FALSE)
    total.model_inter = lmer(mean_ca_corr ~ genotype + (1|animal) + (1|animal:recording), data=data_subset, REML=FALSE)
    
    a_model_inter = ((anova(total.full, total.model_inter)))
    print(a_model_inter)
    p_value = a_model_inter["total.model_inter","Pr(>Chisq)"]
    
    
    df_con = data.frame(p,a,p_value)
    df_mean_pot_LMM = rbind(df_mean_pot_LMM, df_con)
    
    plot(fitted(total.model_inter),
         resid(total.model_inter, type = "pearson"), main = paste(toString(p),  "mM K+,", toString(a), "pM Ang II", sep = " "))
    abline(0, 0, col="red")
  }
}


# Find the correct underlying probability distribution
p = 5
a = 100
data_subset = subset(data_table, potassium == p& angiotensin==a)
ggplot(data_subset, aes(mean_ca_corr, fill = genotype)) + geom_histogram()
m1.norm <- glmmTMB(mean_ca_corr ~ genotype + (1|recording), 
                   data = data_subset)
m1.lnorm <- glmmTMB(mean_ca_corr ~ genotype + (1|recording), 
                    family = gaussian(link="log"), data = data_table)
m1.gamma <- glmmTMB(mean_ca_corr ~ genotype + (1|recording), 
                    family = Gamma(link="log"), data = data_table)
m1.tweedie <- glmmTMB(mean_ca_corr ~ genotype + (1|recording), 
                      family = tweedie(link="log"), data = data_table)

# Selcet one of the distributions
res <- simulateResiduals(m1.norm, plot = T)
res <- simulateResiduals(m1.lnorm, plot = T)
res <- simulateResiduals(m1.gamma, plot = T)
res <- simulateResiduals(m1.tweedie, plot = T)


#select normal
for (a in ang_list){
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
    df_mean_pot_LMM_TMB = rbind(df_mean_pot_LMM_TMB, df_con)
    
    
  }
}
# date
today = Sys.Date()
print(today)
today = format(today, format="%y%m%d")

#export anova results
write.csv(df_mean_pot_LMM_TMB,paste("anova_result_meanCa_lmmTMB_pot_",today,".csv",sep=""))
write.csv(df_mean_pot_LMM,paste("anova_result_meanCa_LMM_pot_",today,".csv",sep=""))

sessionInfo()
