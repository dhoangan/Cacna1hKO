library(lme4)
library(ggplot2)
library(lmerTest)
library(MASS)
library(car)
library(MCMCglmm)
library(DHARMa)
library(glmmTMB)
library(ggbeeswarm)

data_table = read.csv("data_filtered_Calbryte_activity_burstparam_230220.csv", sep=",")

# response variable: intraburst_freq
# fixed factor: genotype
# randomfactor: animal, recording
data_table = data_table[!is.na(data_table$intraburst_freq),]
data_table = subset(data_table, f_count>5)


boxplot(intraburst_freq ~ genotype, data_table)
ggplot(aes(genotype, intraburst_freq), data = data_table) +  geom_point()  

# Find the correct underlying probability distribution
hist(data_table$intraburst_freq, breaks=50)
ggplot(aes(genotype, intraburst_freq), data = data_table) + geom_beeswarm(dodge.width = 1)

mean(data_table[data_table$genotype == 'wt', 'intraburst_freq'])
mean(data_table[data_table$genotype == 'KO', 'intraburst_freq'])

m1.norm <- glmmTMB(intraburst_freq ~ genotype + (1|recording), 
                    data = data_table)
m1.lnorm <- glmmTMB(intraburst_freq ~ genotype + (1|recording), 
                    family = gaussian(link="log"), data = data_table)
m1.gamma <- glmmTMB(intraburst_freq ~ genotype + (1|recording), 
                    family = Gamma(link="log"), data = data_table)
#m1.tweedie <- glmmTMB(intraburst_freq ~ genotype + (1|recording), 
#          family = tweedie(link="log"), data = data_table)

# Selcet one of the distributions
res <- simulateResiduals(m1.norm, plot = T)
res <- simulateResiduals(m1.lnorm, plot = T)
res <- simulateResiduals(m1.gamma, plot = T)
#res <- simulateResiduals(m1.tweedie, plot = T)


# Create models
# use normal distribution
# reduced base model
model.base <- glmmTMB(intraburst_freq ~ 1 + (1|recording), family = Gamma(link="log"),
                       data = data_table)

# full model
model <- glmmTMB(intraburst_freq ~ genotype + (1|recording), family = Gamma(link="log"),
                  data = data_table)

summary(model)

simulationOutput <- simulateResiduals(model, 1000)
plot(simulationOutput)

plot(fitted(model), residuals(model), xlab = "Fitted Values", ylab = "Residuals")
abline(h = 0, lty = 2)
lines(mean(fitted(model)), mean(residuals(model)))

# compare models
anova_result = anova(model.base, model)
print(anova_result)

# date
today = Sys.Date()
print(today)
today = format(today, format="%y%m%d")

#export anova results
write.csv(anova_result,paste("anova_result_intraburst_freq_7p5_cal_",today,".csv",sep=""))

sessionInfo()

