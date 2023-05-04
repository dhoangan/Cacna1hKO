library(lme4)
library(ggplot2)
library(lmerTest)
library(MASS)
library(car)
library(MCMCglmm)
library(DHARMa)
library(glmmTMB)
library(ggbeeswarm)



data = read.csv("ttap2_freq_WT_pos2_positive.csv", sep=",")

data_table <- subset(data, select=c( "recording", 
                                    "genotype_3", 
                                    "cell", 
                                    "freq"))

## make variables as factor
data_table$recording <- as.factor(data_table$recording)
data_table$cell <- as.factor(data_table$cell)
data_table$genotype_3 <- as.factor(data_table$genotype_3)
data_table = data_table[!is.na(data_table$freq),]


# response variable: intraburst_freq
# fixed factor: genotype_2
# randomfactor: recording


################################

################################

###### Plots 

ggplot(aes(y=freq, x=genotype_3), data = data_table) + geom_boxplot()

ggplot(aes(genotype_3, freq), data = data_table) + geom_beeswarm(dodge.width = 1)
ggplot(aes(genotype_3, freq), data = data_table) + geom_quasirandom(dodge.width=.8)

# Check data
unique(data_table["genotype_3"])[[1]]

# Make subsets
data_subset_1 = subset(data_table, (genotype_3 == "wt_1") | (genotype_3 == "wt_ctrl_1"))
data_subset_2 = subset(data_table, (genotype_3 == "wt_2") | (genotype_3 == "wt_ctrl_2"))

ggplot(aes(y=freq, x=genotype_3), data = data_subset_2) + geom_boxplot()
ggplot(data_subset_2, aes(freq, fill=genotype_3)) + geom_histogram()

mean(subset(data_subset_1, genotype_3 == "wt_1")$freq)
mean(subset(data_subset_2, genotype_3 == "wt_2")$freq)
mean(subset(data_subset_1, genotype_3 == "wt_ctrl_1")$freq)
mean(subset(data_subset_2, genotype_3 == "wt_ctrl_2")$freq)

unique(data_subset_2["genotype_3"])[[1]]


################################

################################

# Find the correct underlying probability distribution

data_subset_A = data_subset_1

ggplot(data_subset_A, aes(freq, fill=genotype_3)) + geom_histogram()

m1.norm <- glmmTMB(freq ~ genotype_3 + (1|recording), 
                 data = data_subset_A)
m1.lnorm <- glmmTMB(freq ~ genotype_3 + (1|recording), 
                  family = gaussian(link="log"), data = data_subset_A)
m1.gamma <- glmmTMB(freq ~ genotype_3 + (1|recording), 
                  family = Gamma(link="log"), data = data_subset_A)
#m1.tweedie <- glmmTMB(freq ~ genotype_3 +(1|recording), 
#                    family = tweedie(link="log"), data = data_subset_A)


# Select one of the distributions
res <- simulateResiduals(m1.norm, plot = T)
res <- simulateResiduals(m1.lnorm, plot = T)
res <- simulateResiduals(m1.gamma, plot = T)
#res <- simulateResiduals(m1.tweedie, plot = T)

################################

################################

# use gamma distribution

##############################

##### Create models for pos1 ####

# reduced base model
model_1.base <- glmmTMB(freq ~ 1 + (1|recording),family = Gamma(link="log"),
                      data = data_subset_1)

# full model
model_1 <- glmmTMB(freq ~ genotype_3 + (1|recording),family = Gamma(link="log"),
                 data = data_subset_1)
summary(model_1)

simulationOutput <- simulateResiduals(model_1, 1000)
plot(simulationOutput)

plot(fitted(model_1), residuals(model_1), xlab = "Fitted Values", ylab = "Residuals")
abline(h = 0, lty = 2)
lines(mean(fitted(model_1)), mean(residuals(model_1)))

# compare models
anova_result_1 = anova(model_1.base, model_1)
print(anova_result_1)

##############################

##############################

#### Create models for pos2 ####

# reduced base model
model_2.base <- glmmTMB(freq ~ 1 + (1|recording),family = Gamma(link="log"),
                         data = data_subset_2)

# full model
model_2 <- glmmTMB(freq ~ genotype_3 + (1|recording),family = Gamma(link="log"),
                    data = data_subset_2)
summary(model_2)

simulationOutput <- simulateResiduals(model_2, 1000)
plot(simulationOutput)

plot(fitted(model_2), residuals(model_2), xlab = "Fitted Values", ylab = "Residuals")
abline(h = 0, lty = 2)
lines(mean(fitted(model_2)), mean(residuals(model_2)))

# compare models
anova_result_2 = anova(model_2.base, model_2)
print(anova_result_2)

################################

################################

################################


