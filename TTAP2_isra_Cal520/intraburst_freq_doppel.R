library(lme4)
library(ggplot2)
library(lmerTest)
library(MASS)
library(car)
library(MCMCglmm)
library(DHARMa)
library(glmmTMB)
library(ggbeeswarm)



data = read.csv("doppel_perf_rel_activity_freq.csv", sep=",")

data_table <- subset(data, select=c("animal", 
                                    "recording", 
                                    "genotype_2", 
                                    "cell", 
                                    "intraburst_freq", 
                                    "count",
                                    "pos", 
                                    "control"))

## make variables as factor
data_table$pos <- as.factor(data_table$pos)
data_table$animal <- as.factor(data_table$animal)
data_table$recording <- as.factor(data_table$recording)
data_table$cell <- as.factor(data_table$cell)
data_table$genotype_2 <- as.factor(data_table$genotype_2)
data_table = subset(data_table, count > 5)
data_table = data_table[!is.na(data_table$intraburst_freq),]


# response variable: intraburst_freq
# fixed factor: genotype_2
# randomfactor: recording


################################

################################

###### Plots 

#boxplot(intraburst_freq ~ genotype_2, data_table)
#ggplot(aes(pos, intraburst_freq, color=genotype_2), data = data_table) +  geom_point()
ggplot(aes(y=intraburst_freq, x=pos, fill=genotype_2), data = data_table) + geom_boxplot()

ggplot(aes(pos, intraburst_freq, color=genotype_2), data = data_table) + geom_beeswarm(dodge.width = 1)
ggplot(aes(pos, intraburst_freq, color=genotype_2), data = data_table) + geom_quasirandom(dodge.width=.8)

# Check data
unique(data_table["genotype_2"])[[1]]
unique(data_table["pos"])[[1]]
unique(data_table["control"])[[1]]

# Make subsets
data_subset_1 = subset(data_table, (pos == 1))
data_subset_2 = subset(data_table, pos==2)

unique(data_subset_1["genotype_2"])[[1]]
unique(data_subset_1["pos"])[[1]]
unique(data_subset_1["control"])[[1]]

ggplot(aes(y=intraburst_freq, x=genotype_2), data = data_subset_1) + geom_boxplot()
ggplot(data_subset_1, aes(intraburst_freq, fill=genotype_2)) + geom_histogram()

mean(subset(data_subset_1, genotype_2 == "wt_ctrl")$intraburst_freq)
mean(subset(data_subset_1, genotype_2 == "wt")$intraburst_freq)

################################

################################

# Find the correct underlying probability distribution

data_subset_A = data_subset_1

ggplot(data_subset_A, aes(intraburst_freq, fill=genotype_2)) + geom_histogram()

m1.norm <- glmmTMB(intraburst_freq ~ genotype_2 + (1|recording), 
                 data = data_subset_A)
m1.lnorm <- glmmTMB(intraburst_freq ~ genotype_2 + (1|recording), 
                  family = gaussian(link="log"), data = data_subset_A)
m1.gamma <- glmmTMB(intraburst_freq ~ genotype_2 + (1|recording), 
                  family = Gamma(link="log"), data = data_subset_A)
#m1.tweedie <- glmmTMB(intraburst_freq ~ genotype_2 +(1|recording), 
#                   family = tweedie(link="log"), data = data_subset_A)


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
model_pos1.base <- glmmTMB(intraburst_freq ~ 1 + (1|recording),family = Gamma(link="log"),
                      data = data_subset_1)

# full model
model_pos1 <- glmmTMB(intraburst_freq ~ genotype_2 + (1|recording),family = Gamma(link="log"),
                 data = data_subset_1)
summary(model_pos1)

simulationOutput <- simulateResiduals(model_pos1, 1000)
plot(simulationOutput)

plot(fitted(model_pos1), residuals(model_pos1), xlab = "Fitted Values", ylab = "Residuals")
abline(h = 0, lty = 2)
lines(mean(fitted(model_pos1)), mean(residuals(model_pos1)))

# compare models
anova_result_pos1 = anova(model_pos1.base, model_pos1)
print(anova_result_pos1)

##############################

##############################

##### Create models for pos2 ####

# reduced base model
model_pos2.base <- glmmTMB(intraburst_freq ~ 1 + (1|recording),family = Gamma(link="log"),
                           data = data_subset_2)

# full model
model_pos2 <- glmmTMB(intraburst_freq ~ genotype_2 + (1|recording),family = Gamma(link="log"),
                      data = data_subset_2)
summary(model_pos2)

simulationOutput <- simulateResiduals(model_pos2, 1000)
plot(simulationOutput)

plot(fitted(model_pos2), residuals(model_pos2), xlab = "Fitted Values", ylab = "Residuals")
abline(h = 0, lty = 2)
lines(mean(fitted(model_pos2)), mean(residuals(model_pos2)))

# compare models
anova_result_pos2 = anova(model_pos2.base, model_pos2)
print(anova_result_pos2)

##############################

##############################

################################

################################

################################


