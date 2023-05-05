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
                                    "activity", 
                                    "pos", 
                                    "control"))

## make variables as factor
data_table$pos <- as.factor(data_table$pos)
data_table$animal <- as.factor(data_table$animal)
data_table$recording <- as.factor(data_table$recording)
data_table$cell <- as.factor(data_table$cell)
data_table$genotype_2 <- as.factor(data_table$genotype_2)



# response variable: activity
# fixed factor: genotype_2
# random factor: recording


################################

################################

###### Plots 

#boxplot(activity ~ genotype_2, data_table)
#ggplot(aes(pos, activity, color=genotype_2), data = data_table) +  geom_point()
ggplot(aes(y=activity, x=pos, fill=genotype_2), data = data_table) + geom_boxplot()

ggplot(aes(pos, activity, color=genotype_2), data = data_table) + geom_beeswarm(dodge.width = 1)
ggplot(aes(pos, activity, color=genotype_2), data = data_table) + geom_quasirandom(dodge.width=.8)

# Check data
unique(data_table["genotype_2"])[[1]]
unique(data_table["pos"])[[1]]
unique(data_table["control"])[[1]]

# Make subsets
data_subset_1 = subset(data_table, (pos == 2))
data_subset_2 = subset(data_table, (pos==3))

ggplot(aes(y=activity, x=pos, fill=genotype_2), data = data_subset_1) + geom_boxplot()
ggplot(data_subset_1, aes(activity, fill=genotype_2)) + geom_histogram()

mean(subset(data_subset_1, (genotype_2 == "wt"))$activity)
mean(subset(data_subset_1, (genotype_2 == "wt_ctrl"))$activity)


################################

################################

# Find the correct underlying probability distribution

data_subset_A = data_subset_1

ggplot(data_subset_A, aes(activity, fill=genotype_2)) + geom_histogram()

#m1.norm <- glmmTMB(activity ~ genotype_2 + (1|recording), 
#                 data = data_subset_A)
#m1.lnorm <- glmmTMB(activity ~ genotype_2 + (1|recording), 
#                  family = gaussian(link="log"), data = data_subset_A)
#m1.gamma <- glmmTMB(activity ~ genotype_2 + (1|recording), 
#                  family = Gamma(link="log"), data = data_subset_A)
m1.tweedie <- glmmTMB(activity ~ genotype_2 ++ (1|recording), 
                  family = tweedie(link="log"), data = data_subset_A)


# Select one of the distributions
#res <- simulateResiduals(m1.norm, plot = T)
#res <- simulateResiduals(m1.lnorm, plot = T)
#res <- simulateResiduals(m1.gamma, plot = T)
res <- simulateResiduals(m1.tweedie, plot = T)

################################

################################

# use  tweedie for distribution

##############################

##### Create models for wt ####

# reduced base model
model_1.base <- glmmTMB(activity ~ 1 + (1|recording),family = tweedie(link="log"),
                      data = data_subset_1)

# full model
model_1 <- glmmTMB(activity ~ genotype_2 + (1|recording),family = tweedie(link="log"),
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

#### Create models for KO ####

# reduced base model
model_2.base <- glmmTMB(activity ~ 1 + (1|recording),family = tweedie(link="log"),
                         data = data_subset_2)

# full model
model_2 <- glmmTMB(activity ~ genotype_2 + (1|recording),family = tweedie(link="log"),
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


