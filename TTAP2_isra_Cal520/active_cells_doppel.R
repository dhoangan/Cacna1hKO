library(lme4)
library(ggplot2)
library(lmerTest)
library(MASS)
library(car)
library(MCMCglmm)
library(DHARMa)
library(glmmTMB)
library(ggbeeswarm)



data = read.csv("doppel_perf_active_cells_per_slice.csv", sep=",")

data_table <- subset(data, select=c("animal", 
                                    "recording", 
                                    "genotype_2", 
                                    "n_percentage_active_cells_pos2",
                                    "n_percentage_active_cells_pos3"))

## make variables as factor
data_table$animal <- as.factor(data_table$animal)
data_table$recording <- as.factor(data_table$recording)
data_table$genotype_2 <- as.factor(data_table$genotype_2)


# response variable: n_percentage_active_cells
# fixed factor: genotype_2
# randomfactor: recording


################################

################################

###### Plots 

ggplot(aes(y=n_percentage_active_cells_pos2, x=genotype_2), data = data_table) + geom_boxplot()
ggplot(aes(y=n_percentage_active_cells_pos3, x=genotype_2), data = data_table) + geom_boxplot()

ggplot(aes(genotype_2, n_percentage_active_cells_pos2), data = data_table) + geom_beeswarm(dodge.width = 1)
ggplot(aes(genotype_2, n_percentage_active_cells_pos3), data = data_table) + geom_quasirandom(dodge.width=.8)

# Check data
unique(data_table["genotype_2"])[[1]]

# Make subsets
data_subset_1 = subset(data_table, select=c("recording", 
                                       "genotype_2", 
                                       "n_percentage_active_cells_pos2"))
data_subset_2 = subset(data_table, select=c(
                                            "recording", 
                                            "genotype_2", 
                                             "n_percentage_active_cells_pos3"))


################################

################################

# Find the correct underlying probability distribution

data_subset_A = data_subset_1

ggplot(data_subset_A, aes(n_percentage_active_cells_pos2, fill=genotype_2)) + geom_histogram()

#m1.norm <- glmmTMB(n_percentage_active_cells_pos2 ~ genotype_2 + (1|recording), 
#                data = data_subset_A)
#m1.lnorm <- glmmTMB(n_percentage_active_cells_pos2 ~ genotype_2  + (1|recording), 
#                 family = gaussian(link="log"), data = data_subset_A)
#m1.gamma <- glmmTMB(n_percentage_active_cells_pos2 ~ genotype_2  + (1|recording), 
#                  family = Gamma(link="log"), data = data_subset_A)
m1.tweedie <- glmmTMB(n_percentage_active_cells_pos2 ~ genotype_2 +  + (1|recording), 
                    family = tweedie(link="log"), data = data_subset_A)


# Select one of the distributions
#res <- simulateResiduals(m1.norm, plot = T)
#res <- simulateResiduals(m1.lnorm, plot = T)
#res <- simulateResiduals(m1.gamma, plot = T)
res <- simulateResiduals(m1.tweedie, plot = T)

################################

################################

# use tweedie distribution

#### Create models for pos2 ####

# reduced base model
model_pos2.base <- glmmTMB(n_percentage_active_cells_pos2 ~ 1 + (1|recording),family = tweedie(link="log"),
                           data = data_subset_1)

# full model
model_pos2 <- glmmTMB(n_percentage_active_cells_pos2 ~ genotype_2 + (1|recording),family = tweedie(link="log"),
                      data = data_subset_1)
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

##### Create models for pos3 ####

# reduced base model
model_pos3.base <- glmmTMB(n_percentage_active_cells_pos3 ~ 1 + (1|recording),family = tweedie(link="log"),
                           data = data_subset_2)

# full model
model_pos3 <- glmmTMB(n_percentage_active_cells_pos3 ~ genotype_2 + (1|recording),family = tweedie(link="log"),
                      data = data_subset_2)
summary(model_pos3)

simulationOutput <- simulateResiduals(model_pos3, 1000)
plot(simulationOutput)

plot(fitted(model_pos3), residuals(model_pos3), xlab = "Fitted Values", ylab = "Residuals")
abline(h = 0, lty = 2)
lines(mean(fitted(model_pos3)), mean(residuals(model_pos3)))

# compare models
anova_result_pos3 = anova(model_pos3.base, model_pos3)
print(anova_result_pos3)

################################

################################

################################
