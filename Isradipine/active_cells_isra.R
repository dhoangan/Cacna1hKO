library(lme4)
library(ggplot2)
library(lmerTest)
library(MASS)
library(car)
library(MCMCglmm)
library(DHARMa)
library(glmmTMB)
library(ggbeeswarm)



data = read.csv("isra_active_cells_per_slice.csv", sep=",")

data_table <- subset(data, select=c("animal", 
                                    "recording", 
                                    "genotype_2", 
                                    "n_percentage_active_cells"))

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

#boxplot(n_percentage_active_cells ~ genotype_2, data_table)
#ggplot(aes(pos, n_percentage_active_cells, color=genotype_2), data = data_table) +  geom_point()
ggplot(aes(y=n_percentage_active_cells, x=genotype_2), data = data_table) + geom_boxplot()

ggplot(aes(genotype_2, n_percentage_active_cells), data = data_table) + geom_beeswarm(dodge.width = 1)
ggplot(aes(genotype_2, n_percentage_active_cells), data = data_table) + geom_quasirandom(dodge.width=.8)

# Check data
unique(data_table["genotype_2"])[[1]]

# Make subsets
data_subset_wt = subset(data_table, genotype_2 == "wt" | genotype_2 == "wt_ctrl")
data_subset_KO = subset(data_table, genotype_2 == "KO" | genotype_2 == "KO_ctrl")


################################

################################

# Find the correct underlying probability distribution

data_subset_A = data_subset_KO

ggplot(data_subset_A, aes(n_percentage_active_cells, fill=genotype_2)) + geom_histogram()

m1.norm <- glmmTMB(n_percentage_active_cells ~ genotype_2 + (1|recording), 
                data = data_subset_A)
m1.lnorm <- glmmTMB(n_percentage_active_cells ~ genotype_2  + (1|recording), 
                 family = gaussian(link="log"), data = data_subset_A)
#m1.gamma <- glmmTMB(n_percentage_active_cells ~ genotype_2  + (1|recording), 
 #                 family = Gamma(link="log"), data = data_subset_A)
m1.tweedie <- glmmTMB(n_percentage_active_cells ~ genotype_2 +  + (1|recording), 
                    family = tweedie(link="log"), data = data_subset_A)


# Select one of the distributions
res <- simulateResiduals(m1.norm, plot = T)
res <- simulateResiduals(m1.lnorm, plot = T)
#res <- simulateResiduals(m1.gamma, plot = T)
res <- simulateResiduals(m1.tweedie, plot = T)

################################

################################

# use normal distribution

#### Create models for wt vs wt_ctrl ####

# reduced base model
model_wt.base <- glmmTMB(n_percentage_active_cells ~ 1 + (1|recording),
                         data = data_subset_wt)

# full model
model_wt <- glmmTMB(n_percentage_active_cells ~ genotype_2 + (1|recording),
                    data = data_subset_wt)
summary(model_wt)

simulationOutput <- simulateResiduals(model_wt, 1000)
plot(simulationOutput)

plot(fitted(model_wt), residuals(model_wt), xlab = "Fitted Values", ylab = "Residuals")
abline(h = 0, lty = 2)
lines(mean(fitted(model_wt)), mean(residuals(model_wt)))

# compare models
anova_result_wt = anova(model_wt.base, model_wt)
print(anova_result_wt)

##############################

##############################
#### Create models for KO vs KO_ctrl ####

# reduced base model
model_KO.base <- glmmTMB(n_percentage_active_cells ~ 1 ,
                         data = data_subset_KO)

# full model
model_KO <- glmmTMB(n_percentage_active_cells ~ genotype_2 ,
                    data = data_subset_KO)
summary(model_KO)

simulationOutput <- simulateResiduals(model_KO, 1000)
plot(simulationOutput)

plot(fitted(model_KO), residuals(model_KO), xlab = "Fitted Values", ylab = "Residuals")
abline(h = 0, lty = 2)
lines(mean(fitted(model_KO)), mean(residuals(model_KO)))

# compare models
anova_result_KO = anova(model_KO.base, model_KO)
print(anova_result_KO)

################################

################################

################################
