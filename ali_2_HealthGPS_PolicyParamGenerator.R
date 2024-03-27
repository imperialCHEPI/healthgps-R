################################################################################
######################## Libraries and Packages Initialisation##################

## Package setup

# Package install
package_list <- c("stargazer", "ggplot2", "dplyr", "nnet", "zoo", "moments", 
    "MASS", "tidyr", "corrplot")
new_packages <- package_list[!(package_list %in% 
    installed.packages()[, "Package"])]
if(length(new_packages)) install.packages(new_packages)

# Package load
lapply(package_list, require, character.only = TRUE)

################################################################################
######################## Parameters ############################################

# File to use
# Not at all clear which file I should be pointing to here.
# needs to be something with 'ind_ic_pr_sodium_mg_'
file_name <- "data/maxime-180324/Health_GPS_ind_baseline.csv"

# Filter threshold 
filter <- 0.005 

# Scenario identifier
scenario_number="1"

#######################################################################################################################
######################## Data Loading and Variable Extraction ##########################################################

# Load the CSV data into a data frame
data <- read.csv(file_name)

# Attach the data frame for easier variable referencing
attach(data)

# Construct variable names based on the scenario
scenario = paste("sc",scenario_number,sep="")

carb_variable = paste0("ind_ic_pr_carb_g_", scenario)
fat_variable = paste0("ind_ic_pr_tot_fat_g_", scenario)
protein_variable = paste0("ind_ic_pr_protein_g_", scenario)
sodium_variable = paste0("ind_ic_pr_sodium_mg_", scenario)

# Access the variables in your dataframe
carb_effect = 100 * data[[carb_variable]]
fat_effect = 100 * data[[fat_variable]]
protein_effect = 100 * data[[protein_variable]]
sodium_effect = 100 * data[[sodium_variable]]

# Extracting variables from the 'data' data frame
sex <- data$ind_gender
age <- data$ind_age
age1 <- data$ind_age
age2 <- data$ind_age * data$ind_age
age3 <- data$ind_age * data$ind_age * data$ind_age
inc <- data$hh_tercile
protein <- data$ind_protein_g
carb <- data$ind_carb_g
fat <- data$ind_tot_fat_g
sodium <- data$ind_sodium_mg
sector <- data$hh_sector
energy_original <- data$ind_energy_kcal

# Calculate 'energy' based on the standard formula
energy <- 4 * carb + 9 * fat + 4 * protein + 0 * sodium

############################################################################################
########################Filtering##########################################################

# Create 'subdata' dataframe
subdata <- data.frame(sex, age, age1, age2, age3, inc, sector, carb, fat, protein, sodium, energy, energy_original,carb_effect,fat_effect,protein_effect,sodium_effect)

# Remove rows with non-numeric NaN or empty values
subdata <- subdata[complete.cases(subdata), ]

# Set upper and lower quantiles
upper_q <-1-filter
lower_q <- filter 

# Filter 'subdata' based on conditions
df <- subdata %>%
  filter(
    age <= 100,  # Age less than 100
    carb > quantile(carb, lower_q) & carb < quantile(carb, upper_q),  # Carb within quantiles
    fat > quantile(fat, lower_q) & fat < quantile(fat, upper_q),  # Fat within quantiles
    protein > quantile(protein, lower_q) & protein < quantile(protein, upper_q),  # Protein within quantiles
    sodium > quantile(sodium, lower_q) & sodium < quantile(sodium, upper_q),  # Sodium within quantiles
    energy > quantile(energy, lower_q) & energy < quantile(energy, upper_q),  # Energy within quantiles
    carb_effect > quantile(carb_effect, lower_q) & carb_effect < quantile(carb_effect, upper_q),  # Carb effect within quantiles
    fat_effect > quantile(fat_effect, lower_q) & fat_effect < quantile(fat_effect, upper_q),  # Fat effect within quantiles
    protein_effect > quantile(protein_effect, lower_q) & protein_effect < quantile(protein_effect, upper_q),  # Protein effect within quantiles
    sodium_effect > quantile(sodium_effect, lower_q) & sodium_effect < quantile(sodium_effect, upper_q)  # Sodium effect within quantiles
  )

############################################################################################
######################## Data Visualization and Analysis ##########################################################

# Plot the density of carb and carb_effect with different colors for comparison
plot(density(df$carb), col="blue", lwd=3, main = "Carb Density")
plot(density(df$carb_effect), col="red", lwd=3, main = "Carb Effect Density")

# Plot the density of fat and fat_effect with different colors for comparison
plot(density(df$fat), col="blue", lwd=3, main = "Fat Density")
plot(density(df$fat_effect), col="red", lwd=3, main = "Fat Effect Density")

# Plot the density of protein and protein_effect with different colors for comparison
plot(density(df$protein), col="blue", lwd=3, main = "Protein Density")
plot(density(df$protein_effect), col="red", lwd=3, main = "Protein Effect Density")

# Plot the density of sodium and sodium_effect with different colors for comparison
plot(density(df$sodium), col="blue", lwd=3, main = "Sodium Density")
plot(density(df$sodium_effect), col="red", lwd=3, main = "Sodium Effect Density")

# Calculate the percentage of remaining data after filtering
percentage = dim(df)[1] / dim(data)[1]
percentage 

# Create dataframes for nutrients and effects
nutrients = data.frame(df$carb, df$fat, df$protein, df$sodium)
effects = data.frame(df$carb_effect, df$fat_effect, df$protein_effect, df$sodium_effect)

# Display summary statistics for nutrients and effects
summary(nutrients)
summary(effects)

# Create the heatmap for nutrient correlations
nutrients_cor = cor(nutrients)
effects_cor = cor(effects)

# Set row and column names for better interpretation in the heatmap
rownames(nutrients_cor) <- c('carb', 'fat', 'protein', 'sodium')
colnames(nutrients_cor) <- c('carb', 'fat', 'protein', 'sodium')

rownames(effects_cor) <- c('carb_effect', 'fat_effect', 'protein_effect', 'sodium_effect')
colnames(effects_cor) <- c('carb_effect', 'fat_effect', 'protein_effect', 'sodium_effect')

# Plot nutrient correlation heatmap
main_title <- "Nutrients: Correlation Plot"
corrplot(nutrients_cor, type="upper", title=main_title)

# Plot nutrient effects correlation heatmap
main_title <- "Nutrient Effects: Correlation Plot"
corrplot(effects_cor, type="upper", title=main_title)

############################################################################################
######################## Policy Effect: Modelling ###################################

# Transforming nutrient variables with logarithm
log_carb = log(df$carb)
log_fat = log(df$fat)
log_protein = log(df$protein)
log_sodium = log(df$sodium)

# Regression Models for Policy Effects
# Model for carb_effect
reg_carb = lm(carb_effect ~ sex + age + sector + inc + log_carb + log_fat + log_protein + log_sodium, data = df)
summary(reg_carb)

# Model for fat_effect
reg_fat = lm(fat_effect ~ sex + age + sector + inc + log_carb + log_fat + log_protein + log_sodium, data = df)
summary(reg_fat)

# Model for protein_effect
reg_protein = lm(protein_effect ~ sex + age + sector + inc + log_carb + log_fat + log_protein + log_sodium, data = df)
summary(reg_protein)

# Model for sodium_effect
reg_sodium = lm(sodium_effect ~ sex + age + sector + inc + log_carb + log_fat + log_protein + log_sodium, data = df)
summary(reg_sodium)

# Residuals Analysis
residuals_matrix = data.frame(reg_carb$residuals, reg_fat$residuals, reg_protein$residuals, reg_sodium$residuals)
names(residuals_matrix) = c("carb", "fat", "protein", "sodium")

residuals_correlation = round(cor(residuals_matrix), 2)
residuals_covariance = round(cov(residuals_matrix), 2)

# Write model coefficients and residuals covariance to CSV files
# Construct variable names based on the scenario
scenario = paste("scenario",scenario_number,sep="")
write.csv(reg_carb$coefficients, paste(scenario, sep = "_", "carb_policyeffect_model.csv"))
write.csv(reg_fat$coefficients, paste(scenario, sep = "_", "fat_policyeffect_model.csv"))
write.csv(reg_protein$coefficients, paste(scenario, sep = "_", "protein_policyeffect_model.csv"))
write.csv(reg_sodium$coefficients, paste(scenario, sep = "_", "sodium_policyeffect_model.csv"))
write.csv(residuals_covariance, paste(scenario, sep = "_", "residuals_policyeffect_covariance.csv"))

# Calculate min and max values for nutrients and policy effects
min_nutrients = c(min(df$carb), min(df$fat), min(df$protein), min(df$sodium))
max_nutrients = c(max(df$carb), max(df$fat), max(df$protein), max(df$sodium))

min_policy_effect = c(min(df$carb_effect), min(df$fat_effect), min(df$protein_effect), min(df$sodium_effect))
max_policy_effect = c(max(df$carb_effect), max(df$fat_effect), max(df$protein_effect), max(df$sodium_effect))

# Create dataframes for boundaries of effects and nutrients
boundaries_effect = data.frame(min_policy_effect, max_policy_effect)
rownames(boundaries_effect) = c("carb", "fat", "protein", "sodium")
write.csv(boundaries_effect, paste(scenario, sep = "_", "boundaries_effect.csv"))

boundaries_nutrients = data.frame(min_nutrients, max_nutrients)
rownames(boundaries_nutrients) = c("carb", "fat", "protein", "sodium")
write.csv(boundaries_nutrients, paste(scenario, sep = "_", "boundaries_nutrients.csv"))

#####################################################################################
#####################################################################################
