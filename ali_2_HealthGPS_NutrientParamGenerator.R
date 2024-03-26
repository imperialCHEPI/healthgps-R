#######################################################################################################################
######################## Libraries and Packages Initialisation##########################################################

library('stargazer')
library(ggplot2)
library(dplyr) 
library(nnet)
library(zoo)
library(moments)
library(MASS)
library(tidyr)
library(corrplot)

#######################################################################################################################
######################## Data Loading and Variable Extraction ##########################################################

# Load the CSV data into a data frame
file_name <- "Health_GPS_ind_20231012.csv"
data <- read.csv(file_name)

# Attach the data frame for easier variable referencing
attach(data)

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
subdata <- data.frame(sex, age, age1, age2, age3, inc, sector, carb, fat, protein, sodium, energy, energy_original)

# Remove rows with non-numeric NaN or empty values
subdata <- subdata[complete.cases(subdata), ]

# Set lower and upper quantiles
upper_q <- 0.99
lower_q <- 1-upper_q 

# Filter 'subdata' based on conditions
df <- subdata %>%
  filter(
    age < 100,  # Age less than 100
    carb > quantile(carb, lower_q) & carb < quantile(carb, upper_q),  # Carb within quantiles
    fat > quantile(fat, lower_q) & fat < quantile(fat, upper_q),  # Fat within quantiles
    protein > quantile(protein, lower_q) & protein < quantile(protein, upper_q),  # Protein within quantiles
    sodium > quantile(sodium, lower_q) & sodium < quantile(sodium, upper_q),  # Sodium within quantiles
    energy > quantile(energy, lower_q) & energy < quantile(energy, upper_q)  # Energy within quantiles
  )

############################################################################################
########################Summary Statistics by Age and Sex#######################################

# Calculate mean values for 'carb', 'fat', 'protein', 'sodium', and 'energy' based on 'age' and 'sex'
result <- df %>%
  group_by(age, sex) %>%
  summarize(
    carb_mean = mean(carb),
    fat_mean = mean(fat),
    protein_mean = mean(protein),
    sodium_mean = mean(sodium),
    energy_mean = mean(energy)
  )

# Merge the calculated mean values back to the original dataframe ('df')
merged_df <- merge(df, result, by = c("sex", "age"), all.x = TRUE)

############################################################################################
########################BOX-COX Transformation: Carb#######################################

# Calculate Box-Cox transformation for 'carb' variable
x <- merged_df$carb / merged_df$carb_mean
boxcox_results <- boxcox(lm(x ~ 1))
lambda_carb <- boxcox_results$x[which.max(boxcox_results$y)]

# Apply Box-Cox transformation to 'carb' and create a new variable 'new_x_exact'
merged_df$new_x_exact <- (x^lambda_carb - 1) / lambda_carb

# Fit a linear regression model for the transformed variable
reg_carb <- lm(new_x_exact ~ sex + age1 + age2 + age3 + sector + inc, data = merged_df)

# Display summary statistics of the regression model
summary(reg_carb)

# Calculate the standard deviation of the residuals
sd_carb <- sd(reg_carb$residuals)

# Plot the density distribution of 'carb' variable
plot(density(carb), col="red", lwd=3, main="Carbs: Original Distribution")

# Plot the density distribution of transformed 'carb' variable
plot(density(merged_df$new_x_exact), col="red", lwd=3, main="Carbs: Transformed Distribution")

# Plot the density distribution of standardized residuals for 'carb' regression
plot(density(scale(reg_carb$residuals)), col="red", main="Carb Residuals", lwd=3)

# Standardize 'carb' residuals
carb_transformed <- scale(reg_carb$residuals)

# Adjust the number of random values for comparison
num_values <- 1000000  
random_values <- rnorm(num_values)

# Plot the density distribution of standardized 'carb' residuals and random values from a normal distribution
main_title <- "Carbs: Transformed Distribution vs. Random Normal Distribution"
plot(density(carb_transformed), col="red", lwd=3, main = main_title)
lines(density(random_values), col="blue", lwd=3)

############################################################################################
########################BOX-COX Transformation: Fat#########################################

# Calculate Box-Cox transformation for 'fat' variable
x <- merged_df$fat / merged_df$fat_mean
boxcox_results <- boxcox(lm(x ~ 1))
lambda_fat <- boxcox_results$x[which.max(boxcox_results$y)]

# Apply Box-Cox transformation to 'fat' and create a new variable 'new_x_exact'
merged_df$new_x_exact <- (x^lambda_fat - 1) / lambda_fat

# Plot the density of the transformed variable
plot(density(merged_df$new_x_exact))

# Fit a linear regression model for the transformed variable
reg_fat <- lm(new_x_exact ~ sex + age1 + age2 + age3 + sector + inc, data = merged_df)

# Display summary statistics of the regression model
summary(reg_fat)

# Calculate the standard deviation of the residuals
sd_fat <- sd(reg_fat$residuals)

# Plot the density distribution of 'fat' variable
plot(density(fat), col="red", lwd=3, main="Fat: Original Distribution")

# Plot the density distribution of transformed 'fat' variable
plot(density(merged_df$new_x_exact), col="red", lwd=3, main="Fat: Transformed Distribution")

# Plot the density distribution of standardized residuals for 'fat' regression
plot(density(scale(reg_fat$residuals)), col="red", main="Fat Residuals",lwd=3)

# Standardize 'fat' residuals
fat_transformed <- scale(reg_fat$residuals)

# Adjust the number of random values for comparison
num_values <- 1000000  
random_values <- rnorm(num_values)

# Plot the density distribution of standardized 'fat' residuals and random values from a normal distribution
main_title <- "Fat: Transformed Distribution vs. Random Normal Distribution"
plot(density(fat_transformed), col="red", lwd=3, main = main_title)
lines(density(random_values), col="blue", lwd=3)

############################################################################################
########################BOX-COX Transformation: Protein####################################

# Calculate Box-Cox transformation for 'protein' variable
x <- merged_df$protein / merged_df$protein_mean
boxcox_results <- boxcox(lm(x ~ 1))
lambda_protein <- boxcox_results$x[which.max(boxcox_results$y)]

# Apply Box-Cox transformation to 'protein' and create a new variable 'new_x_exact'
merged_df$new_x_exact <- (x^lambda_protein - 1) / lambda_protein

# Fit a linear regression model for the transformed variable
reg_protein <- lm(new_x_exact ~ sex + age1 + age2 + age3 + sector + inc, data = merged_df)

# Display summary statistics of the regression model
summary(reg_protein)

# Calculate the standard deviation of the residuals
sd_protein <- sd(reg_protein$residuals)

# Plot the density distribution of 'protein' variable
plot(density(protein), col="red", lwd=3, main="Protein: Original Distribution")

# Plot the density distribution of transformed 'protein' variable
plot(density(merged_df$new_x_exact), col="red", lwd=3, main="Protein: Transformed Distribution")

# Plot the density distribution of standardized residuals for 'protein' regression
plot(density(scale(reg_protein$residuals)), col="red", main="Protein Residuals",lwd=3)

# Standardize 'protein' residuals
protein_transformed <- scale(reg_protein$residuals)

# Adjust the number of random values for comparison
num_values <- 1000000  
random_values <- rnorm(num_values)

# Plot the density distribution of standardized 'protein' residuals and random values from a normal distribution
main_title <- "Protein: Transformed Distribution vs. Random Normal Distribution"
plot(density(protein_transformed), col="red", lwd=3, main = main_title)
lines(density(random_values), col="blue", lwd=3)

############################################################################################
########################BOX-COX Transformation: Sodium#######################################################

# Calculate Box-Cox transformation for 'sodium' variable
x <- merged_df$sodium / merged_df$sodium_mean
boxcox_results <- boxcox(lm(x ~ 1))
lambda_sodium <- boxcox_results$x[which.max(boxcox_results$y)]

# Apply Box-Cox transformation to 'sodium' and create a new variable 'new_x_exact'
merged_df$new_x_exact <- (x^lambda_sodium - 1) / lambda_sodium

# Fit a linear regression model for the transformed variable
reg_sodium <- lm(new_x_exact ~ sex + age1 + age2 + age3 + sector + inc, data = merged_df)

# Display summary statistics of the regression model
summary(reg_sodium)

# Calculate the standard deviation of the residuals
sd_sodium <- sd(reg_sodium$residuals)

# Plot the density distribution of 'sodium' variable
plot(density(sodium), col="red", lwd=3, main="Sodium: Original Distribution")

# Plot the density distribution of transformed 'sodium' variable
plot(density(merged_df$new_x_exact), col="red", lwd=3, main="Sodium: Transformed Distribution")

# Plot the density distribution of standardized residuals for 'sodium' regression
plot(density(scale(reg_sodium$residuals)), col="red", main="Sodium Residuals",lwd=3)

# Standardize 'sodium' residuals
sodium_transformed <- scale(reg_sodium$residuals)

# Adjust the number of random values for comparison
num_values <- 1000000  
random_values <- rnorm(num_values)

# Plot the density distribution of standardized 'sodium' residuals and random values from a normal distribution
main_title <- "Sodium: Transformed Distribution vs. Random Normal Distribution"
plot(density(sodium_transformed), col="red", lwd=3, main = main_title)
lines(density(random_values), col="blue", lwd=3)

############################################################################################
########################Box-Cox Transformation Parameters#######################################################

lambda = c(lambda_carb,lambda_fat,lambda_protein,lambda_sodium)
sd = c(sd_carb,sd_fat,sd_protein,sd_sodium)

boxcoxparameters = data.frame(lambda,sd)
rownames(boxcoxparameters) = c("carb","fat","protein","sodium")

#########################################################################################
##############################Nutrient Regression Residuals & Correlations################################

# Create a dataframe 'nutrients_residuals' containing residuals from regression models for nutrients
nutrients_residuals <- data.frame(
  carb_residuals = reg_carb$residuals,
  fat_residuals = reg_fat$residuals,
  protein_residuals = reg_protein$residuals,
  sodium_residuals = reg_sodium$residuals
)

# Calculate the correlation matrix for residuals
correlation_matrix <- cor(nutrients_residuals)

# Plot the correlation matrix
rownames(correlation_matrix) <- c('carb','fat','protein','sodium')
colnames(correlation_matrix) <- c('carb','fat','protein','sodium')
main_title <- "Residuals: Correlation Plot"
corrplot(correlation_matrix, type="upper", title=main_title)

######################################################################################
#######################################Output Files###################################

# Save the coefficients from regression models for transformed variables ('carb', 'fat', 'protein', 'sodium') using Box-Cox transformation
write.csv(reg_carb$coefficients, "boxcox_carb_coefficients.csv")
write.csv(reg_fat$coefficients, "boxcox_fat_coefficients.csv")
write.csv(reg_protein$coefficients, "boxcox_protein_coefficients.csv")
write.csv(reg_sodium$coefficients, "boxcox_sodium_coefficients.csv")

# Save Box-Cox transformation parameters to a CSV file
write.csv(boxcoxparameters, "boxcox_parameters.csv")

# Save correlation matrix to a CSV file
write.csv(cor, "boxcox_correlationmatrix.csv")

#####################################################################################
#####################################################################################