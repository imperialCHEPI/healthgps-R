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
######################## Parameters ##########################################################

#File to use
file_name <- "Health_GPS_final_main_20231109.csv"


#Filter threshold 
rate_to_filter <- 0.01 

#Physical Activity Parameters: mean and pa 
pa_mean = 1.6
pa_sd = 0.06

#number of EPA quantiles to sample
num_epa_quantiles_to_sample <- 10000

#Weight Parameters 
num_weight_quantiles_to_sample = 10000
bounds = c(15.5,18.5,20,25,30,35,40,45)

#BMI Prevalence by range from NCD RisC (e.g. the following is India 2022)
prevalence_male = c(0.1247,0.1159,0.4763,0.2291,0.0439,0.0068,0.0029)
prevalence_female = c(0.1369,0.1144,0.4211,0.2296,0.0744,0.0177,0.0055)

#######################################################################################################################
######################## Data Loading and Variable Extraction ##########################################################

# Load the CSV data into a data frame
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
######################## Filtering ##########################################################

# Create 'subdata' dataframe
subdata <- data.frame(sex, age, age1, age2, age3, inc, sector, carb, fat, protein, sodium, energy, energy_original)

# Remove rows with non-numeric NaN or empty values
subdata <- subdata[complete.cases(subdata), ]

# Set upper and lower quantiles
lower_q <- rate_to_filter 
upper_q <-1-lower_q

# Filter 'subdata' based on conditions
df <- subdata %>%
  filter(
    age <= 100,  # Age less than 100
    carb > quantile(carb, lower_q) & carb < quantile(carb, upper_q),  # Carb within quantiles
    fat > quantile(fat, lower_q) & fat < quantile(fat, upper_q),  # Fat within quantiles
    protein > quantile(protein, lower_q) & protein < quantile(protein, upper_q),  # Protein within quantiles
    sodium > quantile(sodium, lower_q) & sodium < quantile(sodium, upper_q),  # Sodium within quantiles
    energy > quantile(energy, lower_q) & energy < quantile(energy, upper_q)  # Energy within quantiles
  )

############################################################################################
###############################Energy to Physical Activity Quantiles###########################

# Generate normally distributed numbers with mean = 0 and std = 1
number_samples = length(df$energy)

random_numbers <- rnorm(number_samples, mean = 0, sd = 1)

# Physical activity follows a log normal distribution with mean = pa_mean and sd = pa_sd
physicalactivity = pa_mean*exp(pa_sd*random_numbers-0.5*pa_sd*pa_sd)
plot(density(physicalactivity), col="red", main="Physical Activity Density", lwd=3)

df$physicalactivity = physicalactivity
df$energy_to_pa_ratio = df$energy/df$physicalactivity

result <- df %>%
  group_by(age, sex) %>%
  summarize(
    energy_to_pa_ratio_mean = mean(energy_to_pa_ratio)
  )

merged_df <- merge(df, result, by = c("sex", "age"), all.x = TRUE)

# Calculate normalised energy ratio by dividing energy_to_pa_ratio by energy_to_pa_ratio_mean
normalised_energy_to_pa_ratio <- with(merged_df, energy_to_pa_ratio / energy_to_pa_ratio_mean)

# Set a threshold (cap) at the 99th percentile of the normalised energy ratio
cap <- quantile(normalised_energy_to_pa_ratio , 0.99)

# Keep only values below the cap to filter extreme outliers
normalised_energy_to_pa_ratio  <- normalised_energy_to_pa_ratio [normalised_energy_to_pa_ratio  < cap]

#plot distribution of normalised energy to pa ratio
main_title <- "Distribution of Normalised Energy to Physical Activity Ratio"
plot(density(normalised_energy_to_pa_ratio ), lwd=3, col="red", main = main_title)

# Sample random elements from normalised_energy_ratio
sampled_elements <- sample(normalised_energy_to_pa_ratio , num_epa_quantiles_to_sample)

# Order the sampled elements
sampled_elements = sort(sampled_elements)

# Write the sampled elements to a CSV file named "energy_physicalactivity_quantiles.csv"
write.csv(sampled_elements, "energy_physicalactivity_quantiles.csv")

############################################################################################
###############################Weight Quantiles using NCD Risk###########################

# Initialize empty vectors to store sampled BMI values for male and female
sampled_bmi_male <- c()
sampled_bmi_female <- c()

# Get the number of lists from the 'prevalence' vector
num_lists <- length(prevalence_male)

# Use a loop to generate and append BMI values for both male and female
for (i in 1:num_lists) {
  # Generate and append BMI values for male
  current_list <- runif(prevalence_male[i] * num_weight_quantiles_to_sample, min = bounds[i], max = bounds[i + 1])
  sampled_bmi_male <- c(sampled_bmi_male, current_list)
  
  # Generate and append BMI values for female
  current_list <- runif(prevalence_female[i] * num_weight_quantiles_to_sample, min = bounds[i], max = bounds[i + 1])
  sampled_bmi_female <- c(sampled_bmi_female, current_list)
}

# Plot density of BMI distribution by sex
plot(density(sampled_bmi_male), col = "blue", lwd = 3, main = "BMI distribution by sex")
lines(density(sampled_bmi_female), col = "red", lwd = 3)

# Adding legend
legend("topright", legend = c("Male", "Female"), col = c("blue", "red"), lty = 1, lwd = 2)

# Normalize BMI values by dividing each value by the mean of the respective sex
normalised_bmi_male <- sampled_bmi_male / mean(sampled_bmi_male)
normalised_bmi_female <- sampled_bmi_female / mean(sampled_bmi_female)

# Write normalised BMI values to CSV files
write.csv(normalised_bmi_male, "weight_quantiles_NCDRisk_male.csv")
write.csv(normalised_bmi_female, "weight_quantiles_NCDRisk_Female.csv")

#####################################################################################
#####################################################################################