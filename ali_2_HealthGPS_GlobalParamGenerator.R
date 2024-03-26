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

#File to use
file_name <- "data/maxime-180324/Health_GPS_sc_effect_main.csv"

# (e.g., 0.5 for half the rows)
country = "India"
rate_to_print <- 0.2

#Filter threshold 
rate_to_filter <- 0.01 

#Smoothing parameters
smoothing = TRUE
repeated_smoothing = 150

################################################################################
######################## Data Loading and Variable Extraction ##################

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

################################################################################
########################Filtering###############################################

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
    carb > quantile(carb, lower_q) & carb < quantile(carb, upper_q),              # Carb within quantiles
    fat > quantile(fat, lower_q) & fat < quantile(fat, upper_q),                  # Fat within quantiles
    protein > quantile(protein, lower_q) & protein < quantile(protein, upper_q),  # Protein within quantiles
    sodium > quantile(sodium, lower_q) & sodium < quantile(sodium, upper_q),      # Sodium within quantiles
    energy > quantile(energy, lower_q) & energy < quantile(energy, upper_q)       # Energy within quantiles
  )

################################################################################
########################ToolBox#################################################

get_smooth_factor <- function(initial_factor, times) {
  if(times<1)
    return(initial_factor)
  res <- initial_factor
  for (j in 1:times) {
    tmp <- res
    for (p in 1:length(res)) {
      tmp[p] <- res[p]
    }
    sum <- 0
    for (i in 1:length(res)) {
      if (i == 1) {
        res[i] <- (2 * tmp[i] + tmp[i + 1]) / 3
      } else if (i == length(tmp)) {
        res[i] <- (tmp[i - 1] + 2 * tmp[i]) / 3
      } else {
        res[i] <- (tmp[i - 1] + tmp[i] + tmp[i + 1]) / 3
      }   
      sum <- sum + res[i]
    }
  }  
  return(res)
}

################################################################################
########################Calculate Means and Create CSV Files####################

# Calculate mean values for 'carb', 'fat', 'protein', 'sodium', and 'energy' 
# based on 'age' and 'sex'

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
result$sex = ifelse(result$sex=="1","Male","Female")

# Reorder the rows based on the 'sex' column
result <- result[order(result$sex), ]

# Remove the "_mean" suffix from column names
colnames(result) <- sub("_mean$", "", colnames(result))

# Change the "energy" to "energyintake"
colnames(result) <- sub("energy", "energyintake", colnames(result))

# Split the data frame based on 'sex'
male_data <- result[result$sex == "Male", ]
female_data <- result[result$sex == "Female", ]

# List of columns to smooth
columns_to_smooth <- c("carb", "fat", "protein", "sodium","energyintake")

# Apply the smoothing function to selected columns
male_data[, columns_to_smooth] <- lapply(male_data[, columns_to_smooth], get_smooth_factor, repeated_smoothing)
female_data[, columns_to_smooth] <- lapply(female_data[, columns_to_smooth], get_smooth_factor, repeated_smoothing)

# Variables to plot
variables_to_plot <- c("carb", "fat", "protein", "sodium", "energyintake")

# Loop through each variable and create a separate plot
for (variable in variables_to_plot) {
  # Plotting male_data
  plot(male_data$age, male_data[[variable]], 
    type="l", 
    lwd=3, 
    col="blue", 
    main=paste("Average ", variable, " by sex and age"), 
    xlab="Age", 
    ylab=variable,
    ylim = c(1000,2500))
  
  # Plotting female_data
  lines(female_data$age, female_data[[variable]], col="red", lwd=3)
  
  # Adding legend
  legend("topright", legend=c("Male", "Female"), col=c("blue", "red"), lwd=3)
}

# Write CSV files for Male and Female
write.csv(male_data, "out/Male.FactorMeans.csv", row.names = FALSE)
write.csv(female_data, "out/Female.FactorMeans.csv", row.names = FALSE)

################################################################################
######################################Rural Sector: Prevalence##################

# Subset the dataframe into two age groups (under 18 and above 18)
df_under_18 <- df[df$age < 18, ]
df_above_18 <- df[df$age >= 18, ]

# Calculate the sector prevalence for males and females in each age group
prevalence_under_18 <- aggregate(sector ~ sex, data = df_under_18, FUN = function(x) sum(x == 1) / length(x))
prevalence_above_18 <- aggregate(sector ~ sex, data = df_above_18, FUN = function(x) sum(x == 1) / length(x))

# Create bar plots for each age group and gender
par(mfrow = c(2, 1))  # Create a 2x1 grid for plots

# Bar plot for individuals under 18
barplot(prevalence_under_18$sector, 
        beside = TRUE,
        main = "Rural Prevalence for Individuals Under 18",
        xlab = "Gender",
        ylab = "Prevalence",
        col = c("pink", "blue"),
        names.arg = c("Female", "Male"))

# Bar plot for individuals above 18
barplot(prevalence_above_18$sector, 
        beside = TRUE,
        main = "Rural Prevalence for Individuals Above 18",
        xlab = "Gender",
        ylab = "Prevalence",
        col = c("pink", "blue"),
        names.arg = c("Female", "Male"))

# Reset the plot layout
par(mfrow = c(1, 1))

# Combine the results into a single data frame
prevalence_under_18$age = c("Under 18", "Under 18")
prevalence_above_18$age = c("Above 18", "Above 18")
result <- rbind(prevalence_under_18, prevalence_above_18)
result$sex = ifelse(result$sex=="1","Male","Female")

# Reorder columns (sex, age, sector)
result <- result[, c("sex", "age", "sector")]

# Write the result to a CSV file
write.csv(result, "out/rural_prevalence.csv", row.names = FALSE)

################################################################################
####################################Income######################################

# Calculate the prevalence by sex and income group
prevalence_table <- table(df$sex, df$inc) / rowSums(table(df$sex, df$inc))

# Convert the table to a data frame
prevalence_df <- as.data.frame.matrix(prevalence_table)

# Create a bar plot
barplot(as.matrix(prevalence_df),
        beside = TRUE,
        main = "Prevalence of Income Groups by Sex",
        xlab = "Sex (0 = Female, 1 = Male)",
        ylab = "Prevalence",
        col = c("pink", "blue"),
        legend.text = TRUE,
        names.arg = c("Tercile 1", "Tercile 2", "Tercile 3"))


# Subset the data for cases where sector is equal to 1
data_subset <- df[df$sector == 0, ]

# Calculate the prevalence by sex and income group within sector 1
prevalence_table <- table(data_subset$sex, data_subset$inc) / rowSums(table(data_subset$sex, data_subset$inc))

# Convert the table to a data frame
prevalence_df <- as.data.frame.matrix(prevalence_table)

# Create a bar plot
barplot(as.matrix(prevalence_df),
        beside = TRUE,
        main = "Prevalence of Income Groups by Sex (Urban)",
        xlab = "Sex (0 = Female, 1 = Male)",
        ylab = "Prevalence",
        col = c("pink", "blue"),
        legend.text = TRUE,
        names.arg = c("Tercile 1", "Tercile 2", "Tercile 3"))


df$over18 <- ifelse(df$age < 18, 0, 1)

# Subset the data for cases where sector is equal to 1
data_subset <- df

# Calculate the prevalence by sex and income group within sector 1
prevalence_table <- table(data_subset$over18, data_subset$inc) / rowSums(table(data_subset$over18, data_subset$inc))

# Convert the table to a data frame
prevalence_df <- as.data.frame.matrix(prevalence_table)

# Create a bar plot
barplot(as.matrix(prevalence_df),
        beside = TRUE,
        main = "Prevalence of Income Groups by Age Groups",
        xlab = "AgeGroup (0 = Under 18, 1 = Over 18)",
        ylab = "Prevalence",
        col = c("green", "red"),
        legend.text = TRUE,
        names.arg = c("Tercile 1", "Tercile 2", "Tercile 3"))

# Fit a multinomial logistic regression model
model <- multinom(inc ~ sex + over18 + sector, data = df)
summary(model)

# 26/03/24 – 17:09
# This section appears to throw an error
# Error in model$fitted : 
#   (converted from warning) partial match of 'fitted' to 'fitted.values'

model

coefficients <- summary(model)$coefficients

# Create a new row with all values as 0
new_row <- c(0, 0, 0, 0)

# Add the new row to your existing data frame
coefficients <- rbind(coefficients, new_row)
rownames(coefficients)[rownames(coefficients) == "new_row"] <- "1"
coefficients <- coefficients[order(rownames(coefficients)), ]

write.csv(coefficients,"income_model.csv")

# Generate predicted probabilities for each class
predicted_probs <- predict(model, type = "probs")

# Sample from the predicted probabilities
n_samples <- 1  # Replace with the desired number of samples
samples <- apply(predicted_probs, 1, function(probs) sample(1:length(probs), size = n_samples, prob = probs))

df$samples = samples

data_subset <- df[df$sex == 0, ]

# Calculate the prevalence of each category for 'samples' and 'inc'
prevalence_samples <- table(data_subset$samples) / length(data_subset$samples)
prevalence_inc <- table(data_subset$inc) / length(data_subset$inc)

# Convert prevalence values to numeric
prevalence_samples <- as.numeric(prevalence_samples)
prevalence_inc <- as.numeric(prevalence_inc)

# Combine the prevalence data into a single dataframe
prevalence_data <- data.frame(Category = c("Tercile 1", "Tercile 2", "Tercile 3"),
                              Samples = prevalence_samples,
                              Income = prevalence_inc)

# Create a bar plot
barplot(
  t(as.matrix(prevalence_data[, 2:3])),
  beside = TRUE,
  main = "Female Income: Prevalence Comparison of 'model' and 'data'",
  xlab = "Tercile",
  ylab = "Prevalence",
  col = c("pink", "blue"),
  legend.text = TRUE,
  names.arg = prevalence_data$Category
)

######################################################################################
####################################Data File Creation###########################################

output_age = df$age
output_sex = ifelse(df$sex=="1","Male","Female")
output_income  = df$inc
output_sector  = ifelse(df$sector=="0","Rural","Urban")
output_carb = df$carb
output_fat = df$fat
output_protein = df$protein
output_sodium = df$sodium
output_energyintake = df$energy

output_df = data.frame(output_sex, output_age,output_sector,output_income,output_carb,output_fat,output_protein,output_sodium,output_energyintake)
colnames(output_df) = c("sex","age","sector","income","carb","fat","protein","sodium","energyintake")

# Calculate the number of rows to print
rows_to_print <- round(nrow(df) * (rate_to_print))

# Sample the rows to print
sampled_rows <- sample(1:nrow(df), size = rows_to_print)

# Create output dataframe with sampled rows
sampled_df <- output_df[sampled_rows, ]

# Order the dataframe by age and sex
sampled_df <- sampled_df[order(sampled_df$age, sampled_df$sex), ]

# Set row names from 1 to the last row
rownames(sampled_df) <- 1:nrow(sampled_df)

# Write the output dataframe to a CSV file
write.csv(sampled_df, file = paste0(country, ".DataFile.Incomplete.csv"))

#####################################################################################
#####################################################################################