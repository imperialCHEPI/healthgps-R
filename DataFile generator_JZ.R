## Package setup

# Package install
package_list <- c("MASS", "dplyr", "DescTools")
new_packages <- package_list[!(package_list %in% 
    installed.packages()[, "Package"])]
if(length(new_packages)) install.packages(new_packages)

# Package load
lapply(package_list, require, character.only = TRUE)



setwd("C:/Users/jzhu5/OneDrive - Imperial College London/Project/RSTL_India/BMI data")
prevalence <- read.csv("prevalence_data.csv")
mean_sd <- read.csv("mean_sd.csv")

## This should read the 'incomplete' data file created by Ali's scripts
setwd("C:/Users/jzhu5/OneDrive - Imperial College London/Project/RSTL_India/Input data")
data <- read.csv("India.DataFile.JZ.csv")

#########################################
# Generate mean BMI for each age#gender
#########################################
# Weight Parameters 
num_weight_quantiles_to_sample = 10000
bounds = c(15.5,18.5,20,25,30,35,40,45)


for (sex in 0:1) {
  for (age in 18:100) {
    prevalence_values <-  prevalence %>%
      filter(prevalence$Sex==sex & prevalence$Age==age)
    prevalence_values <- as.vector(prevalence_values[,5:11])
    prevalence_values <- as.numeric(prevalence_values)
    # Initialize empty vectors to store sampled BMI values
    sampled_bmi <- c()
    # Get the number of lists from the 'prevalence' vector
    num_lists <- length(prevalence_values)
    # Use a loop to generate and append BMI values
    for (i in 1:num_lists) {
      # Generate and append BMI values
      current_list <- runif(prevalence_values[i] * num_weight_quantiles_to_sample, min = bounds[i], max = bounds[i + 1])
      sampled_bmi <- c(sampled_bmi, current_list)
    }
    # Fit data to log-normal distribution
    fit <- fitdistr(sampled_bmi, densfun = "lognormal")
    # Generate some example data (replace this with your own vector of data)
    sampled_data <- rlnorm(100000, meanlog = fit$estimate[1], sdlog = fit$estimate[2])
    # Filter data to keep only values between min and max
    sampled_data <- sampled_data[sampled_data >= min(bounds) & sampled_data <= max(bounds)]
    mean_bmi = exp(fit$estimate[1])
    sd_bmi = exp(fit$estimate[2])
    print(sd_bmi)
  }
}


########################################
# Generate random physical activity
########################################
pa_normal <- rnorm(86602, mean = 1.6, sd = 0.1)
data$PhysicalActivity <- pa_normal

summary(data$PhysicalActivity)
sd(data$PhysicalActivity)


##################################################
# Generate random height from normal distribution
#################################################

for (sex in 0:1) {
  for (age in 19:100) {
    subdata <-  filter(data, data$Gender==sex & data$Age==age)
    size <- nrow(subdata)
    mean_sd$height_sd[which(mean_sd$Age==age & mean_sd$Gender==sex)] <- mean_sd$height_se[which(mean_sd$Age==age & mean_sd$Gender==sex)] * sqrt(size)
    }
}

data <- data %>%
  left_join(mean_sd, by=c("Age","Gender"))


df <- data.frame()
for (sex in 0:1) {
  for (age in 0:100) {
    subdata <-  filter(data, data$Gender==sex & data$Age==age)
    size <- nrow(subdata)
    sample_height <- rnorm(size, mean = subdata$height_mean, sd = subdata$height_sd)
    subdata$Height <- sample_height
    df <- rbind(df,subdata)
  }
}
summary(df$Height)

write.csv(df,"df_pa_height.csv")

##########################################################
# Generate random BMI from lognormal distribution
# Option 1: Adjuste mean and sd to get obesity prevalence
###########################################################

dff_female <- data.frame()
for (age in 0:90) {
    subdata <-  filter(df, df$Gender==0 & df$Age==age)
    size <- nrow(subdata)
    means <- c(mean(subdata$EnergyIntake),mean(subdata$bmi_mean+1))
    sd_ei <- sd(subdata$EnergyIntake)
    sd_bmi <- mean(subdata$bmi_sd+3.9)
    cov_matrix <- matrix(c(sd_ei^2, 0.004*sd_ei*sd_bmi,
                           0.004*sd_ei*sd_bmi, sd_bmi^2),
                         ncol=2)
    sample_bmi <- mvrnorm(size, mu = means, Sigma = cov_matrix, empirical = TRUE)
    sample_bmi <- as.data.frame(sample_bmi)
    colnames(sample_bmi) <- c("V1","bmi")
    sample_bmi$bmi <- Winsorize(sample_bmi$bmi, probs = c(0.05,1),na.rm = TRUE)
    print(mean(sample_bmi$bmi>=30))
    subdata$order <- order(subdata$EnergyIntake)
    sample_bmi$order <- order(sample_bmi$bmi)
    merged <- subdata %>%
      left_join(sample_bmi, by="order")
    dff_female <- rbind(dff_female,merged)
}

summary(dff_female$bmi)
summary(dff_female$bmi>=30)


dff_male <- data.frame()
  for (age in 0:90) {
    subdata <-  filter(df, df$Gender==1 & df$Age==age)
    size <- nrow(subdata)
    means <- c(mean(subdata$EnergyIntake),mean(subdata$bmi_mean+1))
    sd_ei <- sd(subdata$EnergyIntake)
    sd_bmi <- mean(subdata$bmi_sd+3)
    cov_matrix <- matrix(c(sd_ei^2, 0.004*sd_ei*sd_bmi,
                           0.004*sd_ei*sd_bmi, sd_bmi^2),
                         ncol=2)
    sample_bmi <- mvrnorm(size, mu = means, Sigma = cov_matrix, empirical = TRUE)
    sample_bmi <- as.data.frame(sample_bmi)
    colnames(sample_bmi) <- c("V1","bmi")
    sample_bmi$bmi <- Winsorize(sample_bmi$bmi, probs = c(0.05,1), na.rm = TRUE)
    print(mean(sample_bmi$bmi>=30))
    subdata$order <- order(subdata$EnergyIntake)
    sample_bmi$order <- order(sample_bmi$bmi)
    merged <- subdata %>%
      left_join(sample_bmi, by="order")
    dff_male <- rbind(dff_male,merged)
  }

summary(dff_male$bmi)
summary(dff_male$bmi>=30)

dff_merged <- rbind(dff_female,dff_male)

moredata <-  filter(df, df$Age>=91)
moredata$order <- "NA"
moredata$V1 <- "NA"
moredata$bmi <- moredata$bmi_mean

dff_merged <- rbind(dff_merged,moredata)

dff_merged$Weight <- dff_merged$bmi * (dff_merged$Height/100)^2
summary(dff_merged$Weight)
write.csv(dff_merged,"India.DataFile.obesity.csv")
saveRDS(dff_merged,"India.DataFile.obesity")




##################################################
# Generate random BMI from lognormal distribution
# Option 2: Ignore the prevalence of obesity
#################################################
dff <- data.frame()
for (sex in 0:1) {
  for (age in 0:90) {
    subdata <-  filter(df, df$Gender==sex & df$Age==age)
    size <- nrow(subdata)
    means <- c(mean(subdata$EnergyIntake),mean(subdata$bmi_mean))
    sd_ei <- sd(subdata$EnergyIntake)
    sd_bmi <- mean(subdata$bmi_sd)
    cov_matrix <- matrix(c(sd_ei^2, 0.004*sd_ei*sd_bmi,
                           0.004*sd_ei*sd_bmi, sd_bmi^2),
                         ncol=2)
    sample_bmi <- mvrnorm(size, mu = means, Sigma = cov_matrix, empirical = TRUE)
    sample_bmi <- as.data.frame(sample_bmi)
    colnames(sample_bmi) <- c("V1","bmi")
    print(mean(sample_bmi$bmi>=30))
    subdata$order <- order(subdata$EnergyIntake)
    sample_bmi$order <- order(sample_bmi$bmi)
    merged <- subdata %>%
      left_join(sample_bmi, by="order")
    dff <- rbind(dff,merged)
  }
}
summary(dff$bmi)
summary(dff$bmi>=30)


moredata <-  filter(df, df$Age>=91)
moredata$order <- "NA"
moredata$V1 <- "NA"
moredata$bmi <- moredata$bmi_mean
dff <- rbind(dff,moredata)

dff$Weight <- dff$bmi * (dff$Height/100)^2
summary(dff$Weight)

write.csv(dff,"India.DataFile.noobesity.csv")
saveRDS(dff, "India.DataFile.noobesity")
