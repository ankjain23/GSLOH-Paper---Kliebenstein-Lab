library(tidyverse)
library(dplyr)
library(tidyr)
library(MASS)
library(lme4)
library(car)
library(emmeans)
library(iheatmapr)
library(readxl)
library(ggplot2)
#set working directory to folder with xl sheet
setwd("~/Documents/URC Data Tables/GSLOH_Enantiomers")
dataset <- read_excel("enantiomers_gsloh2.xlsx")

# Filter the dataset for the specified values of X
filtered_data <- dataset %>% filter(X %in% c(16, 4, 8, 24, 25, 26, 22, 28))

# create VIOLIN PLOT
ggplot(filtered_data, aes(x = factor(Event), y = conversion.rate)) +
  geom_violin(trim = FALSE, fill = "lightblue", alpha = 0.7) +  # Violin plot
  geom_jitter(width = 0.2, alpha = 0.5) +  # Adds individual data points
  theme_minimal() + 
  labs(x = "Event", y = "Conversion Rate", title = "Violin Plot of Conversion Rate by Event") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Rotate x-axis labels for readability

# create BOX PLOT for filtered data with onyl certain events
ggplot(filtered_data, aes(x = factor(Event), y = conversion.rate)) +
  geom_boxplot(fill = "lightblue", alpha = 0.7, outlier.color = "red") +  # Boxplot with colored outliers
  geom_jitter(width = 0.2, alpha = 0.5, color = "black") +  # Adds individual data points for visibility
  theme_minimal() + 
  labs(x = "Event", y = "Conversion Rate", title = "Boxplot of Conversion Rate by Event") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for readability
# create BOX PLOT for ALL data with X = X and Y = conversion.rate
ggplot(dataset, aes(x = factor(X), y = conversion.rate)) +
  geom_boxplot(fill = "lightblue", alpha = 0.7, outlier.color = "red") +  # Boxplot
  theme_minimal() + 
  labs(x = "X", y = "Conversion Rate", title = "Boxplot of Conversion Rate by X") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
## color coded by preference (>0.5)
ggplot(dataset, aes(x = factor(X), y = conversion.rate, fill = preference > 0.5)) +
  geom_boxplot(alpha = 0.7, outlier.color = "red") +  # Boxplot
  theme_minimal() + 
  scale_fill_manual(values = c("TRUE" = "blue", "FALSE" = "orange"),
                    name = "Preference > 0.5",
                    labels = c("S", "R")) +  # Custom colors
  labs(x = "X", y = "Conversion Rate", title = "Boxplot of Conversion Rate by X") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
# Scatterplot color coded by preference (>0.5)
ggplot(dataset, aes(x = factor(X), y = conversion.rate, color = preference > 0.5)) +
  geom_jitter(alpha = 0.7, width = 0.2, size = 2) +  # Jittered scatterplot to avoid overlap
  theme_minimal() + 
  scale_color_manual(values = c("TRUE" = "blue", "FALSE" = "orange"),
                     name = "Preference > 0.5",
                     labels = c("S", "R")) +  # Custom colors
  labs(x = "X", y = "Conversion Rate", title = "Scatterplot of Conversion Rate by X") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# create BOXPLOT for ALL data with X = X and Y = preference
ggplot(dataset, aes(x = factor(X), y = preference)) +
  geom_boxplot(fill = "lightblue", alpha = 0.7, outlier.color = "red") +  # Boxplot
  theme_minimal() + 
  labs(x = "X", y = "Conversion Rate", title = "Boxplot of Preference by X") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# average preference (R/S ratio) vs average conversion rate
# Compute the average preference and conversion rate for each X category
avg_data <- dataset %>%
  group_by(X) %>%
  summarise(mean_preference = mean(preference, na.rm = TRUE),
            mean_conversion_rate = mean(conversion.rate, na.rm = TRUE))

# Determine color coding based on whether mean preference is > 0.5
avg_data <- avg_data %>%
  mutate(preference_group = mean_preference > 0.5)

ggplot(avg_data, aes(x = mean_preference, y = mean_conversion_rate, color = preference_group)) +
  geom_point(alpha = 0.7, size = 3) +  # Scatterplot
  geom_smooth(method = "lm", se = FALSE, color = "black") +  # Line of best fit
  scale_color_manual(values = c("TRUE" = "blue", "FALSE" = "orange"),
                     name = "Preference > 0.5",
                     labels = c("S", "R")) +  
  labs(x = "Average Preference", y = "Average Conversion Rate", title = "Average Preference against Average Conversion Rate") +
  theme_minimal()
        