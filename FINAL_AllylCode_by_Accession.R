#loading libraries
library(tidyverse)
library(dplyr)
library(tidyr)
library(MASS)
library(lme4)
library(car)
library(emmeans)
library(iheatmapr)
library(readxl)
#read dataset, store it as Hollyhock
Hollyhock <- read_excel("Allyl28xAllyl54_Cropped_Accessions.xlsx")
#subset data based on your needs
#CVI Only
selected_rows <- Hollyhock[1:205, ]
#Sha-1 Only
selected_rows <- Hollyhock[206:564, ]

#Plot density plot of the Lesion.mm2 variable in the Hollyhock dataset, with x-axis limits set from 0 to 8.
ggplot(selected_rows, aes(Lesion.mm2)) +geom_density() + xlim(0,8)

# the first peak is right before 10, we can fix the lesion size threshold to 8mm2

#-------------------------------------

### Create a "Key" 
#create new variable Genotransfo by concatenating Isolate and Transtreat from Hollyhock Data
Hollyhock$GenoTransfo <- as.factor(paste(Hollyhock$Isolate, Hollyhock$TransTreat, sep="-"))

selected_rows$GenoTransfo <- as.factor(paste(selected_rows$Isolate, selected_rows$TransTreat, sep="-"))

#summarize lesion area stats (Lesion.mm2) for each GenoTransfo group
Summary_GenoIsoTime <- summarise(group_by(Hollyhock, GenoTransfo), 
                                 mean.sp=mean(Lesion.mm2), 
                                 median.sp=median(Lesion.mm2), 
                                 sd.sp=sd(Lesion.mm2), 
                                 min.sp=min(Lesion.mm2), 
                                 max.sp=max(Lesion.mm2), 
                                 length.sp=length(Lesion.mm2), 
                                 cv.sp=sd.sp/mean.sp)
#selected rows version
Summary_GenoIsoTime <- summarise(group_by(selected_rows, GenoTransfo), 
                                 mean.sp=mean(Lesion.mm2), 
                                 median.sp=median(Lesion.mm2), 
                                 sd.sp=sd(Lesion.mm2), 
                                 min.sp=min(Lesion.mm2), 
                                 max.sp=max(Lesion.mm2), 
                                 length.sp=length(Lesion.mm2), 
                                 cv.sp=sd.sp/mean.sp)
#merge summarized stats back into original data set
#Hollyhock <- merge(Hollyhock, Summary_GenoIsoTime, by="GenoTransfo")
Hollyhock <- merge(selected_rows, Summary_GenoIsoTime, by="GenoTransfo")

#plot density plot of median lesion size (median.sp) from summary data
ggplot(Summary_GenoIsoTime, aes(median.sp)) +geom_density() + xlim(0,20)

# Filter out failed lesions and outliers
splt_Data_cleaning_filterOutliers <- list()
Failed_lesions <- list()
##Filter outliers based on lesion size
Hollyhock_S <- mutate(Hollyhock,
median_logic=median.sp>6,
threshold_logic=Lesion.mm2<6,
TBC_param=as.factor(paste(median_logic, threshold_logic, sep="")),
count_TBC=length(Lesion.mm2[TBC_param=="TRUETRUE"]))


# Create new data set for filtered lesions - selects all data points that aren't 'TRUETRUE'
HollyhockfilterOutliers <- Hollyhock_S %>% filter(TBC_param!="TRUETRUE")

# Create new dataset for Failed Lesions - selects all data points that are 'TRUETRUE' aka removed from the data set
Failed_lesions <- Hollyhock_S %>% filter(TBC_param=="TRUETRUE")


#-------------------------------------
# checking the effect plot everything to see how things look like w/o emmeans

HollyhockfilterOutliers$Transformation <- as.factor(HollyhockfilterOutliers$Transformation )
HollyhockfilterOutliers$Isolate <- as.factor(HollyhockfilterOutliers$Isolate)
HollyhockfilterOutliers$Transgene <- as.factor(HollyhockfilterOutliers$Transgene)
levels(HollyhockfilterOutliers$Transgene)[2] <- c("Allyl54")

###Plot Violin and Boxplots of Lesion Size by Transgene
ggplot(HollyhockfilterOutliers, aes(Transgene, Lesion.mm2)) + 
  geom_violin() + 
  geom_boxplot(width = 0.1) +
  theme_bw()+
  ylab("Lesion area [mm2]")+
  ggtitle("Allyl54 vs Allyl28 vs WT")

#Perform ANOVA on linear model with Lesion.mm2 as response variable and Transgene as predictor
Treat_lm <- lm(Lesion.mm2~Transgene, data=HollyhockfilterOutliers)
anova_result <- anova(Treat_lm)
print(anova_result)

#ANOVA with Lesion.mm2 as response and Transgene, Isolate-Transgene interaction, Image, Plant nested within Transfo, and Trasnfo nested within Transgene as predictor
Treat_lm <- lm(Lesion.mm2~Isolate+Transgene+Isolate*Transgene
               +Image+Transformation/Plant
               +Transgene/Transformation, data=HollyhockfilterOutliers)
anova_result <- anova(Treat_lm)
print(anova_result)

#FILTERING DATA TO ONLY INCLUDE ALLYL54 AND ALLYL28
A54xA28_data <- subset(HollyhockfilterOutliers, Transgene %in% c("Allyl54", "Allyl28"))
Treat_lm <- lm(Lesion.mm2~Isolate+Transgene+Isolate*Transgene
               +Image+(Transformation:Transgene)+(Transformation:Transgene:Plant), data=A54xA28_data)
anova_result <- anova(Treat_lm)
print(anova_result)

#FILTERING DATA TO ONLY INCLUDE ALLYL54 AND WT
A54xWT_data <- subset(HollyhockfilterOutliers, Transgene %in% c("WT", "Allyl54"))
Treat_lm <- lm(Lesion.mm2~Isolate+Transgene+Isolate*Transgene
               +Image+(Transformation:Transgene)+(Transformation:Transgene:Plant), data=A54xWT_data)
anova_result <- anova(Treat_lm)
print(anova_result)

#FILTERING DATA TO ONLY INCLUDE ALLYL28 AND WT
A28xWT_data <- subset(HollyhockfilterOutliers, Transgene %in% c("WT", "Allyl28"))
Treat_lm <- lm(Lesion.mm2~Isolate+Transgene+Isolate*Transgene
               +Image+(Transformation:Transgene)+(Transformation:Transgene:Plant), data=A28xWT_data)
anova_result <- anova(Treat_lm)
print(anova_result)

### Plots box plots of lesion size by Transgene, colored by transformation, and faceted by Isolate.
ggplot(HollyhockfilterOutliers, aes(Transgene, Lesion.mm2, color=Transformation)) + 
  geom_boxplot() + 
  theme_bw()+
  ylab("lesion area [mm2]")+
  facet_grid(Isolate~.)


### Plots box plots of lesion size by Transformation, colored by Transgene, and faceted by Isolate.

ggplot(HollyhockfilterOutliers, aes(Transformation, Lesion.mm2, color=Transgene)) + 
  geom_boxplot() + 
  theme_bw()+
  ylab("lesion area [mm2]")+
  facet_grid(Isolate~.)

### Plots box plots of lesion size by Transgene, colored by transformation
ggplot(HollyhockfilterOutliers, aes(Transgene, Lesion.mm2, color=Transformation)) + 
  geom_boxplot() + 
  theme_bw()+
  ylab("lesion area [mm2]")

### Plots box plots of lesion size by Transgene, colored by Plant, and faceted by Isolate.
ggplot(HollyhockfilterOutliers, aes(Transgene, Lesion.mm2, color=Plant)) + 
  geom_boxplot() + 
  theme_bw()+
  ylab("lesion area [mm2]")+
  facet_grid(Isolate~.)

####################### 


# Modeling with bioreps nested within Transgene
###Performs linear modeling / ANOVA
lm72 <- lm(Lesion.mm2~Isolate+ Transgene/Transformation + Image + Plant, data=HollyhockfilterOutliers)
anova(lm72)

lm72 <- lm(Lesion.mm2~Isolate+ Transgene + Image + Plant, data=HollyhockfilterOutliers)
anova(lm72)


Treat_lm <- lm(Lesion.mm2~Transgene, data=HollyhockfilterOutliers)
anova(Treat_lm)

# Testing the individual lines, splits data by transformation and assigns to different variables

bioreps <- split(HollyhockfilterOutliers, HollyhockfilterOutliers$Transformation)

names(bioreps) <- c("d12", "b1", "e1", "g1")

list2env(bioreps, .GlobalEnv)

lm_d12 <- lm(Lesion.mm2~Isolate*Transgene + Image + Plant, data=d12)
anova(lm_d12)

lm_b1 <- lm(Lesion.mm2~Isolate*Transgene + Image + Plant, data=b1)
anova(lm_b1)

lm_e1 <- lm(Lesion.mm2~Isolate*Transgene + Image + Plant, data=e1)
anova(lm_e1)

lm_g1 <- lm(Lesion.mm2~Isolate*Transgene + Image + Plant, data=g1)
anova(lm_g1)


#above fits linear models to data frames (d12, b1, e1, g1) using Lesion.mm2 as response variable and Isolate

##################################


ggplot(HollyhockfilterOutliers, aes(Transgene, Lesion.mm2, color=Isolate)) + geom_boxplot() + facet_grid(Transformation~.)



ggplot(Hollyhock, aes(Lesion.mm2, color=Isolate_name)) + geom_boxplot() + facet_grid(Timepoints~.)

ggplot(Hollyhock, aes(Lesion.mm2, color=Isolates)) + geom_boxplot() + facet_grid(Timepoints~.)


#--------------------------------------

## starting with linear model- this is to separate time points but we don't have to

HollyhockfilterOutliers$GenoTime <- as.factor(paste(HollyhockfilterOutliers$Isolates,HollyhockfilterOutliers$Timepoint, sep="-"))

HollyhockfilterOutliers_time <- split(HollyhockfilterOutliers, HollyhockfilterOutliers$Timepoint)
HollyhockfilterOutliers_time  <- as.array(HollyhockfilterOutliers_time)

list2env(HollyhockfilterOutliers_time, .GlobalEnv)

#ANOVA Isolate is plant isolate is botrytis

# Step 1: Calculate the average lesion size per isolate and transformation event
library(dplyr)
avg_lesion <- HollyhockfilterOutliers %>%
  group_by(Isolate, Transgene) %>%
  summarise(Avg_Lesion = mean(Lesion.mm2, na.rm = TRUE), .groups = "drop")

control_avg <- avg_lesion %>%
  filter(Transgene == "WT") %>%
  rename(Control_Lesion = Avg_Lesion)

# Step 2: Normalize lesion sizes by dividing by control values
normalized_lesion <- avg_lesion %>%
  left_join(control_avg, by = "Isolate") %>%  # Merge control values
  mutate(Normalized_Lesion = Avg_Lesion / Control_Lesion)
print(normalized_lesion)

# Load necessary libraries
install.packages('DT')
library(DT)
datatable(normalized_lesion, options = list(pageLength = 10), caption = "Normalized Lesion Sizes by Isolate and Transgene")

# Load necessary libraries
library(dplyr)
library(knitr)

# Step 1: Compute the average lesion size per isolate, transformation, and transgene
avg_lesion <- HollyhockfilterOutliers %>%
  group_by(Isolate, Transformation, Transgene) %>%
  summarise(Avg_Lesion = mean(Lesion.mm2, na.rm = TRUE), .groups = "drop")

# Step 2: Extract WT control values for each isolate and transformation event
control_avg <- avg_lesion %>%
  filter(Transgene == "WT") %>%
  rename(Control_Lesion = Avg_Lesion)

# Step 3: Normalize lesion sizes by dividing by the corresponding WT lesion size
normalized_lesion <- avg_lesion %>%
  left_join(control_avg, by = c("Isolate", "Transformation")) %>%  # Match by Isolate & Transformation
  mutate(Percent_Difference = ((Avg_Lesion - Control_Lesion) / Control_Lesion) * 100)

# Step 4: Print as a formatted table
kable(normalized_lesion, digits = 2, caption = "Percent Difference in Lesion Size Compared to WT Control")

library(ggplot2)

ggplot(normalized_lesion, aes(x = Percent_Difference)) +
  geom_histogram(binwidth = 5, fill = "gray", color = "black") +
  theme_bw() +
  labs(title = "Distribution of Percent Differences in Lesion Size",
       x = "Percent Difference from WT Control",
       y = "Frequency")

ggplot(normalized_lesion, aes(x = Percent_Difference)) +
  geom_density(fill = "gray", alpha = 0.5, color = "black") +
  theme_bw() +
  labs(title = "Density Plot of Percent Differences",
       x = "Percent Difference from WT Control",
       y = "Density")


#Violin Plot of Allyl28 and Allyl54 for normalized_lesions
ggplot(normalized_lesion[normalized_lesion$Transgene.x %in% c("Allyl28", "Allyl54"), ], 
       aes(x = Transgene.x, y = Percent_Difference, fill = Transgene.x)) +
  geom_violin(trim = FALSE, color = "black", alpha = 0.5) +  # Create violins
  scale_fill_manual(values = c("gray40", "gray70")) +  # Color for Allyl28 and Allyl54
  theme_bw() +
  labs(title = "Distribution of Percent Differences for Allyl28 and Allyl54",
       x = "Transgene",
       y = "Percent Difference from WT Control") +
  theme(text = element_text(size = 10), legend.position = "none")

ggplot(normalized_lesion[normalized_lesion$Transgene.x %in% c("Allyl28", "Allyl54"), ], 
       aes(x = Transgene.x, y = Percent_Difference, fill = Transgene.x)) +
  geom_violin(trim = FALSE, color = "black", alpha = 0.5) +  # Create violins
  scale_fill_manual(values = c("gray40", "gray70")) +  # Color for Allyl28 and Allyl54
  theme_bw() +
  labs(title = "Distribution of Percent Differences for Allyl28 and Allyl54",
       x = "Transgene",
       y = "Percent Difference from WT Control") +
  theme(text = element_text(size = 10), legend.position = "none")

ggplot(normalized_lesion[normalized_lesion$Transgene.x %in% c("Allyl28", "Allyl54"), ], 
       aes(x = Isolate, y = Percent_Difference, fill = Transgene.x)) +
  geom_violin(trim = FALSE, color = "black", alpha = 0.5) +  # Create violins
  scale_fill_manual(values = c("gray40", "gray70")) +  # Color for Allyl28 and Allyl54
  theme_bw() +
  labs(title = "Distribution of Percent Differences for Allyl28 and Allyl54",
       x = "Transgene",
       y = "Percent Difference from WT Control") +
  theme(text = element_text(size = 10), legend.position = "none")

#library(dplyr)

# Calculate the average lesion size for Allyl28 and Allyl54 per isolate
avg_lesion_size <- normalized_lesion %>%
  filter(Transgene.x %in% c("Allyl28", "Allyl54")) %>%
  group_by(Isolate, Transgene.x) %>%
  summarise(Average_Lesion_Size = mean(Avg_Lesion))

# Spread the data so that Allyl28 and Allyl54 are in separate columns
avg_lesion_size_wide <- avg_lesion_size %>%
  spread(key = Transgene.x, value = Average_Lesion_Size)

# Now calculate the percent difference between Allyl54 and Allyl28
avg_lesion_size_wide <- avg_lesion_size_wide %>%
  mutate(Percent_Difference_Allyl54_Allyl28 = 
           ((`Allyl54` - `Allyl28`) / `Allyl28`) * 100)

# Optionally, you can gather the data back into long format if needed
avg_lesion_size_long <- avg_lesion_size_wide %>%
  gather(key = "Transgene", value = "Average_Lesion_Size", Allyl28, Allyl54)

# View the updated data frame with Percent_Difference_Allyl54_Allyl28
print(avg_lesion_size_wide)

library(dplyr)
library(tidyr)

# Step 1: Calculate average lesion size for Allyl28 and Allyl54 per isolate
avg_lesion_size <- normalized_lesion %>%
  filter(Transgene %in% c("Allyl28", "Allyl54")) %>%
  group_by(Isolate, Transgene) %>%
  summarise(Average_Lesion_Size = mean(Lesion_Size)) %>%
  ungroup()
# Step 2: Spread the data so that Allyl28 and Allyl54 have their own columns
avg_lesion_size_wide <- avg_lesion_size %>%
  spread(key = Transgene, value = Average_Lesion_Size)
# Step 3: Calculate percent difference between Allyl54 and Allyl28
avg_lesion_size_wide <- avg_lesion_size_wide %>%
  mutate(Percent_Difference_Allyl54_Allyl28 = 
           ((`Allyl54` - `Allyl28`) / `Allyl28`) * 100)
# Step 4: Join this data back to the original 'normalized_lesion' table by 'Isolate'
normalized_lesion_with_diff <- normalized_lesion %>%
  left_join(avg_lesion_size_wide)
# View the updated table
print(normalized_lesion_with_diff)

