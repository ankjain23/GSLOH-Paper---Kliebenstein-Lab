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
#read dataset, store it as dataset
dataset <- read_excel("ImageJ_Size_Data_2.xlsx")

# FIRST SUBSET DATA based on gene/wt group of interest
##Subset Data for Allyl 28 vs Allyl 54 Data Only (237-277)
selected_rows <- dataset[237:276, ]

##Subset Data for Allyl 28 vs Control Data Only (1-80)
selected_rows <- dataset[1:80, ]

##Subset Data for Allyl 54 vs Control Only (81-152) (OLD data using Sha-1 AND CVI)
selected_rows <- dataset[81:152, ]
## 54/WT Data for CVI ONLY
selected_rows <- dataset[117:152, ]
## 54/WT Data for Sha-1 ONLY
selected_rows <- dataset[81:116, ]
#ggplot for subset
ggplot(selected_rows, aes(Lesion.mm2)) +geom_density() + xlim(0,8)

###Plot Violin and Boxplots of Lesion Size by Treatment
ggplot(selected_rows, aes(Treatment, Lesion.mm2)) + 
  geom_violin() + 
  geom_boxplot(width = 0.1) +
  theme_bw()+
  ylab("Lesion area [mm2]")+
  ggtitle("")

#Binomial Test
binom.test(x=10, n=10, p=0.5)

###Plot Violin of Percent Eaten by Gene - SKIP ALL THIS DOWN TO AVERAGE DATA PART (line 72)
#Convert Percent Eaten to Continuous Variable
#selected_rows$`Percent Eaten` <- as.numeric(as.character(selected_rows$`Percent Eaten`))

#ggplot(selected_rows, aes(x = Treatment, y = 'Percent Eaten')) + 
  geom_boxplot(width=0.1)+
  theme_bw() + 
  ylab("Percent Eaten") + 
  xlab("Treatment") +
  ggtitle("Line Plot of Percent Eaten by Treatment")

#selected_rows$`Percent Eaten` <- as.numeric(as.character(selected_rows$`Percent Eaten`))

#Boxplot of genes by
ggplot(selected_rows, aes(x = Gene, y = `Percent Eaten`, color = Image)) + 
  geom_boxplot(aes(group=Image)) +
  theme_bw()

#Percent Eaten as a Dotted box plot, Compared across Genes and colored by Tray
ggplot(selected_rows, aes(x = Gene, y = `Percent Eaten`, color = Image)) + 
    geom_point()+
    geom_line(aes(group = (Image)))
  theme_bw()
#Percent Eaten Dotted Plot, colored by Tray (Image) for Allyl(X) x WT 
  ggplot(selected_rows, aes(x = Treatment, y = `Percent Eaten`, color = Image)) + 
    geom_point()+
    geom_line(aes(group = Image))+
    scale_x_discrete(labels = c("control" = "Control", "gene" = "Allyl54"))+
  theme_bw()


# Assuming 'selected_rows' is your data frame, plot AVERAGE Allyl28/54 for each image (FINAL CODE!!!!)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
# First, calculate the average percentage eaten for each tray and gene (Allyl28 and Allyl54)
  avg_data <- selected_rows %>%
    group_by(Image, Treatment) %>%
    summarise(Avg_Percent_Eaten = mean(`Percent Eaten`, na.rm = TRUE))
  
# Now create the plot
  ggplot(avg_data, aes(x =Treatment, y = Avg_Percent_Eaten, color = Image, group = Image)) + 
    geom_line(aes(group = Image), size = 1) +  # Draw lines connecting each pair of points for each tray
    geom_point(size = 2) +  # Add points for clarity
    theme_bw() +
    labs(x = "Treatment", y = "Average Percent Eaten", color="Tray (Image)") +
    ggtitle("Averge Percent Eaten of Allyl54 and WT Across Trays for Sha-1")
    theme(
      legend.position = "bottom",
      legend.title = element_text(size=10),
      legend.text = element_text(size=8),
      legend.key.width= unit(1.5, "lines")) # Optional: Remove legend if not needed
  
#Dots connected by lines plot with Image labels on graph
library(ggrepel)
    
    ggplot(avg_data, aes(x = Treatment, y = Avg_Percent_Eaten, color = Image, group = Image)) + 
      geom_line(size = 1) +  # Draw lines connecting each pair of points for each tray
      geom_point(size = 2) +  # Add points for clarity
      geom_text_repel(aes(label = Image), 
                      nudge_x = 0.5,  # Adjust horizontal spacing of labels
                      size = 3,      # Adjust text size
                      show.legend = FALSE) +  # Hide text labels from the legend
      theme_bw() +
      labs(x = "Treatment", y = "Average Percent Eaten", color = "Tray (Image)") +
      ggtitle("Average Percent Eaten of Allyl28 and WT Across Trays") +
      theme(
        legend.position = "right",  # Optional: Move the legend to the bottom
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8)
      )
##ANOVA, first average data (For Allyl28/54 change "Treatment" to "Gene")
    avg_data <- selected_rows %>%
      group_by(Image, Treatment) %>%
      summarise(Avg_Percent_Eaten = mean(`Percent Eaten`, na.rm = TRUE))
#ANOVA for average data of Treatment-Image interaction on Percent Eaten (For Allyl28/54 change Treatment to Gene)
anova_result <- aov(Avg_Percent_Eaten ~ Treatment + Image, data = avg_data)
summary(anova_result)

###Color Replicates together on Graphs (do this after filtering by gene/wt grouping and creating "avg_data" variable)
## Create a new grouping variable in avg_data
#Allyl28vsWT
avg_data$ImageGroup <- case_when(
  avg_data$Image %in% c("IMG_0196", "IMG_0198", "IMG_0199") ~ "Group1",
  avg_data$Image %in% c("IMG_0200", "IMG_0201", "IMG_0202") ~ "Group2",
  avg_data$Image %in% c("IMG_0203", "IMG_0204", "IMG_0205", "IMG_0206") ~ "Group3",
  avg_data$Image %in% c("IMG_0207", "IMG_0208", "IMG_0209", "IMG_0210") ~ "Group4",
  avg_data$Image %in% c("IMG_0211", "IMG_0212", "IMG_0213") ~ "Group5",
  avg_data$Image %in% c("IMG_0214", "IMG_0215", "IMG_0216") ~ "Group6",
  TRUE ~ "Other"  # For any remaining images not listed
)
#Allyl54vsWT CVI AND SHA1 DATA
avg_data$ImageGroup <- case_when(
  avg_data$Image %in% c("IMG_0217", "IMG_0218", "IMG_0219") ~ "Group7",
  avg_data$Image %in% c("IMG_0220", "IMG_0221", "IMG_0222") ~ "Group8",
  avg_data$Image %in% c("IMG_0223", "IMG_0224", "IMG_0225") ~ "Group9",
  avg_data$Image %in% c("IMG_0226", "IMG_0227", "IMG_0228", "IMG_0210") ~ "Group10",
  avg_data$Image %in% c("IMG_0229", "IMG_0230", "IMG_0231") ~ "Group11",
  avg_data$Image %in% c("IMG_0232", "IMG_0233", "IMG_0234") ~ "Group12",
  TRUE ~ "Other"  # For any remaining images not listed
)
#Allyl54vsWT CVI DATA ONLY
avg_data$ImageGroup <- case_when(
  avg_data$Image %in% c("IMG_0226", "IMG_0227", "IMG_0228", "IMG_0210") ~ "Group10",
  avg_data$Image %in% c("IMG_0229", "IMG_0230", "IMG_0231") ~ "Group11",
  avg_data$Image %in% c("IMG_0232", "IMG_0233", "IMG_0234") ~ "Group12",
  TRUE ~ "Other"  # For any remaining images not listed
)
#Allyl54vsWT SHA1 DATA ONLY
avg_data$ImageGroup <- case_when(
  avg_data$Image %in% c("IMG_0217", "IMG_0218", "IMG_0219") ~ "Group7",
  avg_data$Image %in% c("IMG_0220", "IMG_0221", "IMG_0222") ~ "Group8",
  avg_data$Image %in% c("IMG_0223", "IMG_0224", "IMG_0225") ~ "Group9",
  TRUE ~ "Other"  # For any remaining images not listed
)

#Allyl28vs54
avg_data$ImageGroup <- case_when(
  avg_data$Image %in% c("IMG_0256", "IMG_0257", "IMG_0258", "IMG_0259") ~ "Group13",
  avg_data$Image %in% c("IMG_0260", "IMG_0261", "IMG_0262") ~ "Group14",
  avg_data$Image %in% c("IMG_0263", "IMG_0264", "IMG_0265") ~ "Group15",
  TRUE ~ "Other"  # For any remaining images not listed
)
# Plot with the new grouping variable for color (Make sure to rename Title based on Subset of Data)
## Change "Treatment" to "Gene" for Allyl28 vs Allyl54 Data)
ggplot(avg_data, aes(x = Treatment, y = Avg_Percent_Eaten, color = ImageGroup, group = Image)) + 
  geom_line(size = 1) +  # Draw lines connecting each pair of points for each tray
  geom_point(size = 2) +  # Add points for clarity
  theme_bw() +
  labs(x = "Treatment", y = "Average Percent Eaten", color = "Image Group") +
  ggtitle("Average Percent Eaten of Allyl54 and Wildtype Across Trays") +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    legend.key.width = unit(1.5, "lines")
  )


#Binomial Test - success = more WT than gene eaten (based on avg data)
# Ally28vsWT
binom.test(8, 20, p = 0.5, alternative = "two.sided")
# Allyl54vsWT (CVI and Sha-1 Data)
binom.test(4, 18, p = 0.5, alternative = "two.sided")
#Allyl54vsWT (CVI ONLY)
binom.test(1, 9, p=0.5, alternative = "two.sided")
#Allyl54vsWT (SHA1 ONLY)
binom.test(4, 9, p=0.5, alternative = "two.sided")
# Allyl28vs54 - success = more 54 than 28 eaten (based on avg data)
binom.test(1, 10, p = 0.5, alternative = "two.sided")

#Chi Squared Goodness-of-Fit Test - tests if two categorical variables are independent
#Calculate average gene and WT percent eaten across the data
avg_summary <- tapply(avg_data$Avg_Percent_Eaten, avg_data$Treatment, mean, na.rm = TRUE)
print(avg_summary)
### For Allyl28 vs Allyl54 use "Gene" not "Treatment" as neither are controls
avg_summary <- tapply(avg_data$Avg_Percent_Eaten, avg_data$Gene, mean, na.rm = TRUE)
print(avg_summary)

# Create the contingency table data
data <- matrix(
  c(
    4.310413, 4.697683,  # "28/WT" for "gene" and "control"
    6.720285, 3.426256,  # "54/WT" for "gene" and "control"
    10.551163, 2.885945   # "28/54" for "gene" and "control"
  ),
  nrow = 3,              
  byrow = TRUE           
)
# Assign row and column names, create contingency table
rownames(data) <- c("28/WT", "54/WT", "28/54")  # Row categories
colnames(data) <- c("Gene", "Control")          # Column categories
contingency_table <- as.table(data)
print(contingency_table)
# Perform the chi-squared test
chisq_test <- chisq.test(contingency_table)
# View the results
print(chisq_test)

# Mirror Bar Plot
library(ggplot2)
##Data preparation, for OLD binomial tests using both CVI and Sha-1
data <- data.frame(
  Comparison = c("WT/28", "WT/28", "WT/54", "WT/54", "54/28", "54/28"),
  Percent = c(40, -60, 22, -78, -10, 90),  # Negative values for mirroring
  Group = c("More WT Eaten", "More 28 Eaten", "More WT Eaten", "More 54 Eaten", "More 54 Eaten", "More 28 Eaten")  # Bar categories
)
## Data for NEW binomial tests using ONLY CVI across all data
data <- data.frame(
  Comparison = c("WT/28", "WT/28", "WT/54", "WT/54", "54/28", "54/28"),
  Percent = c(40, -60, 11, -89, -10, 90),  # Negative values for mirroring
  Group = c("More WT Eaten", "More 28 Eaten", "More WT Eaten", "More 54 Eaten", "More 54 Eaten", "More 28 Eaten")  # Bar categories
)
###Plot in color
ggplot(data, aes(x = Comparison, y = Percent, fill = Group)) +
  geom_bar(stat = "identity", position = "identity", width = 0.6) +
  coord_flip() +  # Flip coordinates for horizontal bars
  theme_minimal() +
  labs(
    x = "Comparison",
    y = "Percent of Trials",
    title = "Mirrored Bar Plot of Binomial Test Results"
  ) +
  scale_y_continuous(
    breaks = seq(-100, 100, 20),
    labels = abs(seq(-100, 100, 20))  # Absolute values for y-axis labels
  ) +
  scale_fill_manual(values = c("More WT" = "skyblue", "More 28" = "orange", "More 54" = "purple")) +
  theme(
    legend.title = element_blank(),  # Hide legend title
    legend.position = "bottom"       # Place legend at the bottom
  )

###Plot in Black and White
ggplot(data, aes(x = Comparison, y = Percent, fill = Group)) +
  geom_bar(stat = "identity", position = "identity", width = 0.6, color = "black") +  # Black outline
  coord_flip() +  # Flip coordinates for horizontal bars
  theme_minimal() +
  labs(
    x = "Comparison",
    y = "Percent of Trials",
    title = "Mirrored Bar Plot of Binomial Test Results"
  ) +
  scale_y_continuous(
    breaks = seq(-100, 100, 20),
    labels = abs(seq(-100, 100, 20))  # Absolute values for y-axis labels
  ) +
  scale_fill_manual(
    values = c("More WT Eaten" = "gray80", "More 28 Eaten" = "gray50", "More 54 Eaten" = "gray20")  # Grayscale shades
  ) +
  theme(
    legend.title = element_blank(),  # Hide legend title
    legend.position = "bottom",      # Place legend at the bottom
    panel.grid.major = element_line(color = "gray90"),  # Light gray grid
    panel.grid.minor = element_blank()                  # Remove minor grid lines
  )

## Flip order of Comparisons to be Wt/28, Wt/54 and 54/28
# Set the desired order for the Comparison factor
data$Comparison <- factor(data$Comparison, levels = c("54/28", "WT/54", "WT/28"))

# Plot
ggplot(data, aes(x = Comparison, y = Percent, fill = Group)) +
  geom_bar(stat = "identity", position = "identity", width = 0.6, color = "black") +  # Black outline
  coord_flip() +  # Flip coordinates for horizontal bars
  theme_minimal() +
  labs(
    x = "Comparison",
    y = "Percent of Trials",
    title = "Mirrored Bar Plot of Binomial Test Results"
  ) +
  scale_y_continuous(
    breaks = seq(-100, 100, 20),
    labels = abs(seq(-100, 100, 20))  # Absolute values for y-axis labels
  ) +
  scale_fill_manual(
    values = c("More WT Eaten" = "gray80", "More 28 Eaten" = "gray50", "More 54 Eaten" = "gray20")  # Grayscale shades
  ) +
  theme(
    legend.title = element_blank(),  # Hide legend title
    legend.position = "bottom",      # Place legend at the bottom
    panel.grid.major = element_line(color = "gray90"),  # Light gray grid
    panel.grid.minor = element_blank()                  # Remove minor grid lines
  )

