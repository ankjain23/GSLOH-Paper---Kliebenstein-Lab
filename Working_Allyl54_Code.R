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
install.packages("ggpubr")
library(ggpubr)
library(ggplot2)
#setworkingdirectory to Allyl54xAllyl28 folder
setwd("~/Documents/URC Data Tables/Allyl54xAllyl28")
#read dataset, store it as Hollyhock
Hollyhock <- read_excel("Allyl28xAllyl54.xlsx")
#remove columns 13 to 148, save as Hollyhock_S
#Hollyhock_S <- Hollyhock[, -c(13:148)]

#renaming the colname to Isolate_Number so that it's a trait shared between datasets

#colnames(Hollyhock)[3] <- c("Isolate_number")
#Hollyhock$Isolate_number <- as.numeric(Hollyhock$Isolate_number)

#Hollyhock <- merge(Hollyhock, Iso72, by="Isolate_number")

#Hollyhock$Genotypes <- as.factor(Hollyhock$Genotypes)
#Hollyhock$Isolate_name <- as.factor(Hollyhock$Isolate_name)

#Plot density plot of the Lesion.mm2 variable in the Hollyhock dataset, with x-axis limits set from 0 to 8.
ggplot(Hollyhock, aes(Lesion.mm2)) +geom_density() + xlim(0,8)

# the first peak is right before 10, we can fix the lesion size threshold to 8mm2

#-------------------------------------

### Create a "Key" 
#create new variable Genotransfo by concatenating Genotype and Transtreat from Hollyhock Data
Hollyhock$GenoTransfo <- as.factor(paste(Hollyhock$Genotype, Hollyhock$TransTreat, sep="-"))

#summarize lesion area stats (Lesion.mm2) for each GenoTransfo group
Summary_GenoIsoTime <- summarise(group_by(Hollyhock, GenoTransfo), 
                                 mean.sp=mean(Lesion.mm2), 
                                 median.sp=median(Lesion.mm2), 
                                 sd.sp=sd(Lesion.mm2), 
                                 min.sp=min(Lesion.mm2), 
                                 max.sp=max(Lesion.mm2), 
                                 length.sp=length(Lesion.mm2), 
                                 cv.sp=sd.sp/mean.sp)

#merge summarized stats back into original data set
Hollyhock <- merge(Hollyhock, Summary_GenoIsoTime, by="GenoTransfo")

#plot density plot of median lesion size (median.sp) from summary data
ggplot(Summary_GenoIsoTime, aes(median.sp)) +geom_density() + xlim(0,20)

# the first peak is right before 10, we can fix the median threshold to 10mm2

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
HollyhockfilterOutliers$Genotype <- as.factor(HollyhockfilterOutliers$Genotype)
HollyhockfilterOutliers$Treatment <- as.factor(HollyhockfilterOutliers$Treatment)
levels(HollyhockfilterOutliers$Treatment)[2] <- c("GSLOH_Allyl54")

###Plot Violin and Boxplots of Lesion Size by Treatment
ggplot(HollyhockfilterOutliers, aes(Treatment, Lesion.mm2)) + 
  geom_violin() + 
  geom_boxplot(width = 0.1) +
  theme_bw()+
  ylab("Lesion area [mm2]")+
  ggtitle("GSLOH Allyl54")

#Perform ANOVA on linear model with Lesion.mm2 as response variable and treatment as predictor
Treat_lm <- lm(Lesion.mm2~Transgene, data=HollyhockfilterOutliers)
anova_result <- anova(Treat_lm)
print(anova_result)

# Convert binary Y/N to factors
HollyhockfilterOutliers$Allyl28 <- as.factor(HollyhockfilterOutliers$Allyl28)
HollyhockfilterOutliers$Allyl54 <- as.factor(HollyhockfilterOutliers$Allyl54)
HollyhockfilterOutliers$WT <- as.factor(HollyhockfilterOutliers$WT)

# Fit the three-way ANOVA model
model <- aov(Lesion.mm2 ~ Allyl28+Allyl54+WT+ Allyl28:Allyl54 + Allyl54:WT + Allyl28:WT, data = HollyhockfilterOutliers)
#View the summary of the ANOVA
summary(model)

# Perform ANOVA
anova_result <- aov(ThreeWay_lm)

# Print the results
print(anova_result)

# Perform three-way ANOVA
anova_result <- anova(ThreeWay_lm)

# Print the result
print(anova_result)

### Plots box plots of lesion size by treatment, colored by transformation, and faceted by Genotype.
ggplot(HollyhockfilterOutliers, aes(Transgene, Lesion.mm2, color=Transformation)) + 
  geom_boxplot() + 
  theme_bw()+
  ylab("lesion area [mm2]")+
  facet_grid(Isolate~.)


### Plots box plots of lesion size by Transformation, colored by treatment, and faceted by Genotype.

ggplot(HollyhockfilterOutliers, aes(Transgene, Lesion.mm2, color=Treatment)) + 
  geom_boxplot() + 
  theme_bw()+
  ylab("lesion area [mm2]")+
  facet_grid(Isolate~.)

### Plots box plots of lesion size by treatment, colored by transformation
ggplot(HollyhockfilterOutliers, aes(Transgene, Lesion.mm2, color=Transgene)) + 
  geom_boxplot() + 
  theme_bw()+
  ylab("lesion area [mm2]")

### Plots box plots of lesion size by treatment, colored by Plant, and faceted by Genotype.
ggplot(HollyhockfilterOutliers, aes(Treatment, Lesion.mm2, color=Plant)) + 
  geom_boxplot() + 
  theme_bw()+
  ylab("lesion area [mm2]")+
  facet_grid(Genotype~.)

####################### 


# Modeling with bioreps nested within treatment
###Performs linear modeling / ANOVA
lm72 <- lm(Lesion.mm2~Genotype+ Treatment/Transformation + Image + Plant, data=HollyhockfilterOutliers)
anova(lm72)

lm72 <- lm(Lesion.mm2~Genotype+ Treatment + Image + Plant, data=HollyhockfilterOutliers)
anova(lm72)


Treat_lm <- lm(Lesion.mm2~Treatment, data=HollyhockfilterOutliers)
anova(Treat_lm)

# Testing the individual lines, splits data by transformation and assigns to different variables

bioreps <- split(HollyhockfilterOutliers, HollyhockfilterOutliers$Transformation)

names(bioreps) <- c("d12", "b1", "e1", "g1")

list2env(bioreps, .GlobalEnv)

lm_d12 <- lm(Lesion.mm2~Genotype*Treatment + Image + Plant, data=d12)
anova(lm_d12)

lm_b1 <- lm(Lesion.mm2~Genotype*Treatment + Image + Plant, data=b1)
anova(lm_b1)

lm_e1 <- lm(Lesion.mm2~Genotype*Treatment + Image + Plant, data=e1)
anova(lm_e1)

lm_g1 <- lm(Lesion.mm2~Genotype*Treatment + Image + Plant, data=g1)
anova(lm_g1)


#above fits linear models to data frames (d12, b1, e1, g1) using Lesion.mm2 as response variable and genotype

##################################


ggplot(HollyhockfilterOutliers, aes(Treatment, Lesion.mm2, color=Genotype)) + geom_boxplot() + facet_grid(Transformation~.)



ggplot(Hollyhock, aes(Lesion.mm2, color=Isolate_name)) + geom_boxplot() + facet_grid(Timepoints~.)

ggplot(Hollyhock, aes(Lesion.mm2, color=Genotypes)) + geom_boxplot() + facet_grid(Timepoints~.)


#--------------------------------------

## starting with linear model- this is to separate time points but we don't have to

HollyhockfilterOutliers$GenoTime <- as.factor(paste(HollyhockfilterOutliers$Genotypes,HollyhockfilterOutliers$Timepoint, sep="-"))

HollyhockfilterOutliers_time <- split(HollyhockfilterOutliers, HollyhockfilterOutliers$Timepoint)
HollyhockfilterOutliers_time  <- as.array(HollyhockfilterOutliers_time)

list2env(HollyhockfilterOutliers_time, .GlobalEnv)

#ANOVA genotype is plant isolate is botrytis

#HPI72$Experiment <- as.factor(HPI72$Experiment)
lm72 <- lm(Lesion.mm2~Genotype+ Treatment/Transformation + Image + Plant, data=HollyhockfilterOutliers)
anova(lm72)

#Including experiment

#lm96_interactions <- lm(Lesion.mm2 ~ Genotypes * Isolate_name + Experiment + Image + Plants, data = HPI96)
#anova(lm96_interactions)

#lm72_interactions <- lm(Lesion.mm2 ~ Genotypes * Isolate_name + Experiment + Image + Plants, data = HPI72)
#anova(lm72_interactions)
#linear mixed model (| is random effect)
lmm72 <- lmer(Lesion.mm2~Genotypes*Isolate_name + (1|Image) + (1|Plants), data=HPI72)
anova(lmm72)
#lmm96 <- lmer(Lesion.mm2~Genotypes*Isolate_name + (1|Experiment) + (1|Image) + (1|Plants), data=HPI96)
#anova(lmm96)
#most likely to use
lmm72 <- lmer(Lesion.mm2~Genotypes*Isolate_name + (1|Image) + (1|Plants), data=HPI72)
anova(lmm72)

#lmm96 <- lmer(Lesion.mm2~Genotypes*Isolate_name +  (1|Image) + (1|Plants), data=HPI96)
#anova(lmm96)

#####################

HollyhockfilterOutliers_Isotime <- split(HollyhockfilterOutliers, HollyhockfilterOutliers$Genotypes)
HollyhockfilterOutliers_Isotime  <- as.array(HollyhockfilterOutliers_Isotime)

Sp_names <- names(HollyhockfilterOutliers_Isotime)

Lesion_LS <- as.data.frame(matrix(ncol=7))
colnames(Lesion_LS) <- c("Isolate_name", "emmean", "SE" , "df", "lower.CL", "upper.CL","Genotypes" )
emmeans_list <- list()

for (i in c(1:4)) {
  Species_model <- lmer(Lesion.mm2 ~ Isolate_name + (1|Image) +(1|Plants) , data=HollyhockfilterOutliers_Isotime[[i]]) 
  tmp <- emmeans(Species_model, ~Isolate_name, lmer.df='satterthwaite')
  tmp <- as.data.frame(print(tmp))
  tmp$Genotypes <-  rep(Sp_names[[i]],dim(tmp)[1])
  emmeans_list[[i]] <- tmp
  Lesion_LS <- rbind(Lesion_LS, tmp)
  
  
}


Lesion_LS <- Lesion_LS[-1,]
write.csv(Lesion_LS, file="Pepperexp_LSmean96_exp.csv")

Lesion_LS <- read.csv(file="Pepperexp_LSmean96_exp.csv")
#----------------------------------------------



#Lesion_LS  <-  separate(data=Lesion_LS, col=Genotype, into=c("Genotypes", "Time"), sep="\\-")


#Lesion_LS$Time <- as.factor(gsub('HPI72', 'LSM72', Lesion_LS$Time, fixed=F))
#Lesion_LS$Time <- as.factor(gsub('HPI96', 'LSM96', Lesion_LS$Time, fixed=F))


#Lesion_LS_sp <- split(Lesion_LS,Lesion_LS$Time)
#list2env(Lesion_LS_sp, .GlobalEnv)


LSM72s <- Lesion_LS[, -c(1,4,5,6,7)]
LSM96s <- Lesion_LS[, -c(1,4,5,6,7)]

Hollyhock_LSM_col <- pivot_wider(LSM72s, names_from = "Genotypes", values_from = "emmean")
Hollyhock_LSM_col <- pivot_wider(LSM96s, names_from = "Genotypes", values_from = "emmean")
write.table(Hollyhock_LSM_col, file="Pepperexp_LSM_col96_exp.txt")
Hollyhock_LSM_col <- read.csv(file="All4_72.csv")

rownames(Hollyhock_LSM_col) <- Hollyhock_LSM_col[,1]

Hollyhock_LSM_col <- Hollyhock_LSM_col[,-1]
Hollyhock_LSM_col <- as.matrix(Hollyhock_LSM_col)
library(heatmaply)
Eudicot_heatmap_stand <- main_heatmap(Hollyhock_LSM_col, 
                                      name = "Standardized lesion size") %>%
  add_col_clustering() %>%
  add_row_clustering() %>%
  add_row_title("Plant Genotypes") %>%
  add_col_title("Botrytis Isolates") %>%
  add_col_labels() %>%
  # add_col_title() %>%
  add_col_summary()

Hollyhock_LSM_col <- as.matrix(Hollyhock_LSM_col)
heatmap(Hollyhock_LSM_col)
# above makes heatmap of lesion size means

library(pheatmap)
#pheatmap(Hollyhock_LSM_col)
pheatmap(Hollyhock_LSM_col, fontsize = 8, fontfamily = "Arial")

# Plot violin plots and boxplots for lesion size by genotype
ggplot(HollyhockfilterOutliers,aes(Isolate, Lesion.mm2, color=Isolate))+
  geom_violin()+
  geom_boxplot(width=0.2)+
  theme_bw()
# Define specific group comparisons
comparisons <- list(
  c("Allyl28", "Allyl54"),
  c("Allyl54", "WT"),
  c("Allyl28", "WT")
)

#more librarys
library(rstatix)
library(dplyr)
#
ggplot(HollyhockfilterOutliers, aes(x = Transgene, y = Lesion.mm2, color = Transgene)) +
  geom_violin() +
  stat_compare_means(method="anova", label = "p.signif") +
  theme_bw()
#violin plot with P values
stat.test <- HollyhockfilterOutliers %>%
  t_test(Lesion.mm2 ~ Transgene) %>%
  add_y_position()
#manually insert p values from ANOVA here, print to check order of groups
stat.test$p <- c("<0.001", 0.92, "<0.001")
print(stat.test, width = Inf)
#print violin plot with p values and boxplots
ggplot(HollyhockfilterOutliers,aes(Transgene, Lesion.mm2))+
  geom_violin()+
  stat_pvalue_manual(stat.test, label = "p = {p}")+
  geom_boxplot(width = 0.2, color = "black", fill = "white", outlier.shape = NA)
  theme_bw()

#plots violin and boxplots with nice colors
ggplot(HollyhockfilterOutliers, aes(x = Var2, y = value)) +
  geom_violin(aes(fill = Var2)) +
  geom_boxplot(width = 0.2, color = "black", fill = "white", outlier.shape = NA) +
  scale_fill_manual(values = c("Bv1" = "#BFD3C1", "Bv2" = "pink", "Bv3" = "#A4DDED", "Bv4" = "#D8BFD8")) +  # Custom colors
  theme_minimal() +
  labs(x = "Genotypes", y = "LsMeans_cm2", fill = "Variables") 
