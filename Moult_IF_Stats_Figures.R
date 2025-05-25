#*********************************************************************************
#Title----
#Immune function varies during adult primary moult but not juvenile body moult in a European passerine 

#Malin V. Klumpp & Arne Hegemann

#Script by M. Klumpp, May 2025
#This script provides all code necessary to reproduce the statistical analyses and figures in the manuscript and supporting information

#*********************************************************************************
rm(list = ls())
setwd() #set working directory

#*********************************************************************************
#Packages ----
#*********************************************************************************
library(car)
library(dplyr)
library(ggplot2)
library(ggeffects)
library(ggpubr)
library(performance)
library(DHARMa)
library(lubridate)
library(betareg)

#*********************************************************************************
#Data  ----

#load data
moult_data = read.csv("blackbird_database_for_submission.csv") #make sure data file is in the same working directory

#Data pruning
moult_data <- moult_data %>% 
  mutate(lysis_score = as.numeric(lysis_score), #reclassify variables
         killingD = as.numeric(killingD),
         AC_hpconc_OLD = as.numeric(AC_hpconc_OLD),
         aglutt_score = as.numeric(aglutt_score),
         moult.score = as.numeric(moult.score),
         moulting = as.factor(moulting),
         sex_conclusion = if_else(sex_conclusion %in% c("F", "M"), sex_conclusion, NA_character_), #sex is "NA" if not specified
         sex = as.factor(sex_conclusion))#rename column

#Create a new column with Julian dates for capture dates
moult_data<- moult_data %>%
mutate(date = dmy(date), #transform date stored as character into date object
       date = format(date, "%Y-%m-%d"), #convert date into year-month-day format
       julian_date = yday(date)) #convert into Julian date

#*********************************************************************************
#1. Adults ----


#create subset 
moult_ad <- moult_data %>%
  filter(ageF == "ad")

#*********************************************************************************
#__models ----


#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
#1.a. Moult (yes/no) ----

#set moulting = no as the reference level
moult_ad$moulting <- relevel(moult_ad$moulting, ref = "no")

#**BKA**

#beta regression

#data transformation following Smithson & Verkuilen, 2006:
n=39 #sample size
moult_ad$killingD_adjusted <- (moult_ad$killingD * (n - 1) + 0.5) / n

#model
betareg_BKA <- betareg(killingD_adjusted ~ moulting + sex + julian_date, data = moult_ad, link = "logit")

#model output
summary(betareg_BKA)

#check model assumptions

#quantile residuals 
ggplot(data.frame(resid = residuals(betareg_BKA, type = "quantile")), aes(sample = resid)) +
  stat_qq() +
  stat_qq_line(color = "red")


#**Lysis**

#model
lm_lysis <- lm(lysis_score ~ moulting + sex + julian_date, data=moult_ad)

#model output
summary(lm_lysis)
Anova(lm_lysis, type='II')

#check model assumptions
simulate_res <- simulateResiduals(fittedModel = lm_lysis, plot = TRUE)


#**Agglutination**

#model
lm_agglut <- lm(aglutt_score ~ moulting + sex + julian_date, data=moult_ad)

#model output
summary(lm_agglut)
Anova(lm_agglut, type='II')

#check model assumptions
simulate_res <- simulateResiduals(fittedModel = lm_agglut, plot = TRUE)


#**Haptoglobin**

#model
lm_Hp <- lm(AC_hpconc_OLD ~ moulting + sex + julian_date + X405.nm, data=moult_ad) 


#model output
summary(lm_Hp)
Anova(lm_Hp, type='II')

#check model assumptions
simulate_res <- simulateResiduals(fittedModel = lm_Hp, plot = TRUE)



#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
#1.b. Moult score  ----

#0 Not started moulting yet but restricted to birds from 2 August onwards
#1 Moulting primaries 1-3 (starting to count on inner primaries)
#2 Moulting primaries 4-6
#3 Moulting primaries 7-10

#set data structure
moult_ad$moult.score <- as.numeric(moult_ad$moult.score) 


#**BKA**

#beta regression 

#data transformation following Smithson & Verkuilen, 2006:
n=39 #sample size
moult_ad$killingD_adjusted <- (moult_ad$killingD * (n - 1) + 0.5) / n

#model
betareg_BKA_nm <- betareg(killingD_adjusted ~ moult.score + sex + julian_date, data = moult_ad, link = "logit")

#model output
summary(betareg_BKA_nm)

#check model assumptions

#quantile residuals 
ggplot(data.frame(resid = residuals(betareg_BKA_nm, type = "quantile")), aes(sample = resid)) +
  stat_qq() +
  stat_qq_line(color = "red")


#**Lysis**

#model
lm_lysis_nm <- lm(lysis_score ~ moult.score + sex + julian_date, data=moult_ad)

#model output
summary(lm_lysis_nm)

#check model assumptions
simulate_res <- simulateResiduals(fittedModel = lm_lysis_nm, plot = TRUE)


#**Agglutination**

#model
lm_agglut_nm <- lm(aglutt_score ~ moult.score + sex + julian_date, data=moult_ad)

#model output
summary(lm_agglut_nm)

#check model assumptions
simulate_res <- simulateResiduals(fittedModel = lm_agglut_nm, plot = TRUE)


#**Haptoglobin**

#model
lm_Hp_nm <- lm(AC_hpconc_OLD ~ moult.score + sex + julian_date + X405.nm, data=moult_ad) 

#model output
summary(lm_Hp_nm)

#check model assumptions
simulate_res <- simulateResiduals(fittedModel = lm_Hp_nm, plot = TRUE)


#*********************************************************************************
#2. Juveniles ----

#Juveniles
moult_juv <- moult_data %>%
  filter(ageF == "juv") 

#*********************************************************************************
#__models ----


#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
#2.a. Moult score ----

#0 Not yet moulting
#1 1-20% Body moult
#2 21-80% body moult
#3 81-99% body moult
#4 Completed moult

#set data structure
moult_juv$moult.score <- as.factor(moult_juv$moult.score)

#set moult score = 0 (non moulting) as the reference level
moult_juv$moult.score <- relevel(moult_juv$moult.score, ref = "0")

#check order of levels
levels(moult_juv$moult.score)


#**BKA**

#beta regression

#data transformation following Smithson & Verkuilen, 2006:
n=29#sample size
moult_juv$killingD_adjusted <- (moult_juv$killingD * (n - 1) + 0.5) / n

#model
betareg_BKA_juv <- betareg(killingD_adjusted ~ moult.score + julian_date, data = moult_juv, link = "logit")

#model output
summary(betareg_BKA_juv)

#check model assumptions

#quantile residuals 
ggplot(data.frame(resid = residuals(betareg_BKA_juv, type = "quantile")), aes(sample = resid)) +
  stat_qq() +
  stat_qq_line(color = "red")


#**Lysis**

#model
lm_lysis_juv <- lm(lysis_score ~ moult.score + julian_date, data=moult_juv)

#model output
summary(lm_lysis_juv)
Anova(lm_lysis_juv, type='II')

#check model assumptions
sim_res <- simulateResiduals(fittedModel = lm_lysis_juv, plot = TRUE) #simulate residuals


#**Agglutination**

#model
lm_agglut_juv <- lm(aglutt_score ~ moult.score + julian_date, data=moult_juv)

#model output
summary(lm_agglut_juv)
Anova(lm_agglut_juv, type='II')


#check model assumptions
simulate_res <- simulateResiduals(fittedModel = lm_agglut_juv, plot = TRUE)


#**Haptoglobin**

#model
lm_Hp_juv <- lm(AC_hpconc_OLD ~ moult.score + julian_date + X405.nm, data=moult_juv) 

#model output
summary(lm_Hp_juv)
Anova(lm_Hp_juv, type='II')

#check model assumptions
simulate_res <- simulateResiduals(fittedModel = lm_Hp_juv, plot = TRUE)


#*********************************************************************************
#Plots for publication ----


#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
#1. Adults ----

#change structure
moult_ad$moult.score <- as.numeric(moult_ad$moult.score)

#Plots 

#Figure 1A
(plot_BKA_mscore <- ggplot(moult_ad, aes(x = moult.score, y = killingD)) +
  geom_jitter(width = 0.1, height = 0, alpha = 1, color = "black", size = 3) +
  labs(x = "", y=expression(paste("BKA (prop. ",italic("E.coli "), "killed)")))+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5),
        axis.title = element_text(size = 25), 
        axis.text = element_text(size = 25)))

#Figure 1B
(plot_lysis_mscore <- ggplot(moult_ad, aes(x = moult.score, y = lysis_score)) +
  geom_jitter(width = 0.1, height = 0, alpha = 1, color = "black", size = 3) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed", linewidth = 0.5) +
  labs(x = "", y="Lysis (titre)")+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5),
        axis.title = element_text(size = 25), 
        axis.text = element_text(size = 25)))

#Figure 1C
(plot_agglut_mscore <- ggplot(moult_ad, aes(x = moult.score, y = aglutt_score)) +
  geom_jitter(width = 0.1, height = 0, alpha = 1, color = "black", size = 3) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "solid", linewidth = 0.5) +
  labs(x = "Moult score", y="Agglutination (titre)")+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5),
        axis.title = element_text(size = 25), 
        axis.text = element_text(size = 25),
        axis.title.y = element_text(margin = margin(r = 35))))

#Figure 1D
(plot_Hp_mscore <- ggplot(moult_ad, aes(x = moult.score, y = AC_hpconc_OLD)) +
  geom_jitter(width = 0.1, height = 0, alpha = 1, color = "black", size = 3) +
  labs(x = "Moult score", y="Haptoglobin (mg/ml)")+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5),
        axis.title = element_text(size = 25), 
        axis.text = element_text(size = 25)))


#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
#2. Juveniles ----

#change structure
moult_juv$moult.score <- as.factor(moult_juv$moult.score)

#plot only means and standard errors:
tapply(moult_juv$killingD, moult_juv$moult.score, summary)
tapply(moult_juv$lysis_score, moult_juv$moult.score, summary)
tapply(moult_juv$aglutt_score, moult_juv$moult.score, summary)
tapply(moult_juv$AC_hpconc_OLD, moult_juv$moult.score, summary)

#first calculate means and SEs for each moult score
moult_summary_juv <- moult_juv %>%
  group_by(moult.score) %>%
  summarise(
    mean_killingD = mean(killingD),  
    mean_lysis = mean(lysis_score),  
    mean_agglut = mean(aglutt_score),  
    mean_Hp = mean(AC_hpconc_OLD), 
    se_killingD = sd(killingD) / sqrt(n()),
    se_lysis = sd(lysis_score) / sqrt(n()),
    se_agglut= sd(aglutt_score) / sqrt(n()),
    se_Hp = sd(AC_hpconc_OLD) / sqrt(n()))

# Plotting means and SEs together with raw values

(plot_BKA_means_juv <- ggplot(moult_summary_juv, aes(x = moult.score, y = mean_killingD)) +
  geom_jitter(data = moult_juv, aes(x = moult.score, y = killingD), 
              width = 0.1, height = 0, alpha = 0.5, color = "grey30", size = 3) +
  geom_point(size = 4, shape = 15, color = "black") +  
  geom_errorbar(aes(ymin = mean_killingD - se_killingD, ymax = mean_killingD + se_killingD), 
                width = 0.1) +  
  geom_line(aes(group = 1), linetype = "dashed", color = "black") + 
  labs(x = "", y = expression(paste("BKA (prop. ", italic("E.coli"), " killed)"))) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.title = element_text(size = 25), 
    axis.text = element_text(size = 25)))


(plot_lysis_means_juv <- ggplot(moult_summary_juv, aes(x = moult.score, y = mean_lysis)) +
  geom_jitter(data = moult_juv, aes(x = moult.score, y = lysis_score), 
              width = 0.1, height = 0, alpha = 0.5, color = "grey30", size = 3) +
  geom_point(size = 4, shape = 15, color = "black") +  
  geom_errorbar(aes(ymin = mean_lysis - se_lysis, ymax = mean_lysis + se_lysis), width = 0.1) +
  geom_line(aes(group = 1), linetype = "dashed", color = "black") + 
  labs(x = "", y = "Lysis (titre)") +
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5),
        axis.title = element_text(size = 25), 
        axis.text = element_text(size = 25)))


(plot_agglut_means_juv <- ggplot(moult_summary_juv, aes(x = moult.score, y = mean_agglut)) +
  geom_jitter(data = moult_juv, aes(x = moult.score, y = aglutt_score), 
              width = 0.1, height = 0, alpha = 0.5, color = "grey30", size = 3) +
  geom_point(size = 4, shape = 15, color = "black") +  
  geom_errorbar(aes(ymin = mean_agglut - se_agglut, ymax = mean_agglut + se_agglut), width = 0.1) +  
  geom_line(aes(group = 1), linetype = "dashed", color = "black") + 
  labs(x = "Moult Score", y = "Agglutination (titre)") +
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5),
        axis.title = element_text(size = 25), 
        axis.text = element_text(size = 25),
        axis.title.y = element_text(margin = margin(r = 36))))


(plot_Hp_means_juv <- ggplot(moult_summary_juv, aes(x = moult.score, y = mean_Hp)) +
  geom_jitter(data = moult_juv, aes(x = moult.score, y = AC_hpconc_OLD), 
              width = 0.1, height = 0, alpha = 0.5, color = "grey30", size = 3) +
  geom_point(size = 4, shape = 15, color = "black") +  
  geom_errorbar(aes(ymin = mean_Hp - se_Hp, ymax = mean_Hp + se_Hp), width = 0.1) + 
  geom_line(aes(group = 1), linetype = "dashed", color = "black") + 
  labs(x = "Moult Score", y = "Haptoglobin (mg/ml)") +
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5),
        axis.title = element_text(size = 25), 
        axis.text = element_text(size = 25)))


#*********************************************************************************
#Supporting Information ----


#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
#Figure S1 A-D ----

(plot_BKA <- ggplot(moult_ad, aes(x = moulting, y = killingD)) +
  geom_boxplot(width = 0.5, fill="lightgrey")+
  geom_jitter(width = 0.1, size = 2, colour="grey56") +
  labs(y=expression(paste("BKA (prop. ",italic("E.coli "), "killed)")), x="", title = "")+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5),
        axis.title = element_text(size = 28), 
        axis.text = element_text(size = 28)))


(plot_lysis <- ggplot(moult_ad, aes(x = moulting, y = lysis_score)) +
  geom_boxplot(width = 0.5, fill="lightgrey")+
  geom_jitter(width = 0.1, size = 2, colour="grey56") +
  labs(y="Lysis (titre)", x="", title = "")+
  theme_classic()+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5),
        axis.title = element_text(size = 28), 
        axis.text = element_text(size = 28)))


(plot_agglut <- ggplot(moult_ad, aes(x = moulting, y = aglutt_score)) +
  geom_boxplot(width = 0.5, fill="lightgrey")+
  geom_jitter(width = 0.1, size = 2, colour="grey56") +
  labs(y="Agglutination (titre)", x="Moulting", title = "")+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5),
        axis.title = element_text(size = 28), 
        axis.text = element_text(size = 28),
        axis.title.y = element_text(margin = margin(r = 40))))


(plot_Hp <- ggplot(moult_ad, aes(x = moulting, y = AC_hpconc_OLD)) +
  geom_boxplot(width = 0.5, fill="lightgrey")+
  geom_jitter(width = 0.1, size = 2, colour="grey56") +
  labs(y="Haptoglobin (mg/ml)", x="Moulting", title = "")+
  theme_classic()+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5),
        axis.title = element_text(size = 28), 
        axis.text = element_text(size = 28)))

#***********************************************************************************
#END ####
#***********************************************************************************
































