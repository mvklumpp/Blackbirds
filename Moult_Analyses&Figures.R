#*********************************************************************************
#Title----
#Stable immune function during moult regardless of age-specific moulting strategy in a European passerine

#Malin V. Klumpp & Arne Hegemann

#Script by M. V. Klumpp, January 2026
#This script provides all code necessary to reproduce the statistical analyses and figures in the manuscript and supporting information

#*********************************************************************************
rm(list = ls())
setwd("/Users/User/Desktop/Moult_BB_revision_Nov25/Dryad:GitHub/Publication") #set working directory
#*********************************************************************************
#Packages ----
#*********************************************************************************
library(car)
library(dplyr)
library(ggplot2)
library(betareg)

#*********************************************************************************
#Data  ----

#load data

moult_data = read.csv("Moult_immune_function_data.csv")

#*********************************************************************************
#1. Adults ----

#create subset 
moult_ad <- moult_data %>%
  filter(age == "ad")

#*********************************************************************************
#__models ----


#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
#1.a. Moult (yes/no) ----

#prune subset 
moult_ad <- moult_ad %>%
  mutate(moult_status = as.factor(moult_status)) %>%
  filter(!(moult_status == "completed")) #remove the single individual with completed moult as we are only looking at immune function prior to moult and during moult 

#set moult_status = no as the reference level
moult_ad$moult_status <- relevel(moult_ad$moult_status, ref = "no")

#**BKA**

#beta regression

#data transformation following Smithson & Verkuilen, 2006:
n=45 #sample size
moult_ad$BKA_prop_adjusted <- (moult_ad$BKA_prop * (n - 1) + 0.5) / n

#model
betareg_BKA <- betareg(BKA_prop_adjusted ~ moult_status + sex + julian_date, data = moult_ad, link = "logit")

#model output
summary(betareg_BKA)

#check model assumptions
par(mfrow=c(2,2))
plot(betareg_BKA)


#**Lysis**

#model
lm_lysis <- lm(lysis_score ~ moult_status + sex + julian_date, data=moult_ad)

#model output
summary(lm_lysis)
Anova(lm_lysis, type='III')

#check model assumptions

plot(lm_lysis)


#**Agglutination**

#model
lm_agglut <- lm(agglut_score ~ moult_status + sex + julian_date, data=moult_ad)

#model output
summary(lm_agglut)
Anova(lm_agglut, type='III')

#check model assumptions
plot(lm_agglut)


#**Haptoglobin**

#model
lm_Hp <- lm(hp_conc ~ moult_status + sex + julian_date + X450.nm, data=moult_ad) 

#model output
summary(lm_Hp)
Anova(lm_Hp, type='III')

#check model assumptions
plot(lm_Hp)


#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
#1.b. Moult score  ----

#prune subset 
moult_ad <- moult_data %>% #includes the single individual with completed moult again,as we are interested in how immune function varies throughout the moulting period
  filter(age == "ad") %>%
  filter(!is.na(moult_score)) #remove NAs

#set data structure as numeric (nm)
moult_ad$moult_score <- as.numeric(moult_ad$moult_score) 


#**BKA**

#beta regression 

#data transformation following Smithson & Verkuilen, 2006:
n=37 #sample size
moult_ad$BKA_prop_adjusted <- (moult_ad$BKA_prop * (n - 1) + 0.5) / n

#model
betareg_BKA_nm <- betareg(BKA_prop_adjusted ~ moult_score + sex + julian_date, data = moult_ad, link = "logit")

#model output
summary(betareg_BKA_nm)
#check model assumptions
par(mfrow=c(2,2))
plot(betareg_BKA_nm)

#**Lysis**

#model
lm_lysis_nm <- lm(lysis_score ~ moult_score + sex + julian_date, data=moult_ad)

#model output
summary(lm_lysis_nm)

#check model assumptions
plot(lm_lysis_nm)


#**Agglutination**

#model
lm_agglut_nm <- lm(agglut_score ~ moult_score + sex + julian_date, data=moult_ad)

#model output
summary(lm_agglut_nm)


#check model assumptions
plot(lm_agglut_nm)


#**Haptoglobin**

#model
lm_Hp_nm <- lm(hp_conc ~ moult_score + sex + julian_date + X450.nm, data=moult_ad) 

#model output
summary(lm_Hp_nm)

#check model assumptions
plot(lm_Hp_nm)


#*********************************************************************************
#2. Juveniles ----

#create subset 
moult_juv <- moult_data %>%
filter(age == "juv")

#*********************************************************************************
#__models ----


#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
#2.a. Moult percent ----

#set data structure
moult_juv$percent_moult <- as.numeric(moult_juv$percent_moult)

#**BKA**

#beta regression

#data transformation following Smithson & Verkuilen, 2006:
n=26#sample size
moult_juv$BKA_prop_adjusted <- (moult_juv$BKA_prop * (n - 1) + 0.5) / n

#model
betareg_BKA_juv <- betareg(BKA_prop_adjusted ~ percent_moult + julian_date, data = moult_juv, link = "logit")

#model output
summary(betareg_BKA_juv)

#check model assumptions
plot(betareg_BKA_juv)


#**Lysis**

#model
lm_lysis_juv <- lm(lysis_score ~ percent_moult + julian_date, data=moult_juv)

#model output
summary(lm_lysis_juv)

#check model assumptions
plot(lm_lysis_juv)


#**Agglutination**

#model
lm_agglut_juv <- lm(agglut_score ~ percent_moult + julian_date, data=moult_juv)

#model output
summary(lm_agglut_juv)

#check model assumptions
plot(lm_agglut_juv)


#**Haptoglobin**

#model
lm_Hp_juv <- lm(hp_conc ~ percent_moult + julian_date + X450.nm, data=moult_juv) 

#model output
summary(lm_Hp_juv)

#check model assumptions
plot(lm_Hp_juv)


#*********************************************************************************
#Plots for publication ----


#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
#1. Adults ----

#change structure
moult_ad$moult_score <- as.numeric(moult_ad$moult_score)

#Plots 

#__Figure 1 ----

#Figure 1A 
(plot_BKA_mscore <- ggplot(moult_ad, aes(x = moult_score, y = BKA_prop)) +
  geom_jitter(width = 0.1, height = 0, alpha = 1, color = "black", size = 3) +
  labs(x = "", y=expression(paste("BKA (prop. ",italic("E.coli "), "killed)")))+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5),
        axis.title = element_text(size = 25), 
        axis.text = element_text(size = 25))+
    coord_cartesian(xlim = c(0, 20)))

#Figure 1B
(plot_lysis_mscore <- ggplot(moult_ad, aes(x = moult_score, y = lysis_score)) +
  geom_jitter(width = 0.1, height = 0, alpha = 1, color = "black", size = 3) +
  labs(x = "", y="Lysis (titre)")+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5),
        axis.title = element_text(size = 25), 
        axis.text = element_text(size = 25))+
    coord_cartesian(xlim = c(0, 20)))

#Figure 1C
(plot_agglut_mscore <- ggplot(moult_ad, aes(x = moult_score, y = agglut_score)) +
  geom_jitter(width = 0.1, height = 0, alpha = 1, color = "black", size = 3) +
  labs(x = "Moult score", y="Agglutination (titre)")+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5),
        axis.title = element_text(size = 25), 
        axis.text = element_text(size = 25),
        axis.title.y = element_text(margin = margin(r = 38)))+
    coord_cartesian(xlim = c(0, 20)))

#Figure 1D
(plot_Hp_mscore <- ggplot(moult_ad, aes(x = moult_score, y = hp_conc)) +
  geom_jitter(width = 0.1, height = 0, alpha = 1, color = "black", size = 3) +
  labs(x = "Moult score", y="Haptoglobin (mg/ml)")+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5),
        axis.title = element_text(size = 25), 
        axis.text = element_text(size = 25))+
    coord_cartesian(xlim = c(0, 20)))


#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
#2. Juveniles ----

#__Figure 2 ----

(plot_BKA_percent_moult_juv <- ggplot(moult_juv, aes(x = percent_moult, y = BKA_prop)) +
   geom_jitter(width = 0.1, height = 0, alpha = 1, color = "black", size = 3) +
  labs(x = "", y = expression(paste("BKA (prop. ", italic("E.coli"), " killed)"))) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.title = element_text(size = 25), 
    axis.text = element_text(size = 25)))

(plot_lysis_percent_moult_juv <- ggplot(moult_juv, aes(x = percent_moult, y = lysis_score)) +
    geom_jitter(width = 0.1, height = 0, alpha = 1, color = "black", size = 3) +
    labs(x = "", y = "Lysis (titre)") +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.title = element_text(size = 25), 
      axis.text = element_text(size = 25)))

(plot_agglut_percent_moult_juv <- ggplot(moult_juv, aes(x = percent_moult, y = agglut_score)) +
    geom_jitter(width = 0.1, height = 0, alpha = 1, color = "black", size = 3) +
    labs(x = "Body moult (%)", y = "Agglutination (titre)") +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.title = element_text(size = 25), 
      axis.text = element_text(size = 25),
      axis.title.y = element_text(margin = margin(r = 38))))

(plot_Hp_percent_moult_juv <- ggplot(moult_juv, aes(x = percent_moult, y = hp_conc)) +
    geom_jitter(width = 0.1, height = 0, alpha = 1, color = "black", size = 3) +
    labs(x = "Body moult (%)", y = "Haptoglobin (mg/ml)") +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.title = element_text(size = 25), 
      axis.text = element_text(size = 25)))

#*********************************************************************************
#Supporting Information ----


#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
#Figure S1 A-D ----

#prune subset 
moult_ad <- moult_ad %>%
  mutate(moult_status = as.factor(moult_status)) %>%
  filter(!(moult_status == "completed")) #remove the single individual with completed moult as we are only looking at immune function prior to moult and during moult 

#set moulting = no as the reference level
moult_ad$moult_status <- relevel(moult_ad$moult_status, ref = "no")

(plot_BKA <- ggplot(moult_ad, aes(x = moult_status, y = BKA_prop)) +
  geom_boxplot(width = 0.5, fill="lightgrey")+
  geom_jitter(width = 0.1, size = 2, colour="grey56") +
  labs(y=expression(paste("BKA (prop. ",italic("E.coli "), "killed)")), x="", title = "")+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5),
        axis.title = element_text(size = 28), 
        axis.text = element_text(size = 28)))

(plot_lysis <- ggplot(moult_ad, aes(x = moult_status, y = lysis_score)) +
    geom_boxplot(width = 0.5, fill="lightgrey")+
    geom_jitter(width = 0.1, size = 2, colour="grey56") +
    labs(y="Lysis (titre)", x="", title = "") + 
    theme_classic()+
    theme(plot.title = element_text(hjust = 0.5),
          axis.title = element_text(size = 28), 
          axis.text = element_text(size = 28)))

(plot_agglut <- ggplot(moult_ad, aes(x = moult_status, y = agglut_score)) +
  geom_boxplot(width = 0.5, fill="lightgrey")+
  geom_jitter(width = 0.1, size = 2, colour="grey56") +
  labs(y="Agglutination (titre)", x="Moulting", title = "")+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5),
        axis.title = element_text(size = 28), 
        axis.text = element_text(size = 28),
        axis.title.y = element_text(margin = margin(r = 42))))

(plot_Hp <- ggplot(moult_ad, aes(x = moult_status, y = hp_conc)) +
  geom_boxplot(width = 0.5, fill="lightgrey")+
  geom_jitter(width = 0.1, size = 2, colour="grey56") +
  labs(y="Haptoglobin (mg/ml)", x="Moulting", title = "")+
  theme_classic()+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5),
        axis.title = element_text(size = 28), 
        axis.text = element_text(size = 28)))


#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
#Figure S2 A-D ----

par(mfrow = c(4,4), mar = c(2,2,2,1))

plot(betareg_BKA)
plot(lm_lysis)
plot(lm_agglut)
plot(lm_Hp)

#Figure S3 A-D ----

plot(betareg_BKA_nm)
plot(lm_lysis_nm)
plot(lm_agglut_nm)
plot(lm_Hp_nm)

#Figure S4 A-D ----

plot(betareg_BKA_juv)
plot(lm_lysis_juv)
plot(lm_agglut_juv)
plot(lm_Hp_juv)


#***********************************************************************************
#END ####
#***********************************************************************************
