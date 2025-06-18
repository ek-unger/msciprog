#looking at seeds, different types (green/brown, empty/full). going to look both with raintreat and smsum#
#we are assuming that we harvested the green seeds before their growing season was done, and that is why they were green (not because they were never going to develop), but it would still be interesting to see where the different types of seeds are, and if this later drying off effect was only caused by the abiotic/biotic treatment or by the rain treatments as well
library(dplyr)
library(tidyverse)
library(car)
library(performance)
library(glmmTMB)
library(visreg)
library(DHARMa)
library(emmeans)
library(ggeffects)
setwd("C:/Users/emmau/OneDrive/Documents/MSc/field data/datasheets")
seed<-read.csv("PlantGermSeedNTalllabled.csv")
seed%>%count(tubepres == "N") #in this dataset, 90 tubes (15%) are not found. These are not data points for the following analysis and do not count towards non-germination, so they will now be removed from the dataset.
seed<-seed%>%filter(tubepres=="Y")
seed$brc<-as.factor(seed$brc)
seed$brc<-factor(seed$brc, levels=c("1CA","1CB","1DA","1DB","1WA","1WB","2CA","2CB","2DA","2DB","2WA","2WB","3CA","3CB","3DA","3DB","3WA","3WB","4CA","4CB","4DA","4DB","4WA","4WB","5CA","5CB","5DA","5DB","5WA","5WB","6CA","6CB","6DA","6DB","6WA","6WB","7CA","7CB","7DA","7DB","7WA","7WB","8CA","8CB","8DA","8DB","8WA","8WB","9CA","9CB","9DA","9DB","9WA","9WB","10CA","10CB","10DA","10DB","10WA","10WB"))
seed$raintreat<-as.factor(seed$raintreat)
seed$comp<-as.factor(seed$comp)
seed$br<-as.factor(seed$br)
seed$br<-factor(seed$br, levels=c("1C","1D","1W","2C","2D","2W","3C","3D","3W","4C","4D","4W","5C","5D","5W","6C","6D","6W","7C","7D","7W","8C","8D","8W","9C","9D","9W","10C","10D","10W"))
#need to turn germination into binary response where Y = 1 and N = 0
seed <- seed %>% mutate(germ = case_when(germ == "Y" ~ 1, germ == "N" ~ 0))
seed<-seed%>%mutate(green=case_when(green == "N" ~ "No", green == "Y" ~ "Yes"))
seed<-seed%>%mutate(full=case_when(full == "E" ~ "Mushy", full == "R" ~ "Full"))
seed$green<-as.factor(seed$green)
seed$full<-as.factor(seed$full)
seed<-seed%>%mutate(seedtype = paste0(green, full))
seed<-seed%>%mutate(seedtype=case_when(seedtype == "YesFull" ~ "Green Full", seedtype == "NoFull" ~ "Ripe", seedtype == "YesMushy" ~ "Green Empty"))
seed$seedtype<-as.factor(seed$seedtype)

#clean for most confident dataset
seed<-seed%>%filter(!grazed%in%"Y")
seed<-seed%>%filter(!outoftube%in%"Y")
#plants 19, 182, 220, 593 were missing values, so removed
seed<-seed%>%filter(!tubeid%in%c(19,182,220,593))

#plot the raw data
ggplot(seed, aes(x=comp, y=seed, color=seedtype))+
  geom_boxplot()+
  theme_classic()

ggplot(seed, aes(x=comp, y=seed, color=green))+
  geom_boxplot()+
  theme_classic()

ggplot(seed, aes(x=seedtype, y=seed, color=comp))+
  geom_boxplot()+
  theme_classic()

ggplot(seed, aes(x=green, y=seed, color=comp))+
  geom_boxplot()+
  theme_classic()

ggplot(seed, aes(x=raintreat, y=seed, color=seedtype))+
  geom_boxplot()+
  theme_classic()

ggplot(seed, aes(x=raintreat, y=seed, color=green))+
  geom_boxplot()+
  theme_classic()

ggplot(seed, aes(x=seedtype, y=seed, color=raintreat))+
  geom_boxplot()+
  theme_classic()

ggplot(seed, aes(x=green, y=seed, color=raintreat))+
  geom_boxplot()+
  theme_classic()
#there might be some interesting pattens here? but I genuinely cannot tell from these graphs