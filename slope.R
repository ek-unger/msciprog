###does slope vary significantly between plots? what is the general slope of the hill?###
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
slope<-read.csv("Plot_PARSlope.csv")
View(slope)
slope
s<- na.omit(slope[, c("Plot", "RainTreat", "PR", "Slope")])
View(s)
s$Plot<-as.factor(s$Plot)
s$Plot<-factor(s$Plot, levels=c("1","2","3","4","5","6","7","8","9","10"))
s$RainTreat<-as.factor(s$RainTreat)
s$PR<-as.factor(s$PR)

lm<-glmmTMB(Slope~Plot, data=s)
lmnot<-glmmTMB(Slope~1, data=s)
anova(lmnot, lm) #the Plot (block) does not have a significant effect on slope (p~0.6007, and th model is slightly better fit without Plot) (good, this is what we want)

aveslope<-round(mean(s$Slope, ma.rm=TRUE), 7)
aveslope #the average slope is 5 degrees, not bad!
visreg(lm)
summary(lm)
