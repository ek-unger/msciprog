###did the rain treatment work?###
library(lme4)
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
sm<-read.csv("Plot_SoilMoistureCSV.csv")
View(sm)
sm$RainTreat<-as.factor(sm$RainTreat)

#turn SoilMoisture into a persentage
sm <- sm %>% mutate(SoilMoisture = SoilMoisture / 100)

sm<-sm%>%mutate(pr = paste0(Plot, RainTreat))

#visualize the data
hist(sm$SoilMoisture, breaks=20) #slight right tail, fairly normal distribution
ggplot(aes(x=RainTreat, y=SoilMoisture, alpha=0.05), data=sm)+
  geom_point()+
  geom_boxplot(fill=NA)+
  theme_classic()

#create lm
lmsm0<- glmmTMB(SoilMoisture ~ 1 + (1|Plot/pr), beta_family(), data=sm)
lmsm1 <- glmmTMB(SoilMoisture ~ RainTreat+Comp + (1|Plot/pr), beta_family(), data=sm)
lmsm2 <- glmmTMB(SoilMoisture ~ RainTreat + (1|Plot/pr), beta_family(), data=sm)


anova(lmsm0, lmsm2) #raintreat is sig (0.000758)
emm<-emmeans(lmsm2, ~RainTreat, type="response")
emm #C=0.0894, D=0.707, W=0.1100
contrast(emm, method = "trt.vs.ctrl", ref = "C", adjust = "bonferroni")
#both groups sig dif, D (0.0417), W (0.0448). so slightly. 
contrast(emm, method= "trt.vs.ctrl", ref="C", adjust = "tukey")
#runs the error Note: adjust = "tukey" was changed to "sidak" because "tukey" is only appropriate for one set of pairwise comparisons
TukeyHSD(emm, conf.level=.95)#cannot run this function on a glmmtmb object
pairs(emm)
emmeans(lmsm2, pairwise~RainTreat, adjust="tukey") #these seem to do basically the same things, but I trust this line more. Here C/D = 0.0543, C/W = 0.0581, and D/W < 0.0001. 
visreg(lmsm2, scale="response")
summary(lmsm2)

Anova(lmsm2,type=3)
