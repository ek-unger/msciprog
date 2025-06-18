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
ta<-read.csv("testAbiotic.csv")
View(ta)
ggplot(data=ta, aes(x=comp, y=densp))+
  geom_boxplot()+
  geom_point(color="red", alpha=0.5)+
  theme_classic()
#from this plot ot looks like the treatments worked, but there is one major outlier in the biotic treatment (probably due to gophers). Unclear if keeping the categories of treatment would be better or the gradients of density at this point. 

lm<-glmmTMB(densp~comp + (1|block), data=ta)
visreg(lm)
lm_drop<-glmmTMB(densp~1+ (1|block), data=ta)
anova(lm_drop, lm)
#adding or dropping /raintreat makes basically no difference. lm is the preferred model. adding the comp term is significant (p<2.2e-16). So, the abiotic/biotic treatments worked. could probably either run assessment using catigories or individual densities (which I think would be better).
#3db is the main outlier, which has been flagged in the data as gopher disturbed. 