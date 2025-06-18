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
par<-read.csv("Plot_PARslope.csv")
par<-par%>%select(!slope)

#filter out all rows with the ambient par
par<-par%>%filter(raintreat!="amb")
View(par)
par$block<-as.factor(par$block)
par$raintreat<-as.factor(par$raintreat)
xo<-glmmTMB(par~1 + (1|block), data=par)
x1<-glmmTMB(par~raintreat + (1|block), data=par)
anova(xo, x1)
#significant impact of the treatment on PAR (0.002017)

summary(x1) #looks like only the dry treatment differs significantly, not the wet treatment (confirms our hypothesis)

#try to do this better
emm<-emmeans(x1, ~raintreat) #model based estimates of means
emm
pairs<- contrast(emm, method = "trt.vs.ctrl", ref = "C", adjust = "bonferroni")
summary(pairs) #this confirms that the difference between C and W is not sig for par but D and C is (0.0035)
visreg(x1)


#now check how much par under the tarps differed from ambient PAR
par2<-read.csv("Plot_PARslope.csv")
par2<-par2%>%select(!slope)
View(par2)
par2$block<-as.factor(par2$block)
par2$raintreat<-as.factor(par2$raintreat)

x00<-glmmTMB(par~1 + (1|block), data=par2)
x11<-glmmTMB(par~raintreat + (1|block), data=par2)
anova(x00, x11)

emm2<-emmeans(x11, ~raintreat, type="response") #model based estimates of means
emm2

pairs(emm2, adjust="tukey", type="response")
pairs2<- contrast(emm2, method = "trt.vs.ctrl", ref = "amb", adjust = "bonferroni")
summary(pairs2) #this confirms that the difference between C and W is not sig for par but D and C is (0.0035)
visreg(x1)