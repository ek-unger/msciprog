#check gopher impacts (rough)
library(dplyr)
library(tidyverse)
library(car)
library(performance)
library(glmmTMB)
library(visreg)
library(emmeans)
library(ggeffects)
library(DHARMa)
goph<-read.csv("data/plotsmdenspgopherCSV.csv")
goph$block<-as.factor(goph$block)
goph$raintreat<-as.factor(goph$raintreat)
goph$comp<-as.factor(goph$comp)
goph$gopher<-as.factor(goph$gopher)
View(goph)

#plots
#basic boxplot
ggplot(data=goph, aes(x=comp, y=densp))+
  geom_boxplot(fill=NA, color="black")+
  labs(x = "Neighbourhood treatment",y = "Neighbourhood density (individuals)")+ 
  theme_classic()

#add data points coded by disturbance level
ggplot(data = goph, aes(x = comp, y = densp)) +
  geom_boxplot(fill = NA, color = "black") +
  geom_point(aes(color = percov)) +
  scale_color_gradient(low = "cyan3", high = "darkred") +
  labs(x = "Neighbourhood treatment",y = "Neighbourhood density (individuals)",color = "Gopher disturbance\n(% cover)")+ 
  theme_classic()

ggsave(filename = "gopher_disturbance.tiff", width = 7, height = 5, device='tiff', dpi=700)

#jitter data points
ggplot(data = goph, aes(x = comp, y = densp)) +
  geom_boxplot(fill = NA, color = "black") +
  geom_point(aes(color = percov), position = position_jitter(width = 0.07, height = 0)) +
  scale_color_gradient(low = "cyan3", high = "darkred") +
  labs(x = "Neighbourhood treatment",y = "Neighbourhood density (individuals)",color = "Gopher disturbance\n(% cover)")+ 
  theme_classic()



#run some lms. my data is not normal, but is positive and skewed, so using gamma with a log link
hist(goph$densp, breaks = 20)
dg<-glmmTMB(densp~percov + (1|block/comp), family=Gamma(link = log), data=goph)
dg0<-glmmTMB(densp~1 + (1|block/comp), family=Gamma(link = log), data=goph)
anova(dg0, dg, test="Chisq")

#check residuals
simres <- simulateResiduals(dg)
plot(simres)
visreg(dg)
#ok running into some issues here, probably need a zero inflation or something? but this is not a major point here. I have my graphs. leave for now. 