###logger data create summed continuius variable###
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
rms<-read.csv("RainfallMasterSum.csv")
rms$block<-as.factor(rms$block)
rms$raintreat<-as.factor(rms$raintreat)
rms$comp<-as.factor(rms$comp)
rms$brc<-as.factor(rms$brc)
sum<- rms %>%
  group_by(block, raintreat, comp, brc) %>%
  summarise(smsum = sum(sm, na.rm = TRUE))

#plot it
ggplot(sum, aes(x=raintreat, y=smsum, color=comp))+
  geom_point()+
  theme_classic()

ggplot(sum, aes(x=comp, y=smsum, color=raintreat))+
  geom_point()+
  theme_classic()

sum2<- sum %>%
  filter(!brc%in% c("USWS", "DSWS"))
sum2<-sum2%>%filter(!brc%in% c("10CA", "8DB"))
sum2$brc<-factor(sum2$brc, levels=c("2WA","1WA","4WB","1CB","10WB","6WA","2CB","1CA","10DB","8WB","6WB","8WA","2DB","10DA","2CA","4DB","8DA","1DB","2WB","1DA","2DA","6DA","4WA","6CA","6DB","4CA","6CB","10CB","1WB","4CB","8CA","8CB","4DA","10WA"))


ggplot(sum2, aes(y = smsum, x = brc, color = raintreat)) +
  geom_point() +
  scale_color_manual(values = c("C" = "grey1", "D" = "red", "W" = "blue")) +
  theme_classic()

ggplot(sum2, aes(y = smsum, x = brc, color = comp)) +
  geom_point() +
  scale_color_manual(values = c("A" = "blue", "B" = "red")) +
  theme_classic()

#there is something very strange happening with 10WA, it is a huge outlier
summary(rms %>% filter(brc == "10WA") %>% pull(sm))
summary(rms$sm)
#i wonder if there was another poor connection with the soil? what are the ethics of performing the analysis with and without the outlier? 
write.csv(sum2, "newsmvalues.csv", row.names = FALSE)