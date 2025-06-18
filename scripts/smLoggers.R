###did the rain treatment work? but with logger data###
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
rm<-read.csv("RainfallMaster.csv")
View(rm)
#properly set up df
colnames(rm)<-(c("block", "raintreat", "comp", "brc", "datetime", "watcont"))
rm$watcont<-as.numeric(rm$watcont) #keep in mind this now includes NA values
rm$block<-as.factor(rm$block)
rm$block<-factor(rm$block, levels=c("1","2","4","6","8","10"))
rm$raintreat<-as.factor(rm$raintreat)
structure(rm$raintreat)
rm$comp<-as.factor(rm$comp)
rm$datetime<-as.POSIXct(rm$datetime, format = "%m-%d-%Y %H:%M:%S")
rm$brc<-as.factor(rm$brc)
rm$brc<-factor(rm$brc, levels=c("1CA","1CB","1DA","1DB","1WA","1WB","2CA","2CB","2DA","2DB","2WA","2WB","4CA","4CB","4DA","4DB","4WA","4WB","6CA","6CB","6DA","6DB","6WA","6WB","8CA","8CB","8DA","8DB","8WA","8WB","10CA","10CB","10DA","10DB","10WA","10WB","UPWS","DSWS"))

#noticing some negative values and very low values in 10CA and 8DB. Could be due to poor contact with the soil, potentially due to gopher disturbance? May need to remove these logger's data. look into this.
#create new column w/o problematic loggers
rm2 <- rm %>% select(brc, watcont)
View(rm2)
colnames(rm2)<-c("brc", "watcont2")
rm2<-rm2%>%filter(!(brc %in% c("10CA", "8DB")))
rm3 <- rm %>%
  left_join(rm2, by = "brc")


#ok for some reason that did not work.
rm <- rm %>% mutate(watcont2 = ifelse(!(brc %in% c("10CA", "8DB")), watcont2, NA))

#now, rm contains all brcs, and rm2 does not (10CA and 8DB are gone. watcont2 and watcont are the same)
rm2 <- rm %>%
  mutate(watcont2 = ifelse(brc %in% c("10CA", "8DB"), NA, watcont))


ggplot(data=rm, aes(x=datetime, y=watcont, color=brc))+
  geom_abline() #did not work because of NA values

rm2 %>%
  filter(!is.na(watcont2)) %>%
  ggplot(aes(x = datetime, y = watcont, color = brc)) +
  geom_line()+
  theme_classic()

rm2 %>%
  filter(!is.na(watcont2)) %>%
  ggplot(aes(x = datetime, y = watcont, color = raintreat)) +
  geom_point()+
  theme_classic()


#removing 10CA and 8DB removes all neg values so now we can try the lm
lm<-glmmTMB(watcont2~raintreat + (1|block), beta_family(), data=rm2)
#apparently I cannot run this on my computer because I do not have enough RAM. Either try to run it on a dif computer with more RAM, or I need to summarize my values better. 

library(lubridate)

rm2_daily <- rm2 %>%
  mutate(date = as.Date(datetime)) %>%
  group_by(block, raintreat, brc, date) %>%
  summarise(daily_avg = mean(watcont2, na.rm = TRUE))
View(rm2_daily)

lm<-glmmTMB(daily_avg~raintreat + (1|date/block/brc), beta_family(), data=rm2_daily)

lm_simple<-glmmTMB(daily_avg~raintreat + (1|block/brc), beta_family(), data=rm2_daily)

lm_drop<-glmmTMB(daily_avg~1 +(1|date/block/brc), beta_family(), data=rm2_daily)
anova(lm_drop,lm)
#this model shows raintreat as significant p<2.2e-16
visreg(lm)
check_model(lm)
#ok this is bad, i think I need to use an additive model (maybe a gam?) because it is a time series. https://a-little-book-of-r-for-time-series.readthedocs.io/en/latest/src/timeseries.html

#need to remove the WS for the model
rm2_daily <- rm2_daily %>%
    filter(!brc%in% c("USWS", "DSWS"))
#there are some NA and NaN that could be messing with my model, remove them:
rm2_daily <- rm2_daily %>%
  filter(!is.na(daily_avg), !is.na(raintreat), !is.na(date), !is.na(block), !is.na(brc))


#try a gam
library(mgcv)
gam_rm <- gam(
  daily_avg ~ raintreat + s(date, bs = "cs") + 
    s(block, bs = "re") + 
    s(brc, bs = "re"),
  data = rm2_daily,
  family = betar(link = "logit"),
  method = "REML"
)

#Still running an error, remove random effects?
gam_test <- gam(
  daily_avg ~ raintreat + s(date, bs = "cs"),
  data = rm2_daily,
  family = betar(link = "logit")
)
#nope
#try gaussian instead?
gam_gaussian <- gam(
  daily_avg ~ raintreat + s(date, bs = "cs"),
  data = rm2_daily,
  family = gaussian(),
  method = "REML"
)
#nope
class(rm2_daily$date)
#ok it is a date class and still having issues

gam_test <- gam(
  daily_avg ~ raintreat + s(as.numeric(date)),
  data = rm2_daily,
  family = betar(link = "logit"),
  method = "REML"
)
#ok that ran
gam_test <- gam(
  daily_avg ~ raintreat + s(as.numeric(date), bs="cs") +
    s(block, bs = "re") + 
    s(brc, bs = "re"),
  data = rm2_daily,
  family = betar(link = "logit"),
  method = "ML"
)

gam_drop <- gam(
  daily_avg ~ 1 + s(as.numeric(date), bs="cs") +
    s(block, bs = "re") + 
    s(brc, bs = "re"),
  data = rm2_daily,
  family = betar(link = "logit"),
  method = "ML"
)
#poked around online, aparently running date data as numeric is actually fine. so this gam is ok I think? betar is for beta family, block and brc are input as random effects. also you cannot use REML to compare a fixed effect so switched to ML
anova(gam_drop, gam_test, test = "Chisq")
#alright apparently anova does not work for betar
install.packages("itsadug")
library(itsadug)
compareML(gam_drop, gam_test)
#ok so with this model raintreat is not significant. (p=0.397)

AIC(lm, gam_test)
anova(lm_simple, lm) #lm is a much better fit when it includes the crude time piece as a re. However, the gam is a better fit and this means, apparently, that raintreat is not significant in this data. 
visreg(gam_test)


