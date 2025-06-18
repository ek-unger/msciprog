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
rm$brc<-factor(rm$brc, levels=c("1CA","1CB","1DA","1DB","1WA","1WB","2CA","2CB","2DA","2DB","2WA","2WB","4CA","4CB","4DA","4DB","4WA","4WB","6CA","6CB","6DA","6DB","6WA","6WB","8CA","8CB","8DA","8DB","8WA","8WB","10CA","10CB","10DA","10DB","10WA","10WB","DSWS", "USWS"))
rm$brc
View(rm)
setdiff(unique(rm$brc), levels(rm$brc))


#noticing some negative values and very low values in 10CA and 8DB. Could be due to poor contact with the soil, potentially due to gopher disturbance? May need to remove these logger's data. look into this.
#now, rm contains all brcs, and rm2 does not (10CA and 8DB are gone. watcont2 and watcont are the same)
rm2 <- rm %>%
  mutate(watcont2 = ifelse(brc %in% c("10CA", "8DB"), NA, watcont))
View(rm2)

#two less than amazing graphs
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

rm2_daily <- rm2 %>%
  mutate(date = as.Date(datetime)) %>%
  group_by(block, raintreat, brc, date) %>%
  summarise(daily_avg = mean(watcont2, na.rm = TRUE))
View(rm2_daily)

lm<-glmmTMB(daily_avg~raintreat + (1|date/block/brc), beta_family(), data=rm2_daily)
#do add time as re but not as nested
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

library(mgcv)
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

library(itsadug)
compareML(gam_drop, gam_test) #using this function bc apparently you cannot use anova to compare gams
#ok so with this model raintreat is not significant. (p=0.397)

AIC(lm, gam_test)
visreg(gam_test)
plot(gam_test)
summary(gam_test)
#try next: dates as julien date? plot gam_test so I can see the cubic spine against date? Add explanatory variables for rainfall, days since rainfall, etc? days since large rainfall as a decay function
#https://m-clark.github.io/generalized-additive-models/application.html


ggplot(data=rm2_daily, aes(x=raintreat, y=daily_avg))+
  geom_boxplot()+
  geom_point(color="red", alpha=0.1)+
  theme_classic()
#ok I think maybe there is just no difference? The last question that I am currious to answer is if the timing in the rain cycle matters, which would mean fitting the days since large rain event as well. there are 5 main rain events if we ignore the january data.

rl<-read.csv("rainmaclag.csv")

rm2_daily <- rm2 %>%
  mutate(date = as.Date(datetime)) %>%
  group_by(block, raintreat, brc, date) %>%
  summarise(daily_avg = mean(watcont2, na.rm = TRUE))

rl2<-rl%>%mutate(date=as.Date(time))%>%
  group_by(date)%>%
  summarise(dayrain=sum(rainfall, na.rm=TRUE))
View(rl2)
#ok now can i create a number for days since rainfall (dsr) with a weighting function for how much it rained?

#create vector of days with rain
rain_days <- rl2$dayrain > 0

#calc days since last rainfall
rl2$dsr <- NA
counter <- 0
for (i in seq_along(rain_days)) {
  if (rain_days[i]) {
    counter <- 0
  } else {
    counter <- counter + 1
  }
  rl2$dsr[i] <- counter
}


# Sum of rain over past 7 days (rolling sum, including today)
rl2$rain7day <- zoo::rollapply(rl2$dayrain, width = 7, align = "right", fill = NA, partial = TRUE, FUN = sum)

#and 5 days
rl2$rain5day <- zoo::rollapply(rl2$dayrain, width = 5, align = "right", fill = NA, partial = TRUE, FUN = sum)

# Weighted dryness index
rl2$rain_weighted <- rl2$dsr / (1+rl2$rain5day)
#bad weighting?

#now join both together,need to also make the key match (1-13 to 5-6)
rl2_del <- rl2[-c(1:12, 127), ]
View(rl2_del)

master<-rm2_daily%>%left_join(rl2_del, by="date")
master<- na.omit(master[, c("block", "raintreat", "brc", "date","daily_avg","dayrain","dsr","rain7day","rain5day", "rain_weighted")])


View(master)

#remove 8DB and 10CA again
master <- master[!master$brc %in% c("8DB", "10CA"), ]


master<-master[-c(3878:3990),]
# when i delete the last row of NA data if creates issues further up (maybe because it is running formulas?) so i am keeping it for now and will exclude it in the models

#check to make sure all my brc are still good
ggplot(data=master, aes(x=brc, y=daily_avg))+
  geom_boxplot()+
  geom_point(color="red", alpha=0.1)+
  theme_classic()

#try some models

library(MuMIn)
options(na.action = "na.fail")
full<-glmmTMB(daily_avg~(rain_weighted+dsr+raintreat+dayrain)^2 +(1|block/brc), data=master)

lm0<-glmmTMB(daily_avg~dayrain+dsr+raintreat +(1|block/brc), data=master)
lm00<-glmmTMB(daily_avg~dayrain+dsr +(1|block/brc), data=master)
lm1<-glmmTMB(daily_avg~rain_weighted*raintreat +(1|block/brc), data=master)
lm2<-glmmTMB(daily_avg~dayrain+dsr*raintreat +(1|block/brc), data=master)
lm5<-glmmTMB(daily_avg~dayrain*dsr*raintreat +(1|block/brc), data=master)
lm3<-glmmTMB(daily_avg~(rain_weighted+raintreat) +(1|block/brc), data=master)
lm4<-glmmTMB(daily_avg~(dsr+raintreat) +(1|block/brc), data=master)
AIC(lm4, lm3) #the weighted rain treatment seems a better fit for the data
AIC(lm1, lm3) #the model with the interaction term seems a better fit. I believe that the treatments effect how the days since rain impact the daily average, so this fits my hypothesis.
AIC(lm2, lm1) #Including dayrain and dsr separatly does a better job at modeling than my janky weighting system
AIC(lm2, lm3) #lm2 is still favoured

AIC(gam_test, lm2) #the gam test still is a better fit than lm2. could I incorperate the dsr and dayrain into the gam? not sure if i can even compare these two since they are using difference df at this point

visreg(lm2)
summary(lm2)
anova(lm0, lm2) #adding the interaction term to raintreat makes a significant difference
anova(lm00, lm2)
anova(lm00, lm0) #raintreat with the interaction with dsr is significant, without the interaction is not. So, it might be that the data is running as insignificant because so many days in the dataset have very low sm, or the sm is not impacted by the rainfall, because it has been too long since it last rained. Overall based on this and the data from deeper in the soil, I believe the rain treatments worked. 

#saving master as a csv
write.csv(master, "loggersmReworked.csv", row.names = FALSE)

#ok lets do this one more time: create a re for plot - which needs to be block*raintreat, not sure why that was not included. double check using the right distibution. add time as a random effect. keep days since rain. also keep dayrain. 
m<-read.csv("loggersmReworked.csv")
View(m)

#check: on days when it rains, is there an impact? because this is when the rain shelters should be doing their work. filter for dsr==0
m2<-m%>%filter(dsr==0)
View(m2)
m2<-m2%>%mutate(br = paste0(block, raintreat))
hist(m2$daily_avg, breaks=30)
lm5<-glmmTMB(daily_avg~raintreat +(1|block/br/brc) +(1|date), beta_family(), data=m2)
lmnot<-glmmTMB(daily_avg~1 +(1|block/br/brc) +(1|date), beta_family(), data=m2)

anova(lmnot, lm5, test="Chisq") #yes it matters! (p~0.0002125)
emm2<-emmeans(lm5, ~raintreat, type="response")
emm2 #C=0.229, D=0.194, W=0.246
contrast(emm2, method = "trt.vs.ctrl", ref = "C")
#D sig dif than C (0.0006), W NOT sig dif than C (.1853)
visreg(lm5)
summary(lm5)
Anova(lm5, type=3)

#now test how many days after, how do things change
m3<-m%>%filter(dsr==c(0,1))
m3<-m3%>%mutate(br = paste0(block, raintreat))
lm6<-glmmTMB(daily_avg~raintreat +(1|block/br/brc) +(1|date), beta_family(), data=m3)
lmnot3<-glmmTMB(daily_avg~1 +(1|block/br/brc) +(1|date), beta_family(), data=m3)
lm<-glmmTMB(daily_avg~raintreat+dsr +(1|block/br/brc) +(1|date), beta_family(), data=m3)
visreg(lm6)
visreg(lm)
anova(lmnot3, lm6, test="Chisq") #yes it matters! (p~0.000008406)
emm3<-emmeans(lm6, ~raintreat, type="response")
emm3 #C=0.233, D=0.183, W=0.254
contrast(emm3, method = "trt.vs.ctrl", ref = "C", adjust = "bonferroni")
#D sig <0.0001, W Not sig but more so 0.1080
View(m3)

#days 0-2
m4<-m%>%filter(dsr==c(0,1,2))
m4<-m4%>%mutate(br = paste0(block, raintreat))
lm7<-glmmTMB(daily_avg~raintreat +(1|block/br/brc) +(1|date), beta_family(), data=m4)
lmnot4<-glmmTMB(daily_avg~1 +(1|block/br/brc) +(1|date), beta_family(), data=m4)

anova(lmnot4, lm7, test="Chisq") #yes it matters! (p~0.002717)
emm4<-emmeans(lm7, ~raintreat, type="response")
emm4
contrast(emm4, method = "trt.vs.ctrl", ref = "C", adjust = "bonferroni")
#D is decreasing (0.0686), W is decreasing (.1497)

#I wonder if I move the window if it will show that the wet treatment holds onto water longer
m5<-m%>%filter(dsr==c(1,2,3))
m5<-m5%>%mutate(br = paste0(block, raintreat))
lm8<-glmmTMB(daily_avg~raintreat +(1|block/br/brc) +(1|date), beta_family(), data=m5)
lmnot5<-glmmTMB(daily_avg~1 +(1|block/br/brc) +(1|date), beta_family(), data=m5)

anova(lmnot5, lm8, test="Chisq") #yes it matters! (p~0.003041)
emm5<-emmeans(lm8, ~raintreat, type="response")
emm5 
contrast(emm5, method = "trt.vs.ctrl", ref = "C", adjust = "bonferroni") #now its D 0.0275, W 0.1501

#0:5
m6<-m%>%filter(dsr==c(0:5))
m6<-m6%>%mutate(br = paste0(block, raintreat))
lm9<-glmmTMB(daily_avg~raintreat +(1|block/br/brc) +(1|date), beta_family(), data=m6)
lmnot6<-glmmTMB(daily_avg~1 +(1|block/br/brc) +(1|date), beta_family(), data=m6)

anova(lmnot6, lm9, test="Chisq") #yes it matters! (p~0.000179)
emm6<-emmeans(lm9, ~raintreat, type="response")
emm6 
contrast(emm6, method = "trt.vs.ctrl", ref = "C", adjust = "bonferroni") #D 0.0083, W 0.0770
#0:6, #0:4,  does not have sig results. This seems to be the closest. 

#Maybe W holds water better but does not sig get more water:
m7<-m%>%filter(dsr==c(4:5))
m7<-m7%>%mutate(br = paste0(block, raintreat))
lm10<-glmmTMB(daily_avg~raintreat +(1|block/br/brc) +(1|date), beta_family(), data=m7)
lmnot7<-glmmTMB(daily_avg~1 +(1|block/br/brc) +(1|date), beta_family(), data=m7)

anova(lmnot7, lm10, test="Chisq") #yes it matters! (p~0.001421)
emm7<-emmeans(lm10, ~raintreat, type="response")
emm7#C = 0.151, D=0.134, W=0.217
contrast(emm7, method = "trt.vs.ctrl", ref = "C", adjust = "bonferroni") #D 0.6618, W 0.0029
