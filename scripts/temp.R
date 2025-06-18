###check for effect of treatments on temp###
library(dplyr)
library(tidyverse)
library(car)
library(performance)
library(glmmTMB)
library(visreg)
library(emmeans)
library(ggeffects)
library(DHARMa)
setwd("C:/Users/emmau/OneDrive/Documents/MSc/field data/datasheets")
t<-read.csv("TempMaster.csv")
View(t)

tsum <- t %>%
  mutate(time = ceiling(row_number() / 12)) %>%
  group_by(block, raintreat, comp, brc, time) %>%
  summarise(temp_avg = mean(temp, na.rm = TRUE))
tsum

tsum <- t %>%
  mutate(time = ceiling(row_number() / 12)) %>%
  group_by(block, raintreat, comp, brc, time) %>%
  summarise(temp_avg = mean(temp, na.rm = TRUE)) %>%
  mutate(datetime = seq(
    from = as.POSIXct("2025-01-13 00:00:00"),
    by = "6 hours",
    length.out = n()
  )) %>%
  select(block, raintreat, comp, brc, datetime, temp_avg)
View(tsum) #it shifts by a hour halfway through which is weird but does not actually matter for analysis
tsum$block<-as.factor(tsum$block)
tsum$raintreat<-as.factor(tsum$raintreat)
tsum$raintreat<factor(tsum$raintreat, levels = c("AMB", "C", "D", "W"))
tsum$comp<-as.factor(tsum$comp)
tsum$brc<-as.factor(tsum$brc)
tsum$temp_avg<-as.numeric(tsum$temp_avg)
class(tsum$datetime)

#look at temp over time
ggplot(data=tsum, aes(x=datetime, y=temp_avg, color=brc))+
  geom_line()+
  theme_classic()

ggplot(data=tsum, aes(x=datetime, y=temp_avg, color = raintreat))+
  geom_line()+
  theme_classic()

ggplot(data=tsum, aes(x=datetime, y=temp_avg, color=raintreat))+
  facet_wrap( ~raintreat) +
  geom_line()+
  theme_classic()


#get the averages for all values
at1<-glmmTMB(temp_avg~raintreat +(1|datetime) + (1|block/brc), data=tsum)
emmeans(at1, ~raintreat)
emmeans(at1, pairwise~raintreat, adjust = "tukey")
summary(tsum)
#make new df with no ambients
tfree<-tsum[!tsum$brc%in% c("USWS", "DSWS"),]
View(tfree)

#looks like temp is rising in a more or less linear fashion. going to try to treat temp as a linear variable and run a lm. 
l0<-glmmTMB(temp_avg~1 + (1|block/brc), data=tfree)
l1<-glmmTMB(temp_avg~datetime + (1|block/brc), data=tfree)
anova(l0, l1, type="Chisq") #adding time is significant
l2<-glmmTMB(temp_avg~datetime+raintreat + (1|block/brc), data=tfree)
anova(l1, l2, type="Chisq") #adding raintreat without interaction is not significant. however, I expect that the time of year may impact how much the temp of the shelters would be impacted by the treatment
l3<-glmmTMB(temp_avg~(as.numeric(datetime))*raintreat + (1|block), data=tfree) #ok this does not run, I've tried pairing it down (this is the paired down version but i think it is just too much) The graphs do look pretty much the same? I could also filter the data for daylight hours and see if then their is a sig difference, as we are assuming the difference is mostly coming from shading and traping heat. However, perhaps the combo of shade and heat traped is evening out to make no sig difference. 

#ok so maybe the treatments only impact temp when the ambient temp is greater or the sun is highest in the sky. so, select for hours where most likely to find impact of treatment on temp, if no impact here, than no impact anywhere. try to select for hours 11,12,13,14.

allowed_times <- c("10:30:00", "11:00:00", "11:30:00", "12:00:00", "12:30:00","13:00:00", "13:30:00", "14:00:00")

#filter out for only times we want
t<-read.csv("TempMaster.csv")
t$block<-as.factor(t$block)
t$raintreat<-as.factor(t$raintreat)
t$comp<-as.factor(t$comp)
t$brc<-as.factor(t$brc)
t<-t[-c(1:20),]
write.csv(t, "t.csv", row.names = FALSE)

# Total number of rows
n <- nrow(t2)
n
# Create a repeating pattern: keep 8, skip 40
pattern <- rep(c(rep(TRUE, 8), rep(FALSE, 40)), length.out = n)
t2<-read.csv("t2.csv")
# Subset your data
View(t2)
t2$datetime<-as.POSIXct(t2$datetime, format = "%m-%d-%Y %H:%M:%S")
t2$datetime<-as.Date(t2$datetime)

  group_by(block, raintreat, comp, brc, time) %>%
  summarise(temp_avg = mean(temp, na.rm = TRUE)) %>%
  mutate(datetime = seq(
    from = as.Date("2025-01-13"),
    by = "1 day",
    length.out = 114
  )) %>%
  select(block, raintreat, comp, brc, datetime, temp_avg)

tsum2 <- t2 %>%
  group_by(block, raintreat, comp, brc, datetime) %>%
  summarise(temp_avg = mean(temp, na.rm = TRUE))

t2<-t2%>%group_by(datetime, brc)%>%mutate(temp_avg=mean(temp))
tsum2<-t2%>%select(!temp)
tsum2<-unique(tsum2)

view(tsum2)

tsum3<-tsum2[!tsum2$brc %in% c("USWS", "DSWS"),]

#create models for subset data
a0<-glmmTMB(temp_avg~1 +(1|datetime) + (1|block/brc), data=tsum3)
a1<-glmmTMB(temp_avg~raintreat +(1|datetime) + (1|block/brc), data=tsum3)
anova(a0, a1, type="Chisq") #still no sig effect (p=.253)
summary(a1)
emm<-emmeans(a1, ~raintreat)
#check if my residuals are normal to see if I am using the right family ie if I have to use a log link
hist(tsum3$temp_avg, breaks=20) #slight right skew, but not enough that i feel I can't use Gaussian with identity link
visreg(a1)
simres<-simulateResiduals(a1)
plot(simres)
#hmm these don't look great. will try the gamma but I assume its just because there is a lot of variance in the data due to time. ok nevermind I can't. I think this is just what it is. 
#final checks:
aveC <- mean(tsum3$temp_avg[tsum3$raintreat == "C"])
aveD <- mean(tsum3$temp_avg[tsum3$raintreat == "D"])
aveW <- mean(tsum3$temp_avg[tsum3$raintreat == "W"])
#just in case D is sig but W is not.
contrast(emm, method = "trt.vs.ctrl", ref = "C", adjust = "bonferroni") #d/c = 0.9673, w/c = 1
emmeans(a1, pairwise~raintreat, adjust = "tukey")

#repeat this for df including AMB values
a31<-glmmTMB(temp_avg~raintreat +(1|datetime) + (1|block/brc), data=tsum2)
emmeans(a31, pairwise~raintreat, adjust = "tukey") #significantly different from ambient values, not sig dif from each other
