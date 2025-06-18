###did the rain treatment work? but with logger data###
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
#properly set up df
colnames(rm)<-(c("block", "raintreat", "comp", "brc", "datetime", "watcont"))
rm$watcont<-as.numeric(rm$watcont) #keep in mind this now includes NA values
rm$block<-as.factor(rm$block)
rm$block<-factor(rm$block, levels=c("1","2","4","6","8","10"))
rm$raintreat<-as.factor(rm$raintreat)
rm$comp<-as.factor(rm$comp)
rm$datetime<-as.POSIXct(rm$datetime, format = "%m-%d-%Y %H:%M:%S")
rm$brc<-as.factor(rm$brc)
rm$brc<-factor(rm$brc, levels=c("1CA","1CB","1DA","1DB","1WA","1WB","2CA","2CB","2DA","2DB","2WA","2WB","3CA","3CB","3DA","3DB","3WA","3WB", "4CA","4CB","4DA","4DB","4WA","4WB","5CA","5CB","5DA","5DB","5WA","5WB","6CA","6CB","6DA","6DB","6WA","6WB","7CA","7CB","7DA","7DB","7WA","7WB","8CA","8CB","8DA","8DB","8WA","8WB","9CA","9CB","9DA","9DB","9WA","9WB","10CA","10CB","10DA","10DB","10WA","10WB","DSWS", "USWS"))
setdiff(unique(rm$brc), levels(rm$brc))


#noticing some negative values and very low values in 10CA and 8DB. Could be due to poor contact with the soil, potentially due to gopher disturbance? May need to remove these logger's data. look into this.
#now, rm contains all brcs, and rm2 does not 
rm2<-rm%>%filter(!brc%in%c("10CA", "8DB", "DSWS", "USWS"))
rm2<-na.omit(rm2)

#rough graphs of raw data
ggplot(rm2, aes(x = datetime, y = watcont, color = brc)) +
  geom_line()+
  theme_classic()

ggplot(rm2, aes(x = datetime, y = watcont, color = raintreat)) +
  geom_point(alpha=0.5)+
  theme_classic()

#condense data into daily values
rm2_daily <- rm2 %>%
  mutate(date = as.Date(datetime)) %>%
  group_by(block, raintreat, brc, date) %>%
  summarise(daily_avg = mean(watcont, na.rm = TRUE))

ggplot(rm2_daily, aes(x = date, y = daily_avg, color = raintreat)) +
  geom_smooth(se=TRUE)+
  theme_classic()

#plot raw data from rm2_daily
ggplot(data=rm2_daily, aes(x=raintreat, y=daily_avg))+
  geom_boxplot()+
  geom_point(color="red", alpha=0.1)+
  theme_classic()
#ok I think maybe there is just no difference? The last question that I am curious to answer is if the timing in the rain cycle matters, which would mean fitting the days since large rain event as well. there are 5 main rain events if we ignore the january data.

rl<-read.csv("rainmaclag.csv")

rm2_daily <- rm2 %>%
  mutate(date = as.Date(datetime)) %>%
  group_by(block, raintreat, brc, date) %>%
  summarise(daily_avg = mean(watcont, na.rm = TRUE))

rl2<-rl%>%mutate(date=as.Date(time))%>%
  group_by(date)%>%
  summarise(dayrain=sum(rainfall, na.rm=TRUE))

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


#saving master as a csv
write.csv(master, "loggersmReworked.csv", row.names = FALSE)

#ok lets do this one more time: create a re for plot - which needs to be block*raintreat, not sure why that was not included. double check using the right distibution. add time as a random effect. keep days since rain. also keep dayrain. 
m<-read.csv("loggersmReworked.csv")
m1<-m%>%mutate(br=paste0(block,raintreat))
View(m1)
m1$block<-as.factor(m1$block)
m1$block<-factor(m1$block, levels=c("1","2","4","6","8","10"))
m1$raintreat<-as.factor(m1$raintreat)
m1$date<-as.POSIXct(m1$date)
m1$brc<-as.factor(m1$brc)
m1$brc<-factor(m1$brc, levels=c("1CA","1CB","1DA","1DB","1WA","1WB","2CA","2CB","2DA","2DB","2WA","2WB","4CA","4CB","4DA","4DB","4WA","4WB","6CA","6CB","6DA","6DB","6WA","6WB","8CA","8CB","8DA","8WA","8WB","10CB","10DA","10DB","10WA","10WB"))
m1$br<-as.factor(m1$br)

#create models for m1 data. note: when running the models with beta_family, i received an error: Warning message:
#In (function (start, objective, gradient = NULL, hessian = NULL,  :
#NA/NaN function evaluation

#it seems there was over paramatization. checking the daily_avg data shows it is normally distributed (mostly), and I have no values at 0 and 1 so it should be mostly safe to run the models with gaussian and an identity link. 
hist(m1$daily_avg, breaks = 30)

#models
lm0<-glmmTMB(daily_avg~1 + (1|block/br/brc) + (1|date), data=m1)
lm1<-glmmTMB(daily_avg~raintreat + (1|block/br/brc) + (1|date), data=m1)
lm2<-glmmTMB(daily_avg~raintreat+dsr + (1|block/br/brc) + (1|date), data=m1)
lm3<-glmmTMB(daily_avg~raintreat*dsr + (1|block/br/brc) + (1|date), data=m1) #recieve warning 
lm3$fit$convergence #value of 0, so the model did converge
lm4<-glmmTMB(daily_avg~raintreat*dsr+dayrain + (1|block/br/brc) + (1|date), data=m1)
lm04<-glmmTMB(daily_avg~raintreat*dsr+dayrain +(1|block/brc)+(1|date), beta_family(), data=m1)
#beta family only runs if I take out the br random effect
Anova(lm4, type=3) #everything is significant. in the simple model, raintreat ns, but when accounting for amount of rain per day, days since rain, and the interaction between dsr and the raintreatment, raintreat is sig (0.0003137). I am assuming this is because the model fits better and some of the other noise is cleared so the patterns can be viewed better. the interaction between raintreat and dsr is highly significant (<2.2e-16), showing that in different parts of the rain cycle, it impacts the treatments differently(which matches with my next models).
Anova(lm04, type=3) #honestly this doesn't change things much so I'll stick with my original model
anova(lm3, lm4)
anova(lm2, lm4)
anova(lm0,lm1)
anova(lm1, lm2)
visreg(lm4)
simres<-simulateResiduals(lm4)
plot(simres)
summary(lm4)
emmeans(lm4, pairwise~raintreat, adjust="tukey") #by itself, the raintreat contrasts are not significantly different
plot(residuals(lm4) ~ fitted(lm4))  #this check of the residuals shows no major problematic patterns, and is centered around 0, so looks like it should be ok to use gaussian
qqnorm(residuals(lm4)); qqline(residuals(lm4))  #cerntainly a bit of devation from the normal at the ends of the plot, but I think I will say it is good enough. 

#visualize the full model
#Create a grid of values for prediction
emm <- emmeans(lm4, ~ raintreat | dsr, 
               at = list(dsr = seq(min(m1$dsr, na.rm = TRUE), 
                                   max(m1$dsr, na.rm = TRUE), 
                                   by = 1)))

#Convert to data.frame for ggplot
emm_df <- as.data.frame(emm)

#Plot
ggplot(emm_df, aes(x = dsr, y = emmean, color = raintreat)) +
  geom_line(size = 1.2) +
  geom_ribbon(aes(ymin = lower.CL, ymax = upper.CL, fill = raintreat), alpha = 0.2, color = NA) +
  labs(x = "Days Since Rain (dsr)", y = "Predicted Daily Avg", color = "Rain Treatment", fill = "Rain Treatment") +
  theme_minimal()

##looks like there might be a difference in variance in treatment. try to calculate the CV for the different treatments so I can show this##
cv_summary <- m1 %>%
  group_by(raintreat) %>%
  summarise(
    mean_val = mean(daily_avg, na.rm = TRUE),
    sd_val = sd(daily_avg, na.rm = TRUE),
    cv = sd_val / mean_val *100
  )
#CVs: C = 44.3, D = 35.5, W = 41.8
#testing the sig is hard bc they are ratios. so will use bootstrapping to help with this
library(boot)

cv_func <- function(data, indices) {
  d <- data[indices, ]
  mean_val <- mean(d$daily_avg, na.rm = TRUE)
  sd_val <- sd(d$daily_avg, na.rm = TRUE)
  return(sd_val / mean_val)
}

#run it
boot_c_list <- m1 %>%
  group_by(raintreat) %>%
  group_map(~ boot::boot(.x, statistic = cv_func, R = 1000), .keep = TRUE)

#name the results
names(boot_c_list) <- unique(m1$raintreat)

cv_c <- boot_c_list[["C"]]$t
cv_d <- boot_c_list[["D"]]$t
cv_w <- boot_c_list[["W"]]$t

diff_cd <- cv_c - cv_d
diff_cw <- cv_c - cv_w
diff_dw <- cv_d - cv_w

quantile(diff_cd, probs = c(0.025, 0.975)) #none of these values are zero so apparently this means they are sig different? but idk will have to look into this
quantile(diff_cw, probs = c(0.025, 0.975))
quantile(diff_dw, probs = c(0.025, 0.975))

mean(abs(diff_cd) >= abs(mean(diff_cd))) #emperical p-value, 0.499
mean(abs(diff_cw) >= abs(mean(diff_cw))) #0.501
mean(abs(diff_dw) >= abs(mean(diff_dw))) #0.5


#hist to see impact of different dsr. do dif values carry dif weight?
hist(m1$dsr, breaks = seq(floor(min(0)), ceiling(max(25)), by = 1), 
     right = FALSE,  # Optional: left-closed intervals [a, b)
     main = "Histogram with Integer Bins", xlab = "Value")

emm_filt<- emm_df %>%
  select(dsr, raintreat, emmean) %>%
  tidyr::pivot_wider(names_from = raintreat, values_from = emmean)



#quantify how much time spent in dif dsr
sum(m1$dsr==0) #884 22.81%
sum(m1$dsr==1) #306 7.90%
sum(m1$dsr==2) #203 5.24%
sum(m1$dsr==3) #170 4.39%
sum(m1$dsr==4) #170 4.39%
sum(m1$dsr==5) #170 4.39%
sum(m1$dsr==6) #102 2.63%
sum(m1$dsr==7) #136 3.51%
sum(m1$dsr==8) #136 3.51%
sum(m1$dsr==9) #102 2.63%
sum(m1$dsr==10) #102 2.63%
sum(m1$dsr==11) #102 2.63%
sum(m1$dsr==12) #102 2.63%
sum(m1$dsr==13) #102 2.63%
sum(m1$dsr==14) #102 2.63%
sum(m1$dsr==15) #102 2.63%
sum(m1$dsr==16) #102 2.63%
sum(m1$dsr==17) #102 2.63%
sum(m1$dsr==18) #102 2.63%
sum(m1$dsr==19) #102 2.63%
sum(m1$dsr==20) #102 2.63%
sum(m1$dsr==21) #102 2.63%
sum(m1$dsr==22) #102 2.63%
sum(m1$dsr==23) #102 2.63%
sum(m1$dsr==24) #68 1.75%
#total = 3875

# Count how many times C > D, W > D, etc
plot(emm_df) #this plot helps count
sum(emm_filt$C > emm_filt$D, na.rm = TRUE) #D is dryer than C for the first 8 days
sum(emm_filt$W>emm_filt$D, na.rm = TRUE) #D is dryer than W for the first 13 days
sum(emm_filt$W>emm_filt$C, na.rm = TRUE) #W wetter than C all 25 days
#need sum of time 0-7, 8-24, 0-12, 13-24
#time with D<C = 55.26%
#time with C<D = 44.74%
#time with W>D = 69.29%
#time with D>W = 30.71%

#now this is not testing if all these differences are sig, to do that we have to run a model for each indiviual day and see how many there are
library(glmmTMB)
library(emmeans)
library(dplyr)

# Storage for emmeans and pairwise contrasts
emm_list <- list()
contrast_list <- list()

for (i in 0:25) {
  m_sub <- m1 %>% filter(dsr == i)
  
  # Skip if not enough data
  if (nrow(m_sub) < 10) next
  
  # Fit model
  model <- try(glmmTMB(daily_avg ~ raintreat + (1|block/br/brc) + (1|date),
                       data = m_sub,
                       family = beta_family()), silent = TRUE)
  if (inherits(model, "try-error")) next
  
  # Get emmeans
  emm <- emmeans(model, ~ raintreat)
  emm_df <- as.data.frame(emm) %>% mutate(dsr = i)
  emm_list[[as.character(i)]] <- emm_df
  
  # Tukey-adjusted pairwise comparisons
  cont <- pairs(emm, adjust = "tukey")
  cont_df <- as.data.frame(cont) %>% mutate(dsr = i)
  contrast_list[[as.character(i)]] <- cont_df
}

# Combine into single dataframes
all_emm <- bind_rows(emm_list)
all_contrasts <- bind_rows(contrast_list)
all_contrasts
all_emm

#next steps: find which ones are significant, and what to do with that information. This represents sig without the interaction term which we found to be highly significant.
#days where D is sig drier than C: 0,1: 30.71%
#days where D is sig drier than W: 0,1,2,3: 40.34%
#days where D is sig wetter than C: 16,17,18,19,20,21,22,23,24: 22.79%
#days where D is sig wetter than W: N/A 0%
#days where W is sig wetter than C: N/A 0%
#days where W is sig drier than C: N/A 0%

##filtering dataset for specific days##
#check: on days when it rains, is there an impact? because this is when the rain shelters should be doing their work. filter for dsr==0
m2<-m1%>%filter(dsr==0)
View(m2)
hist(m2$daily_avg, breaks=30) #normal distribution
lm5<-glmmTMB(daily_avg~raintreat +(1|block/br/brc) +(1|date), beta_family(), data=m2)
lmnot<-glmmTMB(daily_avg~1 +(1|block/br/brc) +(1|date), beta_family(), data=m2)

anova(lmnot, lm5, test="Chisq") #yes it matters! (p~0.0002125)
emm2<-emmeans(lm5, ~raintreat, type="response")
emm2 #C=0.229, D=0.194, W=0.246
emmeans(lm5, pairwise~raintreat, adjust="tukey") #c/d 0.0009, c/w 0.2123, d/w <0.0001
visreg(lm5)
summary(lm5)
Anova(lm5, type=3)

#now test how many days after, how do things change
m3<-m1%>%filter(dsr==c(0,1))
m3<-m3%>%mutate(br = paste0(block, raintreat))
lm6<-glmmTMB(daily_avg~raintreat +(1|block/br/brc) +(1|date), beta_family(), data=m3)
lmnot3<-glmmTMB(daily_avg~1 +(1|block/br/brc) +(1|date), beta_family(), data=m3)
lm<-glmmTMB(daily_avg~raintreat*dsr +(1|block/br/brc) +(1|date), beta_family(), data=m3)
visreg(lm6)
visreg(lm)
anova(lmnot3, lm6, test="Chisq") #yes it matters! (p~0.000008406)
anova(lm6, lm, test="Chisq") #the interaction also still matters (0.002288)
emm3<-emmeans(lm6, ~raintreat, type="response")
emm3 #C=0.233, D=0.183, W=0.254
contrast(emm3, method = "trt.vs.ctrl", ref = "C", adjust = "bonferroni")
#D sig <0.0001, W Not sig but more so 0.1080
emmeans(lm6, pairwise~raintreat, adjust="tukey") #c/d <0.0001, c/w 0.1312, d/w <0.0001
Anova(lm, type=3) #dsr is ns (its only two days so it makes sense) but the interaction is still sig, and raintreat is very sig. so which days are included is very important

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
emmeans(lm7, pairwise~raintreat, adjust="tukey") #c/d 0.0864, c/w 0.1758, d/w 0.0002

#I wonder if I move the window if it will show that the wet treatment holds onto water longer
m5<-m%>%filter(dsr==c(1,2,3))
m5<-m5%>%mutate(br = paste0(block, raintreat))
lm8<-glmmTMB(daily_avg~raintreat +(1|block/br/brc) +(1|date), beta_family(), data=m5)
lmnot5<-glmmTMB(daily_avg~1 +(1|block/br/brc) +(1|date), beta_family(), data=m5)

anova(lmnot5, lm8, test="Chisq") #yes it matters! (p~0.003041)
emm5<-emmeans(lm8, ~raintreat, type="response")
emm5 
contrast(emm5, method = "trt.vs.ctrl", ref = "C", adjust = "bonferroni") #now its D 0.0275, W 0.1501
emmeans(lm8, pairwise~raintreat, adjust="tukey") #c/d 0.0366, c/w 0.1762, D/W <0.0001

#0:5
m6<-m%>%filter(dsr==c(0:5))
m6<-m6%>%mutate(br = paste0(block, raintreat))
lm9<-glmmTMB(daily_avg~raintreat +(1|block/br/brc) +(1|date), beta_family(), data=m6)
lmnot6<-glmmTMB(daily_avg~1 +(1|block/br/brc) +(1|date), beta_family(), data=m6)

anova(lmnot6, lm9, test="Chisq") #yes it matters! (p~0.000179)
emm6<-emmeans(lm9, ~raintreat, type="response")
emm6 
contrast(emm6, method = "trt.vs.ctrl", ref = "C", adjust = "bonferroni") #D 0.0083, W 0.0770
emmeans(lm9, pairwise~raintreat, adjust="tukey") #c/d 0.0115 c/w 0.0962 d/w <0.0001
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
emmeans(lm10, pairwise~raintreat, adjust="tukey") #c/d 0.5944 #c/w 0.0041 #d/w 0.0001

#compare 2" probe to 4.5" probe for more information: input new datasheet that includes the info for the same day/time (as close as possible)
depth<-read.csv("2vs4_5smDepthComp.csv")
depth$length<-as.factor(depth$length)
depth$block<-as.factor(depth$block)
depth$block<-factor(depth$block, levels=c("1","2","3","4","5","6","7","8","9","10"))
depth$raintreat<-as.factor(depth$raintreat)
depth$brc<-as.factor(depth$brc)
depth$br<-as.factor(depth$br)
depth$comp<-as.factor(depth$comp)

ggplot(depth, aes(x=length, y=sm, color=raintreat))+
  geom_boxplot()+
  theme_classic()
#clearly there is a difference, maybe leave it here unless RG has other thoughts/ideas on what I should asses 