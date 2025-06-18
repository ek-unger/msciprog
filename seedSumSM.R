###running basic seed lm###
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
seed<-read.csv("PlantGermSeedNT.csv")
seed%>%count(tubepres == "N") #in this dataset, 90 tubes (15%) are not found. These are not data points for the following analysis and do not count towards non-germination, so they will now be removed from the dataset.
seed<-seed%>%filter(tubepres=="Y")
seed$brc<-as.factor(seed$brc)
seed$brc<-factor(seed$brc, levels=c("1CA","1CB","1DA","1DB","1WA","1WB","2CA","2CB","2DA","2DB","2WA","2WB","3CA","3CB","3DA","3DB","3WA","3WB","4CA","4CB","4DA","4DB","4WA","4WB","5CA","5CB","5DA","5DB","5WA","5WB","6CA","6CB","6DA","6DB","6WA","6WB","7CA","7CB","7DA","7DB","7WA","7WB","8CA","8CB","8DA","8DB","8WA","8WB","9CA","9CB","9DA","9DB","9WA","9WB","10CA","10CB","10DA","10DB","10WA","10WB"))
seed$raintreat<-as.factor(seed$raintreat)
seed$comp<-as.factor(seed$comp)
seed$br<-as.factor(seed$br)
seed$br<-factor(seed$br, levels=c("1C","1D","1W","2C","2D","2W","3C","3D","3W","4C","4D","4W","5C","5D","5W","6C","6D","6W","7C","7D","7W","8C","8D","8W","9C","9D","9W","10C","10D","10W"))
#need to turn germination into binary response where Y = 1 and N = 0
seed <- seed %>% mutate(germ = case_when(germ == "Y" ~ 1, germ == "N" ~ 0))

#load in sm data
sumsm<-read.csv("newsmvalues.csv")
sumsm<-sumsm%>%select(c(brc, smsum))

#join with key
seed2 <- left_join(sumsm, seed, by = "brc")
seed2<-seed2%>%filter(!brc %in% c("8DB", "10CA"))

#the problem with using the smsum values is that we only have values for 6/10 blocks, 36 total subplots, but actually 34, because removal of 10CA and 8DB.
##germination rate##
#check if different treatments had impact on germination rate. germ is binary so use family = binomial(link = "logit")
seedgerm2<-seed2%>%select(!c(seed, green, full, outoftube, grazed))
seedgerm2<-seedgerm%>%group_by(block, raintreat, comp, br, brc, tubeid)%>%distinct()

lm0<-glmmTMB(germ~1 + (1/block/br/brc), family = binomial(link = "logit"), data=seedgerm2)
lm1<-glmmTMB(germ~smsum + (1/block/br/brc), family = binomial(link = "logit"), data=seedgerm2)
anova(lm0,lm1, test="Chisq") #ns, p=0.0.641
visreg(lm1) #I wonder if there is a better way to visualize this?
lm2<-glmmTMB(germ~comp + (1/block/br/brc), family = binomial(link = "logit"), data=seedgerm2)
anova(lm0,lm2, teset="Chisq") #ns, p=0.9751
lm3<-glmmTMB(germ~comp*smsum + (1/block/br/brc), family = binomial(link = "logit"), data=seedgerm2)
anova(lm0,lm3, test="Chisq") #ns, p=0.9518
visreg(lm3, scale="response")

lm4<-glmmTMB(germ~comp+smsum + (1/block/br/brc), family = binomial(link = "logit"), data=seedgerm2)
anova(lm0, lm4) #ns, 0.8955
#there was no significant effect of any of the treatments on the germination rate using smsum


##seed count##
seedcomb2<-seed2%>%group_by(block, raintreat, comp, br, brc, smsum, tubeid)%>%summarise(allseed= sum(seed))
seedcomb2$block<-as.factor(seedcomb2$block)
seedcomb2$brc<-as.factor(seedcomb2$brc)
seedcomb2$brc<-factor(seedcomb2$brc, levels=c("1CA","1CB","1DA","1DB","1WA","1WB","2CA","2CB","2DA","2DB","2WA","2WB","3CA","3CB","3DA","3DB","3WA","3WB","4CA","4CB","4DA","4DB","4WA","4WB","5CA","5CB","5DA","5DB","5WA","5WB","6CA","6CB","6DA","6DB","6WA","6WB","7CA","7CB","7DA","7DB","7WA","7WB","8CA","8CB","8DA","8WA","8WB","9CA","9CB","9DA","9DB","9WA","9WB","10CB","10DA","10DB","10WA","10WB"))
seedcomb2$comp<-factor(seedcomb2$comp, levels=c("B", "A"))


#make some models. origonally this first set did not work, hence why there are so many other sets after this. but after I remembered to removed the two subplots with bad data (10CA and 8DB), the model ran. 
hist(seedcomb2$allseed, breaks=30) #defiantly need zero-inflation

s0<-glmmTMB(allseed~1+(1|block/br/brc),ziformula = ~., family=nbinom2, control=glmmTMBControl(), data=seedcomb2)
s1<-glmmTMB(allseed~smsum + (1|block/br/brc), ziformula = ~., family=nbinom2(), data=seedcomb2, na.action = na.omit)
resid_dev <- residuals(s1, type = "pearson")
dispersion <- sum(resid_dev^2) / df.residual(s1)
dispersion  #1.9

resid_dev <- residuals(s0, type = "pearson")
dispersion <- sum(resid_dev^2) / df.residual(s0)
dispersion #1.9, not good for s0 with poisson but also its the not-model

s2<-glmmTMB(allseed~comp+(1|block/br/brc),ziformula=~., family=nbinom2, data=seedcomb2)
s3<-glmmTMB(allseed~smsum*comp+(1|block/br/brc),ziformula=~., family=nbinom2, data=seedcomb2)
s4<-glmmTMB(allseed~smsum+comp+(1|block/br/brc),ziformula=~., family=nbinom2, data=seedcomb2)

#test models
anova(s0, s1, test="Chisq") #ns at 0.1828
anova(s0, s2, test="Chisq") #very sig at 4.276e-10
anova(s2, s3, test="Chisq") #ns at 0.6461
anova(s2, s4, test="Chisq") #ns at 0.4067
anova(s3, s4) #technically s3 is the best model, but its not significant
summary(s3) #bunch of NaN values again, have to try to reduce the random effects

#its just too much for the model to handle with all the nested random effects plus the complex ziformula. It only allows me to run zi~1, probably since smsum is continuous now and not catigorical (I'm assuming it was too much to handle). 
s0<-glmmTMB(allseed~1+(1|block/br/brc),ziformula = ~1, family=nbinom2, control=glmmTMBControl(), data=seedcomb2) #with less complex ziformula, no error, and dispersion at 1.9
s1<-glmmTMB(allseed~smsum + (1|block/br/brc), ziformula = ~1, family=nbinom2, control=glmmTMBControl(), data=seedcomb2) #also 1.9 dispersion, much better
resid_dev <- residuals(s0, type = "pearson")
dispersion <- sum(resid_dev^2) / df.residual(s0)
dispersion 

resid_dev <- residuals(s1, type = "pearson")
dispersion <- sum(resid_dev^2) / df.residual(s1)
dispersion

s2<-glmmTMB(allseed~comp+(1|block/br/brc),ziformula=~1, family=nbinom2, data=seedcomb2)
#still overdispersed at 2.04 but better

resid_dev <- residuals(s2, type = "pearson")
dispersion <- sum(resid_dev^2) / df.residual(s2)
dispersion

s3<-glmmTMB(allseed~smsum*comp+(1|block/br/brc),ziformula=~1, family=nbinom2, data=seedcomb2)

resid_dev <- residuals(s3, type = "pearson")
dispersion <- sum(resid_dev^2) / df.residual(s3)
dispersion #2.06

s4<-glmmTMB(allseed~smsum+comp+(1|block/br/brc),ziformula=~1, family=nbinom2, data=seedcomb2)

#test models
anova(s0, s1, test="Chisq") #ns at 0.3432
anova(s0, s2, test="Chisq") #very sig at 8.009e-11
anova(s2, s3, test="Chisq") #ns at 0.2897
anova(s2, s4, test="Chisq") #ns at 0.1811
#however, technically s4 is the best model so I will keep it
summary(s4)
summary(s2)
emm<-emmeans(s2, ~ comp, type="response") #predicted means and CIs by comp
pairs(emm, adjust="tukey") #p < 0.0001

#removing br because it's too complex and I'm not convinced it matters anyway
s00<-glmmTMB(allseed~1+(1|block/brc),ziformula = ~1, family=nbinom2, control=glmmTMBControl(), data=seedcomb2)
s01<-glmmTMB(allseed~smsum + (1|block/brc), ziformula = ~1, family=nbinom2, control=glmmTMBControl(), data=seedcomb2)
s02<-glmmTMB(allseed~comp+(1|block/brc),ziformula=~1, family=nbinom2, data=seedcomb2)
s03<-glmmTMB(allseed~smsum*comp+(1|block/brc),ziformula=~1, family=nbinom2, data=seedcomb2)
s04<-glmmTMB(allseed~smsum+comp+(1|block/brc),ziformula=~1, family=nbinom2, data=seedcomb2)


#test models
anova(s00, s01, test="Chisq") #ns at 0.3432
anova(s00, s02, test="Chisq") #very sig at 8.009e-11
anova(s02, s03, test="Chisq") #ns at 0.2897
anova(s02, s04, test="Chisq") #ns at 0.1811
#however, technically s04 is the best model so I will keep it. Revoming br did not change any of the results and allowed all the models to run. 

summary(s04)
summary(s02)
Anova(s04, type=2) #comp V SIG at <2e-16. smsum at 0.1906, not sig but worth keeping in model
Anova(s02, type=2)
emm1<-emmeans(s04, ~ comp + smsum, type="response") #predicted means and CIs by comp
pairs(emm1, adjust="tukey")
#ok need to find out if snsum has a sig impact on only, say, the A treatment in a different way (subset the data and run an lm?)
visreg(s04)

#try to better visualize
vis<-ggpredict(s04, terms=c("smsum","comp"),type="fe.zi")
plot(vis)
vis2<-ggpredict(s04, terms=c("smsum","comp"),type="zi_prob")
plot(vis2)
vis3<-ggpredict(s04, terms=c("smsum","comp"),type="fe")
plot(vis3)

#one last time, to incorperate back the ~. zi
s000<-glmmTMB(allseed~1+(1|block/brc),ziformula = ~., family=nbinom2, control=glmmTMBControl(), data=seedcomb2)
s001<-glmmTMB(allseed~smsum + (1|block/brc), ziformula = ~., family=nbinom2, control=glmmTMBControl(), data=seedcomb2)
s002<-glmmTMB(allseed~comp+(1|block/brc),ziformula=~., family=nbinom2, data=seedcomb2)
s003<-glmmTMB(allseed~smsum*comp+(1|block/brc),ziformula=~., family=nbinom2, data=seedcomb2)
s004<-glmmTMB(allseed~smsum+comp+(1|block/brc),ziformula=~., family=nbinom2, data=seedcomb2)


#test models
anova(s000, s001, test="Chisq") #ns at 0.6254
anova(s000, s002, test="Chisq") #very sig at 4.62e-10
anova(s002, s003, test="Chisq") #ns at 0.6357
anova(s002, s004, test="Chisq") #ns at 0.4022
#however, technically s04 is the best model so I will keep it. 

summary(s004)
summary(s002)
Anova(s004, type=2) #comp V SIG at <2e-16. smsum at 0.1912, not sig but worth keeping in model
Anova(s002, type=2)
emm1<-emmeans(s004, ~ comp + smsum, type="response") #predicted means and CIs by comp
pairs(emm1, adjust="tukey")
#ok need to find out if snsum has a sig impact on only, say, the A treatment in a different way (subset the data and run an lm?)
visreg(s004)

#try to better visualize
vis<-ggpredict(s004, terms=c("smsum","comp"),type="fe.zi")
plot(vis)
anova(s004,s04)
sr04<-simulateResiduals(s04)
plot(sr04) #no problems at all!
plot(residuals(s04), fitted(s04)) #looks kinda weird

sr004<-simulateResiduals(s004)
plot(sr004) #no significant problems but the lines def look worse
plot(residuals(s004), fitted(s004)) #also looks kinda weird
#feel better about using the ~1 zi bc it fits better and also Im not conviced the treatments impacted the 0 - they did not impact germination, which is probably a large portion of the zeros... #s04 selected as final base model. 

#I do not think this will matter, but for due diligence, a new clean df with no grazed or out of tube plants
seed3<-seed2%>%filter(!grazed%in%"Y")
seed3<-seed3%>%filter(!outoftube%in%"Y")
seedcomb3<-seed3%>%group_by(block, raintreat, comp, br, brc, smsum, tubeid)%>%summarise(allseed= sum(seed))
seedcomb3$block<-as.factor(seedcomb3$block)
seedcomb3$brc<-as.factor(seedcomb3$brc)
seedcomb3$brc<-factor(seedcomb3$brc, levels=c("1CA","1CB","1DA","1DB","1WA","1WB","2CA","2CB","2DA","2DB","2WA","2WB","3CA","3CB","3DA","3DB","3WA","3WB","4CA","4CB","4DA","4DB","4WA","4WB","5CA","5CB","5DA","5DB","5WA","5WB","6CA","6CB","6DA","6DB","6WA","6WB","7CA","7CB","7DA","7DB","7WA","7WB","8CA","8CB","8DA","8WA","8WB","9CA","9CB","9DA","9DB","9WA","9WB","10CB","10DA","10DB","10WA","10WB"))
seedcomb3$comp<-factor(seedcomb3$comp, levels=c("B", "A"))

s30<-glmmTMB(allseed~1+(1|block/brc),ziformula=~1, family=nbinom2, data=seedcomb3)
s34<-glmmTMB(allseed~smsum+comp+(1|block/brc),ziformula=~1, family=nbinom2, data=seedcomb3)
s32<-glmmTMB(allseed~comp+(1|block/brc),ziformula=~1, family=nbinom2, data=seedcomb3)
anova(s30, s34) #p=2.287e-10
anova(s32, s34) #p = 0.1598
Anova(s34, type=2) #comp <2e-16, smsum 0.1572
sr34<-simulateResiduals(s34)
plot(sr34) #no problems at all!
plot(residuals(s34), fitted(s34)) #looks a bit weird
#anyways, it does not matter. the rain treatment did not significantly impact the total seed count, but the competition had a very significant impact on the total seed count. 


#it looked like there might be a trend for just the A
seedA<-seedcomb3%>%filter(comp=="A")
sA0<-glmmTMB(allseed~1+(1|block/brc),ziformula = ~1, family=nbinom2, control=glmmTMBControl(), data=seedA)
sA1<-glmmTMB(allseed~smsum+(1|block/brc),ziformula = ~1, family=nbinom2, control=glmmTMBControl(), data=seedA)
anova(sA0, sA1)
Anova(sA1, type=2)
visreg(sA1)
vis<-ggpredict(sA1, terms=c("smsum"),type="fe.zi")
plot(vis) # there is a trend, but it is NS

#exclude 10WA bc it is a huge outlier
seedcomb3<-seedcomb3%>%filter(!brc%in%"10WA")
s30<-glmmTMB(allseed~1+(1|block/brc),ziformula=~1, family=nbinom2, data=seedcomb3)
s34<-glmmTMB(allseed~smsum+comp+(1|block/brc),ziformula=~1, family=nbinom2, data=seedcomb3)
s32<-glmmTMB(allseed~comp+(1|block/brc),ziformula=~1, family=nbinom2, data=seedcomb3)
s33<-glmmTMB(allseed~smsum+(1|block/brc),ziformula=~1, family=nbinom2, data=seedcomb3)
anova(s30, s33) #NS
anova(s30, s34) #p=9.149e-10
anova(s32, s34) #p = 0.1956
Anova(s34, type=2) #comp <2e-16, smsum 0.1784
#the outlier does not make a difference, it can stay in

#next steps: split up into dif seed types, look for different ways to quanitify soil moisture as a continous variable, look into the beverton-holt seed model (or whatever it is) and figure out how to do that in R. look into seasonality of rainfall. is there a period of time when the sm would matter the most. should we focus on sm from that time period? 
