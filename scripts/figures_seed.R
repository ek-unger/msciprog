#figures from seed data
###running basic seed lm###
#this script checks how germination and seed count changes with rain and competition#
library(dplyr)
library(tidyverse)
library(car)
library(performance)
library(glmmTMB)
library(visreg)
library(DHARMa)
library(emmeans)
library(ggeffects)
seed<-read.csv("data_no_touch/PlantGermSeedNT.csv")
View(seed)
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


###start of Rachel's tests

#definitely some bonkers outliers
ggplot(data=seed,aes(x=raintreat, y=seed, color=comp)) +
  geom_boxplot() +
  theme_classic() +
  facet_wrap(~comp)

ggplot(data=subset(seed,germ==1),aes(x=raintreat, y=seed, color=comp)) +
  geom_boxplot() +
  theme_classic() +
  facet_wrap(~comp,scales="free")

#this is the main base model - needs ziformula or else massively fails ZI formula test, which also makes biological sense given low germ rates
lm<-glmmTMB(seed ~ raintreat * comp + (1|block/br/brc), data=seed, family="poisson", ziformula=~.)

testDispersion(lm)
testUniformity(lm)
res <- simulateResiduals(lm)
plot(res)
testZeroInflation(lm)
testOutliers(lm)

#KS test failed even though dispersion is okay. The model thinks the biotic treatment is significantly less variable than the abiotic treatment. This is making me think we need nbinom2
lm<-glmmTMB(seed ~ raintreat * comp + (1|block/br/brc), data=seed, family="nbinom1",ziformula=~., dispformula =~raintreat * comp)
summary(lm) #note sig disp

#with nbinom2, passes all tests
#here checking that we haven't overparameterized the zi portion. Although lm3 is the "best", the results are qualitatively identical across models, which is nice
lm1<-glmmTMB(seed ~ raintreat * comp + (1|block/br/brc), data=seed, family="nbinom2",ziformula=~1)
lm2<-glmmTMB(seed ~ raintreat * comp + (1|block/br/brc), data=seed, family="nbinom2",ziformula=~comp)
lm3<-glmmTMB(seed ~ raintreat * comp + (1|block/br/brc), data=seed, family="nbinom2",ziformula=~raintreat)
lm4<-glmmTMB(seed ~ raintreat * comp + (1|block/br/brc), data=seed, family="nbinom2",ziformula=~raintreat + comp)
lm5<-glmmTMB(seed ~ raintreat * comp + (1|block/br/brc), data=seed, family="nbinom2",ziformula=~raintreat * comp)
lm6<-glmmTMB(seed ~ raintreat * comp + (1|block/br/brc), data=seed, family="nbinom2",ziformula=~.)

anova(lm1,lm2,lm3,lm4,lm5,lm6)

lm<-lm3

#all looks good
testDispersion(lm)
testUniformity(lm)
res <- simulateResiduals(lm)
plot(res)
testZeroInflation(lm)
testOutliers(lm)

summary(lm)
Anova(lm,type=3)
Anova(lm,type=2)

#quick outlier sensitivity check
outliers<-outliers(lm, lowerQuantile = 0, upperQuantile = 1, return = c("index", "logical"))

lm_outliers<-glmmTMB(seed ~ raintreat * comp + (1|block/br/brc), data=seed[-outliers,], family="nbinom2",ziformula=~raintreat)
summary(lm_outliers)
Anova(lm_outliers,type=3)
Anova(lm_outliers,type=2)

#count only
vis <- ggpredict(lm,terms=c("raintreat","comp"), type="count")
plot(vis)

#zi only
vis <- ggpredict(lm,terms=c("raintreat"), type="zi_prob")
plot(vis)

##end of Rachel

##EK repeat Rachel with summed seed count data
##seed count##
seedcomb<-seed%>%group_by(block, raintreat, comp, br, brc, germ, tubeid)%>%summarise(allseed= sum(seed))
seedcomb$block<-as.factor(seedcomb$block)
seedcomb$comp<-factor(seedcomb$comp, levels=c("B", "A"))

#definitely some bonkers outliers
ggplot(data=seedcomb,aes(x=raintreat, y=allseed, color=comp)) +
  geom_boxplot() +
  theme_classic() +
  facet_wrap(~comp)

ggplot(data=subset(seedcomb,germ==1),aes(x=raintreat, y=allseed, color=comp)) +
  geom_boxplot() +
  theme_classic() +
  facet_wrap(~comp,scales="free")

#this is the main base model - needs ziformula or else massively fails ZI formula test, which also makes biological sense given low germ rates
lm<-glmmTMB(allseed ~ raintreat * comp + (1|block/br/brc), data=seedcomb, family="poisson", ziformula=~.)

testDispersion(lm)
testUniformity(lm)
res <- simulateResiduals(lm)
plot(res)
testZeroInflation(lm)
testOutliers(lm)
summary(lm)
#This is all this true: KS test failed even though dispersion is okay. The model thinks the biotic treatment is significantly less variable than the abiotic treatment. This is making me think we need nbinom2
lm<-glmmTMB(allseed ~ raintreat * comp + (1|block/br/brc), data=seedcomb, family="nbinom1",ziformula=~., dispformula =~raintreat * comp)
summary(lm) #note sig disp

#still true
#with nbinom2, passes all tests
#here checking that we haven't overparameterized the zi portion. Although lm3 is the "best", the results are qualitatively identical across models, which is nice
lm1<-glmmTMB(allseed ~ raintreat * comp + (1|block/br/brc), data=seedcomb, family="nbinom2",ziformula=~1)
lm2<-glmmTMB(allseed ~ raintreat * comp + (1|block/br/brc), data=seedcomb, family="nbinom2",ziformula=~comp)
lm3<-glmmTMB(allseed ~ raintreat * comp + (1|block/br/brc), data=seedcomb, family="nbinom2",ziformula=~raintreat)
lm4<-glmmTMB(allseed ~ raintreat * comp + (1|block/br/brc), data=seedcomb, family="nbinom2",ziformula=~raintreat + comp)
lm5<-glmmTMB(allseed ~ raintreat * comp + (1|block/br/brc), data=seedcomb, family="nbinom2",ziformula=~raintreat * comp)
lm6<-glmmTMB(allseed ~ raintreat * comp + (1|block/br/brc), data=seedcomb, family="nbinom2",ziformula=~.)

anova(lm1,lm2,lm3,lm4,lm5,lm6)

lmSum<-lm3

#all looks good, still true
testDispersion(lmSum)
testUniformity(lmSum)
res <- simulateResiduals(lmSum)
plot(res)
testZeroInflation(lmSum)
testOutliers(lmSum)

summary(lmSum) #looks like there are still more zeros in the D treatment compared to C
Anova(lmSum,type=3) #interaction ns
Anova(lmSum,type=2) #raintreat ns at 0.05799

#calculate treatment specific means
emm<-emmeans(lmSum, ~raintreat * comp, type="response")
emm
pairs(emm, adjust="tukey") #all pairs differing in comp treatment sig, all treatments with same comp treatment ns. 

vis<-ggpredict(lmSum, terms=c("raintreat","comp"),type="fe.zi")
plot(vis)
vis<-ggpredict(lmSum, terms=c("comp", "raintreat"),type="fe.zi")
plot(vis)

#quick outlier sensitivity check
outliers<-outliers(lm, lowerQuantile = 0, upperQuantile = 1, return = c("index", "logical"))

lm_outliers<-glmmTMB(allseed ~ raintreat * comp + (1|block/br/brc), data=seedcomb[-outliers,], family="nbinom2",ziformula=~raintreat) #this doesn't run for me but maybe it does for you?
summary(lm_outliers)
Anova(lm_outliers,type=3)
Anova(lm_outliers,type=2)

#count only
vis <- ggpredict(lmSum,terms=c("raintreat","comp"), type="count")
plot(vis)

#zi only
vis <- ggpredict(lmSum,terms=c("raintreat"), type="zi_prob")
plot(vis)

#RG final model
lm3<-glmmTMB(allseed ~ raintreat * comp + (1|block/br/brc), data=seedcomb, family="nbinom2",ziformula=~raintreat)
#EK final model
s3<-glmmTMB(allseed ~ raintreat * comp + (1|block/br/brc), ziformula=~., family="nbinom2", data=seedcomb)
#so it's the same I just did not check for zi's other than . and 1
Anova(s3, type=2) #the results are essentially the same, great. comp v sig, rain ns at 0.05762
##end of RG/EK code