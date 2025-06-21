###running basic seed lm###
#this script checks how germination and seed count changes with rain and competition. also includes germination figures#
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


##germination rate##
#check if different treatments had impact on germination rate. germ is binary so use family = binomial(link = "logit")
seedgerm<-seed%>%select(!c(seed, green, full, outoftube, grazed))
seedgerm<-seedgerm%>%group_by(block, raintreat, comp, br, brc, tubeid)%>%distinct()

lm0<-glmmTMB(germ~1 + (1/block/br/brc), family = binomial(link = "logit"), data=seedgerm)
lm1<-glmmTMB(germ~raintreat + (1/block/br/brc), family = binomial(link = "logit"), data=seedgerm)
anova(lm0,lm1, test="Chisq") #ns, p=0.3029
visreg(lm1) #I wonder if there is a better way to visualize this?
lm2<-glmmTMB(germ~comp + (1/block/br/brc), family = binomial(link = "logit"), data=seedgerm)
anova(lm0,lm2, teset="Chisq") #ns, p=0.8443
lm3<-glmmTMB(germ~comp*raintreat + (1/block/br/brc), family = binomial(link = "logit"), data=seedgerm)
anova(lm0,lm3, test="Chisq") #ns, p=0.2575
emmeans(lm3, ~comp*raintreat, type="response")
Anova(lm3, type=3)
Anova(lm3, type=2)
visreg(lm3, scale="response")
vis<-ggpredict(lm3, terms=c("raintreat","comp"),type="fe")
plot(vis)
vis2<-ggpredict(lm3, terms=c("comp","raintreat"),type="fe")
plot(vis2) + labs(title= "Predicticted probability of germination", y="Germination rate", x="Competition treatment", color="Rain treatment")
ggsave(filename = "germination_predict.tiff", width = 6, height = 4, device='tiff', dpi=300)

#visualize raw data

summary_df <- seedgerm %>%
  group_by(comp, raintreat) %>%
  summarise(
    mean_val = mean(germ, na.rm = TRUE),
    se = sd(germ, na.rm = TRUE) / sqrt(n()),
    n = n(),
    .groups = "drop"
  ) %>%
  mutate(
    ci95 = se * qt(0.975, df = n - 1)  # 95% CI
  )

summary_df <- summary_df %>%
  mutate(raintreat = case_when(
    raintreat == "C" ~ "Control",
    raintreat == "D" ~ "Dry",
    raintreat == "W" ~ "Wet",
    TRUE ~ raintreat
  ))

summary_df<-summary_df%>%mutate(comp=case_when(comp=="A"~"Abiotic", comp=="B"~"Biotic"))

# Plot with error bars
ggplot(summary_df, aes(x = comp, y = mean_val, color = raintreat, group = raintreat)) +
  geom_point(position = position_dodge(width = 0.3), size = 3) +
  geom_errorbar(
    aes(ymin = mean_val - ci95, ymax = mean_val + ci95),
    position = position_dodge(width = 0.3),
    width = 0.2
  ) +
  scale_color_manual(values=c("Control"="brown4", "Dry"="darkorange2", "Wet"="blue2"))+
  labs(
    x = "Competition treatment",
    y = "Germination rate",
    color = "Rain treatment"
  ) +
  theme_classic()
ggsave(filename = "germination_raw.tiff", width = 6, height = 4, device='tiff', dpi=300)


summary(lm3)

emm_options(rg.limit = 10800)  # increase temporarily
emm<-emmeans(lm3, ~ comp * raintreat, type="response") #predicted means and CIs by comp and raintreat
emm #gives probabilities for each combination of treatments
pairs(emm, adjust="tukey") #none of the contrasts are sig, as expected
#in summary, the treatments did not sig impact germination rate, although we can see some interesting non-significant patterns (particularly the difference in wet treatment germination in A vs B comp)




#make some models. here, poisson had some convergence issues in some cases with the full set of nested random effects, so I removed br to make it run because I thought it overlaped with the rain treatment and would be the best re to drop. I later switched to a nbinom because of overdispersion, which solves this issue
s0<-glmmTMB(allseed~1+(1|block/br/brc),ziformula = ~., family="poisson", data=seedcomb)
s1<-glmmTMB(allseed~raintreat+(1|block/br/brc),ziformula=~., family="poisson", data=seedcomb)
s2<-glmmTMB(allseed~comp+(1|block/br/brc),ziformula=~., family="poisson", data=seedcomb)
s00<-glmmTMB(allseed~1+(1|block/br),ziformula = ~., family="poisson", data=seedcomb)
s02<-glmmTMB(allseed~comp+(1|block/br),ziformula=~., family="poisson", data=seedcomb)
s3<-glmmTMB(allseed~raintreat*comp+(1|block/br/brc),ziformula=~., family="poisson", data=seedcomb)
s03<-glmmTMB(allseed~raintreat*comp+(1|block/br),ziformula=~., family="poisson", data=seedcomb)
s4<-glmmTMB(allseed~raintreat+comp+(1|block/br/brc),ziformula=~., family="poisson", data=seedcomb)
s30<-glmmTMB(allseed~raintreat*comp+(1|block/br/brc), family="poisson", data=seedcomb)
s04<-glmmTMB(allseed~raintreat+comp+(1|block/br),ziformula=~., family="poisson", data=seedcomb)



#test models
anova(s0, s1, test="Chisq") #ns, 0.3719 including rain only does NOT matter
anova(s00, s02, test="Chisq") #VERY SIG <2.2e-16 including comps alone does matter
anova(s0, s3, test="Chisq") #sig, 1.355e-12 adding rain and comps with interaction does matter
anova(s02, s03, test="Chisq") #VERY SIG <2.2e-16 adding rain with interaction does matter
anova(s4, s3, test="Chisq") #ns 0.1453 the interaction does NOT matter
anova(s02, s04, test="Chisq") #ns 0.05278, both no brc, then adding rain no interaction does NOT matter 
anova(s30, s3, test="Chisq") #very sig, zero-inflation is important to include for this model, which makes sense for the data type and aligns with our hypotheses
summary(s3)
Anova(s3, type=3)
visreg(s3)

#try to better visualize
vis<-ggpredict(s3, terms=c("raintreat","comp"),type="fe.zi")
plot(vis)
vis1<-ggpredict(s3, terms=c("comp","raintreat"),type="fe.zi")
plot(vis1)
vis2<-ggpredict(s3, terms=c("raintreat","comp"),type="zi_prob")
plot(vis2)
vis3<-ggpredict(s3, terms=c("raintreat","comp"),type="fe")
plot(vis3)
vis4<-ggpredict(s3, terms=c("comp","raintreat"),type="fe")
plot(vis4)


#need to check for overdispersion and move to nbinom2 if needed
resid_dev <- residuals(s1, type = "pearson")
dispersion <- sum(resid_dev^2) / df.residual(s1)
dispersion #yes, there is overdispersion, switch to nbinom2

s0<-glmmTMB(allseed~1+(1|block/br/brc),ziformula = ~., family="nbinom2", data=seedcomb)
s1<-glmmTMB(allseed~raintreat+(1|block/br/brc),ziformula=~., family="nbinom2", data=seedcomb)
s2<-glmmTMB(allseed~comp+(1|block/br/brc),ziformula=~., family="nbinom2", data=seedcomb)
s3<-glmmTMB(allseed~raintreat*comp+(1|block/br/brc),ziformula=~., family="nbinom2", data=seedcomb)
s4<-glmmTMB(allseed~raintreat+comp+(1|block/br/brc),ziformula=~., family="nbinom2", data=seedcomb)



#test models
anova(s0, s1, test="Chisq") #ns, 0.3265 including rain only does NOT matter
anova(s0, s3, test="Chisq") #sig, 1.048e-12 adding rain and comps with interaction does matter
anova(s0, s2, test="Chisq") #sig, 9.319e-15 adding comps alone does matter
anova(s2, s3, test="Chisq") #ns, 0.0884 adding rain to comps with interaction does NOT matter
anova(s4, s3, test="Chisq") #ns 0.2421 the interaction does NOT matter
anova(s2, s4, test="Chisq") #ns
summary(s3) 
Anova(s3, type=3)
visreg(s3)

emm_options(rg.limit = 10800)  # increase temporarily
emm<-emmeans(s3, ~ comp * raintreat, type="response") #predicted means and CIs by comp and raintreat
pairs(emm, adjust="tukey") 
#from the s3 model, the treatments that are sig different (p=<0.0001) are:
#B C / A C
#B C / A D
#B C / A W
#A C / B D
#A C / B W
#B D / A D
#B D / A W
#A D / B W
#B W / A W
#Every single raincomp is sig different than all other raincomps with the opposite comp treatment. there is no sig difference between any raincomps of the same rain treatment
plot(emm)
vis<-ggpredict(s3, terms=c("raintreat","comp"),type="fe.zi")
plot(vis)
vis1<-ggpredict(s3, terms=c("comp","raintreat"),type="fe.zi")
plot(vis1)
vis2<-ggpredict(s3, terms=c("raintreat","comp"),type="zi_prob")
plot(vis2)
vis3<-ggpredict(s3, terms=c("raintreat","comp"),type="fe") #conditional model, not including the zero inflation portion of the model
plot(vis3)
vis4<-ggpredict(s3, terms=c("comp","raintreat"),type="fe")
plot(vis4)
