###estimate lambdas and growth rates for different treatments###
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
seedgerm<-seed%>%select(!c(seed, green, full, outoftube, grazed))
seedgerm<-seedgerm%>%group_by(block, raintreat, comp, br, brc, tubeid)%>%distinct()

#calculate germination rates for each subplot
germrate<-seedgerm%>%group_by(brc)%>%mutate(germrate=mean(germ))
germrate<-germrate%>%select(!c(tubeid, tubepres, germ))%>%distinct()

#calculate number of total seeds per subplot, then growth rate (base pop is number of tubes that survived)
seed<-seed%>%filter(!outoftube%in%"Y")
seed<-seed%>%filter(!grazed%in%"Y")

seedcomb<-seed%>%group_by(block, raintreat, comp, br, brc, tubeid)%>%summarise(allseed= sum(seed))
seedcomb$block<-as.factor(seedcomb$block)
seedcomb$brc<-as.factor(seedcomb$brc)
seedcomb$comp<-factor(seedcomb$comp, levels=c("B", "A"))
seedsum<-seedcomb%>%group_by(brc)%>%mutate(seedsum=sum(allseed))
seedsum<-seedsum%>%group_by(brc)%>%mutate(count=n())
seedsum<-seedsum%>%select(!c(tubeid, allseed))%>%distinct()
seedsum<-seedsum%>%group_by(brc)%>%mutate(gr=(seedsum/count))
bhprep<-left_join(germrate, seedsum, by = "brc")
bhprep<-bhprep%>%select(!c(block.y, raintreat.y, comp.y, br.y))
bhprep<-bhprep%>%rename(block=block.x, raintreat=raintreat.x, comp=comp.x, br=br.x)
write.csv(bhprep, "bhprep.csv", row.names = FALSE)
bhprep$group<-as.numeric(factor(bhprep$comp))
#try to run beverton holt? rjags requires that you pass the model to it in a string. 
data_jags<-list(N=nrow(bhprep), Ni=bhprep$count, Ni_plus1=bhprep$seedsum, group=bhprep$group, n_group = length(unique(bhprep$comp)))

model_string <- "
model {
  for (i in 1:N) {
    Ni_plus1[i] ~ dnorm(mu[i], tau) #likelihood
    mu[i] <- (r[group[i]] * Ni[i]) / (1 + alpha[group[i]] * Ni[i]) #beverton-holt formula
  }

  # Priors for treatment groups
  for (j in 1:n_group) {
    r[j] ~ dunif(0, 10) #intrinsic growth rate
    alpha[j] ~ dunif(0, 1) #competition strength, set between 0 and 1
  }
  #priors for observation error
  tau <- pow(sigma, -2)
  sigma ~ dunif(0, 10)
}
"

install.packages("rjags")
library(rjags)

model <- rjags::jags.model(
  textConnection(model_string),
  data = data_jags,
  n.chains = 3,
  n.adapt = 1000
)

update(model, 1000)  # Burn-in

samples <- coda.samples(
  model,
  variable.names = c("r", "alpha", "sigma"),
  n.iter = 5000
)

summary(samples) #alpha here is interspecific competition, so how strongly bromus is competing with everyone else (inter+intra) in the [1] and just intra in the [2]. and the r values are the intrinsic growth rate or the growth rate under low comp/ideal conditions, aka the lambda. ideally would repeat this for each rain treatment and see how they differ
plot(samples)

#need to change the dnorm in the model to a dnbinom but it threw an error... problem solve?

#try with negative binome (dnegbin in rjags. r_orbs is the overdispersion, higher = less overdispersion) i still need to figure out what all the numbers mean and if i need to change them. and add random effects
model_string <- "
model {
  for (i in 1:N) {
    # NEGATIVE BINOMIAL: overdispersed count data
    Ni_plus1[i] ~ dnegbin(p[i], r_obs)
    
    # Convert mu to probability p
    mu[i] <- (r[group[i]] * Ni[i]) / (1 + alpha[group[i]] * Ni[i])
    p[i] <- r_obs / (r_obs + mu[i])
  }

  # PRIORS for treatment-specific parameters
  for (j in 1:n_group) {
    r[j] ~ dunif(0, 20)         # Intrinsic growth rate
    alpha[j] ~ dunif(0, 1)      # Competition strength
  }

  # PRIOR for overdispersion (size parameter in negbin)
  r_obs ~ dunif(0.01, 100)
}
"
model2 <- rjags::jags.model(
  textConnection(model_string),
  data = data_jags,
  n.chains = 3,
  n.adapt = 1000
)

update(model2, 1000)  #burn-in

samples2 <- coda.samples(
  model2,
  variable.names = c("r", "alpha", "r_obs"),
  n.iter = 5000
)
summary(samples2)

#add in random intercepts (1|block/br)
model_string <- "
model {
  for (i in 1:N) {
    # Likelihood
    Ni_plus1[i] ~ dnegbin(p[i], r_obs)

    # Beverton-Holt growth with random intercepts
    mu_raw[i] <- (r[group[i]] * Ni[i]) / (1 + alpha[group[i]] * Ni[i])
    mu[i] <- mu_raw[i] * exp(u_block[block[i]] + u_br[br[i]])  # random intercepts

    # Negative binomial parametrization
    p[i] <- r_obs / (r_obs + mu[i])
  }

  # Fixed effect priors
  for (j in 1:n_group) {
    r[j] ~ dunif(0, 20)
    alpha[j] ~ dunif(0, 1)
  }

  # Overdispersion prior
  r_obs ~ dunif(0.01, 100)

  # Random effect priors
  for (b in 1:n_block) {
    u_block[b] ~ dnorm(0, tau_block)
  }
  for (s in 1:n_br) {
    u_br[s] ~ dnorm(0, tau_br)
  }

  # Hyperpriors for random effect SDs
  sigma_block ~ dunif(0, 10)
  sigma_br ~ dunif(0, 10)

  tau_block <- pow(sigma_block, -2)
  tau_br <- pow(sigma_br, -2)
}
"
data_jags2 <- list(
  N = nrow(bhprep),
  Ni = bhprep$count,
  Ni_plus1 = bhprep$seedsum,
  group = as.numeric(factor(bhprep$comp)),   # 1 or 2
  block = as.numeric(factor(bhprep$block)),       # numeric block index
  br = as.numeric(factor(bhprep$br)),             # numeric sub-block (nested in block)
  n_group = length(unique(bhprep$comp)),
  n_block = length(unique(bhprep$block)),
  n_br = length(unique(bhprep$br))
)


model3 <- jags.model(textConnection(model_string), data = data_jags2, n.chains = 3)
update(model3, 1000)  # Burn-in

samples3 <- coda.samples(model3, variable.names = c("r", "alpha", "r_obs", "sigma_block", "sigma_br"), n.iter = 5000)
summary(samples3)
plot(samples3)
png("traceplot.png", width = 800, height = 600)
plot(samples3)
dev.off()
