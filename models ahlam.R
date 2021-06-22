
library(lme4)
library(lmerTest)
library(ggplot2)
library(dplyr)
library(tidyr)
library(MASS)
library(moments)
library(LambertW)

surveys <- read.table("data_glm.txt",header=T) %>%
  mutate(Abundance_100 = AB_100 + 37) %>%
  mutate(Abundance_75 = AB_75 + 33) %>%
  mutate(Richness_100 = SR_100 + 3) %>%
  mutate(Richness_75 = SR_75 + 3) %>%
  mutate(AB_75gaussianized = AB_100gaussianized)

names(surveys)

library(LambertW)
AB_75gaussianized <- Gaussianize(data = surveys$AB_75)

hist(AB_100gaussianized)

hist(((surveys$AB_100))) # this is normal enough
hist(surveys$AB_75)
table(surveys$AB_100)

hist(log(surveys$Abundance_100) )
hist(surveys$SR_100) table(surveys$SR_100)
hist(abs(surveys$SR_100))
hist((surveys$SR_75))

surveys_richness <- surveys %>%
  mutate(Abundance_100 = AB_100+37) %>%

  hist(sqrt(surveys_richness$Abundance_100))
  # dplyr::filter(SR_75<5)%>%
  # dplyr::filter(SR_100<12)

######### considering all the crops together, to see the global trends
ab100 <- lmer((AB_75gaussianized)~mean_phylogenetic_distance+mean_T+prt+(1|crop)+(1|year/region), data=surveys) # I included the v
car::vif(ab100) # values have to be lower than 4
summary(ab100)
hist(resid(ab100)) # normal enough
qqnorm(residuals(ab100))

print(kurtosis(residuals(ab100))) 


ab100 <- lmer((AB_75))~mean_phylogenetic_distance+mean_T+prt+(1|crop)+(1|year/region), data=surveys) # I included the v
car::vif(ab100) # values have to be lower than 4
summary(ab100)
hist(resid(ab100)) # normal enough
qqnorm(residuals(ab100))

ab100 <- lmer(((AB_100))~mean_phylogenetic_distance+mean_T+prt+(1|crop)+(1|year/region), data=surveys) # I included the v
car::vif(ab100) # values have to be lower than 4
summary(ab100)
hist(resid(ab100)) # normal enough


ggplot(surveys, aes(x=mean_phylogenetic_distance,y=AB_100)) +
  geom_point()+
  stat_smooth(method = "lm", col = "red")+#automatic, not exactly your values
  geom_abline(aes(intercept=-35, slope=0.36),colour="black")+#exact, but ugly
  theme_classic()+
  labs(x = "Phylogenetic distance", y = "AB_100")


########### considering the crops independently, to see if the pattern varies between species
ab100 <- lmer(AB_100~mean_phylogenetic_distance+mean_T+prt+(1|region)+(1|year)+mean_phylogenetic_distance*crop, data=surveys)
car::vif(ab100) # values have to be lower than 4
hist(resid(ab100)) # totally normal
summary(ab100)
## there is a problem about the model, I need to check it. 

ggplot(surveys, aes(x=mean_phylogenetic_distance)) +
  geom_point(aes(y=AB_100))+
  geom_abline(aes(intercept=53, slope=-0.61, species="overall"),colour="black")+ #cucurbita maxima
  geom_abline(aes(intercept=-46, slope=0.28, species="c.pepo"),colour="blue")+ # cucurbita pepo
  geom_abline(aes(intercept=-115, slope=0.97, species="m.domestica"),colour="red")+ # malus domestica
  geom_abline(aes(intercept=-2907, slope=26, species="s.lycopersicom"),colour="yellow")+ # solanum lycopersicom
  geom_abline(aes(intercept=-307, slope=3.2, species="s.melogena"),colour="orange")+ # solanum melogena
  geom_abline(aes(intercept=-51, slope=0.51, species="v.faba"),colour="green")+ # vicia faba
  theme_classic()+
  labs(x = "Phylogenetic distance", y = "AB_100")

