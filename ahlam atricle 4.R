setwd("C:/Users/saret/Desktop/Students")


library(dplyr)
library(tidyr)
library(ggpubr)
library(cowplot)
library(gridExtra)
library(ggExtra)
library(ggplot2)
library(grid)


dataset <- read.table("Article 4_database.txt",header=T) %>%
  filter(Identified == "Yes") %>%
  group_by(crop_trial, trial, type_plant) %>% 
  dplyr::summarize(Abundance=sum(Abundance),Pollinator_richness=n_distinct(po_species))



ggplot(data = dataset, aes(x=type_plant, y=Abundance,fill=type_plant)) + 
  geom_boxplot(aes(fill=type_plant))+
  theme_classic()+
  facet_grid(. ~ crop_trial)+
  theme(legend.position = "none") +
  labs(y = "Abundance per field")


ggplot(data = dataset, aes(x=type_plant, y=Pollinator_richness,fill=type_plant)) + 
  geom_boxplot(aes(fill=type_plant))+
  theme_classic()+
  facet_grid(. ~ crop_trial)+
  theme(legend.position = "none") +
  labs(y = "Pollinator richness per field")




library(VennDiagram)

crop <- read.table("Article 4_database.txt",header=T) %>%
  filter(type_plant == "Crop") 
MHEP <- read.table("Article 4_database.txt",header=T) %>%
  filter(type_plant == "MHEP") 
wild <- read.table("Article 4_database.txt",header=T) %>%
  filter(type_plant == "Wild") 

venn.diagram(
  x = list(crop$po_species, MHEP$po_species, wild$po_species),
  category.names = c("crop" , "MHEP " , "wild"),
  filename = '#14_venn_diagramm.png',
  output=TRUE
)






library(bipartite)

dataset <- read.table("Article 4_database.txt",header=T) %>%
  filter(Identified == "Yes") %>%
  filter(crop_trial == "Cucumis_melo") 


datasettypeplants <- read.table("Article 4_database.txt",header=T) %>%
  filter(Identified == "Yes") %>%
  filter(crop_trial == "Cucumis_melo") %>%
  group_by(crop_trial, po_species, type_plant) %>% 
  dplyr::summarize(Abundance=sum(Abundance))%>%
  filter(Abundance >= 2) 

w1<-data.frame(frame2webs(dataset, varnames= c("type_plant", "po_species","crop_trial", "Abundance"), type.out="list", emptylist="F"))

plotweb(w1, method="cca")

specieslevel(w1, index="d", level="higher", logbase=exp(1), low.abun=NULL, 
             high.abun=NULL, PDI.normalise=TRUE, PSI.beta=c(1,0), nested.method="NODF", 
             nested.normalised=TRUE, nested.weighted=TRUE, empty.web=TRUE)


