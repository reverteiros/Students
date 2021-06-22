

library(tidyr)
library(dplyr)
library(ggpubr)
library(cowplot)
library(gridExtra)
library(ggExtra)
library(ggplot2)
library(grid)
library(vegan)
library(SpatialTools)
library(betapart)
library(PERMANOVA)
library(VennDiagram)
library(ape)
library(rgdal)
library(ggmap)
library(gstat)


setwd("C:/Users/saret/Desktop/Students")

# as I use commonly the same names of dataframes, I recommend to run each time every chunk independently, to re-run the dataand overrun the other dataframes with the same names





################### Descriptive analyses and plots between community and altitude

database1<-read.table("dataMira.txt", header=T)%>%
  dplyr::select(Alti_Category, Taxon) 
  
# this is doing a venn diagram. It is generating a png file, check it in the folder you have the working directory
venn.diagram(
  x = list(
    database1 %>% filter(Alti_Category=="High") %>% dplyr::select(Taxon) %>% unlist() , 
    database1 %>% filter(Alti_Category=="Mid") %>% dplyr::select(Taxon) %>% unlist() , 
    database1 %>% filter(Alti_Category=="Low") %>% dplyr::select(Taxon) %>% unlist()
  ),
  category.names = c("High" , "Mid" , "Low"),
  filename = 'venn.png',
  output = TRUE ,
  imagetype="png" ,
  height = 480 , 
  width = 480 , 
  resolution = 300,
  compression = "lzw",
  lwd = 1,
  col=c("#440154ff", '#21908dff', '#fde725ff'),
  fill = c(alpha("#440154ff",0.3), alpha('#21908dff',0.3), alpha('#fde725ff',0.3)),
  cex = 0.5,
  fontfamily = "sans",
  cat.cex = 0.3,
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  cat.col = c("#440154ff", '#21908dff', '#fde725ff'),
  rotation = 1
)


# graphs with aubndance and richness
database1<-read.table("dataMira.txt", header=T)%>%
  mutate(Abundance=1)%>%
  dplyr::group_by(Alti_Category, Join_ID) %>% 
  dplyr::summarize(Abundance=sum(Abundance),Richness=n_distinct(Taxon))

ggplot(data = database1, aes(x=Alti_Category, y=Abundance)) + 
  geom_boxplot(aes(fill=Alti_Category))+
  theme_classic()+
  theme(legend.position = "none") +
  labs(y = "Abundance per site")


ggplot(data = database1, aes(x=Alti_Category, y=Richness)) + 
  geom_boxplot(aes(fill=Alti_Category))+
  theme_classic()+
  theme(legend.position = "none") +
  labs(y = "Pollinator richness per site")




##########################SPATIAL AUTOCORRELATION

# to represent the variables in space

database1_spatial <-read.table("dataMira.txt", header=T)%>%
  mutate(Abundance=1) %>%
  group_by(Join_ID) %>% 
  dplyr::summarize(Latitude=mean(Latitude),Longitude=mean(Longituse),Abundance=sum(Abundance),Richness=n_distinct(Taxon)) %>% 
  mutate(X = Latitude) %>%
  mutate(Y = Longitude) %>%
  dplyr::select(-c(Join_ID,Latitude,Longitude))%>% 
  ungroup()


coords <- SpatialPoints(database1_spatial[, c("X", "Y")], proj4string = CRS("+proj=longlat"))
plots <- SpatialPointsDataFrame(coords, database1_spatial)
ddll <- spTransform(plots, CRS("+proj=longlat"))
pts <- as.data.frame(coordinates(ddll))
names(pts) <- c("lon", "lat")

print(bubble(plots, "Abundance", maxsize = 5,key.entries = 10*(1:20),col="blue"))
print(bubble(plots, "Richness", maxsize = 5,key.entries = 4*(1:20),col="blue"))




##### moran's I, for abundance and richness

database1_spatial_aurocorrelaiton <-read.table("dataMira.txt", header=T)%>%
  mutate(Abundance=1) %>%
  group_by(Join_ID) %>% 
  dplyr::summarize(Latitude=mean(Latitude),Longitude=mean(Longituse),Abundance=sum(Abundance),Richness=n_distinct(Taxon)) %>% 
  mutate(X = Latitude) %>%
  mutate(Y = Longitude) %>%
  dplyr::select(-c(Join_ID,Latitude,Longitude))%>% 
  ungroup()
  

zone.dists <- as.matrix(dist(cbind(database1_spatial_aurocorrelaiton$X, database1_spatial_aurocorrelaiton$Y)))
zone.dists.inv <- 1/zone.dists
diag(zone.dists.inv) <- 0

Moran.I(database1_spatial_aurocorrelaiton$Abundance, zone.dists.inv) #no pattern
Moran.I(database1_spatial_aurocorrelaiton$Richness, zone.dists.inv) # a bit of pattern. check


##### mantel correlogram. Spatial variation between sites

# here only 2018. you can do the same for 2019 or put them together (remove line with filter for year for both together)
database1_betadiv_with_sites <-read.table("dataMira.txt", header=T)%>%
  mutate(Abundance=1) %>%
  dplyr::filter(Year=="Year_2018") %>% #this line
  group_by(Join_ID, Taxon) %>% 
  dplyr::summarize(Abundance=sum(Abundance)) %>% 
  tidyr::spread(Taxon, Abundance) %>%
  ungroup()

database1_betadiv_with_sites[is.na(database1_betadiv_with_sites)] <- 0 # the above generates NA, transform to 0

database1_betadiv <- database1_betadiv_with_sites %>%
  dplyr::select(-Join_ID) # remove the site column, as every row is a site is enough


#calculate beta-diversity. Quantitative with bray-curtis index
bray_curtis_distance<-bray.part(database1_betadiv) 

# and calculate geographical distance
coordinates <- read.table("dataMira.txt", header=T)%>%
  dplyr::filter(Year=="Year_2018") %>% #this line
  group_by(Join_ID) %>% 
  dplyr::summarize(Latitude=mean(Latitude),Longituse=mean(Longituse))  %>%
  mutate(X = Latitude) %>%
  mutate(Y = Longituse) %>%
  dplyr::select(-c(Join_ID,Latitude,Longituse))%>% 
  ungroup()
  
coordinates_matrix <- as.matrix(coordinates)
f<-dist1(coordinates_matrix)
d.dist<-as.dist(f*100)

##### do the mantel correlogram
mite.correlog <- mantel.correlog(bray_curtis_distance$bray, D.geo=d.dist, nperm=999)
mite.correlog  
plot(mite.correlog)





############################ How elevation is affecting the beta-diversity patterns

####### calculate the qualitative beta diversity (presence/absence), where you can calculate the difference due to subsets of species or differences because of different species are present

database1_betadiv_with_sites <-read.table("dataMira.txt", header=T)%>%
  mutate(Abundance=1) %>%
  dplyr::filter(Year=="Year_2018") %>% #this line
  group_by(Join_ID, Taxon) %>% 
  dplyr::summarize(Abundance=sum(Abundance)) %>% 
  tidyr::spread(Taxon, Abundance) %>%
  ungroup()

database1_betadiv_with_sites[is.na(database1_betadiv_with_sites)] <- 0 # the above generates NA, transform to 0

database1_betadiv <- database1_betadiv_with_sites %>%
  dplyr::select(-Join_ID) 

database1_betadiv_matrix <- as.matrix(database1_betadiv)
database1_betadiv_matrix[database1_betadiv_matrix>0] <- 1 #convert the quantitative matrix into zeros and ones

# run the beta diversity
qualitative_betadiv<-beta.pair(database1_betadiv_matrix, index.family="sor") # I selected the sorensen distance, is the most commonly used for presence/absence (=qualitative) community data



altitude_vector <- read.table("dataMira.txt", header=T) %>%
  dplyr::filter(Year=="Year_2018") %>% #to make it coherent with the previous matrix
  group_by(Join_ID) %>% 
  dplyr::summarize(Altitude=mean(Altitude)) %>% 
  ungroup()

altdist <- dist(altitude_vector$Altitude) #calculate the distance in altitude in metres between sites


########## turnover of species
mantel(qualitative_betadiv$beta.sim, altdist, method = "pearson", permutations = 9999, na.rm = FALSE)

########## nestedness of species
mantel(qualitative_betadiv$beta.sne, altdist, method = "pearson", permutations = 9999, na.rm = FALSE)

# for details on what is nestedness and turnover, check in the function beta.pair
# https://cran.r-project.org/web/packages/betapart/betapart.pdf

# a mantel test is a test to do a correlation between two distance matrices. everytime you work with beta diversity almost always you have to deal with it because beta diversity is always in the form of distance matrices. since the turnover is significant but the nestednes is not, looks like the difference between altitides is due to change in the species present but not due to the communities of higher altitudes being subsets of the ones found at lower altitudes. worth exploring further






##################### permANOVA


####################### NOt working for the moment####################

database1_betadiv_with_sites_andclasses <-read.table("dataMira.txt", header=T)%>%
  mutate(Abundance=1) %>%
  dplyr::filter(Year=="Year_2018") %>% #this line
  group_by(Join_ID, Taxon,Alti_Category) %>% 
  dplyr::summarize(Abundance=sum(Abundance)) %>% 
  tidyr::spread(Taxon, Abundance) %>%
  ungroup()

database1_betadiv_with_sites_andclasses[is.na(database1_betadiv_with_sites_andclasses)] <- 0 # the above generates NA, transform to 0


database1_betadiv <- database1_betadiv_with_sites_andclasses %>%
  dplyr::select(-c(Join_ID,Alti_Category)) # remove the site column, as every row is a site is enough


#calculate beta-diversity. Quantitative with bray-curtis index
bray_curtis_distance<-bray.part(database1_betadiv) 

bray_curtis_distance$bray


PERMANOVA(Distance, group, C = NULL, Effects = NULL, nperm = 1000, seed = NULL,
          CoordPrinc = FALSE, dimens = 2, PCoA = "Standard", ProjectInd = TRUE, tol = 1e-04,
          DatosIni = TRUE)


factorvector <- as.factor(database1_betadiv_with_sites_andclasses$Alti_Category)


data(wine)
X = wine[,4:21]
X=IniTransform(database1_betadiv)
D = DistContinuous (X)
perwine=PERMANOVA(D, factorvector)
perwine



C = matrix(c(1, 1, -1, 1, 1, -1, 1, 1, 1, -1, -1, 1), nrow=3, byrow=TRUE)
rownames(C)=c("C1", "C2", "C3")
colnames(C)=levels(wine$Group)

effects=factor(c(1,2,3))
levels(effects)=c("Origin", "Year", "Interaction")
perwine2=PERMANOVA(D, wine$Group, C=C, Effects=effects, CoordPrinc = TRUE)
summary(perwine2)



###### I didn't finish it, I will send you the final script when I figure out what's wrong












