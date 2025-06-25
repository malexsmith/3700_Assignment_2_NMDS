# ZOO*3700 Assignment 2 : Community richness and similarity between deep-sea vents around the world. 

# Welcome to your second invertebrateR assignments where we continue our exploration of the deep-sea thermal vent communities.  This time, you will map locations, estimate alpha and beta-diversity and ask whether these communities are different and isolated and what effect mining in international waters will have on this biodiversity. 

# As in Assignment 1, to complete this assignment you will need to be running the most up to date version of R in Windows.  There are some new packages that we will be installing to complete this work. 

# As before, this Assignment is designed to be followed along in steps from start to finish.  Each step creates an output that subsequent steps depend on.  Remember to follow the code blocks in order!

# This first block of code will clear your working environment in case you've been using R for something else in the past.

rm(list=ls())

# The next block of code will tell you where your working directory is (i.e. any figures you generate you can find in that folder).  The 'wd' is also where to put your input files.

getwd()

# If you want to change the 'wd' - use the optional setwd command below, and replace the ____ between the "_____" with where you'd like to save your material (and then remove the '#')

# setwd("_____")
setwd("C:\\Users\\malex\\Sync\\R")

data2 <- read.csv(file = "3700 bigger test nmds community.csv",head=TRUE,row.names = 1, sep=",")
metadata2 <- read.csv(file = "3700 bigger test nmds metadata.csv",head=TRUE,row.names = 1, sep=",")

# This next block of code will install a package called BioManager that we will use to install 10 further packages (if you have not already installed them). If you have already installed them, you can skip to the library() commands below which will open the packages you need to complete this assignment. 

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager");

BiocManager::install("vegan");
BiocManager::install("ggplot2");
BiocManager::install("bipartite");
BiocManager::install("ggordiplots");
BiocManager::install("sf");
BiocManager::install("mapview");
BiocManager::install("viridis");
BiocManager::install("webshot");
BiocManager::install("pagedown");
BiocManager::install("reshape2");


# The next block are the library() commands which will 10 the packages you need to complete this assignment. 

library(vegan); 
library(ggplot2);
library(bipartite);
library(ggordiplots)
library(sf)
library(mapview)
library(viridis)
library(webshot)
library(pagedown)
library(reshape2)

# this next block of code will upload incidence matrices (presence, or 1's, and absence, or 0's) for the thermal vents were we have data.  In addition, the metadata files will add details about the sites (including their location, region, area, chemistry etc.)

data <- read.csv(file = "3700 test nmds community.csv",head=TRUE,row.names = 1, sep=",")
metadata <- read.csv(file = "3700 test nmds metadata.csv",head=TRUE,row.names = 1, sep=",")
taxon <- read.csv(file = "3700 test nmds taxon metadata.csv",head=TRUE,row.names = 1, sep=",")
metadata2 <- read.csv(file = "3700 bigger test nmds metadata.csv",head=TRUE,row.names = 1, sep=",")

## First, lets make a map of the locations using the package mapview. 
## this map is dynamic and you can move around and plot against several background layers
## https://map-rfun.library.duke.edu/01_georeference.html

colors = viridis(7)

ThermalVents <- st_as_sf(metadata, coords = c("Long", "Lat"),  crs = 4326)
mapview(ThermalVents, zcol = "Region")
mapview(ThermalVents, zcol = "Region", map.types = "Esri.WorldImagery")
m=mapview(ThermalVents, zcol = "Region", map.types = "Esri.WorldImagery")

## Your first map will be against a standard cartographic background of shorelines
ThermalVents2 <- st_as_sf(metadata2, coords = c("Long", "Lat"),  crs = 4326)
m=mapview(ThermalVents2, zcol = "Region", burst = TRUE, map.types = "CartoDB.Positron", col.regions = colors)
m

## Your second map will be against a ESRI World Imagery background - which shows elevations and depths. 
## You should be able to zoom in and see your sites siting on large geological ridges. 

mapview(ThermalVents2, zcol = "Region", burst = TRUE, col.regions = brewer.pal(7, "Paired"), map.types = "Esri.WorldImagery")

## Now, we're going to follow several steps to export this map to an image so that you can include it in your final pdf working document. 

# First, export your map to html in your working directory using the package webshot. 
webshot::install_phantomjs()

mapshot(m, url = paste0(getwd(), "/map.html"))

## Next we'll use the package pagedown to export this map to a png file in your working directory
chrome_print(
  "map.html",
  output = "nmds_map.png",
  format = "png")

## Last maping step is to add this file to your R environment so that it can be included in your final print to pdf. 
map_static<-load.image("nmds_map.png")
plot(map_static,axes=FALSE, main = "Deep-sea vents where spp. were sampled")

## Now, lets explore the alpha and beta-diversity within and between these plots.  
## First off, we're going to use the vegan package to estimate the species richness of the vents in each region. 

richness = specnumber(data2)
df = as.data.frame(richness)

## Add your specific Site and the Region that site is within to that richness calculation. 
df$Region = metadata2$Region
df$Site = metadata2$Site

## Now we'll plot this as a box-plot using ggplot2
g=ggplot(df, aes(x=Region, y=richness, fill=Region)) + 
  geom_boxplot(alpha=1, show.legend = T) 

alpha = g+scale_fill_viridis_d()+labs(title = "Box-plot of deep-sea vent community alpha-diversity (Species Richness)")

alpha
## Which regions are the most diverse? 

## Next, we're going to use the vegan package to estimate similar or different the community of species living at vents in each region are. 
## We're going to be doing this using the analysis called Non-Metric multiDimentional Scaling (or NMDS), in the vegan and ggordiplot packages. 

set.seed(123)
ord2 <- metaMDS(data2)

## The next block creates this NMDS map and then we'll use the ouput to make a good plot in ggplot. 

my.plot = gg_ordiplot(ord2, groups=metadata2$Region, kind="se", conf = 0.99)

#  this extracts the data from the nmds and then plots the species as points that I've coloured by taxon (class).

a.plot <- my.plot$plot
a.plot + labs(color = "Region", x = "NMDS1", y = "NMDS2", title = "NMDS") +
  theme(plot.title = element_text(hjust = 0.5)) +scale_color_viridis(discrete = TRUE)

beta=a.plot + labs(color = "Region", x = "NMDS1", y = "NMDS2", title = "NMDS plot of deep-sea vent community beta-diversity (Bray-Curtis Distance)") +
  scale_color_viridis(discrete = TRUE)
beta
## Are the regions characterised by the same thermal-vent species? 
 
## ONe way to determine how different they are is to use ANOSIM (Clarke 1993) or the PERMANOVA 

# The ANalysis Of SIMilarity (ANOSIM) test has some similarity to an ANOVA-like hypothesis test, however, it is used to evaluate a dissimilarity matrix rather than raw data (Clarke, 1993). The PERmutational Multivariate ANalysis of VAriance (PERMANOVA) compares the variation between groups to the variation within groups

## Essentially, if two groups of sampling units are really different in their species composition, then compositional dissimilarities between the groups ought to be greater than those within the groups

comm.bc.dist <- vegdist(data2, method = "bray") 
attach(metadata2)

deep.sea.anosim <- anosim(comm.bc.dist , Region)
summary(deep.sea.anosim)
plot(deep.sea.anosim)

## or using adonis2 from the vegan package 
stats = adonis2(data2 ~ Region, data = metadata2,permutations = 999,
                method = "bray")
stats
## so does the region affect the community of species living at these vents? 
##https://uw.pressbooks.pub/appliedmultivariatestatistics/chapter/permanova/


## The next block will create a three page pdf of your map, alpha- and beta-diversity analyses.  Print these off and use them 


pdf("3700_NMDS_deep_sea_thermal_vents_250625.pdf", width = 12, height = 8) # Open a new pdf file
plot(map_static,axes=FALSE, main = "Deep-sea vents where spp. were sampled")
alpha
beta
dev.off() # Close the file


# So - hats off to you!! You've mapped the locations of many of the thermal-vents around the world and then used species incidence data from these sites to compare alpha- and beta-diversity. Now, print your pdf, examine the map, box-plot and NMDS, and prepare to speak for three minutes (!!without notes!!) about the conclusions you might make regarding the diversity, ecological similarity and vulnerability of these thermal-vents. 
#  
# How species rich are these vents
# Are the vent sites distinct? 
# Which is the most diverse? 
# Is that diversity nested?  
# How vulnerable are these sites and regions and what consequences would mining have on species living at and around deep-sea vents? 

# good luck!



