# ZOO*3700 Assignment 2 : Community richness and similarity between deep-sea vents around the world. 

# Welcome to your second invertebrateR assignments where we continue our exploration of the deep-sea thermal vent communities.  This time, you will map locations, estimate alpha and beta-diversity and ask whether these communities are different and isolated and what effect mining in international waters will have on this biodiversity. 

# As in Assignment 1, to complete this assignment you will need to be running the most up to date version of R in Windows.  There are some new packages that we will be installing to complete this work. 

# As before, in this Assignment, do not simply run the entire code at once here.  This assignment is designed to be followed along in steps (code blocks) from start to finish.  Each step will create an output that subsequent steps depend on.  So remember to follow the code blocks in order and run them one at a time!

# This first block of code will clear your working environment in case you've been using R for something else in the past.

rm(list=ls())

# The next block of code will tell you where your working directory is (i.e. any figures you generate you can find in that folder).  The 'wd' is also where to put your input files.

getwd()

# This next block of code will install a package called BioManager that we will use to install 11 further packages (if you have not already installed them - at least two of them (ggplot2 and imager) you would have installed for the first assignment). For any of the packages that you might have already installed, you can skip down to the library() commands below which will open the packages needed to complete this assignment. 

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager");

# used in assignment 1
BiocManager::install("ggplot2");
BiocManager::install("imager");

# new for this assignment 
BiocManager::install("vegan");
BiocManager::install("bipartite");
BiocManager::install("ggordiplots");
BiocManager::install("sf");
BiocManager::install("mapview");
BiocManager::install("viridis");
BiocManager::install("webshot");
BiocManager::install("pagedown");
BiocManager::install("betapart");

# The next block are the library() commands which will 11 the packages you need to complete this assignment. 

library(vegan); 
library(ggplot2);
library(bipartite);
library(ggordiplots)
library(sf)
library(mapview)
library(viridis)
library(webshot)
library(pagedown)
library(imager)
library(betapart)

# this next block of code will upload incidence matrices (presence, or 1's, and absence, or 0's) for the thermal vents were we have data.  In addition, the metadata files will add details about the sites (including their location, region, area, chemistry etc.)

data2 <- read.csv(file = "3700 bigger test nmds community.csv",head=TRUE,row.names = 1, sep=",");

metadata2 <- read.csv(file = "3700 bigger test nmds metadata.csv",head=TRUE,row.names = 1, sep=",")

## First, lets make a map of the locations using the package mapview. 
## this map is dynamic and you can move around and plot against several background layers

colors = viridis(7)

## Your first map will be against a standard cartographic background of shorelines
ThermalVents2 <- st_as_sf(metadata2, coords = c("Long", "Lat"),  crs = 4326)

mapview(ThermalVents2, zcol = "Region")

m=mapview(ThermalVents2, zcol = "Region", burst = TRUE, map.types = "CartoDB.Positron", col.regions = colors)


## Your second map will be against a ESRI World Imagery background - which shows elevations and depths. 
## You should be able to zoom in and see your sites siting on large geological ridges. 

mapview(ThermalVents2, zcol = "Region", map.types = "Esri.WorldImagery")


## Now, we're going to follow several steps to export this map to an image so that you can include it in your final pdf working document. 

# First, export your map to html in your working directory using the package webshot which needs to install phantomjs.

webshot::install_phantomjs()

mapshot(m, url = paste0(getwd(), "/map.html"))

## Next we'll use the package pagedown to export this html map to a png file in your working directory

chrome_print(
  "map.html",
  output = "nmds_map.png",
  format = "png")

## Last maping step is to add this file to your R environment so that it can be included in your final print to pdf. (Make sure to add you name to the title! Find main and insert your name inbetween the quotes). 


map_static<-load.image("nmds_map.png")
plot(map_static,axes=FALSE, main = "Deep-sea vents where spp. were sampled")

## Now, lets explore the alpha and beta-diversity within and between these deep-sea communities.  
## First off, we're going to use the vegan package to estimate the species richness of the vents in each region. 

richness = specnumber(data2)
df = as.data.frame(richness)

## Add your specific Site and the Region that site is within to that richness calculation. 
df$Region = metadata2$Region
df$Site = metadata2$Site

## Now we'll plot this as a box-plot using ggplot2. (Make sure to add you name to the title! Find title below and insert your name inbetween the quotes). 

g=ggplot(df, aes(x=Region, y=richness, fill=Region)) + 
  geom_boxplot(alpha=1, show.legend = T) 

alpha = g+scale_fill_viridis_d()+labs(title = "Box-plot of deep-sea vent community alpha-diversity (Species Richness)")

alpha
## Which regions are the most diverse? Are the regions significantly different in terms of total richness? One simple way of looking at this is with an ANOVA. 

alpha.anova <- aov(richness ~ Region, data=df)
summary(alpha.anova)

## Next, we're going to use the vegan package to estimate similar or different the community of species living at vents in each region are. 
## We're going to be doing this using the analysis called Non-Metric multiDimentional Scaling (or NMDS), in the vegan and ggordiplot packages. 

set.seed(123)
ord2 <- metaMDS(data2)

## The next block creates this NMDS map and then we'll use the ouput to make a good plot in ggplot. 

my.plot = gg_ordiplot(ord2, groups=metadata2$Region, kind="se", conf = 0.99)

#  The above code for 'my.plot' is your NMDS map, but the colours don't align with your earlier plots - so let's extracts the data from the nmds and then plots the species as points that I've coloured by taxon (class).

a.plot <- my.plot$plot
a.plot + labs(color = "Region", x = "NMDS1", y = "NMDS2", title = "NMDS") +
  theme(plot.title = element_text(hjust = 0.5)) +scale_color_viridis(discrete = TRUE)

# Now to prepare your final NMDS plot (Make sure to add you name to the title! Find title below and insert your name inbetween the quotes).

beta=a.plot + labs(color = "Region", x = "NMDS1", y = "NMDS2", title = "NMDS plot of deep-sea vent community beta-diversity (Bray-Curtis Distance)") +
  scale_color_viridis(discrete = TRUE)
beta

## Are the regions characterised by the same thermal-vent species? 
 
## One way to determine how different they are is to use ANOSIM (Clarke 1993) or the PERMANOVA 

# The ANalysis Of SIMilarity (ANOSIM) test has some similarity to an ANOVA-like hypothesis test, however, it is used to evaluate a dissimilarity matrix rather than raw data (Clarke, 1993). The PERmutational Multivariate ANalysis of VAriance (PERMANOVA) compares the variation between groups to the variation within groups

## Essentially, if two groups of sampling units are really different in their species composition, then compositional dissimilarities between the groups ought to be greater than those within the groups

comm.bc.dist <- vegdist(data2, method = "bray") 
attach(metadata2)

# This next block will run the ANOSIM to compare the variation between groups to the variation within groups
deep.sea.anosim <- anosim(comm.bc.dist , Region)
summary(deep.sea.anosim)

## This next block will run the PERMANOVA to compare the variation between groups to the variation within groups
stats = adonis2(data2 ~ Region, data = metadata2,permutations = 999,
                method = "bray")
stats
## so does the region affect the community of species living at these vents? Compare the variation between groups to the variation within groups.

# What component of this betadiversity turns over across space, and what component is nested, one within another? The next code block uses the package betapart to differentiate the importance of these components.
# beta.JTU is the  value of the turnover component, measured as Simpson dissimilarity.  
# beta.JNE is the value of the nestedness component, measured as nestedness-resultant fraction of Sorensen dissimilarity
# beta.JAC is the value of the overall beta diversity, measured as Sorensen dissimilarity


turnover_or_nestedness = beta.multi(data2, index.family="jaccard")
turnover_or_nestedness=as.data.frame(turnover_or_nestedness)
turnover_or_nestedness


## Finally - the next block will create a three page pdf of your map, alpha- and beta-diversity analyses.  Print these off and use them 


pdf("3700_Assignment_2_NMDS_deep_sea_thermal_vents_250625.pdf", width = 12, height = 8) # Open a new pdf file
plot(map_static,axes=FALSE, main = "Deep-sea vents where spp. were sampled")
alpha
beta
dev.off() # Close the file


# You made it - Amazing!! You've mapped the locations of many of the thermal-vents around the world and then used species incidence data from these sites to compare patterns of alpha- and beta-diversity. Now, print your pdf, examine the map, box-plot and NMDS, and prepare to speak for three minutes (!!without notes!!) about the conclusions you might make regarding the diversity, ecological similarity and vulnerability of these thermal-vents. 

  
# How species rich are these vents? Which is the most diverse? Is alpha diversity different between regions?

# Are the vent Regions distinct? How does variation between regions compare to the variation within regions?

# Is that diversity nested within the diversity of another?  

# How vulnerable are these sites and regions? What consequences would mining have on the species living at and around deep-sea vents? 

# good luck!
#


