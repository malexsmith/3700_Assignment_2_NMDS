

# ZOO*3700 Assignment 2 : Community richness and similarity between deep-sea vents around the world. 
# 250624
# https://doi.org/10.1098/rsos.200110op 


# First, remove any objects that may still be in your work space (can alter results if you don't and they're there)
rm(list=ls())

# First, you have to install and load the packages used here: 
#install.packages("vegan") 
library(vegan)

#Next, I will assume that you have the data (the .csv files) stored in your working directory:
setwd("C:\\Users\\malex\\Sync\\R")
# You  change the line above to where ever you saved the data files.

# This next block of code will install a package called BioManager that we will use to install 9 further packages (if you have not already installed them). If you have already installed them, you can skip to the library() commands below which will open the packages you need to complete this assignment. 

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager");

BiocManager::install("vegan");
BiocManager::install("ggplot2");
BiocManager::install("bipartite");
BiocManager::install("ggordiplots");

# The next block are the library() commands which will 9 the packages you need to complete this assignment. 

library(vegan); 
library(ggplot2);
library(bipartite);
library(ggordiplots)



# next, upload your community matrix (species incidence or abundance against site or collection)
# and the metadata (the things you know about the sites)
data <- read.csv(file = "3700 test nmds community.csv",head=TRUE,row.names = 1, sep=",")
metadata <- read.csv(file = "3700 test nmds metadata.csv",head=TRUE,row.names = 1, sep=",")
taxon <- read.csv(file = "3700 test nmds taxon metadata.csv",head=TRUE,row.names = 1, sep=",")
metadata2 <- read.csv(file = "3700 bigger test nmds metadata.csv",head=TRUE,row.names = 1, sep=",")

## First, lets make a map
## https://map-rfun.library.duke.edu/01_georeference.html

library(tidyverse)
library(sf)
library(mapview)
library(RColorBrewer)

library(viridis)
colors = viridis(7)

ThermalVents <- st_as_sf(metadata, coords = c("Long", "Lat"),  crs = 4326)
mapview(ThermalVents, zcol = "Region")
mapview(ThermalVents, zcol = "Region", map.types = "Esri.WorldImagery")
m=mapview(ThermalVents, zcol = "Region", map.types = "Esri.WorldImagery")

ThermalVents2 <- st_as_sf(metadata2, coords = c("Long", "Lat"),  crs = 4326)
m=mapview(ThermalVents2, zcol = "Region", burst = TRUE, map.types = "CartoDB.Positron", col.regions = colors)
mapview(ThermalVents2, zcol = "Region", burst = TRUE, col.regions = brewer.pal(7, "Paired"), map.types = "Esri.WorldImagery")

## exporting your map to html
library(webshot)
webshot::install_phantomjs()

mapshot(m, url = paste0(getwd(), "/map.html"))

library(pagedown)
#chrome_print("map.html", output = "nmds_map.pdf")

chrome_print(
  "map.html",
  output = "nmds_map.png",
  format = "png")


map_static<-load.image("nmds_map.png")
plot(map_static,axes=FALSE, main = "Deep-sea vents where spp. were sampled")



data_bipartite <- read.csv(file = "3700 test nmds community.csv",head=TRUE,row.names = 1, sep=",")


plotweb(data_bipartite, method="cca", text.rot = 90)

pdf("3700_bipartite_nmds_250624.pdf", width = 8, height = 12) # Open a new pdf file
plotweb(data_bipartite, method="cca", text.rot = 90)
dev.off() # Close the file
# changed interaction and forest/field colours by hand in illustrator


region=metadata[,1]
environ=metadata[,3:16]

# note that the community and metatdata file must be arrayed in the same order


# 250508  this to colour species by class (or whatever taxon you want)
#environ3 <- read.csv(file = "anibal_meta_column_250508.csv",head=TRUE,row.names = 1, sep=",")


##-----------------------------------------------

# first, let's calculate dissimilarity measures between sites and then cluster them to see what is the most similar. 

comm2 = decostand(data, method="hellinger")
# Anibal - read about what decostand does in the r package vegan for the nmds
# https://r.qcbs.ca/workshop09/book-en/transformations.html

# calculate Bray-curtis distance among samples (or Cao)
comm.bc.dist <- vegdist(comm2, method = "bray")

#cluster communities using average-linkage algorithm.  (You'll use this output "comm.bc.clust" again later on in your biplot)
comm.bc.clust <- hclust(comm.bc.dist, method = "average")

# Now, plot the cluster diagram
plot(comm.bc.clust, ylab = "Bray")

# is there the clustering you expected?  -

var.hel <- decostand(data, method="hellinger")
nmds1 <- metaMDS(var.hel, autrotransform=FALSE)


# ORDINATION

# The metaMDS function automatically transforms data and checks solution robustness

comm.bc.mds <- metaMDS(comm2, dist = "bray") # dist = cao is better for many non-overlapping sites. 

# for anibal's dairy bush data, bray seems like a better plot

# You can adjust the distance measure to others including euclidean, jaccard and manhattan
# and then compare the stress values between your observed and ordination to see which is the best measure.  
# Bray-Curtis ("bray" above) is pretty resilient. 

# comm.bc.mds2 <- metaMDS(comm, dist = "manhattan") #euclidean, jaccard too

# Next, make the stress plot that allows you to make a Shepard diagram that 
# shows you a goodness of fit measure for points in nonmetric multidimensional scaling
stressplot(comm.bc.mds)
# the actual stress value calculated for each distance measure can be seen using this command
comm.bc.mds$stress
# lower values of stress (<0.2 or 0.3) are useful




##=============================
# simplify and then build up to the Figure 
# build up the figure in parts


# First set up the plotting
#area but don't plot anything yet
mds.fig <- ordiplot(comm.bc.mds, type = "none")

# now, let's plot just the sites, colour by habitat, (btw - pch=19 means plot a circle)

#now let's add in the forest sites
points(mds.fig, "sites", pch = 16, col = "darkgoldenrod", select =
         metadata$Region ==
         "SWIR")
points(mds.fig, "sites", pch = 16, col = "darkgreen", select =
         metadata$Region ==
         "CIR")
# let's add confidence ellipses around habitat types
ordiellipse(comm.bc.mds, metadata$Region, conf = 0.95, label =
              TRUE)
## Now, lets ADD ENVIRONMENTAL DATA TO ORDINATION

# the projection of the points here have maximum correlation with the corresponding environmental variables. 
env.z <- decostand(metadata[, 3:16], method = "standardize")
# this uses decostand to create z-scores for enviromental variables with different units. 
plot(envfit(comm.bc.mds, metadata[, 3:16]))
# right now there is an error because i haven't put in na.rm = TRUE somewhere to control for missing dadta. 

## Now, lets ADD the taxon data that you loaded in environ3
# to show how species within each taxa are split across the ecotone
# 

points(mds.fig, "species", pch = 18, col = "cornflowerblue", select =
         taxon$Phylum ==
         "Cnidaria")
points(mds.fig, "species", pch = 18, col = "red", select =
         taxon$Phylum  ==
         "Annelida")
points(mds.fig, "species", pch = 18, col = "darkorchid1", select =
         taxon$Phylum  ==
         "Mollusca")
points(mds.fig, "species", pch = 18, col = "orange", select =
         taxon$Phylum  ==
         "Arthropoda")
points(mds.fig, "species", pch = 18, col = "black", select =
         taxon$Phylum  ==
         "Echinodermata")


## 250508 - now work on doing ggplot2!
## here's one way of doing the  NMDS PLOT IN GGPLOT2888888888888888888888888888888888888888888888888888888888888888888888888888888888888


set.seed(123)
ord <- metaMDS(data)
gg_ordiplot(ord, groups=metadata$Region, kind="se", conf = 0.99)
io=gg_ordiplot(ord, groups=metadata$Region, kind="se", conf = 0.99)


b.plot <- io$plot
b.plot + labs(color = "Region", x = "NMDS1", y = "NMDS2", title = "NMDS") +
  theme(plot.title = element_text(hjust = 0.5)) +scale_color_viridis(discrete = TRUE)

species_plot <- ord$species
species_plot

#  this extracts the data from the nmds and then plots the species as points that I've coloured by taxon (class). 

environ2 = metadata[, 3:6]

gg_envfit(ord=ord, env=environ2, groups=metadata$Region)
two = gg_envfit(ord = ord, env = environ2, perm = 9999, pt.size = 2, alpha = 0.2, groups=metadata$Region)

c.plot <- two$plot
c.plot + labs(color = "Region", x = "NMDS1", y = "NMDS2", title = "NMDS") +
  theme(plot.title = element_text(hjust = 0.5)) +scale_color_viridis(discrete = TRUE)

specnumber(data2)
specnumber(data2, metadata2$Region)


gg_ordibubble(ord, env.var=environ2$Temperature.C., var.label="Temp")
gg_ordibubble(ord, env.var=environ2$Depth.m., var.label="Depth")

##https://github.com/SusanneWalden/OnthePhenologyofProtists/blob/main/03_Alpha_Diversity/AlphaBoxplotGrouped.md

library(reshape2)


richness = specnumber(data2)
df = as.data.frame(richness)

df$Region = metadata2$Region
df$Site = metadata2$Site
df_melted = melt(df)

g = ggplot(df_melted, aes(x = Region, y = value, fill = Region)) + 
  stat_boxplot(geom = "errorbar", width = 0.1, show.legend = F) +
  geom_boxplot(show.legend = T) 

g=ggplot(df, aes(x=Region, y=richness, fill=Region)) + 
  geom_boxplot(alpha=1, show.legend = T) 

alpha = g+scale_fill_viridis_d()+labs(title = "Box-plot of deep-sea vent community alpha-diversity (Species Richness)")
alpha

data2 <- read.csv(file = "3700 bigger test nmds community.csv",head=TRUE,row.names = 1, sep=",")
metadata2 <- read.csv(file = "3700 bigger test nmds metadata.csv",head=TRUE,row.names = 1, sep=",")
data3 <- read.csv(file = "3700 bigger test nmds community_phylum.csv",head=TRUE,row.names = 1, sep=",")

plotweb(data3, method="cca", text.rot = 90)
plotweb(data3,text.rot=90, col.high = c("darkseagreen2","cadetblue2","tan2","lightpink","peachpuff","plum1", "black"))
# this doesn't work - colours are all over the place - even on a phylum dataset. 

set.seed(123)
ord2 <- metaMDS(data2)
my.plot = gg_ordiplot(ord2, groups=metadata2$Region, kind="se", conf = 0.99)


a.plot <- my.plot$plot
a.plot + labs(color = "Region", x = "NMDS1", y = "NMDS2", title = "NMDS") +
  theme(plot.title = element_text(hjust = 0.5)) +scale_color_viridis(discrete = TRUE)

beta=a.plot + labs(color = "Region", x = "NMDS1", y = "NMDS2", title = "NMDS plot of deep-sea vent community beta-diversity (Bray-Curtis Distance)") +
  scale_color_viridis(discrete = TRUE)
beta

ggplot(shapes, aes(x, y)) +
  geom_point(aes(shape = shape_names), fill = "red", size = 5) +
  geom_text(aes(label = shape_names), nudge_y = -0.3, size = 3.5) 

# centers main title, ggplot2 version 2.2+
ord.data <- my.plot$df_ord
head(ord.data)
p= ggplot(data = ord.data, aes(x = x, y = y, color = metadata2$Region)) + geom_point(size = 3) +
  xlab("PCA 1") + ylab("PCA 2")
p+scale_color_viridis(discrete = TRUE)

## ANOSIM (Clarke 1993) or the PERMANOVA 

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




sortweb(data2, sort.order="dec")
plotweb(data2, method="cca", text.rot = 90)


pdf("3700_NMDS_deep_sea_thermal_vents_250625.pdf", width = 12, height = 8) # Open a new pdf file
plot(map_static,axes=FALSE, main = "Deep-sea vents where spp. were sampled")
alpha
beta
dev.off() # Close the file





pdf("3700_NMDS_deep_sea_thermal_vents_250624.pdf", width = 12, height = 8) # Open a new pdf file
plotweb(data_bipartite, method="cca", text.rot = 90)
gg_ordiplot(ord, groups=metadata$Region, kind="se", conf = 0.99)
gg_ordiplot(ord2, groups=metadata2$Region, kind="se", conf = 0.99)
plotweb(data2, method="cca", text.rot = 90)
dev.off() # Close the file


