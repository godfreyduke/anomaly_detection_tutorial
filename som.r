# Updated August 12, 2017
# set working directory
setwd("put R working directory name here")
# modify this to the path for your directory
# Download kohonen package documentation from:
#https://cran.r-project.org/web/packages/kohonen/kohonen.pdf


# load packages
library(rmarkdown)
library(knitr)
library(dplyr)
library(kohonen)
library(dummies)
library(ggplot2)
library(sp)
library(reshape2)
library(RColorBrewer)
library(magrittr)

#read data
SR1 <- read.csv(file = "SR1-SOM-example.csv", head=TRUE, sep =",")

#exploratory analysis
summary(SR1)
SR1_SOM <- SR1[, -c(1)]  #remove first column
summary(SR1_SOM)
str(SR1_SOM)
names(SR1_SOM)

hist(SR1_SOM$EventDur)
select(SR1_SOM, EventDur) %>% filter(EventDur <= 1000) -> d  #select values less tha 1000
hist(d$EventDur, breaks = 15)

#create frequency tables of error codes by minute
EC1_table <- table(SR1_SOM$Time, SR1_SOM$EC1)
margin.table(EC1_table, 1)
margin.table(EC1_table, 2)
plot(margin.table(EC1_table, 1), 
     col = c("darkblue", "maroon"), xlab = "Time", ylab = "Counts")  #summed over Time
margin.table(EC1_table, 1)  #summed over EC1
EC2_table <- table(SR1_SOM$Time, SR1_SOM$EC2)
margin.table(EC2_table, 2)
EC2_counts <- table(SR1_SOM$EC2, SR1_SOM$Time)
barplot(EC2_counts, main="Distribution of EC2 by Minute",
        xlab="Time", col = rainbow(7),
        legend = rownames(EC2_counts), beside=TRUE)
EC3_table <- table(SR1_SOM$Time, SR1_SOM$EC3)


# Colour palette definition
display.brewer.all()
#blue color for large cluster
colors <- brewer.pal(10, "Paired") #first 2 colors are blue
pal <- colorRampPalette(colors)
my_palette <- c(pal(10))

# Palette defined by kohonen package
coolBlueHotRed <- function(n, alpha = 1) {
  rainbow(n, end=4/6, alpha=alpha)[n:1]
}

# SOM Model #################################################
# scale the data
#center is TRUE then centering done by subtracting the column means (omitting NAs) of x from their 
#corresponding columns
#scale is TRUE then scaling done by dividing the (centered) columns of x by their standard deviations
#if center is TRUE, and the root mean square otherwise.
SR1_SOM.sc <- scale(SR1_SOM, center = TRUE, scale = TRUE)
summary(SR1_SOM.sc)

set.seed(3)
# you can experiment with other grid sizes
som_grid <- somgrid(xdim = 15, ydim=10, topo="hexagonal")
som_model <- som(SR1_SOM.sc, grid = som_grid, rlen = 100) 
#rlen is the number of times the complete data set will be presented to the network 
summary(som_model)
print(som_model)

# Changes by iteration (specified with rlen, try different values)
plot(som_model, type = "changes", main = "SR1: SOM")

#plot SOMs
#counts per node - empty nodes shown in gray
plot(som_model, type = "counts", main="SR1: Node Counts")

#shows the sum of the distances to all immediate neighbours.  
#also known as a U-matrix plot.
#Units near a class boundary likely to have higher average distances to their neighbours
plot(som_model, type="dist.neighbours", main = "SR1: SOM neighbour distances", palette.name=grey.colors)

#code spread
plot(som_model, type = "codes", main = "SR1: Codebook Vectors")

#shows the mean distance of objects mapped to a unit to the codebook vector of that unit.
#The smaller the distances, the better the objects are represented by the codebook vectors
plot(som_model, type = "quality", main="SR1: Node Quality/Distance")

# Plot the original scale heatmap for all variables 
# (it will be from training set if a training dataset was created; 
# we did not do that, used all the variables for the model)
var <- 1  #column number 1 is EventDuration
plot(som_model, type = "property", property = as.data.frame(som_model$codes)[,var], 
     main=names(som_model$SR1_SOM)[var],palette.name=coolBlueHotRed )

# Plot the original scale heatmap for all variables
par(mfrow=c(2,3))  #6 plots per page
for (i in 1:6) {
  var <- i #define the variable to plot
  var_unscaled <- aggregate(as.numeric(SR1_SOM[,var]), 
                            by=list(som_model$unit.classif), FUN=mean, simplify=TRUE)[,2]
  plot(som_model, type = "property", property=var_unscaled, main=names(SR1_SOM)[var],
       palette.name=coolBlueHotRed)
#rm(var_unscaled, var)
}
par(mfrow=c(1,1))

#Create plots by minute
SR1_SOM_4 <- filter(SR1_SOM, Time == 4)
#Create SOM models by minute
#Experiment with different grid size
#
#


# ------------------ Clustering SOM results -------------------

# Show the WCSS (within cluster sum of squares) metric for kmeans 
# for different clustering sizes.
# Can be used as a "rough" indicator of the ideal number of clusters
# --> have to convert som_model$codes from list to dataframe, else gives error
mySR1_SOM <- as.data.frame((som_model$codes))
wss <- ((nrow(mySR1_SOM))-1)*sum(apply(mySR1_SOM,2,var))
for (i in 2:15) wss[i] <- sum(kmeans(mySR1_SOM,
                                     centers=i)$withinss)
par(mar=c(5.1,4.1,4.1,2.1))
plot(1:15, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares", main="Within cluster sum of squares (WCSS)")

# Form clusters on grid
## use hierarchical clustering to cluster the codebook vectors; 
## last parameter is number of clusters
som_cluster <- cutree(hclust(dist(as.data.frame((som_model$codes)))), 8)
som_cluster

# Show the map with different colours for every cluster  					  
plot(som_model, type="mapping", 
#pch = ".",
     labels = as.integer(SR1_SOM$Time), col = as.integer(SR1_SOM$Time),
     bgcol = my_palette[som_cluster], main = "SR1 - 8 Clusters")
add.cluster.boundaries(som_model, som_cluster, lwd = 5, col = "maroon")
#identify.kohonen

#show the same plot with the codes instead of colours and points
plot(som_model, type="codes", codeRendering = "segments", bgcol = my_palette[som_cluster], 
     main = "SR1 - 8 Clusters")
add.cluster.boundaries(som_model, som_cluster, lwd = 3, col = "brown")
identify(som_model, som_cluster)

#######################################################################