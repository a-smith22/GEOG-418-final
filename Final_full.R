#Geog 418/518 Final Project
install.packages("plyr")
install.packages("dplyr")
install.packages("spdep")
install.packages("GISTools")
install.packages("raster")
install.packages("maptools")
install.packages("spatstat")
install.packages("sp")
install.packages("spgwr")
install.packages("rgdal")
install.packages("tmap")
install.packages("BAMMtools")
install.packages("shinyjs")
install.packages("gridExtra")
install.packages("grid")
install.packages("gtable")
install.packages("gstat")

library(tmap)
library(BAMMtools)
library(shinyjs)
library(gridExtra)
library(grid)
library(gtable)
library(gstat)
library(dplyr)
library(spdep)
library(GISTools)
library(raster)
library(maptools)
library(rgdal)
library(spatstat)
library(sp)
library(spgwr)

#Set working directory
setwd("D:Final/working")

png("sutdyarea.png")
tm_shape(income.tracts) + 
  tm_polygons() +
  tm_layout(title = "Greater Vancouver, BC", title.position = c("centre", "TOP")) +
  tm_compass(position=c("left", "bottom"))
dev.off()


png("map_studysitecalifornia.png")
map_studysitecalifornia
dev.off()


#Reading in particulate matter dataset
pm25 <- read.csv("PM25.csv") #Read in PM2.5 data
#Select only columns 1 and 2
pm25 <- pm25[,1:2]
#Change the column names 
colnames(pm25) <- c("POSTALCODE", "PM25")
#Reading in postal code shapefile
postalcodes <- readOGR("BC_Postal_Codes") #Read in related postal code data


#Reading in the income dataset
income <- read.csv("Income.csv") #Read in census income data  
#Change the column names
colnames(income) <- c("DAUID", "Income") #Select only ID and Income columns
#Read in the dissemination tract shapefile
census.tracts <- readOGR("BC_DA") 
#Merge the income dataset and the DA shapefile
income.tracts <- merge(census.tracts,income, by = "DAUID") 
#Remove any NA's from the merged dataset
income.tracts <- income.tracts[!is.na(income.tracts$Income),]


#Create choropleth map of income
png("income_choro.png")
tm_shape(income.tracts) + 
  tm_layout(legend.position = c("left","BOTTOM")) +
  tm_polygons(col = "Income", 
              title = "Median Income", 
              style = "jenks", 
              palette = "Greens", n = 6)
dev.off()

png("income_choro_noborders.png")
tm_shape(income.tracts, simplify =  0.5) + 
  tm_layout(legend.position = c("left","BOTTOM")) +
  tm_fill(col = "Income", 
          title = "Median Income", 
          style = "jenks", 
          palette = "Greens", n = 6,
          lwd = 0)
dev.off()

#Select postal codes that fall within dissemination tracts)
postalcodes <- intersect(postalcodes,income.tracts)
plot(postalcodes) #See what the data looks like spatially
head(postalcodes) #See what the data looks like in tabular form


#Join PM2.5 data with postal code data
pm25.spatial <- merge(postalcodes,pm25,by = "POSTALCODE")
pm25.spatial <- pm25.spatial[!is.na(pm25.spatial$PM25),]
#Aggregate the PM2.5 values in each DA in order to have a single value per DA. Here we aggregate based on the mean.
pm25.aggregate <- aggregate((as.numeric(pm25.spatial$PM25)/10)~pm25.spatial$DAUID,FUN=max)

#Re-join aggregated data to the income.tracts layer.
colnames(pm25.aggregate) <- c("DAUID", "PM25AGG") #Select only ID and Income columns
income.pm25 <- merge(income.tracts,pm25.aggregate, by = "DAUID") #Merge income and dissemination data

#Re-join aggregated data to the pm25.spatial points layer.
pm25.points.aggregate <- merge(pm25.spatial, pm25.aggregate, by = "DAUID")
pm25.points.aggregate <- pm25.points.aggregate[!is.na(pm25.points.aggregate$PM25AGG),]

#Create a subsample of the datapoints provided in the PM2.5 dataset using the sample n provided on CourseSpaces
sampleSize=310
spSample <- pm25.points.aggregate[sample(1:length(pm25.points.aggregate),sampleSize),]

####ASSESS PATTERN OF SAMPLED POINTS
##Convert the nearestNeighbor object into a dataframe.
nearestNeighbour <- nndist(spSample$PM25AGG)
nearestNeighbour=as.data.frame(as.numeric(nearestNeighbour))
##Change the column name to "Distance"
colnames(nearestNeighbour) = "Distance"
##Calculate the nearest neighbor statistic to test for a random spatial distribution.
income.tracts <- spTransform(income.tracts, CRS("+init=epsg:26911"))
studyArea <- gArea(income.tracts, byid= FALSE)

nnd = (sum(nearestNeighbour$Distance))/310 #MEAN NND
pointDensity <- (310/studyArea) 
r.nnd = 1/(2*sqrt(pointDensity)) #RANDOM NND
d.nnd = 1.07453/sqrt(pointDensity) #DISPERSED NND
R = nnd/r.nnd #STANDARDIZED MEAN NND
SE.NND <- 0.26136/sqrt(310*pointDensity) #STANDARD ERROR OF MEAN NND
z = (nnd-r.nnd)/SE.NND #Z-SCORE

########################################
#####Weighted Matrix
crd.nb <- poly2nb(income.tracts)
crd.net <- nb2lines(crd.nb,coords=coordinates(income.tracts))##Converts neighbour matrix into plot

png("poly_neighbours.png")
tm_shape(income.tracts) + tm_borders(col='lightgrey') + 
  tm_shape(crd.net) + tm_lines(col='red')
dev.off()
##This graphic shows us the neighbours of the polygons. We are using queen weight in red line (queen = default)

crd.nb2 <- poly2nb(income.tracts, queen = FALSE)#we are using the rook weight for the yellow line one
crd.net2 <- nb2lines(crd.nb2,coords=coordinates(income.tracts))

png("queenandrook.png")
tm_shape(income.tracts) + tm_borders(col='lightgrey') + 
  tm_shape(crd.net) + tm_lines(col='blue', lwd = 2) + #blue lines are corner connections with polygons
  tm_shape(crd.net2) + tm_lines(col='yellow', lwd = 2)
dev.off()
########################

crd.lw <- nb2listw(crd.nb, zero.policy = TRUE, style = "W")
print.listw(crd.lw, zero.policy = TRUE) #These tell us our list matrix

########################  

lisa.test <- localmoran(income.tracts$Income, crd.lw) ##Very similar equation to moran's i

income.tracts$Ii <- lisa.test[,1] ##Lisa i value
income.tracts$E.Ii<- lisa.test[,2] ##Expected lisa i value
income.tracts$Var.Ii<- lisa.test[,3] ##Lisa variance value
income.tracts$Z.Ii<- lisa.test[,4] ##lisa provides a z value, we don't have to calculate for each polygon
income.tracts$P<- lisa.test[,5]  ##Will give us a bunch of descreptive statistics for each polygon


png("income_lisa.png")
map_LISA <- tm_shape(income.tracts) + ##make a map out of our crd data, moran's i, and lisa
              tm_polygons(col = "Ii", 
                  title = "Local Moran's I", 
                  style = "fisher", 
                  palette = "viridis", n = 6)
map_LISA
dev.off()

#####Moran's I
mi <- moran.test(income.tracts$Income, crd.lw, zero.policy = TRUE) ##This is our whole moran's i test.
mi ##shows our results
moran.range <- function(lw) {
  wmat <- listw2mat(lw)
  return(range(eigen((wmat + t(wmat))/2)$values)) ##shows us what the full range of mi test can be for this data
  
}

moran.range(crd.lw)

mI <- mi$estimate[[1]] ##moran's i value
eI <- mi$estimate[[2]] ##expected i value
var <- mi$estimate[[3]]
z <- (mI - eI) / (sqrt(var))

mI.r <- round(mI, digits = 4)
eI.r <- round(eI, digits = 4)
var.r <- round(var, digits = 4)
z.r <- round(z, digits = 4)

data.for.table2 = data.frame(mI.r, eI.r, var.r, z.r)

png("MoranI_income.png")
table2 <- tableGrob(data.for.table2, cols = c("Moran's I","Expected I", "Variance", "Z-Score"), rows = " ")  #make a table
t2Caption <- textGrob("Table 2: Global Moran's I Values for Income", gp = gpar(fontsize = 08)) 
padding <- unit(5, "mm") 

table2 <- gtable_add_rows(table2, 
                          heights = grobHeight(t2Caption) + padding, 
                          pos = 0)

table2 <- gtable_add_grob(table2,
                          t2Caption, t = 1, l = 2, r = ncol(data.for.table2) + 1)
grid.arrange(table2, newpage = TRUE)
dev.off()

#######Kriging
#Create a grid called grd to use in your interpolation
# Create an empty grid where n is the total number of cells
grd <- as.data.frame(spsample(spSample, "regular", n = 50000))
names(grd)       <- c("X", "Y")
coordinates(grd) <- c("X", "Y")
gridded(grd)     <- TRUE  # Create SpatialPixel object
fullgrid(grd)    <- TRUE  # Create SpatialGrid object
proj4string(grd) <- proj4string(spSample)
head(spSample)

#find the right order to define the trend 
spSample$X <- coordinates(spSample)[,1] #taking the coordinates and storing in columns labeled x,y
spSample$Y <- coordinates(spSample)[,2]

f.1 <- as.formula(value ~ X + Y) 
f.2 <- as.formula(PM25AGG ~ X + Y + I(X*X)+I(Y*Y) + I(X*Y))

var.smpl <- variogram(f.2, spSample, cloud = FALSE) #, cutoff=1000000, width=89900)
dat.fit  <- fit.variogram(var.smpl, fit.ranges = FALSE, fit.sills = FALSE,
                          vgm(model="Sph"))
png("trend_plot.png")
plot(var.smpl, dat.fit)
dev.off()

#do the kriging
dat.krg <- krige( f.2, spSample, grd, dat.fit) ##This is what interpolates

# Convert kriged surface to a raster object for clipping
r <- raster(dat.krg)
r.m <- mask(r, income.tracts)

# Plot the map
png("kriging_pm25.png")
tm_shape(r.m) + 
  tm_raster(n=10, palette="-RdBu",  
            title="Predicted PM2.5 \n(ug/m3)") +
  tm_shape(spSample) + tm_dots(size=0.2) +
  tm_legend(legend.outside=TRUE)
dev.off()

#Plot the Variance
r   <- raster(dat.krg, layer="var1.var")
r.m <- mask(r, income.tracts)

png("variance_krig.png")
tm_shape(r.m) + 
  tm_raster(n=7, palette ="-RdBu",
            title="Variance map \n(in squared ug/m3)") +tm_shape(spSample) + tm_dots(size=0.2) +
  tm_legend(legend.outside=TRUE)
dev.off()

#Plot the Confidence Interval
r   <- sqrt(raster(dat.krg, layer="var1.var")) * 1.96
r.m <- mask(r, income.tracts)

png("CI_krig.png")
tm_shape(r.m) + 
  tm_raster(n=7, palette ="-RdBu",
            title="95% CI map \n(ug/m3)") +tm_shape(spSample) + tm_dots(size=0.2) +
  tm_legend(legend.outside=TRUE)
dev.off()


###join income and pm25

pm.income.poly <- extract(r, income.tracts, fun = mean, sp = TRUE) #r is the raster used for the kriging earlire (pm25)
colnames(pm.income.poly@data)[25] <- "PM25"
pm.income.poly <- pm.income.poly[!is.na(pm.income.poly$PM25),]

######Linear Regression##########
#Plot income and PM2.5 from the pm.income.poly dataset you created
plot(pm.income.poly$PM25~pm.income.poly$Income)
#Notice that there are a lot of 0's in this dataset. If you decide to remove them, use the following line:
pm.income.poly <- pm.income.poly[(pm.income.poly$PM25 >= 0.1), ]
#Now plot the data again
plot(pm.income.poly$PM25~pm.income.poly$Income)

#Perform a linear regression on the two variables. You should decide which one is dependent.
lm.model <- lm(pm.income.poly$PM25~pm.income.poly$Income)
#Add the regression model to the plot you created
abline(lm.model)
#Get the summary of the results
summary(lm.model)

#You want to determine if the model residuals are spatially clustered. 
#First obtain the residuals from the model
model.resids <- as.data.frame(residuals.lm(lm.model))
#Then add the residuals to your spatialpolygon dataframe
pm.income.poly$residuals <- residuals.lm(lm.model)
#Observe the result to make sure it looks correct
head(pm.income.poly)

#Now, create choropleth map of residuals
png("residuals.png")
tm_shape(pm.income.poly) + 
  tm_layout(title = "Linear Regression Residuals", title.position = c("center","TOP")) +
  tm_layout(legend.position = c("left","BOTTOM")) +
  tm_polygons(col = "residuals", 
              title = "Residuals", 
              style = "fisher", 
              palette = "-PuBu", n = 6)
dev.off()

##MoransI2
pm.income.poly <- spTransform(pm.income.poly, CRS("+init=epsg:26911"))
View(pm.income.poly@data)
#####Weighted Matrix
crd.nb.reg <- poly2nb(pm.income.poly)
crd.net.reg <- nb2lines(crd.nb.reg,coords=coordinates(pm.income.poly))##Converts neighbour matrix into plot

png("inc_pm25.png")
tm_shape(pm.income.poly) + tm_borders(col='lightgrey') + 
  tm_shape(crd.net.reg) + tm_lines(col='red')
dev.off()
##This graphic shows us the neighbours of the polygons. We are using queen weight in red line (queen = default)

crd.nb2.reg <- poly2nb(pm.income.poly, queen = FALSE)#we are using the rook weight for the yellow line one
crd.net2.reg <- nb2lines(crd.nb2.reg,coords=coordinates(pm.income.poly))

plot(pm.income.poly)
plot(step.5a)

png("pm.in.weight.png")
tm_shape(pm.income.poly) + tm_borders(col='lightgrey') + 
  tm_shape(crd.net.reg) + tm_lines(col='blue', lwd = 2) + #blue lines are corner connections with polygons
  tm_shape(crd.net2.reg) + tm_lines(col='yellow', lwd = 2)
dev.off()
########################

crd.lw.reg <- nb2listw(crd.nb.reg, zero.policy = TRUE, style = "W")
print.listw(crd.lw.reg, zero.policy = TRUE) #These tell us our list matrix

#####Moran's I
mi.reg <- moran.test(pm.income.poly$residuals, crd.lw.reg, zero.policy = TRUE) ##This is our whole moran's i test.
mi.reg ##shows our results

moran.range <- function(lw) {
  wmat <- listw2mat(lw)
  return(range(eigen((wmat + t(wmat))/2)$values)) ##shows us what the full range of mi test can be for this data
  
}

moran.range(crd.lw.reg)


mI.reg <- mi.reg$estimate[[1]] ##moran's i value
eI.reg <- mi.reg$estimate[[2]] ##expected i value
var.reg <- mi.reg$estimate[[3]]
z.reg <- (mI.reg - eI.reg) / (sqrt(var))

mI.reg <- round(mI.reg, digits = 4)
eI.reg <- round(eI.reg, digits = 4)
var.reg <- round(var.reg, digits = 4)
z.reg <- round(z.reg, digits = 4)

data.for.table3 = data.frame(mI.reg, eI.reg, var.reg, z.reg)

png("pm.in.moran.stat.table.png")
table3 <- tableGrob(data.for.table3, cols = c("Moran's I","Expected I", "Variance", "Z-Score"), rows = " ")  #make a table
t3Caption <- textGrob("Table 3: Global Moran's I Values for Regression Model", gp = gpar(fontsize = 08)) 
padding <- unit(5, "mm") 

table3 <- gtable_add_rows(table3, 
                          heights = grobHeight(t3Caption) + padding, 
                          pos = 0)

table3 <- gtable_add_grob(table3,
                          t3Caption, t = 1, l = 2, r = ncol(data.for.table3) + 1)
grid.arrange(table3, newpage = TRUE)
dev.off()


########################  

lisa.test <- localmoran(pm.income.poly$residuals, crd.lw.reg) ##Very similar equation to moran's i

pm.income.poly$Ii <- lisa.test[,1] ##Lisa i value
pm.income.poly$E.Ii<- lisa.test[,2] ##Expected lisa i value
pm.income.poly$Var.Ii<- lisa.test[,3] ##Lisa variance value
pm.income.poly$Z.Ii<- lisa.test[,4] ##lisa provides a z value, we don't have to calculate for each polygon
pm.income.poly$P<- lisa.test[,5]  ##Will give us a bunch of descreptive statistics for each polygon


png("lisa_residuals.png")
map_LISA <- tm_shape(pm.income.poly) + ##make a map out of our crd data, moran's i, and lisa
  tm_polygons(col = "Ii", 
              title = "Local Moran's I", 
              style = "fisher", 
              palette = "viridis", n = 6) 
map_LISA
dev.off()


####Geographically Weighted Regression

#The first thing you need to do is to add the polygon coordinates to the spatialpolygondataframe.
#You can obtain the coordinates using the "coordinates" function from the sp library
pm.income.poly.coords <- sp::coordinates(pm.income.poly)
#Observe the result
head(pm.income.poly.coords)
#Now add the coordinates back to the spatialpolygondataframe
pm.income.poly$X <- pm.income.poly.coords[,1]
pm.income.poly$Y <- pm.income.poly.coords[,2]
head(pm.income.poly)

pm.income.poly <- pm.income.poly[!is.na(pm.income.poly$PM25),] ##remove NA values
View(pm.income.poly@data) ##View data of this spatial data frame


###Determine the bandwidth for GWR: this will take a while
GWRbandwidth <- gwr.sel(pm.income.poly$PM25~pm.income.poly$Income, 
                        data=pm.income.poly, coords=cbind(pm.income.poly$X,pm.income.poly$Y),adapt=T) 

###Perform GWR on the two variables with the bandwidth determined above
###This will take a looooooong while
gwr.model = gwr(pm.income.poly$PM25~pm.income.poly$Income, ###anywhere there is PM25, change to layer
                data=pm.income.poly, coords=cbind(pm.income.poly$X,pm.income.poly$Y), 
                adapt=GWRbandwidth, hatmatrix=TRUE, se.fit=TRUE) 

#Print the results of the model
gwr.model

#Look at the results in detail
results<-as.data.frame(gwr.model$SDF)
head(results)

#Now for the magic. Let's add our local r-square values to the map
pm.income.poly$localr <- results$localR2

#map the GWR R2 values
png("gwr_r2.png")
tm_shape(pm.income.poly) + 
  tm_layout(title = "R-square Values", title.position = c("left","BOTTOM")) +
  tm_polygons(col = "localr", 
              title = "R-square Values", 
              style = "fisher", 
              palette =  "Spectral", n = 6, midpoint = NA)
dev.off()
png("gwr_r2_noline.png")
tm_shape(pm.income.poly) + 
  tm_layout(title = "R-square Values", title.position = c("left","BOTTOM")) +
  tm_fill(col = "localr", 
              title = "R-square Values", 
              style = "fisher", 
              palette = "-RdBu", n = 6,
              midpoint = NA,
              lwd = 0)
dev.off()
#Time for more magic. Let's map the coefficients
pm.income.poly$coeff <- results$pm.income.poly.Income

#Create choropleth map of the coefficients
png("coeff.png")
tm_shape(pm.income.poly) + 
  tm_layout(title = "Coefficients", title.position = c("center","TOP")) +
  tm_polygons(col = "coeff", 
              title = "Coefficients", 
              style = "fisher", 
              palette = "-RdBu", n = 6)

