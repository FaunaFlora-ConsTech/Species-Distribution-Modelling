
# 23 Apr, 2013 by Jeremy Yoder
#https://www.molecularecologist.com/2013/04/23/species-distribution-models-in-r/

# Load packages
require(maptools)
require(maps)
require(mapdata)
require(dismo)# dismo has the SDM analyses we"ll need

# Load the table of latitude and longitude coordinates
locs <- read.csv(file = "JoTrPresence02202008_dryad.txt",
                 header = TRUE,
                 sep = "\t")
str(locs)

# Load the database with the U.S. state borders
data(stateMapEnv)

# Plot these points to check them
# We're limiting the extent of the map to focus on the Mojave Desert region
plot(x = c(-119, -113),
     y = c(33.5, 38),
     mar = par("mar"),
     xlab = "longitude",
     ylab = "latitude",
     xaxt = "n",
     yaxt = "n",
     type = "n",
     main = "Joshua tree presence data")
rect(xleft = par("usr")[1],
     ybottom = par("usr")[3],
     xright = par("usr")[2],
     ytop = par("usr")[4],
     col = "lightblue")
map("state",
    xlim = c(-119, -113),
    ylim = c(33.5, 38),
    fill = TRUE,
    col = "cornsilk",
    add = TRUE)

# add some nice state labels ...
text(x = -117.5, y = 35.5, "California", col = "cornsilk3", cex = 3)
text(x = -116, y = 37.5, "Nevada", col = "cornsilk3", cex = 3)
text(x = -113, y = 34.5, "Arizona", col = "cornsilk3", cex = 3)
text(x = -113, y = 37.75, "Utah", col = "cornsilk3", cex = 3)

# plot the points
points(locs$longitude,
       locs$latitude,
       col = "darkolivegreen4",
       pch = 20,
       cex = 0.5)

# add some axes
axis(side = 1, las = 1)
axis(side = 2, las = 1)

# and add the box around the map
box()

# Uneven distribution: Joshua trees aren’t very evenly distributed/sampling
# effort isn’t very evenly distributed
# we can divide the map up into a grid and randomly draw a
# single point from each grid square

# create sequences of latitude and longitude values to define the grid
longrid <- seq(from = -119, to = -113, by = 0.05)
latgrid <- seq(from = 33.5, to = 38, by = 0.05)

# identify points within each grid cell, draw one at random
# original code has dim(gridsq)[1],
# I replaced with seq_len(dim(gridsq)[1]) to obey sequence linter rules
subs <- c()
for (i in 1:(length(longrid) - 1)) {
  for (j in 1:(length(latgrid) - 1)) {
    gridsq <- subset(locs, latitude > latgrid[j] &
                       latitude < latgrid[j + 1] &
                       longitude > longrid[i] &
                       longitude < longrid[i + 1])
    if (dim(gridsq)[1] > 0) {
      subs <- rbind(subs, gridsq[sample(1:seq_len(dim(gridsq)[1]), 1), ])
    }
  }
}
dim(subs) # confirm that you have a smaller dataset than you started with
head(subs)

points(subs$longitude,
       subs$latitude,
       col = "palegreen",
       pch = 16,
       cex = 0.5)

# Create absence points (pseudo-absence points)
# need some points that define regions where Joshua trees aren’t found
# define a “background” region to sample at random,
# which captures environmental conditions that Joshua trees
# could disperse to, but hasn’t
# To do this, we’ll define circular regions,
# each with a radius of 50km, centered on each point in our presence list,
# then draw random points that must fall within those circles.
# Effectively, this draws random points that must be within 50 of a point where
# Joshua trees have been spotted. These will be in regions that have reasonably
# similar environments to the points where we’ve identified Joshua trees.

# define circles with a radius of 50 km around the subsampled points
x <- circles(p = subs[, c("longitude", "latitude")],
             d = 50000,
             lonlat = TRUE)
str(x)

# draw random points that must fall within the circles in object x
bg <- spsample(x = x@polygons,
               n = 1000,
               type = "random",
               iter = 1000)

plot(x,
     col = adjustcolor("black", 0.4),
     border = "blue",
     add = TRUE)
points(bg,
       col = "brown",
       cex = 0.6)

# We have presence points and absence points.

# Get bioclim data (or any other raster layers of your interest)
bio_clim <- getData("worldclim", var = "bio", res = 2.5, path = "")

ybrev_range <- extent(-119.25, -112.75, 33.25, 38.25) # define the extent

# crop the Bioclim data
bio_clim <- crop(x = bio_clim, y = ybrev_range)
writeRaster(bio_clim,
            filename = "YbrevBC_2.5.grd",
            overwrite = TRUE)

# re-load that cropped dataset
bio_clim <- brick("YbrevBC_2.5.grd")
bio_clim

# this format plots the first (of 19) variables stored in BClim;
# change the 1 to 2-19 for the others
plot(bio_clim, 1,
     cex = 0.5,
     legend = TRUE,
     mar = par("mar"),
     xaxt = "n",
     yaxt = "n",
     main = "Annual mean temperature (ºC x 10)")
map("state",
    xlim = c(-119, -113),
    ylim = c(33.5, 38),
    fill = FALSE,
    col = "cornsilk",
    add = TRUE)

# state names
text(x = -117.5, y = 35.5, "California", col = rgb(1, 1, 1, 0.4), cex = 3)
text(x = -116, y = 37.5, "Nevada", col = rgb(1, 1, 1, 0.4), cex = 3)
text(x = -113, y = 34.5, "Arizona", col = rgb(1, 1, 1, 0.4), cex = 3)
text(x = -113, y = 37.75, "Utah", col = rgb(1, 1, 1, 0.4), cex = 3)

# plot the presence points
points(subs$longitude,
       subs$latitude,
       pch = 20,
       cex = 0.5,
       col = "darkgreen")

# and the pseudo-absence points
points(bg,
       cex = 0.5,
       col = "darkorange3")

# add axes
axis(side = 1, las = 1)
axis(side = 2, las = 1)

# restore the box around the map
box()

# use the extract() function

# pulling bioclim values
# for the subsampled presence points
ybrev_bc <- extract(bio_clim, subs[, c("longitude", "latitude")])
bg_bc <- extract(bio_clim, bg) # for the pseudo-absence points

# create dataframes
# for presence points
ybrev_bc <- data.frame(lon = subs$longitude, lat = subs$latitude, ybrev_bc)
str(ybrev_bc)

# for pseudo-absences
bgpoints <- bg@coords
colnames(bgpoints) <- c("lon", "lat")
bg_bc <- data.frame(cbind(bgpoints, bg_bc))
length(which(is.na(bg_bc$bio1))) # double-check for missing data
bg_bc <- bg_bc[!is.na(bg_bc$bio1), ] # and pull out the missing lines

# to check the effectiveness of a predictive model,
# we can use cross validation approach
# vector of group assignments splitting the Ybrev_bc into 5 groups
group_p <- kfold(ybrev_bc, k = 5)
group_a <- kfold(bg_bc, k = 5) # ditto for bg_bc

# Now we can build species distribution model
# From here you can either use packages or standalone software
# dismo calls Maxent but first need to install javascript maxent.jar

# Since we defined 5 groups, let’s pick a number between 1 at 5 to use
# for our validation set-aside
test <- 3

# Then we use test & the kfold groupings to divide the presence & absence points
train_p <- ybrev_bc[group_p != test, c("lon", "lat")]
train_a <- bg_bc[group_a != test, c("lon", "lat")]
test_p <- ybrev_bc[group_p == test, c("lon", "lat")]
test_a <- bg_bc[group_a == test, c("lon", "lat")]

# build a maxent model using the “training” points and the Bioclim data
# install:
# java (https://www.oracle.com/java/technologies/downloads/#jdk21-windows)
# maxent (https://biodiversityinformatics.amnh.org/open_source/maxent/)
me <- maxent(bio_clim, p = train_p, a = train_a)

# validate the model,
e <- evaluate(p = test_p,
              a = test_a,
              model = me,
              x = bio_clim)
e

# visualize the predictions from our model
pred_me <- predict(object = me,
                   x = bio_clim) # generate the predictions

# make a nice plot
plot(x = pred_me,
     y = 1,
     cex = 0.5,
     legend = TRUE,
     mar = par("mar"),
     xaxt = "n",
     yaxt = "n",
     main = "Predicted presence of Joshua trees")
map("state",
    xlim = c(-119, -113),
    ylim = c(33.5, 38),
    fill = FALSE,
    col = "cornsilk",
    add = TRUE)

# state names
text(x = -117.5, y = 35.5, "California", col = rgb(1, 1, 1, 0.4), cex = 3)
text(x = -116, y = 37.5, "Nevada", col = rgb(1, 1, 1, 0.4), cex = 3)
text(x = -113, y = 34.5, "Arizona", col = rgb(1, 1, 1, 0.4), cex = 3)
text(x = -113, y = 37.75, "Utah", col = rgb(1, 1, 1, 0.4), cex = 3)

# plot the presence points
points(subs$longitude,
       subs$latitude,
       pch = 20,
       cex = 0.5,
       col = "dodgerblue")

# and the pseudo-absence points
points(bg, cex = 0.5, col = "brown")

# add axes
axis(1, las = 1)
axis(2, las = 1)

# restore the box around the map
box()

################################### END ########################################
