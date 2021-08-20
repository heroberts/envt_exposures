## GPS DATA PREPROCESSING -----

library(geosphere); library(readr)
library(sp); library(rgdal)
library(sf); library(lwgeom)
library(dplyr); library(ggplot2)
library(lubridate); library(e1071)


setwd("H:/depression/data/temp")

# read in GPS data
# this file contains GPS points of all participants that recorded at least one GPS point (n=629)
points <- st_read("gps_data.csv", stringsAsFactors=FALSE)

# formatting
points$LONGITUDE <- as.numeric(gsub(",", ".", gsub("\\.", "", points$LONGITUDE)))
points$LATITUDE <- as.numeric(gsub(",", ".", gsub("\\.", "", points$LATITUDE)))
points$timestamp <- strptime(x= as.character(points$timestamp), format = "%d-%m-%Y %H:%M:%S")
tz(points$timestamp) <- "UTC"
points$mode <- as.factor(points$mode)
points$SPEED <- as.numeric(gsub(",", ".", gsub("\\.", "", points$SPEED)))
points$DISTANCE <- as.numeric(gsub(",", ".", gsub("\\.", "", points$DISTANCE)))

# order by userid and timestamp
points <- points[with(points, order(points$userid_crypt, points$timestamp)), ]
rownames(points) <- NULL

## STEP ONE: REMOVE PARTICIPANTS WITH FEW DATA POINTS ------
obs.count <- points %>% count(userid_crypt)
obs.count <- obs.count[order(obs.count$n),]
hist(obs.count$n, breaks = 25)
e1071::skewness(obs.count$n) 
e1071::kurtosis(obs.count$n)

# square root transformation and skewness/kurtosis check
obs.count$n.sqrt <- sqrt(obs.count$n)
hist(obs.count$n.sqrt, breaks = 25)
e1071::skewness(obs.count$n.sqrt)
e1071::kurtosis(obs.count$n.sqrt)

# median
median(obs.count$n.sqrt)

# median absolute deviation
mad(obs.count$n.sqrt)

median(obs.count$n.sqrt) - (2.5 * mad(obs.count$n.sqrt))

# removal of participants with fewer data points than 2.5 times the median absolute deviation
obs.count$keep_id <- ifelse(obs.count$n.sqrt <= (median(obs.count$n.sqrt) - (2.5 * mad(obs.count$n.sqrt))), 0, 1)
obs.count <- obs.count[obs.count$keep_id==1, ]

points <- points[points$userid_crypt %in% obs.count$userid_crypt, , drop=FALSE ]

length(unique(points$userid_crypt)) #592 users

rm(obs.count, drop.obs.count, points.dropped, time1, time2, i, gps.duration)

## STEP TWO: REMOVE POINTS WHERE SPEED WAS GREATER THAN 200KM/HR --------
# format GPS coords as sf object
points <- st_as_sf(points, coords = c("LONGITUDE","LATITUDE"), crs=4326)
st_is_longlat(points)

# calculate distance between each GPS point if user_id matches with row above, if not assign distance as 0

for (i in 2:nrow(points)){
  
  points[1, "distance.m"] <- 0
  
  
  if(points[i, "userid_crypt"]== points[i-1, "userid_crypt"]){
    
    points[i, "distance.m"] <- st_distance(points$geometry[i], points$geometry[i-1], by_element = TRUE)
    
  } else {
    
    points[i, "distance.m"] <- 0
    
  }
}

# default units for above are metres, divide by 1000 for km
points$distance.km <- points$distance.m / 1000

# dropping redundant variables
rm(i)


# calculate time interval (in hours) between each point
# for each row, if the observation above matches by user_id, then calculate the difference in timestamps, else assign 0
# note: some people have multiple times that are the same, add a flag for these

for (i in 2:length(points$userid_crypt)){
  
  if(points$userid_crypt[i] == points$userid_crypt[i-1]){
    
    time1 <- points$timestamp[i]
    
    time2 <- points$timestamp[i-1]
    
    points[i, "time.diff.hr"] <- signif(difftime(time1[[1]],time2[[1]], units="hours"), 3)
    
    if(points$timestamp[i] == points$timestamp[i-1]) { 
      
      points$flag[i] <- 1
      
    } else {
      
      points$flag[i] <- 0
      
    }
  } else {
    
    points[i, "time.diff.hr"] <- 0
    points$flag[i] <- 0
    
  }
}

# assign 0 to first rows
points$time.diff.hr[1] <- 0
points$flag[1] <- 0

# remove points flagged with 1 (multiple times the same)
points <- points[!points$flag== 1, ]

points$flag <- NULL

# drop redundant variables
rm(time1, time2, i)

# calculate speed (km per hour)
points$time.diff.hr <- as.numeric(points$time.diff.hr)
points$time.diff.hr <- format(points$time.diff.hr, scientific=F)
points$speed.km.hr <- points$distance.km / points$time.diff.hr
points$speed.km.hr[is.nan(points$speed.km.hr)] <- 0

# write out file with speed calculated to temp folder
st_write(points, "points_speed.csv", layer_options = "GEOMETRY=AS_XY")

# removal of points with speed greater than 200
points <- points[which(points$speed.km.hr <= 200), ]


## STEP THREE: REMOVE POINTS THAT FALL WITHIN 100M ALONG GERMAN/BELGIAN BORDER -----
# read in NL shape
NL.shape <- st_read("H:/shapes", "Gemeentegrenzen_2014_dissolve_simplified") %>%
  st_cast(., "POLYGON") %>%
  st_set_crs(28992) %>%
  slice(27) #use mainland only

# set start and end point of border
start.end.poi <- st_sfc(st_point(x=c(14400, 377600)), 
                        st_point(x=c(276500, 584800)), 
                        crs=28992) %>% 
  st_buffer(., 1000) %>% 
  st_combine()
plot(st_geometry(NL.shape))
plot(st_geometry(start.end.poi), add=T, col="red")

# split line at start/end points
NL.shape.tmp <- st_cast(NL.shape, "MULTILINESTRING", group_or_split=F)
NL.shape.parts <- st_collection_extract(
  st_split(NL.shape.tmp$geometry, start.end.poi), "LINESTRING")
NL.shape.parts <- st_as_sf(data.frame(id=1:length(NL.shape.parts),
                                      geometry=NL.shape.parts))
plot(st_geometry(NL.shape.tmp))
plot(st_geometry(NL.shape.parts[3,]), add=T, col="red")

# remove respondents in 100m buffer 
NL.shape.buf100 <- NL.shape.parts[3,] %>% 
  st_buffer(., 100)

points <- st_transform(points, 28992)

# 1108168 points not in buffer, 878 in buffer
table(st_intersects(NL.shape.buf100, points1, sparse=F)) 

# remove points that are in the 100m buffer
points.buf.rm <- points[!lengths(st_intersects(points, NL.shape.buf100)), ] 

# drop dead variables
rm(start.end.poi, NL.shape.parts, NL.shape.tmp, NL.shape.buf100)

## STEP FOUR: REMOVE PARTICIPANTS THAT LEFT THE NETHERLANDS DURING DATA COLLECTION -----
# points outside NL, including those 100m within DE/BE border
points.out <- points[NL.shape, , op=st_disjoint]

# points that do not have a userid that matches an id that left NL
points.in <- subset(points.buf.rm, !(userid_crypt %in% points.out$userid_crypt))

# check number of participants in NL sample and out NL sample
length(unique(points.in$userid_crypt)) #419 participants, 723217 points, in NL, within buffer, and did not leave NL

# drop dead variables
rm(points, points.buf.rm, NL.shape, points.out)

# write out data
st_write(points.in, "points_in.csv", layer_options="GEOMETRY=AS_XY")


## STEP FIVE: REMOVE POINTS THAT ARE MORE THAN 50METRES AWAY FROM A ROAD/RAILWAY/FOOTPATH/BIKE PATH -----
# read in distance file (computed separately)

gps_distance <- read_csv("H:/depression/data/temp/gps_distance.csv",
                         col_types = cols(timestamp = col_datetime(format = "%Y-%m-%d %H:%M:%S")))

# order by user_id and timestamp
gps_distance <- gps_distance[with(gps_distance, order(gps_distance$userid_crypt, gps_distance$timestamp)), ]

# check for multiple timestamps (as done earlier in preprocessing)
gps_distance <- gps_distance %>% 
  group_by(userid_crypt, timestamp) %>%
  mutate(count = seq(n()))

# remove repeated timestamps
gps_distance <- gps_distance[gps_distance$count < 2, ]
gps_distance$count <- NULL

# participant check
length(unique(gps_distance$userid_crypt)) 

# reformat timestamp because dplyr does not like POSTXlt
points.in$timestamp <- as.POSIXct(points.in$timestamp)

# join distance to travel network calculation to each point
points.in <- left_join(points.in, gps_distance, by=c("WE_ID_crypt" = "we_id_crypt", "persoonsidentificatie_crypt" = "persoonsidentificatie_crypt", 
                                                     "userid_crypt" = "userid_crypt", "timestamp" = "timestamp"))

rm(gps_distance)

# remove points more than 50m from network
points.in <- points.in[which(points.in$distance.to.rail < 50 | points.in$distance.to.bike < 50 | points.in$distance.to.highways < 50 | points.in$distance.to.localroads < 50 | points.in$distance.to.footpath < 50 ), ]


length(unique(points.in$userid_crypt)) # 419 participants remain

# removing unnecessary columns
points.in$distance.to.rail <- NULL
points.in$distance.to.bike <- NULL
points.in$distance.to.highways <- NULL
points.in$distance.to.localroads <- NULL
points.in$distance.to.footpath <- NULL
points.in$SPEED <- NULL
points.in$DISTANCE <- NULL
points.in$distance.m <- NULL
points.in$distance.km <- NULL
points.in$time.diff.hr <- NULL
points.in$speed.km.hr <- NULL
rm(gps_distance)

# Write out final pre-processed data frame ------
st_write(points.in, "points_distance.csv", layer_options="GEOMETRY=AS_XY")
