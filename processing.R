## GPS DATA PROCESSING -----

library(sf); library(sp)
library(doParallel); library(foreach)
library(stringr); library(data.table)
library(magrittr); library(haven)
library(tidyverse); library(raster)
library(foreign); library(readstata13)
library(dplyr); library(diverse)
library(nngeo); library(lubridate)

setwd("H:/depression/data/temp")

rasterOptions(tmpdir="H:/depression/results/trash")

## Load final pre-processed GPS data ------------
points.in <- st_read("points_distance.csv", stringsAsFactors=FALSE, options=c("X_POSSIBLE_NAMES=X", "Y_POSSIBLE_NAMES=Y")) %>%
  st_set_crs(28992)

# formatting df
points.in$WE_ID_crypt <- as.factor(points.in$WE_ID_crypt)
points.in$persoonsidentificatie_crypt <- as.factor(points.in$persoonsidentificatie_crypt)
points.in$userid_crypt <- as.factor(points.in$userid_crypt)
points.in$tripId_crypt <- as.factor(points.in$tripId_crypt)
points.in$mode <- as.factor(points.in$mode)
points.in$timestamp <- as.POSIXct(strptime(points.in$timestamp, "%Y/%m/%d %H:%M:%S"), tz= "UTC")
points.in$timestamp <- with_tz(points.in$timestamp, tzone="Europe/Amsterdam")

points.in$timestamp.char <- as.character(points.in$timestamp)

## Load environmental data --------
# air pollution
airpath <- "H:/buffers/data/pollution/"
airfiles <- dir(path = airpath, pattern = ".tif$")
airnames <- substr(airfiles, 1, 3)
airstack <- stack(paste0(airpath, airfiles))
##airstack <- projectRaster(airstack, crs=CRS(st_crs(muni)$proj4string))

# water
watpath <- "H:/buffers/data/LGN2018blue/"
watfiles <- dir(path = watpath, pattern = "lgn2018_blue_space.tif$")
watnames <- substr(watfiles, 1, 18)
watstack <- raster(paste0(watpath, watfiles))

# ndvi
grepath <- "H:/buffers/data/landsat/"
grefiles <- dir(path = grepath, pattern = "NDVI2018.tif$")
grenames <- substr(grefiles, 1, 8)
grestack <- raster(paste0(grepath, grefiles))
##grestack <- projectRaster(grestack, crs=CRS(st_crs(muni)$proj4string))

# noise
noipath <- "H:/buffers/data/noise/"
noifiles <- dir(path = noipath, pattern = "noise")
noinames <- substr(noifiles, 1, 19)
noistack <- raster(paste0(noipath, noifiles))
#noistack <- projectRaster(noistack, crs=CRS(st_crs(muni)$proj4string))
table(getValues(noistack))
reclass.noise <- c(1, 45,
                   2, 50,
                   3, 55,
                   4, 60,
                   4.5, 65,
                   5, 65,
                   6, 70,
                   7, 75,
                   8, 80)
reclass.noise <- matrix(reclass.noise, ncol = 2, byrow = T)
noistack <- reclassify(noistack, reclass.noise)

## Loop: Environmental exposures per GPS point  -------------

bufw <- 100
#50m

# 419 is number of participants
for(j in levels(points.in$userid_crypt)[1:419]){
  
  ppt.j <- points.in[points.in$userid_crypt == j, ]
  
  # LOOP PER RESPONDENT
  
  cores <- detectCores()
  cl <- makeCluster(cores,useXDR=F)
  registerDoParallel(cl)
  myresults <- foreach(i = 1:nrow(ppt.j),
                       .packages= c("sf", "dplyr", "data.table", "stringr", 
                                    "magrittr", "raster", "vegan"), 
                       .combine = function(...) rbindlist(list(...)), 
                       .multicombine = T) %dopar%
    {
      
      # select respondent i + buffer
      pointsi <- ppt.j[i, ]
      userID_crypt <- pointsi$userid_crypt
      persoonsidentificatie_crypt <- pointsi$persoonsidentificatie_crypt
      WE_ID_crypt <- pointsi$WE_ID_crypt
      timestamp.char <- pointsi$timestamp.char
      bufi <- st_buffer(st_geometry(pointsi), dist = bufw)
      bufisp <- as(st_geometry(bufi), "Spatial")
      
      
      # AIR POLLUTION
      
      air.crop.buf <- crop(airstack, extent(bufisp))
      air.mask.buf <- mask(air.crop.buf, bufisp)
      airlist <- list()
      for(imo in 1:nlayers(air.mask.buf)){
        airlist[[imo]] <- mean(getValues(air.mask.buf[[imo]]), na.rm = T)
      }
      airmat <- matrix(unlist(airlist, use.names = F),
                       ncol = nlayers(air.mask.buf), byrow = F)
      colnames(airmat) <- airnames
      RESair <- as.data.frame(airmat)
      
      # NDVI
      
      gre.crop <- crop(grestack, extent(bufisp))
      gre.mask <- mask(gre.crop, bufisp)
      gre.mask[gre.mask < 0] <- NA
      RESgre <- as.data.frame(mean(getValues(gre.mask), na.rm = T))
      colnames(RESgre) <- grenames
      
      # LGN18 BLUE
      
      wat.crop <- raster::crop(watstack, extent(bufisp))
      wat.crop <- calc(wat.crop, fun = function(x) ifelse(is.na(x), 0, 1))
      wat.mask <- mask(wat.crop, bufisp)
      wat.cells <- sum(!is.na(getValues(rasterize(bufisp, wat.mask))))
      wat.count <- sum(getValues(wat.mask), na.rm = T)
      if (is.na(wat.count) | wat.count == 0) {
        RESwat <- 0
      } else {
        RESwat <- as.vector((wat.count/wat.cells)*100)
      }
      
      # NOISE
      
      noise.crop.buf <- crop(noistack, extent(bufisp))
      noise.mask.buf <- mask(noise.crop.buf, bufisp)
      RESnoi <- table(getValues(noise.mask.buf)) %>%
        as.data.frame() %>%
        na.omit() %>%
        mutate(noi.class = as.numeric(as.character(Var1))) %>%
        mutate(., cellprop = (Freq / sum(Freq))) %>%
        mutate(noisew = cellprop * noi.class) %>%
        summarise(NOIw = sum(noisew))
      
      # collect results
      
      data.table(WE_ID_crypt, persoonsidentificatie_crypt, userID_crypt, timestamp.char, RESgre, RESwat, RESnoi, RESair)
      
    }
  stopCluster(cl)
  
  ## Write results per participant ---------- 
  # note: check folder is correct (50/100)
  
  exppath <- paste0("H:/depression/results/gps/buf100m/", j, "_buf100", ".csv")
  
  fwrite(myresults, exppath, sep=",", dec=".", nThread=getDTthreads())
  
  print(j)
}


## Put files together ---------

getwd()
readpath <- "H:/depression/results/gps/buf100m/"
envfiles <- list.files(readpath, pattern = ".csv", full.names = F)
envfiles <- str_sub(envfiles, end=-11)


comb.env <- NULL
for(i in envfiles){
  print(i)
  tmpimp <- fread(paste0(readpath, i, "buf100.csv", sep = ""), sep = ";", dec = ",", header = T)
  comb.env <- rbindlist(list(comb.env, tmpimp), use.names = T, fill=F) 
}
dim(comb.env)
summary(comb.env)

fwrite(comb.env, "H:/depression/results/gps/buf100m/linked_gps_buf100m.csv", 
       sep = ";", dec=",", nThread = getDTthreads())
