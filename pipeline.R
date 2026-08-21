# Pipeline for Della Penna, Braun, Guinet, O'Callaghan, and Arostegui

#### Libraries & Parameters ####
## Libraries

# require(remotes)
# remotes::install_github("afsc-gap-products/trawllight")

library(moments)
library(maptools)
library(trawllight)
library(data.table)
library(lubridate)
#library(TeleFunc)
library(RNetCDF)
library(dplyr)
library(ggplot2)

## Functions
source("Kprofile_metrics_eTUFF.BGC.R")
source("filter_stepwise_eTUFF.V2.R")
source("calculate_attenuation_eTUFF.BGC.R")
source("Event.R")

## Parameters
# Choose number of days (must be odd #) across which to build light-depth profile composites
composite.length <- 1

# Choose interval (in meters) of depth bins for light-depth/attenuation profiles
bin.size <- 10

# Choose allowable WITHIN PROFILE gap (in meters) for binned light data used to generate attenuation profiles
bin.gap <- 30

# Choose allowable FROM THE SURFACE gap (in meters) for start of binned light data used to generate attenuation profiles
surface.gap <- NULL

# Assign tag manufacturer (determines downstream 'D' value, can be "Wildlife Computers", "Lotek", or "BGC-Argo")
tag.manufacturer <- "Wildlife Computers"

# Specify potential maximum depth (must be >/= max depth of all resultant profiles)
potential.max.depth <- 2000

# Maximum acceptable fraction of points in a profile where the modelled particles component exceeds the micronekton one
th_below<- 0.1

# Threshold on how much the modelled particles component can exceed the micronekton one at max for the profile to pass QC
th_int<-0.05

#### example seal read-in script ####
id <- 'mk10-12A1644-12'
input_dir <- "../Data/L0_meop/"

nc_traj <- open.nc(list.files(input_dir, pattern = paste0(id,"_traj"), full.names = TRUE))

ref_dt <- as.POSIXct('1950-01-01 00:00:00', tz='UTC')
time <- var.get.nc(nc_traj, 'TIME')
time[which(time <= 0)] <- NA
datetime_utc <- as.POSIXct.numeric(time * 24 * 3600, origin=ref_dt, tz='UTC')
light <- var.get.nc(nc_traj, 'LIGHT')

lat <- var.get.nc(nc_traj, 'LATITUDE')
lon <- var.get.nc(nc_traj, 'LONGITUDE')
pres <- try(var.get.nc(nc_traj, 'PRES_ADJUSTED'), TRUE)
temp <- try(var.get.nc(nc_traj, 'TEMP_ADJUSTED'), TRUE)
df <- data.frame(id=rep(id, length.out=length(dt)), datetime_utc, lat, lon, pres, temp, light)
df <- df[which(!is.na(df$light)),]

track <- df %>% mutate(date=as.Date(datetime_utc)) %>%
  group_by(date) %>%
  summarise(mean_lat=mean(lat, na.rm=TRUE), mean_lon = mean(lon, na.rm=TRUE)) %>%
  filter(!is.na(mean_lat)) %>%
  as.data.frame()

df_trim <- df %>% filter(datetime_utc >= as.POSIXct(min(track$date) - 1) &
                           datetime_utc <= as.POSIXct(max(track$date) + 1))

out_dir <- "../Data/L1_csv/"
fwrite(track, paste0(out_dir, id, '_traj-TRACK.csv'))
fwrite(df_trim, paste0(out_dir, id, '_traj-SERIES_trim.csv'))

#p1 <- ggplot(df_trim) + geom_point(aes(x=datetime_utc, y=pres, colour=light)) + ylim(1100,0)
#ggsave(paste0("~/Dropbox/WHOI/NZA/R/DAB/Figures/", id, '_traj-SERIES-PLOT.png'), width=20, height=10, p1)

print(paste0(id, ' complete.'))
rm(nc_traj); rm(track); rm(df); rm(light); gc()

#### ATTENUATION PROFILING ####
dir <- '../Data/L1_csv/'

files <- list.files(dir)
seals <- unique(substr(files,1,15))
#seals[which(substr(seals,1,3) == "mk9")] <- substr(seals[which(substr(seals,1,3) == "mk9")],1,14)

for(z in 1:length(seals)) {
  # Identify archival series (vertical movement) and track (horizontal movement) files
  file.series <- paste(dir,paste0(seals[z],'_traj-SERIES_trim.csv'),sep="")
  file.track <- paste(dir,paste0(seals[z],'_traj-TRACK.csv'),sep="")
  
  # Read in archival series field and standardize nomenclature
  ts <- fread(file.series)
  colnames(ts) <- c("id","DateTime","latitude","longitude","depth","temp","light")
  sealID <- unique(ts$id)
  
  ## Conduct dive phase assignment
  # Calculate vertical velocity between each time step
  ts$VV <- c(0, diff(ts$depth))
  
  # Determine temporal resolution of the data (1, 2, or 10 seconds) -- ***CUSTOM TO THESE DATA***
  if (length(which(ts$DateTime == unique(ts$DateTime)[2])) == 84 | length(which(ts$DateTime == unique(ts$DateTime)[2])) == 85) {res = 2}
  if (length(which(ts$DateTime == unique(ts$DateTime)[2])) == 16 | length(which(ts$DateTime == unique(ts$DateTime)[2])) == 17) {res = 10}
  if (length(which(ts$DateTime == unique(ts$DateTime)[2])) == 169) {res = 1}
  
  # Identify descent events in the archival series
  devent = event(ts$VV > .01, duration = c(120 / res, NA), ends=2)
  
  for (i in 1:nrow(devent)) {
    ind <- devent$start[i]:devent$end[i]
    if (i == 1) {
      ind.d <- ind
    } else {
      ind.d <- c(ind.d,ind)
    }
  }
  
  # Subset data to only descent events
  ts <- ts[ind.d,]
  
  #### Generate K profiles ####
  # Remove entries with negative depth (i.e., above surface) readings
  ts <- ts[which(ts$depth>=0),] 
  
  # Read in track file and assign spatial coordinates to archival series
  track <- fread(file.track)
  ts$mean_lat <- NA; ts$mean_lon <- NA
  for (i in 1:nrow(track)){
    ind.date <- which(as.Date(ts$DateTime) == track$date[i])
    ts$mean_lat[ind.date] <- track$mean_lat[i]
    ts$mean_lon[ind.date] <- track$mean_lon[i]
  }
  
  # Tap into code pipeline and filter out missing entries
  df_etuff <- ts
  
  # Remove any rows with NA locations, depth, or light
  df_etuff <- df_etuff[which(is.na(df_etuff$longitude)==F&is.na(df_etuff$latitude)==F&is.na(df_etuff$depth)==F&is.na(df_etuff$light)==F),]
  
  
  ## Assign solar noon
  # Change hms (hours, minutes, seconds) of tag input date/time to seconds [0-86399]
  df_etuff$time.seconds <- period_to_seconds(hms(format(df_etuff$DateTime,"%H:%M:%S")))
  
  # Identify solar noon in UTC  and change its hms to seconds [0-86399]
  df_etuff$solarnoon <- period_to_seconds(hms(format(
    maptools::solarnoon(crds = cbind(df_etuff$mean_lon,df_etuff$mean_lat), 
                        dateTime = df_etuff$DateTime, POSIXct.out=T)[,2],"%H:%M:%S")))
  
  # Identify +/- 1 hr window around solar noon (accounting for 23:59-00:00 transitions)
  df_etuff$solarnoon.lower <- ifelse(df_etuff$solarnoon - 60*60 < 0, df_etuff$solarnoon - 60*60 + 86400, df_etuff$solarnoon - 60*60)
  df_etuff$solarnoon.upper <-ifelse(df_etuff$solarnoon + 60*60 >= 86400, df_etuff$solarnoon + 60*60 - 86400, df_etuff$solarnoon + 60*60)
  
  # Assign data within 1 hr of solar noon
  df_etuff$noonwindow <- ifelse (df_etuff$solarnoon.upper > df_etuff$solarnoon.lower, 
                                 ifelse(df_etuff$time.seconds >= df_etuff$solarnoon.lower & df_etuff$time.seconds <= df_etuff$solarnoon.upper, "noon", "notnoon"), 
                                 ifelse(df_etuff$time.seconds >= df_etuff$solarnoon.lower | df_etuff$time.seconds <= df_etuff$solarnoon.upper, "noon", "notnoon"))
  
  #### Composite Profiling ####
  # Calculate number of days at start/end of track to exclude from being a composite center
  composite.exclude <- composite.length/2 - 0.5
  
  # Create vector of sorted dates (without times)
  unique.days <- sort(unique(as.Date(df_etuff$DateTime)))
  
  # Identify spatial coordinates during specified time window of each day
  latitude <- NA
  longitude <- NA
  df_etuff$longitude[which(df_etuff$longitude<0)] <- df_etuff$longitude[which(df_etuff$longitude<0)] + 360
  for (i in 1:length(unique.days)){
    ind <- which(as.Date(df_etuff$DateTime)==unique.days[i] & df_etuff$noonwindow=="noon")
    latitude[i] <- mean(df_etuff$latitude[ind])
    longitude[i] <- mean(df_etuff$longitude[ind])
  }
  
  # Assign tag-specific light floor with buffer
  light.floor <- min(df_etuff$light,na.rm=T) + 10
  
  # Clear up memory
  rm(ts)
  gc()
  
  
  ## Build daily light attenuation profiles
  for (i in (1+composite.exclude):(length(unique.days)-composite.exclude)){
    
    # Subset time series archive within 1 hr of local noon on dates of composite window 
    for (j in 1:composite.length){
      date <- unique.days[i-composite.exclude+j-1]
      noon.index <- which(as.Date(df_etuff$DateTime)==date & df_etuff$noonwindow=="noon")
      if (j == 1){
        noon.data <- df_etuff[noon.index,]
      } else {
        noon.data <- rbind.data.frame(noon.data,df_etuff[noon.index,])
      }
    }
    
    if (nrow(noon.data) == 0){
      # Do nothing
    } else{
      
      # Bin/filter input raw light data and conduct quality control (pre)
      profile.binned <- filter_stepwise_eTUFF.V2(cast.data = noon.data,
                                                 light.col = "light",
                                                 depth.col = "depth",
                                                 bin.size = bin.size,
                                                 bin.gap = bin.gap,
                                                 surface.gap = surface.gap,
                                                 agg.fun = geometric.mean,
                                                 light.floor = light.floor,
                                                 filter = T)
      
      # Build LOESS light profile, calculate derived K values, and conduct quality control (post)
      output <- calculate_attenuation_eTUFF.BGC(x = profile.binned,
                                            light.col = "light",
                                            depth.col = "depth",
                                            loess.criterion = "aicc",
                                            loess.degree = 1,
                                            kz.binsize = bin.size,
                                            bin.gap = bin.gap,
                                            min.range = 200,
                                            light.predict = T,
                                            tag.manufacturer = tag.manufacturer,
                                            plot=F)
      if (class(output)!="list"){
        # Do nothing
      } else{
        # Calculate K profile metrics and add to output
        Kmetrics <- Kprofile_metrics_eTUFF.BGC(output)
        output[["Kmetrics"]] <- Kmetrics
        
        # Add metadata to output
        output$meta[["date"]] <- unique.days[i]
        output$meta[["latitude"]] <- latitude[i]
        output$meta[["longitude"]] <- longitude[i]
        output$meta[["instrument_name"]] <- sealID
        
        # Save output with attenuation profile, light profile, K metrics, and metadata
        saveRDS(output,
                file = paste("../Data/L2_profiles/",sealID,"_",unique.days[i],".rds",sep=""))
      } 
    }
  }
}


#### Create dataframe of all K profiles, depth matched ####
# Identify directory containing daily profiles
profiles <- list.files("../Data/L2_profiles/")

# Read in profile files and extract requisite information
for (i in 1:length(profiles)){
  profile <- readRDS(paste("../Data/L2_profiles/",profiles[i],sep=""))
  
  # Spatial coordinates
  latitude <- profile$meta$latitude
  longitude <- profile$meta$longitude
  
  # Date
  date <- profile$meta$date
  
  # Unique ID of the tag
  instrument_name <- profile$meta$instrument_name
  
  # Tag Manufacturer
  tag.manufacturer <- as.character(profile$meta$tag.manufacturer)
  
  # Min, max, and mean depth of the profile
  min.depth <- min(profile$attenuation$depth)
  max.depth <- max(profile$attenuation$depth)
  mean.depth <- (min.depth + max.depth) / 2
  
  # Mean absolute error of the profile's LOESS light model
  mae <- profile$meta$mae
  
  # Depth-specific values of light attenuation
  kdz_extract <- profile$attenuation$kdz
  kdz <- rep(NA,potential.max.depth/bin.size)
  kdz[which(sprintf("%sm",seq(bin.size,potential.max.depth,bin.size))==paste(min.depth,"m",sep="")):which(sprintf("%sm",seq(bin.size,potential.max.depth,bin.size))==paste(max.depth,"m",sep=""))] <- 
    kdz_extract
  
  if (i == 1){
    profile_info <- c(instrument_name,tag.manufacturer,date,latitude,longitude,min.depth,max.depth,mean.depth,mae,
                      kdz)
  } else {
    add_info <- c(instrument_name,tag.manufacturer,date,latitude,longitude,min.depth,max.depth,mean.depth,mae,
                  kdz)
    profile_info <- rbind(profile_info,add_info)
  }
}

# Assign data frame and standardize nomenclature
profile_info <- as.data.frame(profile_info,stringsAsFactors = FALSE)
colnames(profile_info) <- c("instrument_name","tag.manufacturer","date","lat","lon","min","max","mean","mae",
                            sprintf("%sm",seq(bin.size,potential.max.depth,bin.size)))

# Assign date format
profile_info$date <- as.Date(as.numeric(profile_info$date),origin='1970-01-01')

# Assign numeric format to coordinates and profile metrics/kdz
profile_info[,4:ncol(profile_info)] <- sapply(profile_info[,4:ncol(profile_info)],as.numeric)

# Filter profiles with high MAE (> 0.1 * units/decade)
I <- profile_info$mae > 2 & profile_info$tag.manufacturer=="Wildlife Computers"|profile_info$mae > 3.2 & profile_info$tag.manufacturer=="Lotek"
profile_info <- profile_info[!I,]

# Specify data frame with depth-matched light attenuation values
profile_kdz <- profile_info; profile_kdz$type <- "kd_total"

#### Create dataframe of all Light (not K) profiles, depth matched ####
# Identify directory containing daily profiles
profiles <- list.files("../Data/L2_profiles/")

# Read in profile files and extract requisite information
for (i in 1:length(profiles)){
  profile <- readRDS(paste("../Data/L2_profiles/",profiles[i],sep=""))
  
  # Spatial coordinates
  latitude <- profile$meta$latitude
  longitude <- profile$meta$longitude
  
  # Date
  date <- profile$meta$date
  
  # Unique ID of the tag
  instrument_name <- profile$meta$instrument_name
  
  # Tag Manufacturer
  tag.manufacturer <- as.character(profile$meta$tag.manufacturer)
  
  # Min, max, and mean depth of the profile
  min.depth <- min(profile$attenuation$depth)
  max.depth <- max(profile$attenuation$depth)
  mean.depth <- (min.depth + max.depth) / 2
  
  # Mean absolute error of the profile's LOESS light model
  mae <- profile$meta$mae
  
  # Depth-specific values of light attenuation
  kdz_extract <- profile$attenuation$predict_light
  kdz <- rep(NA,potential.max.depth/bin.size)
  kdz[which(sprintf("%sm",seq(bin.size,potential.max.depth,bin.size))==paste(min.depth,"m",sep="")):which(sprintf("%sm",seq(bin.size,potential.max.depth,bin.size))==paste(max.depth,"m",sep=""))] <- 
    kdz_extract
  
  if (i == 1){
    profile_info <- c(instrument_name,tag.manufacturer,date,latitude,longitude,min.depth,max.depth,mean.depth,mae,
                      kdz)
  } else {
    add_info <- c(instrument_name,tag.manufacturer,date,latitude,longitude,min.depth,max.depth,mean.depth,mae,
                  kdz)
    profile_info <- rbind(profile_info,add_info)
  }
}

# Assign data frame and standardize nomenclature
profile_info <- as.data.frame(profile_info,stringsAsFactors = FALSE)
colnames(profile_info) <- c("instrument_name","tag.manufacturer","date","lat","lon","min","max","mean","mae",
                            sprintf("%sm",seq(bin.size,potential.max.depth,bin.size)))

# Assign date format
profile_info$date <- as.Date(as.numeric(profile_info$date),origin='1970-01-01')

# Assign numeric format to coordinates and profile metrics/kdz
profile_info[,4:ncol(profile_info)] <- sapply(profile_info[,4:ncol(profile_info)],as.numeric)

# Filter profiles with high MAE (> 0.1 * units/decade)
I <- profile_info$mae > 2 & profile_info$tag.manufacturer=="Wildlife Computers"|profile_info$mae > 3.2 & profile_info$tag.manufacturer=="Lotek"
profile_info <- profile_info[!I,]

# Specify data frame with depth-matched light level values
profile_light <- profile_info; profile_light$type <- "light level"

# Join depth-matched light attenuation and light level data frames
profile_both <- rbind(profile_light,profile_kdz)

# Write out the joint output
write.csv(profile_both,"../Data/L3_output/testsealoutput.csv",row.names=F)


#### Parsing Contributions to Attenuation ####
data <- as.data.frame(fread("../Data/L3_output/testsealoutput.csv"))
data[,10:(ncol(data)-1)] <- as.numeric(unlist(data[,10:(ncol(data)-1)]))

if (tag.manufacturer == "Wildlife Computers") {
  D <- 20
} else if (tag.manufacturer == "Lotek") {
  D <- 32
} else if (tag.manufacturer == "BGC-Argo") {
  D <- 10
} else {
  warning("Your tag manufacturer is not yet recognized by this function, which will need to be updated to include its D value.")
  return(NA)
}

tags <- unique(data$instrument_name) # for each tag
profile <- list() # it creates an empty list
all_good <-list()
qq <- 1 # and initializes a counter to 1

for (tag in tags) { # then for each tag
  print(paste("Processing tag:", tag))
  
  idx <- which(data$instrument_name == tag) # we identify the rows belonging to that tag
  days <- unique(data$date[idx]) # and we list all the days on that list
  
  for (day in days) { # then for each day 
    print(paste("Processing day:", day))
    
    idd <- which(data$date[idx] == day) # we find which of that given tag indices belong to that day
    idxx_light <- which(data$type[idx[idd]] == "light level") # with the lines having the light profile
    idxx_kd <- which(data$type[idx[idd]] == "kd_total") # and the attenuation profile
    
    profile[[qq]] <- list( # we use this to create a profile
      depth = seq(bin.size, potential.max.depth, by = bin.size),
      light = as.matrix(data[idx[idd[idxx_light]], 10:(ncol(data)-1)]),
      kd = as.matrix(data[idx[idd[idxx_kd]], 10:(ncol(data)-1)]),
      lon = data[idx[idd[idxx_kd]], 5],
      lat = data[idx[idd[idxx_kd]], 4]
    )
    
    # Flag unrealistic values as NA
    profile[[qq]]$light[profile[[qq]]$light < -900] <- NA
    profile[[qq]]$kd[profile[[qq]]$kd < -900] <- NA
    
    # Identify the euphotic zone
    id_clear <- which(!is.na(profile[[qq]]$light)) # we locate the valid points along the profile using light
    if (length(id_clear) > 0) { # if there is at least one
      dz <- profile[[qq]]$depth[id_clear[1]] # we define the depth of the first valid observation in the profile
      
      # Interpolate up the water column to identify the surface value
      surface_value <- profile[[qq]]$light[id_clear[1]] - D * log(1 - dz * profile[[qq]]$kd[id_clear[1]]) 
      
      # Find the depth values where we are below 1% of surface light
      id_dark <- which(profile[[qq]]$light < surface_value - 2*D)
     
      # Ff there are no dark values, it assumes that the index of the euphotic depth it will be the first one, otherwise it will pick the last depth before the dark condition is satisfied
      idx_zeu <- ifelse(length(id_dark) > 0, id_dark[1] - 1, 1) 
      
      if (is.na(profile[[qq]]$kd[idx_zeu])) {
        nnan <- which(!is.na(profile[[qq]]$kd))
        idx_zeu <- nnan[1]
      }
      
      # We then create a kd profile where the observations from the euphotic layers are removed
      profile[[qq]]$kd_filt <- profile[[qq]]$kd
      profile[[qq]]$kd_filt[1:idx_zeu] <- NA
      temp <- which(diff(as.numeric(profile[[qq]]$kd_filt)) > 0) # this quite a difference from findpeaks in matlab
      
      # Now we create the estimate particle contribution profile
      first <- idx_zeu
      last <- ifelse(length(temp) > 0, temp[1], which.min(profile[[qq]]$kd)) # here we locate the depth of the local minimum of kd
      
      z1 <- profile[[qq]]$depth[first]
      kdz1 <- profile[[qq]]$kd[first]
      z2 <- profile[[qq]]$depth[last]
      kdz2 <- profile[[qq]]$kd[last]
      
      lambda <- (1 / (z2 - z1)) * log(kdz1 / kdz2)
      profile[[qq]]$retrieved <- kdz1 * exp(- (profile[[qq]]$depth - z1) * lambda)
      
      ## Quality control
      # Compute the differences
      diffs <- profile[[qq]]$kd_filt - profile[[qq]]$retrieved[1:length(profile[[qq]]$kd_filt)]
      
      # Count the negative values below threshold
      num_neg <- length(which(diffs < -th_int))
      
      # Count the total non-NA values
      total_pts <- length(which(!is.na(profile[[qq]]$kd)))
      
      # Condition check
      if (num_neg / total_pts < th_below) {
        profile[[qq]]$pass <- 1
        all_good <- c(all_good, qq)  # Append qq to all_good so that we have a list of how many good profiles we have
      } else {
        profile[[qq]]$pass <- 0
      }
      
      profile[[qq]]$kd_micronekton <- profile[[qq]]$kd_filt - profile[[qq]]$retrieved[1:length(profile[[qq]]$kd_filt)]
      
      print(paste("Profile generated:", qq))
      
      qq <- qq + 1
    }
  }
}

## Join depth-matched kd_micronekton and kd_particles to data frame of kd_total and light level
profile_kdparsed <- data[c(1:(nrow(data)/2)),]; profile_kdparsed[,10:(ncol(data)-1)] <- NA; profile_kdparsed[,"type"] <- "kd_micronekton"

for (i in 1:nrow(profile_kdparsed)){
  profile_kdparsed[i,10:(ncol(profile_kdparsed)-1)] <- as.data.frame(profile[[i]]$kd_micronekton)
  if (profile[[i]]$pass == 0) {profile_kdparsed[i,c(10:(ncol(profile_kdparsed)-1))] <- NA}
}

data_all <- rbind(data,profile_kdparsed)

write.csv(data_all,"../Data/L3_output/final_output.csv",row.names=F)
