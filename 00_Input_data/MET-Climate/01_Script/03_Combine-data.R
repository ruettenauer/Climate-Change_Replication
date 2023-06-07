
#### Prepare MET office climate data ####
#### Tobias Ruettenauer ####
#### 2020/ 08 / 10 ####

rm(list=ls())

### Load packages
library(sf)
library(ncdf4)
library(RCurl)
library(httr)
library(tidyr)




### Working Directory
setwd("C:/work/Forschung/Daten/MET-Climate/02_Data")


########################################
### Download averages over the years ###
########################################

### Set up connection metadata 
user <- "ruettenauer"
pwd <- "XX"


###############################
### Gen spatial data object ###
###############################


### Set up connection metadata 
user <- "ruettenauer"
pwd <- "wasserschaden"

### Tasmax
name <- "tasmax"
# Set up connection
url <- "ftp://ftp.ceda.ac.uk/badc/ukmo-hadobs/data/insitu/MOHC/HadOBS/HadUK-Grid/v1.0.1.0/5km/tasmax/mon-20y/v20190808/"
filenames <- getURL(url, userpwd = paste0(user, ":", pwd),
                    ftp.use.epsv = FALSE, dirlistonly = TRUE) 
files <- strsplit(filenames,"\r\n")

# Choose file
file <-  files[[1]][1]
urlfile <- paste0(url, files[[1]][1])

# Save to drive
res <- GET(urlfile, write_disk(file, overwrite = TRUE), authenticate(user, pwd))
nc_data <- nc_open( res$request$output$path)

# Extract information
lon <- ncvar_get(nc_data, "longitude")
lat <- ncvar_get(nc_data, "latitude", verbose = F)
t <- ncvar_get(nc_data, "month_number")
data.array <- ncvar_get(nc_data, name)
fillvalue <- nc_data$var[[name]]$missval

# Close connection
nc_close(nc_data) 
unlink(file)

# Replace NAs
oo <- which(data.array >= fillvalue)
data.array[oo] <- NA

# Bind to dataset
lon <- c(lon)
lat <- c(lat)
data <- data.frame(matrix(data.array, ncol = dim(data.array)[3], byrow = FALSE))
names(data) <- paste0(name, "_mon_", t)
data <- data.frame(lat, lon, data)


### Generate spatial object
data <- data.frame(id = rownames(data), data)
data <- data[!is.na(data$tasmax_mon_1), ]

haduk.sp <- st_as_sf(x = data[, c("id", "lat", "lon")], coords = c("lon", "lat"))
haduk.sp <- st_set_crs(haduk.sp, 4326)
haduk.sp <- st_transform(haduk.sp, 27700)
haduk.sp <- st_buffer(haduk.sp, dist = 2500, nQuadSegs = 1, endCapStyle = "SQUARE")

# Save
save(haduk.sp, file = "HADUK_sp.RData")

rm(data.array)
Sys.sleep(2); gc(); Sys.sleep(2)


##################################
### Combine different measures ###
##################################

load("HADUK.RData")
load("HADUK2.RData")
tasmin <- weather_uk2.df$tasmin
rm(weather_uk2.df)
Sys.sleep(2); gc(); Sys.sleep(2)
load("HADUK3.RData")
rainfall <- weather_uk3.df$rainfall
rm(weather_uk3.df)
Sys.sleep(2); gc(); Sys.sleep(2)

# Add to df
weather_uk.df$tasmin <- tasmin
rm(tasmin)
Sys.sleep(2); gc(); Sys.sleep(2)
weather_uk.df$rainfall <- rainfall
rm(rainfall)
Sys.sleep(2); gc(); Sys.sleep(2)

weather_uk.df <- data.frame(weather_uk.df)

# Save
save(weather_uk.df, file = "HADUK_all.RData")
# unlink(c("HADUK.RData", "HADUK2.RData", "HADUK3.RData"))
Sys.sleep(2); gc(); Sys.sleep(2)


#############################
### Combine with averages ###
#############################

load("rainfall_avg.RData")
load("tasmax_avg.RData")
load("tasmin_avg.RData")

# Add ids
tasmax_avg.df$id <- as.numeric(rownames(tasmax_avg.df))
tasmin_avg.df$id <- as.numeric(rownames(tasmin_avg.df))
rainfall_avg.df$id <- as.numeric(rownames(rainfall_avg.df))

# Drop NAs
tasmax_avg.df <- tasmax_avg.df[!is.na(tasmax_avg.df$tasmax_mon_1), ]
tasmin_avg.df <- tasmin_avg.df[!is.na(tasmin_avg.df$tasmin_mon_1), ]
rainfall_avg.df <- rainfall_avg.df[!is.na(rainfall_avg.df$rainfall_mon_1), ]

# To long
tasmax_avg.df <- pivot_longer(tasmax_avg.df, cols = tasmax_mon_1:tasmax_mon_12,
                       names_to = c("mon"),
                       names_pattern = paste0("tasmax_mon_(.*)"),
                       values_to = "tasmax_avg")
tasmax_avg.df$mon <- as.numeric(tasmax_avg.df$mon)

tasmin_avg.df <- pivot_longer(tasmin_avg.df, cols = tasmin_mon_1:tasmin_mon_12,
                              names_to = c("mon"),
                              names_pattern = paste0("tasmin_mon_(.*)"),
                              values_to = "tasmin_avg")
tasmin_avg.df$mon <- as.numeric(tasmin_avg.df$mon)

rainfall_avg.df <- pivot_longer(rainfall_avg.df, cols = rainfall_mon_1:rainfall_mon_12,
                              names_to = c("mon"),
                              names_pattern = paste0("rainfall_mon_(.*)"),
                              values_to = "rainfall_avg")
rainfall_avg.df$mon <- as.numeric(rainfall_avg.df$mon)

### Merge averages
avg.df <- tasmax_avg.df
avg.df <- merge(avg.df, tasmin_avg.df, by = c("id", "lon", "lat", "mon"))
avg.df <- merge(avg.df, rainfall_avg.df, by = c("id", "lon", "lat", "mon"))


### Merge averages to full df using dplyr inner join
Sys.sleep(2); gc(); Sys.sleep(2)
avg.df <- avg.df[, which(!names(avg.df) %in% c("lon", "lat"))]
weather_uk.df <- dplyr::inner_join(weather_uk.df, avg.df, by = c("id", "mon"))


# Save
save(weather_uk.df, file = "HADUK_all.RData")



