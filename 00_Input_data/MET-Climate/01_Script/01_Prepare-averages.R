
#### Prepare MET office climate data ####
#### Tobias Ruettenauer ####
#### 2020/ 08 / 10 ####

rm(list=ls())

### Load packages
library(sf)
library(ncdf4)
library(RCurl)
library(httr)




### Working Directory
setwd("C:/work/Forschung/Daten/MET-Climate/02_Data")


########################################
### Download averages over the years ###
########################################

### Set up connection metadata 
user <- "ruettenauer"
pwd <- "XXX"

### Tasmax
name <- "tasmax"
# Set up connection
url <- "ftp://ftp.ceda.ac.uk/badc/ukmo-hadobs/data/insitu/MOHC/HadOBS/HadUK-Grid/v1.0.2.1/5km/tasmax/mon-20y/v20200731/"
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

# save
assign(paste0(name, "_avg.df"), data)
save(list = (paste0(name, "_avg.df")), file = paste0(name, "_avg.RData"))



### Tasmin
name <- "tasmin"
# Set up connection
url <- "ftp://ftp.ceda.ac.uk/badc/ukmo-hadobs/data/insitu/MOHC/HadOBS/HadUK-Grid/v1.0.2.1/5km/tasmin/mon-20y/v20200731/"
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

# save
assign(paste0(name, "_avg.df"), data)
save(list = (paste0(name, "_avg.df")), file = paste0(name, "_avg.RData"))



### Tasmin
name <- "rainfall"
# Set up connection
url <- "ftp://ftp.ceda.ac.uk/badc/ukmo-hadobs/data/insitu/MOHC/HadOBS/HadUK-Grid/v1.0.2.1/5km/rainfall/mon-20y/v20200731/"
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

# save
assign(paste0(name, "_avg.df"), data)
save(list = (paste0(name, "_avg.df")), file = paste0(name, "_avg.RData"))



### Plot example

tasmax.sp <- st_as_sf(x = tasmax_avg.df, coords = c("lon", "lat"))
tasmax.sp <- tasmax.sp[!is.na(tasmax.sp$tasmax_mon_1), ]
tasmax.sp <- st_set_crs(tasmax.sp, 4326)
tasmax.sp <- st_transform(tasmax.sp, 27700)
tasmax.sp <- st_buffer(tasmax.sp, dist = 2500, nQuadSegs = 1, endCapStyle = "SQUARE")
plot(tasmax.sp[,7], border = NA)

