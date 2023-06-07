
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


##############################
### Download daily maximum ###
##############################

### Tasmax
name <- "tasmax"
# Set up connection
url <- "ftp://ftp.ceda.ac.uk/badc/ukmo-hadobs/data/insitu/MOHC/HadOBS/HadUK-Grid/v1.0.1.0/5km/tasmax/day/v20190808/"
filenames <- getURL(url, userpwd = paste0(user, ":", pwd),
                    ftp.use.epsv = FALSE, dirlistonly = TRUE) 
files <- strsplit(filenames,"\r\n")

# Get files since 1990
years <- gsub(paste0(name, "_hadukgrid_uk_5km_day_"), "", files[[1]])
years <- substr(years, 1, 8)
oo <- which(as.numeric(substr(years,1, 4)) >= 1990)
files <- files[[1]][oo]
years <- years[oo]

# Loop through the files
for(i in 1:length(files)){
  file <-  files[i]
  year <-  as.numeric(substr(years[i], 1, 4))
  mon <-  as.numeric(substr(years[i], 5, 6))
  urlfile <- paste0(url, file)
  
  # Save to drive
  res <- GET(urlfile, write_disk(file, overwrite = TRUE), authenticate(user, pwd))
  nc_data <- nc_open(res$request$output$path)

  # Extract information
  lon <- ncvar_get(nc_data, "longitude")
  lat <- ncvar_get(nc_data, "latitude", verbose = F)
  t <- paste("d", mon, c(1:nc_data$dim$time$len), sep = "_")
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
  names(data) <- paste0(name, "_", t)
  data <- data.frame(lat, lon, data)

  # Drop NAs
  data <- data[which(!is.na(data[,3])), ]

  # save
  fn <- paste0(name, "_", year, "_", mon, ".df")
  assign(fn, data)
  save(list = (fn), file = gsub("df", "RData", fn))
  if(i == 1){
    names.list <- fn
  }else{
    names.list <- c(names.list , fn)
  }
}


### Create long format data

df.list <- vector(mode = "list", length = length(names.list))
k <- 1
for(i in names.list){
  df <- get(i)
  df <- data.frame(id = as.numeric(rownames(df)), df)
  df$year <- as.numeric(gsub(".*?_([0-9]+)_.*", "\\1", i))
  dflong <- pivot_longer(df, cols = names(df[grep(name, names(df))]),
                         names_to = c("mon", "day"),
                        names_pattern = paste0(name, "_d_(.*)_(.*)"),
                        values_to = name)
  dflong$mon <- as.numeric(dflong$mon)
  dflong$day <- as.numeric(dflong$day)
  df.list[[k]] <- dflong
  k <- k + 1
}

# Make one df
weather_uk.df <- data.table::rbindlist(df.list)
save(weather_uk.df, file = "HADUK.RData")

# clean up
ls <- ls()
ls <- ls[grep(paste0(name, "_.*"), ls)]
rm(list = ls)
Sys.sleep(2); gc(); Sys.sleep(2)
ls <- gsub(".df", ".RData", ls)
unlink(ls)


# Free memory
rm(weather_uk.df)
Sys.sleep(2); gc(); Sys.sleep(2)


##############################
### Download daily minimum ###
##############################

### Tasmin
name <- "tasmin"
# Set up connection
url <- "ftp://ftp.ceda.ac.uk/badc/ukmo-hadobs/data/insitu/MOHC/HadOBS/HadUK-Grid/v1.0.1.0/5km/tasmin/day/v20190808/"
filenames <- getURL(url, userpwd = paste0(user, ":", pwd),
                    ftp.use.epsv = FALSE, dirlistonly = TRUE) 
files <- strsplit(filenames,"\r\n")

# Get files since 1990
years <- gsub(paste0(name, "_hadukgrid_uk_5km_day_"), "", files[[1]])
years <- substr(years, 1, 8)
oo <- which(as.numeric(substr(years,1, 4)) >= 1990)
files <- files[[1]][oo]
years <- years[oo]

# Loop through the files
for(i in 1:length(files)){
  file <-  files[i]
  year <-  as.numeric(substr(years[i], 1, 4))
  mon <-  as.numeric(substr(years[i], 5, 6))
  urlfile <- paste0(url, file)
  
  # Save to drive
  res <- GET(urlfile, write_disk(file, overwrite = TRUE), authenticate(user, pwd))
  nc_data <- nc_open(res$request$output$path)
  
  # Extract information
  lon <- ncvar_get(nc_data, "longitude")
  lat <- ncvar_get(nc_data, "latitude", verbose = F)
  t <- paste("d", mon, c(1:nc_data$dim$time$len), sep = "_")
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
  names(data) <- paste0(name, "_", t)
  data <- data.frame(lat, lon, data)
  
  # Drop NAs
  data <- data[which(!is.na(data[,3])), ]
  
  # save
  fn <- paste0(name, "_", year, "_", mon, ".df")
  assign(fn, data)
  save(list = (fn), file = gsub("df", "RData", fn))
  if(i == 1){
    names.list <- fn
  }else{
    names.list <- c(names.list , fn)
  }
}


### Create long format data

df.list <- vector(mode = "list", length = length(names.list))
k <- 1
for(i in names.list){
  df <- get(i)
  df <- data.frame(id = as.numeric(rownames(df)), df)
  df$year <- as.numeric(gsub(".*?_([0-9]+)_.*", "\\1", i))
  dflong <- pivot_longer(df, cols = names(df[grep(name, names(df))]),
                         names_to = c("mon", "day"),
                         names_pattern = paste0(name, "_d_(.*)_(.*)"),
                         values_to = name)
  dflong$mon <- as.numeric(dflong$mon)
  dflong$day <- as.numeric(dflong$day)
  df.list[[k]] <- dflong
  k <- k + 1
}


# clean up
ls <- ls()
ls <- ls[grep(paste0(name, "_.*"), ls)]
rm(list = ls)
rm(df)
rm(dflong)
Sys.sleep(2); gc(); Sys.sleep(2)
ls <- gsub(".df", ".RData", ls)
unlink(ls)


# Make one df
weather_uk2.df <- data.table::rbindlist(df.list)
save(weather_uk2.df, file = "HADUK2.RData")

# Free memory
rm(weather_uk2.df)
Sys.sleep(2); gc(); Sys.sleep(2)




####################################
### Download daily Precipitation ###
####################################

### Tasmin
name <- "rainfall"
# Set up connection
url <- "ftp://ftp.ceda.ac.uk/badc/ukmo-hadobs/data/insitu/MOHC/HadOBS/HadUK-Grid/v1.0.1.0/5km/rainfall/day/v20190808/"
filenames <- getURL(url, userpwd = paste0(user, ":", pwd),
                    ftp.use.epsv = FALSE, dirlistonly = TRUE) 
files <- strsplit(filenames,"\r\n")

# Get files since 1990
years <- gsub(paste0(name, "_hadukgrid_uk_5km_day_"), "", files[[1]])
years <- substr(years, 1, 8)
oo <- which(as.numeric(substr(years,1, 4)) >= 1990)
files <- files[[1]][oo]
years <- years[oo]

# Loop through the files
for(i in 1:length(files)){
  file <-  files[i]
  year <-  as.numeric(substr(years[i], 1, 4))
  mon <-  as.numeric(substr(years[i], 5, 6))
  urlfile <- paste0(url, file)
  
  # Save to drive
  res <- GET(urlfile, write_disk(file, overwrite = TRUE), authenticate(user, pwd))
  nc_data <- nc_open(res$request$output$path)
  
  # Extract information
  lon <- ncvar_get(nc_data, "longitude")
  lat <- ncvar_get(nc_data, "latitude", verbose = F)
  t <- paste("d", mon, c(1:nc_data$dim$time$len), sep = "_")
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
  names(data) <- paste0(name, "_", t)
  data <- data.frame(lat, lon, data)
  
  # Drop NAs
  data <- data[which(!is.na(data[,3])), ]
  
  # save
  fn <- paste0(name, "_", year, "_", mon, ".df")
  assign(fn, data)
  save(list = (fn), file = gsub("df", "RData", fn))
  if(i == 1){
    names.list <- fn
  }else{
    names.list <- c(names.list , fn)
  }
}


### Create long format data

df.list <- vector(mode = "list", length = length(names.list))
k <- 1
for(i in names.list){
  df <- get(i)
  df <- data.frame(id = as.numeric(rownames(df)), df)
  df$year <- as.numeric(gsub(".*?_([0-9]+)_.*", "\\1", i))
  dflong <- pivot_longer(df, cols = names(df[grep(name, names(df))]),
                         names_to = c("mon", "day"),
                         names_pattern = paste0(name, "_d_(.*)_(.*)"),
                         values_to = name)
  dflong$mon <- as.numeric(dflong$mon)
  dflong$day <- as.numeric(dflong$day)
  df.list[[k]] <- dflong
  k <- k + 1
}


# clean up
ls <- ls()
ls <- ls[grep(paste0(name, "_.*"), ls)]
rm(list = ls)
rm(df)
rm(dflong)
Sys.sleep(2); gc(); Sys.sleep(2)
ls <- gsub(".df", ".RData", ls)
unlink(ls)


# Make one df
weather_uk3.df <- data.table::rbindlist(df.list)
save(weather_uk3.df, file = "HADUK3.RData")








