
#### Combine Defra Flood Reports and LSOA ####
#### Tobias Ruettenauer ####
#### 2020/ 01 / 31 ####

rm(list=ls())

### Load packages
library(rgdal)
library(sf)
library(spdep)
library(rgeos)
library(doParallel)
library(foreign)
library(GISTools)
library(cleangeo)

library(tmap)
library(tmaptools)
library(extrafont)
loadfonts()


### Working Directory
setwd("C:/work/Forschung/Climate Change_Replication/02_Data")

# Directory for weather data
weada <- "../00_Input_data/MET-Climate/02_Data/"



#########################
### Load weather data ###
#########################

load(paste0(weada, "HADUK_all.RData"))
load(paste0(weada, "HADUK_sp.RData"))

### Add date
weather_uk.df$date <- as.Date(paste(weather_uk.df$year, weather_uk.df$mon, weather_uk.df$day, sep = "-"),
                            format = "%Y-%m-%d")

### Order
weather_uk.df <- weather_uk.df[order(weather_uk.df$id, weather_uk.df$date), ]



#################################
### Compute weather anomalies ###
#################################

# Compute difference to long-term average

weather_uk.df$tasmax_diff <- weather_uk.df$tasmax - weather_uk.df$tasmax_avg
weather_uk.df$tasmin_diff <- weather_uk.df$tasmin - weather_uk.df$tasmin_avg
weather_uk.df$rainfall_diff <- weather_uk.df$rainfall - weather_uk.df$rainfall_avg

hist(weather_uk.df$tasmax_diff)
hist(weather_uk.df$tasmin_diff)
hist(weather_uk.df$rainfall_diff)

summary(weather_uk.df$tasmax_diff)
sd(weather_uk.df$tasmax_diff)
summary(weather_uk.df$tasmin_diff)
sd(weather_uk.df$tasmin_diff)
summary(weather_uk.df$rainfall_diff)
sd(weather_uk.df$rainfall_diff)


hist(weather_uk.df$tasmax)

summary(weather_uk.df$tasmax)
sd(weather_uk.df$tasmax)


### Heatwave, use 3 different rules (29, 31 and 2 sd above long run average)
buf <- c("1", "2", "3")

for(b in buf){
  v1 <- paste0("heatwave", b)
  weather_uk.df[, v1] <- 0
  
  # 3 rules
  if(b == 1){
    oo <- which(weather_uk.df$tasmax >= 29)
  }
  if(b == 2){
    oo <- which(weather_uk.df$tasmax >= 30)
  }
  if(b == 3){
    h2sd <- 2 * sd(weather_uk.df$tasmax)
    oo <- which(weather_uk.df$tasmax_diff >= h2sd)
  }

  weather_uk.df[oo, v1] <- 1
}  

# Save subsample as df
oo <- which(weather_uk.df$heatwave1 == 1 | weather_uk.df$heatwave2 == 1 | weather_uk.df$heatwave3 == 1)
heatwave.df <- weather_uk.df[oo, ]
table(heatwave.df$year)

# Clean up to save memory
N <- nrow(weather_uk.df)
rm(weather_uk.df)
Sys.sleep(2); gc(); Sys.sleep(2)


### Reduce to events with 3 days in a row
for(b in buf){
  v1 <- paste0("heatwave", b)
  
  heatwave.df$tmp1 <- ave(heatwave.df[, v1], heatwave.df$id,
                          FUN = function(x) dplyr::lag(x, default = 0))
  heatwave.df$tmp2 <- ave(heatwave.df[, v1], heatwave.df$id,
                          FUN = function(x) dplyr::lag(x, 2, default = 0))
  # # Also include leads for anticipated heatwave? omit for now!
  # heatwave.df$tmp3 <- ave(heatwave.df[, v1], heatwave.df$id,
  #                         FUN = function(x) dplyr::lead(x, default = 0))
  # heatwave.df$tmp4 <- ave(heatwave.df[, v1], heatwave.df$id,
  #                         FUN = function(x) dplyr::lead(x, 2, default = 0))
  
  
  oo <- which(heatwave.df[, v1] + heatwave.df$tmp1 + heatwave.df$tmp2 < 3)
  heatwave.df[oo, v1] <- 0
  heatwave.df$tmp1 <- NULL; heatwave.df$tmp2 <- NULL
} 

for(b in buf){
  v1 <- paste0("heatwave", b)
  print(table(heatwave.df[, v1]))
  print(table(heatwave.df[, v1])/N * 100)
  hist(heatwave.df$tasmax[heatwave.df[, v1] == 1], main = v1)
  summary(heatwave.df$tasmax[heatwave.df[, v1] == 1])
}


### Reduce to relevant dates
oo <- which(heatwave.df$heatwave1 == 1 | heatwave.df$heatwave2 == 1 | heatwave.df$heatwave3 == 1)
heatwave.df <- heatwave.df[oo, ]



### Save
save(heatwave.df, file = "Heatwave.RData")





# #################################
# ### Look at spatial variation ###
# #################################
# 
# table(heatwave.df$date[heatwave.df$year == 2011])
# 
# # Example day: 2011-10-01
# tmp.df <- weather_uk.df[which(weather_uk.df$year == 2011 & weather_uk.df$mon == 10 & weather_uk.df$day == 3), ]
# 
# # Merge to spatial data
# haduk_2011.sp <- merge(haduk.sp, tmp.df, by = "id")
# 
# par(mfrow = c(1, 2))
# plot(haduk_2011.sp[, "tasmax"], key.pos = NULL, reset = FALSE, border = NA)
# plot(haduk_2011.sp[, "heatwave"], key.pos = NULL, reset = FALSE, border = NA)
# par(mfrow = c(1, 1))



########################################
### Merge the heatwave data to lsoas ###
########################################

### Use the lsoa population averages
# England
centroids_en.spdf <- read_sf(dsn="../00_Input_data/Lower_Layer_Super_Output_Areas_December_2001_Population_Weighted_Centroids", 
                          layer="Lower_Layer_Super_Output_Areas_December_2001_Population_Weighted_Centroids")

names(centroids_en.spdf)[which(names(centroids_en.spdf) == "lsoa01cd")] <- "lsoa01"

# Scotland
centroids_sc.spdf <- read_sf(dsn="../00_Input_data/SG_DataZoneCent_2001", 
                             layer="SG_DataZone_Cent_2001")

names(centroids_sc.spdf)[which(names(centroids_sc.spdf) == "DZ_CODE")] <- "lsoa01"



# NorthernIreland
ni.spdf <- read_sf(dsn="../00_Input_data/soa", 
                   layer = "soa")
projold <- st_crs(ni.spdf)
centroids_ni.spdf <- st_centroid(ni.spdf)
centroids_ni.spdf <- st_transform(centroids_ni.spdf, 27700)
names(centroids_ni.spdf)[which(names(centroids_ni.spdf) == "SOA_CODE")] <- "lsoa01"


### Merge centroids
centroids_en.spdf <- centroids_en.spdf[, c("lsoa01")]
centroids_sc.spdf <- centroids_sc.spdf[, c("lsoa01")]
centroids_ni.spdf <- centroids_ni.spdf[, c("lsoa01")]
centroids.spdf <- rbind(centroids_en.spdf, centroids_sc.spdf, centroids_ni.spdf)



### Point in polygon
centroids.spdf <- st_transform(centroids.spdf, st_crs(haduk.sp))
ID <- st_join(centroids.spdf, haduk.sp, join = st_intersects)
ID <- st_drop_geometry(ID[, c("lsoa01", "id")])
ID$id <- as.numeric(ID$id)

### Join with heatwave data
heatwave_lsoa.df <- dplyr::inner_join(ID, heatwave.df, by = "id")
heatwave_lsoa.df <- data.frame(heatwave_lsoa.df)

# Save
save(heatwave_lsoa.df, file = "heatwave_lsoa.RData")






#####################
### Plot examples ###
#####################


### Gen UK shape
gb.sp <- read_sf(dsn="../00_Input_data/Countries__December_2016__Boundaries_GB_BUC-shp", 
                 layer="Countries__December_2016__Boundaries_GB_BUC")
ni.sp <- read_sf(dsn="../00_Input_data/OSNI_Open_Data_-_Largescale_Boundaries_-_NI_Outline-shp", 
                 layer="OSNI_Open_Data_-_Largescale_Boundaries_-_NI_Outline")
ni.sp <- st_transform(ni.sp, st_crs(gb.sp))
ni.sp <- st_simplify(ni.sp, dTolerance = 500)

gb.sp <- gb.sp[, "OBJECTID"]
ni.sp <- ni.sp[, "OBJECTID"]
uk.sp <- rbind(gb.sp, ni.sp)

uk.sp <- st_transform(uk.sp, st_crs(haduk.sp))


### Load UK LSOA shape
load("lsoa2001_sp.RData")
lsoa2001.spdf <- st_as_sf(lsoa2001.spdf)
lsoa2001.spdf <- st_transform(lsoa2001.spdf, st_crs(haduk.sp))
lsoa2001.spdf <- st_simplify(lsoa2001.spdf, dTolerance = 400)


### Select two dates
tmp.df <- heatwave.df[which(heatwave.df$year == 2009 & heatwave.df$mon == 07 & heatwave.df$day == 1), ]
tmp2.df <- heatwave.df[which(heatwave.df$year == 2013 & heatwave.df$mon == 08 & heatwave.df$day == 1), ]

tmp_lsoa.df <- heatwave_lsoa.df[which(heatwave_lsoa.df$year == 2009 & heatwave_lsoa.df$mon == 07 & heatwave_lsoa.df$day == 1), ]
tmp2_lsoa.df <- heatwave_lsoa.df[which(heatwave_lsoa.df$year == 2013 & heatwave_lsoa.df$mon == 08 & heatwave_lsoa.df$day == 1), ]



# Merge to spatial data
haduk_2011_1.sp <- merge(haduk.sp, tmp.df, by = "id")
haduk_2011_2.sp <- merge(haduk.sp, tmp2.df, by = "id")

lsoa_2011_1.sp <- merge(lsoa2001.spdf, tmp_lsoa.df, by = "lsoa01", all.x = TRUE)
lsoa_2011_2.sp <- merge(lsoa2001.spdf, tmp2_lsoa.df, by = "lsoa01", all.x = TRUE)
lsoa_2011_1.sp$heatwave1[is.na(lsoa_2011_1.sp$heatwave1)] <- 0
lsoa_2011_2.sp$heatwave1[is.na(lsoa_2011_2.sp$heatwave1)] <- 0


### Get weather data for days
load(paste0(weada, "HADUK_all.RData"))

weather.df <- weather_uk.df[which(weather_uk.df$year == 2009 & weather_uk.df$mon == 07 & weather_uk.df$day == 1), ]
weather2.df <- weather_uk.df[which(weather_uk.df$year == 2013 & weather_uk.df$mon == 08 & weather_uk.df$day == 1), ]
rm(weather_uk.df)
weather.df$id <- as.character(weather.df$id)
weather2.df$id <- as.character(weather2.df$id)
weather.spdf <- merge(haduk.sp, weather.df, by = "id")
weather2.spdf <- merge(haduk.sp, weather2.df, by = "id")




### Plot Example 1
mp1 <- tm_shape(weather.spdf) +
  tm_fill(col = "tasmax", style = "cont",
          #breaks = seq(10, 30, by = 5),
          title = "Degrees \nCelsius") +
  tm_shape(uk.sp) +
  tm_borders(col = "black", lwd = 2) +
  tm_layout(legend.frame = TRUE, legend.bg.color = TRUE,
            legend.position = c("left", "bottom"),
            title.position = c("left", "bottom"),
            fontfamily = "CM Roman",
            main.title = "Maximum Temperature 2009/07/01", main.title.position = "center",
            main.title.size = 1.6,
            legend.title.size = 1.4,
            legend.text.size = 1.2,
            title.snap.to.legend = TRUE) 

mp1  

mp2 <- tm_shape(uk.sp) +
  tm_borders(col = "black", lwd = 2) +
  tm_shape(lsoa_2011_1.sp) +
  tm_fill(col = "heatwave1", style = "cat",
          palette = c(ggplot2::alpha("#E8E8E8", 0.5), "#e36809") ,
          title = "Heatwave",
          labels = c("Not affected", "Affected")) +
  tm_shape(uk.sp) +
  tm_borders(col = "black", lwd = 2) +
  tm_layout(legend.frame = TRUE, legend.bg.color = TRUE,
            legend.position = c("left", "bottom"),
            fontfamily = "CM Roman",
            main.title = "Heatwave 2009/07/01", main.title.position = "center",
            main.title.size = 1.6,
            legend.title.size = 1.4,
            legend.text.size = 1.2,) 

mp2  



cairo_pdf(file = paste("../03_Output/Example_Heat1.pdf", sep=""), width = 9, height=8, bg = "white", family="CM Roman")
par(mar=c(0,0,3,0))
par(mfrow=c(1,1),oma=c(0,0,0,0))
tmap_arrange(mp1, mp2, ncol = 2, nrow = 1)
dev.off()

png(file = paste("../03_Output/Example_Heat1_2.png", sep=""), width = 4.5, height=8, bg = "white", family="CM Roman", units = "in", res = 300)
par(mar=c(0,0,3,0))
par(mfrow=c(1,1),oma=c(0,0,0,0))
mp1 + tm_layout(main.title.size = 1.4)
dev.off()




### Plot Example 2
mp1 <- tm_shape(weather2.spdf) +
  tm_fill(col = "tasmax", style = "cont",
          #breaks = seq(10, 30, by = 5),
          title = "Degrees \nCelsius") +
  tm_shape(uk.sp) +
  tm_borders(col = "black", lwd = 2) +
  tm_layout(legend.frame = TRUE, legend.bg.color = TRUE,
            legend.position = c("left", "bottom"),
            title.position = c("left", "bottom"),
            fontfamily = "CM Roman",
            main.title = "Maximum Temperature 2013/08/01", main.title.position = "center",
            main.title.size = 1.6,
            legend.title.size = 1.4,
            legend.text.size = 1.2,
            title.snap.to.legend = TRUE) 

mp1  

mp2 <- tm_shape(uk.sp) +
  tm_borders(col = "black", lwd = 2) +
  tm_shape(lsoa_2011_2.sp) +
  tm_fill(col = "heatwave1", style = "cat",
          palette = c(ggplot2::alpha("#E8E8E8", 0.5), "#e36809") ,
          title = "Heatwave",
          labels = c("Not affected", "Affected")) +
  tm_shape(uk.sp) +
  tm_borders(col = "black", lwd = 2) +
  tm_layout(legend.frame = TRUE, legend.bg.color = TRUE,
            legend.position = c("left", "bottom"),
            fontfamily = "CM Roman",
            main.title = "Heatwave 2013/08/01", main.title.position = "center",
            main.title.size = 1.6,
            legend.title.size = 1.4,
            legend.text.size = 1.2,) 

mp2  



cairo_pdf(file = paste("../03_Output/Example_Heat2.pdf", sep=""), width = 9, height=8, bg = "white", family="CM Roman")
par(mar=c(0,0,3,0))
par(mfrow=c(1,1),oma=c(0,0,0,0))
tmap_arrange(mp1, mp2, ncol = 2, nrow = 1)
dev.off()






