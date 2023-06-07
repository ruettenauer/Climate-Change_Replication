
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
library(future.apply)

library(OpenStreetMap)
library(tmap)
library(tmaptools)
library(extrafont)
loadfonts()


### Working Directory
setwd("C:/work/Forschung/Climate Change_Replication/02_Data")


#################
### Load data ###
#################


flood.spdf <- read_sf(dsn="../00_Input_data/EA_RecordedFloodOutlines_SHP_Full-2020-06-23/data", 
                      layer = "Recorded_Flood_Outlines")


### Start year 

sdate <- strsplit(as.character(flood.spdf$start_date), "-")
sdate <- do.call("rbind", sdate)

edate <- strsplit(as.character(flood.spdf$end_date), "-")
edate <- do.call("rbind", edate)
                  
flood.spdf$start_year <- as.numeric(sdate[,1])
flood.spdf$start_mont <- as.numeric(sdate[,2])
flood.spdf$start_day <- as.numeric(sdate[,3])

flood.spdf$end_year <- as.numeric(edate[,1])
flood.spdf$end_mont <- as.numeric(edate[,2])
flood.spdf$end_day <- as.numeric(edate[,3])

# Set date

flood.spdf$start_date <- as.Date(as.character(flood.spdf$start_date), "%Y-%m-%d")
flood.spdf$end_date <- as.Date(as.character(flood.spdf$end_date), "%Y-%m-%d")

### Set NA

flood.spdf$start_year[flood.spdf$start_year == 2050] <- NA 
flood.spdf$end_year[flood.spdf$end_year == 2050] <- NA 


### Restrict to obs from 1990

flood.spdf <- flood.spdf[which(flood.spdf$start_year >= 1990), ]
table(flood.spdf$start_year)


### Compute duration

flood.spdf$duration <- flood.spdf$end_date - flood.spdf$start_date
table(flood.spdf$duration )

# View(st_drop_geometry(flood.spdf[flood.spdf$duration >= 340, ]))


# XX What about unreastic dates / duration: eg. from 01.01.XX to 12.12.XX ? XX



### Clean / correct geometries

valid <- st_is_valid(flood.spdf)
flood.spdf <- st_make_valid(flood.spdf)


# ### Assign ids
flood.spdf$id <- c(1: nrow(flood.spdf))
# flood.spdf <- spChFIDs(flood.spdf, as.character(flood.spdf$id))

### Save
save(flood.spdf, file = "Floods_sp.RData")



##############################
### Clean Flood data dates ###
##############################

# For many observations: Date range very high
# Use date in name to approx. month
# Rule: Clean all with more than 340 days duration (patter xx.01.01 - xx.12.12)

flood.spdf$date_approx <- 0
flood.spdf$date_approx[flood.spdf$duration >= 340] <- 1


### Use event name to identify month
flood.spdf$new_month <- NA

m <- c("Jan",	"Feb",	"Mar",	"Apr",	"May",	"June",	"July",	"Aug",	"Sept",	"Oct",	"Nov",	"Dec")

k <- 1
for(i in m){
  cont <- which(grepl(i, flood.spdf$name[flood.spdf$date_approx == 1]))
  
  flood.spdf$new_month[flood.spdf$date_approx == 1][cont] <- k
  
  k <- k + 1
  
}

# Correct Easter 1998 flood (April)
flood.spdf$new_month[which(flood.spdf$name == "Easter 1998" & flood.spdf$date_approx ==1)] <- 4

# Autumn 200 approx on November
flood.spdf$new_month[which(flood.spdf$name == "06Autumn2000" & flood.spdf$date_approx ==1)] <- 11

# Event "15/06/2015"
flood.spdf$start_date[which(flood.spdf$name == "15/06/2015" & 
                              flood.spdf$date_approx ==1)] <- as.Date("15/06/2015", "%d/%m/%y" )
flood.spdf$end_date[which(flood.spdf$name == "15/06/2015" & 
                            flood.spdf$date_approx ==1)] <- as.Date("15/06/2015", "%d/%m/%y" )
flood.spdf$date_approx[which(flood.spdf$name == "15/06/2015" & 
                               flood.spdf$date_approx ==1)] <- 0 # to avoid overriding with NA


# Event 	"River Pont 12/01/16 - date wrack recorded"
flood.spdf$start_date[which(flood.spdf$name == "River Pont 12/01/16 - date wrack recorded" & 
                              flood.spdf$date_approx ==1)] <- as.Date("12/01/16", "%d/%m/%y" )
flood.spdf$end_date[which(flood.spdf$name == "River Pont 12/01/16 - date wrack recorded" & 
                            flood.spdf$date_approx ==1)] <- as.Date("12/01/16", "%d/%m/%y" )
flood.spdf$date_approx[which(flood.spdf$name == "River Pont 12/01/16 - date wrack recorded" & 
                               flood.spdf$date_approx ==1)] <- 0 # to avoid overriding with NA


### Replace date with approx month

# Start date (first of month)
flood.spdf$start_date[flood.spdf$date_approx == 1] <- as.Date(paste(flood.spdf$start_year[flood.spdf$date_approx == 1], 
                                                                    flood.spdf$new_month[flood.spdf$date_approx == 1], 
                                                                    "01", sep = "/"), "%Y/%m/%d" )

# End date (last of month)
m1 <- c(1, 3, 5, 7, 8, 10, 12)
flood.spdf$end_date[which(flood.spdf$new_month %in% m1)] <- flood.spdf$start_date[which(flood.spdf$new_month %in% m1)] + 30

m2 <- c(4, 6, 9, 11)
flood.spdf$end_date[which(flood.spdf$new_month %in% m2)] <- flood.spdf$start_date[which(flood.spdf$new_month %in% m2)] + 29

flood.spdf$end_date[which(flood.spdf$new_month == 2)] <- flood.spdf$start_date[which(flood.spdf$new_month == 2)]  + 27 


### Drop remaining instances (anyway only small percentage of lsoa flooded)
flood.spdf <- flood.spdf[-which(is.na(flood.spdf$start_date)), ]

### New duration
flood.spdf$new_duration <- flood.spdf$end_date - flood.spdf$start_date + 1 #(plus 1 for only one day floodings)

table(flood.spdf$new_duration)


###################################################################
### Drop duplicates (different event id but same date and area) ###
###################################################################

# spdupl <- gEquals(flood.spdf, byid = TRUE)
# 
# spdupl2 <- spdupl
# diag(spdupl2) <- FALSE
# 
# spdupl2 <- apply(spdupl2, 1, FUN = function(x) which(x == TRUE))
# spdupl2 <- spdupl2[lapply(spdupl2, length)>0]
# 
# save(spdupl2, file = "Flood_duplicates.RData")
# 
# load("Flood_duplicates.RData")
# 
# ### loop through to assign dupl id
# 
# flood.spdf$duplid <- 1:nrow(flood.spdf)
# 
# for(i in names(spdupl2)){
#   id <- as.numeric(i)
#   dd <- as.numeric(unlist(spdupl2[i]))
#   
#   duplid <- flood.spdf$duplid[id]
#   
#   flood.spdf$duplid[dd] <- duplid
# }


### Equal geoms
equal <- st_equals(flood.spdf)
# Coerce to string
equal <- lapply(equal, FUN = function(x) paste0(x, collapse = "_"))
# Add to df
flood.spdf$duplid <- unlist(equal)

# Delete duplicate later from overlap in long data (as some days overlap)




#################
### Load LSOA ###
#################

### Load centroids only
centroids.spdf <- read_sf(dsn="../00_Input_data/Lower_Layer_Super_Output_Areas_December_2001_Population_Weighted_Centroids", 
                      layer="Lower_Layer_Super_Output_Areas_December_2001_Population_Weighted_Centroids")

names(centroids.spdf)[which(names(centroids.spdf) == "lsoa01cd")] <- "lsoa01"


### LSOA spdf
lsoa2001.spdf <- read_sf(dsn="../00_Input_data/Lower_Layer_Super_Output_Areas__December_2001__Boundaries_EW_BGC-shp", 
                         layer="Lower_Layer_Super_Output_Areas__December_2001__Boundaries_EW_BGC")

names(lsoa2001.spdf)[which(names(lsoa2001.spdf) == "LSOA01CD")] <- "lsoa01"

# Transform
lsoa2001.spdf <- st_transform(lsoa2001.spdf, st_crs(centroids.spdf))
save(lsoa2001.spdf, file = "lsoa2001_sp.RData")


#######################################
### Create buffers around centroids ###
#######################################

# Buffer width in meters
b1 <- 1000
b2 <- 2000
b3 <- 5000


buff1.spdf <- st_buffer(centroids.spdf, units::set_units(b1, m))
buff2.spdf <- st_buffer(centroids.spdf, units::set_units(b2, m))
buff3.spdf <- st_buffer(centroids.spdf, units::set_units(b3, m))


### Plot example

ox.sp <- st_geometry(lsoa2001.spdf[which(grepl("Oxford ", as.character(lsoa2001.spdf$LSOA01NM))), ])

ids <- centroids.spdf[ox.sp, "lsoa01", drop = TRUE]
cent <- centroids.spdf[buff1.spdf$lsoa01 %in% ids, ]
cent$leg <- "Pop-weighted centroid"

buffox.spdf <- buff1.spdf[buff1.spdf$lsoa01 %in% ids, ]

ids2 <- st_intersection(buffox.spdf, st_geometry(flood.spdf[flood.spdf$start_year == 2007, ])[ox.sp])
ids2$area <- as.numeric(st_area(ids2))
ids2 <- aggregate(st_drop_geometry(ids2)[, c("area")],
                  by = list(lsoa01 = ids2$lsoa01), 
                  FUN = function(x) sum(x))
ids2 <- ids2[ids2$area >= 10000, ]

buffox.spdf <- buffox.spdf[buffox.spdf$lsoa01 %in% ids, ]
buffox.spdf$hit <- 0
buffox.spdf$hit[buffox.spdf$lsoa01 %in% ids2$lsoa01] <- 1
buffox.spdf$hit <- factor(buffox.spdf$hit)

# Get OSM data for background
osm_tmp <- read_osm(st_bbox(ox.sp), ext = 1, type = "stamen-toner", minNumTiles = 9) 
ox.sp <- st_transform(ox.sp, st_crs(osm_tmp))
buffox.spdf <- st_transform(buffox.spdf, st_crs(osm_tmp))
cent <- st_transform(cent, st_crs(osm_tmp))
floodox <- st_transform(flood.spdf[flood.spdf$start_year == 2007, ], st_crs(osm_tmp))[ox.sp, ]


# Plot
# mp <- tm_shape(st_as_sfc(st_bbox(ox.sp))) + tm_borders(col = alpha("black", 0)) + # add first layer to get correct scale bar
mp <-  tm_shape(osm_tmp) + tm_rgb() + 
  tm_shape(ox.sp) +
  tm_borders(col = "black") +
  tm_shape(floodox) +
  tm_fill(col = ggplot2::alpha("blue", 0.5)) +
  tm_shape(ox.sp) +
  tm_borders(col = "black") +
  tm_shape(cent) +
  tm_symbols(col = "#01665e", shape = 4, size = 1.3, border.lwd = 1.2) +
  tm_shape(buffox.spdf[buffox.spdf$hit == 0,]) +
  tm_borders(col = ggplot2::alpha("#999999", 0.8), lwd = 1.3) +
  tm_shape(buffox.spdf[buffox.spdf$hit == 1,]) +
  tm_borders(col = "#c51b7d", lwd = 2) +
  tm_add_legend("symbol", shape = c(4, 1, 1), size = 1.3,
                labels = c("Population-weighted centroid",
                           "1km buffer not affected",
                           "1km buffer affected"),
                col = c("#006d2c", ggplot2::alpha("#999999", 0.8), "#c51b7d")) +
  tm_layout(legend.frame = TRUE, legend.bg.color = TRUE,
            legend.position = c("left", "bottom"),
            fontfamily = "CM Roman",
            main.title = "Flood incidences Oxford 2007", main.title.position = "center",
            main.title.size = 1.6,
            legend.title.size = 1.2,
            legend.text.size = 1.2,) 
  # tm_scale_bar(position = c("RIGHT", "TOP"),
  #              width = 0.25, text.size = 1.2, )
            
mp  


cairo_pdf(file = paste("../03_Output/Example_Oxford.pdf", sep=""), width = 9, height=10, bg = "white", family="CM Roman")
par(mar=c(0,0,3,0))
par(mfrow=c(1,1),oma=c(0,0,0,0))
mp
dev.off()




##############################################
### Intersection between buffer and floods ###
##############################################

intersect_par <- function(n, s, x, y){
  n2 <- min(n+s-1, nrow(x))
  pol1 <- x[n:n2, ]
  pol2 <- y[pol1, ]
  res <- st_intersection(pol1, pol2)
  res
}

plan(multisession)

# Select steps
s <- 1000

# Intersection
int1 <- future_lapply(seq(1, nrow(buff1.spdf), s), FUN = function(n) intersect_par(n, s, buff1.spdf[, "lsoa01"], flood.spdf[, "id"]))
int1 <- do.call(rbind, int1)
int1$area1 <- as.numeric(st_area(int1))

# Intersection
int2 <- future_lapply(seq(1, nrow(buff2.spdf), s), FUN = function(n) intersect_par(n, s, buff2.spdf[, "lsoa01"], flood.spdf[, "id"]))
int2 <- do.call(rbind, int2)
int2$area2 <- as.numeric(st_area(int2))

# Intersection
int3 <- future_lapply(seq(1, nrow(buff3.spdf), s), FUN = function(n) intersect_par(n, s, buff3.spdf[, "lsoa01"], flood.spdf[, "id"]))
int3 <- do.call(rbind, int3)
int3$area3 <- as.numeric(st_area(int3))


### Create overlap matrix for all buffers (start with 3 as must include all 2 and 1 combination)
overlap.mat <- merge(st_drop_geometry(int3), st_drop_geometry(int2), by = c("lsoa01", "id"), all.x = TRUE)
overlap.mat <- merge(overlap.mat, st_drop_geometry(int1), by = c("lsoa01", "id"), all.x = TRUE)

# Replace NA by zero
overlap.mat$area2[which(is.na(overlap.mat$area2))] <- 0
overlap.mat$area1[which(is.na(overlap.mat$area1))] <- 0


# Save for reload
save(overlap.mat, file = "Overlap_buffer_flood2.RData")

# ### Reload if saved
# load("Overlap_buffer_flood2.RData")

plan(sequential)





###############################
### Add further information ###
###############################

lsoa_flood.df <- overlap.mat

### Sizes of buffers
rm(pi)
s1 <- b1^2 * pi
s2 <- b2^2 * pi
s3 <- b3^2 * pi

# Proportion overlap
lsoa_flood.df$per_flooded3 <- lsoa_flood.df$area3 / s3
lsoa_flood.df$per_flooded2 <- lsoa_flood.df$area2 / s2
lsoa_flood.df$per_flooded1 <- lsoa_flood.df$area1 / s1

summary(lsoa_flood.df[, c("per_flooded3", "per_flooded2", "per_flooded1")])



### Add flood meta data to each instance
lsoa_flood.df <- plyr::join(lsoa_flood.df, data.frame(flood.spdf), 
                            match = "all", by = "id", type = "left")


# Save
# save(lsoa_flood.df, file = "Lsoa_flood.RData")





######################################
### Generate long format over time ###
######################################

### Number of imprecise /corrected cases per year
table(flood.spdf$start_year)
table(flood.spdf$start_year[flood.spdf$date_approx == 1])


### gen long format for each day

# create id for each row
lsoa_flood.df$runid <- 1:nrow(lsoa_flood.df)

# Get each day between start and end
days <- lsoa_flood.df[, c("runid", "start_date", "end_date")]

gendays <- function(x){
  id <- x[1]
  sd <- x[2]
  ed <- x[3]
  
  ds <- as.Date(seq(as.Date(sd), as.Date(ed), by = "days"))
  
  res <- data.frame(cbind(runid = id, date = ds))
  
  return(res)
  
}

days <- apply(days, 1, FUN =  function(x) gendays(x))

# Reshape to df
days.df <- data.table::rbindlist(days)
days.df$runid <- as.numeric(as.character(days.df$runid))
days.df$date <- as.Date(as.numeric(as.character(days.df$date)), origin = "1970-01-01")


### Add further data

vars <- c("lsoa01", "id", "area3", "area2", "area1", "per_flooded3", "per_flooded2", "per_flooded1", 
          "new_duration", "date_approx", "runid", "duplid")

lsoa_flood_long.df <- plyr::join(days.df, lsoa_flood.df[, vars], 
                                  match = "all", by = "runid", type = "left")

### Remove duplicated geometries in floods (to avoid multiple counting of same instance)
id <- paste(lsoa_flood_long.df$lsoa01, lsoa_flood_long.df$duplid, lsoa_flood_long.df$date, sep = "_")

lsoa_flood_long.df$n <- ave(lsoa_flood_long.df$duplid, id,
                          FUN = function(x) c(1:length(x)))
table(lsoa_flood_long.df$n)

lsoa_flood_long.df <- lsoa_flood_long.df[which(lsoa_flood_long.df$n == 1), ]

### variable by lsoa and date

lsoa_flood_long.df$flood_area3 <- ave(lsoa_flood_long.df$per_flooded3, 
                                     paste(lsoa_flood_long.df$lsoa01, lsoa_flood_long.df$date, sep = "_"),
                                     FUN = function(x) sum(x))
lsoa_flood_long.df$flood_area2 <- ave(lsoa_flood_long.df$per_flooded2, 
                                      paste(lsoa_flood_long.df$lsoa01, lsoa_flood_long.df$date, sep = "_"),
                                      FUN = function(x) sum(x))
lsoa_flood_long.df$flood_area1 <- ave(lsoa_flood_long.df$per_flooded1, 
                                      paste(lsoa_flood_long.df$lsoa01, lsoa_flood_long.df$date, sep = "_"),
                                      FUN = function(x) sum(x))

summary(lsoa_flood_long.df[, c("flood_area3", "flood_area2", "flood_area1")])

test <- lsoa_flood_long.df[which(lsoa_flood_long.df$flood_area3 > 1),]

View(test[order(test$lsoa01, test$date),])


# !! For some dates overlapping events -> overlap percentage counted multiple times !! 



#######################################################
### Correct for spatial intersection between events ###
#######################################################

### Order by lsoa, date, id
lsoa_flood_long.df <- lsoa_flood_long.df[order(lsoa_flood_long.df$lsoa01, lsoa_flood_long.df$date, lsoa_flood_long.df$id), ]

### Count number of ids per lsoa and day
lsoa_flood_long.df$n <- ave(lsoa_flood_long.df$id, 
                                     paste(lsoa_flood_long.df$lsoa01, lsoa_flood_long.df$date, sep = "_"),
                                     FUN = function(x) length(x))

### For each lsoa and day, calculate the intersecting ids
lsoa_flood_long.df$combid <- ave(lsoa_flood_long.df$id, 
                            paste(lsoa_flood_long.df$lsoa01, lsoa_flood_long.df$date, sep = "_"),
                            FUN = function(x) paste(x, collapse = ""))

lsoa_flood_long.df$combid <- paste0(lsoa_flood_long.df$lsoa01, lsoa_flood_long.df$combid)

# unique id combinations with more than one intersection
nids <- by(lsoa_flood_long.df$id[lsoa_flood_long.df$n > 1],
           lsoa_flood_long.df$combid[lsoa_flood_long.df$n > 1],
           FUN = function(x) unique(x))

# intersects within combid?
plan(multisession)
intersecting <- future_lapply(nids, FUN = function(x) lapply(st_intersects(flood.spdf[which(flood.spdf$id %in% x), ]),
                                                      FUN = function(z) length(z)))

save(intersecting, file = "Intersecting_combids.RData")


# load("Intersecting_combids.RData")


intersecting <- lapply(intersecting, FUN = function(x) max(unlist(x)))

plan(sequential)




### Get those which intersect between each other
l <- which(intersecting > 1)

# New nids
nids <- names(intersecting)[l]


### computed combined overlap for each lsoa date bombo
lsoa_flood_long.df$area3_new <- NA
lsoa_flood_long.df$area2_new <- NA
lsoa_flood_long.df$area1_new <- NA


k <- 1
for(i in nids){
  lsoa_id <- unique(lsoa_flood_long.df$lsoa01[lsoa_flood_long.df$combid == i])
  flood_id <- unique(lsoa_flood_long.df$id[lsoa_flood_long.df$combid == i])

  buf_tmp1 <- buff1.spdf[buff1.spdf$lsoa01 %in% lsoa_id, ]
  buf_tmp2 <- buff2.spdf[buff2.spdf$lsoa01 %in% lsoa_id, ]
  buf_tmp3 <- buff3.spdf[buff3.spdf$lsoa01 %in% lsoa_id, ]

  flood_tmp <- flood.spdf[flood.spdf$id %in% flood_id, ]
  flood_tmp <- flood_tmp[buf_tmp3, ] # Reduce to intersecting
  flood_tmp <- st_union(flood_tmp) # Combine to one geometry

  int1 <- st_intersection(st_geometry(buf_tmp1), flood_tmp)
  overalp_area1 <- as.numeric(st_area(int1))
  int2 <- st_intersection(st_geometry(buf_tmp2), flood_tmp)
  overalp_area2 <- as.numeric(st_area(int2))
  int3 <- st_intersection(st_geometry(buf_tmp3), flood_tmp)
  overalp_area3 <- as.numeric(st_area(int3))
  if(length(overalp_area2) == 0 ){overalp_area2 <- 0}
  if(length(overalp_area1) == 0 ){overalp_area1 <- 0}

  lsoa_flood_long.df$area3_new[lsoa_flood_long.df$combid == i] <- overalp_area3
  lsoa_flood_long.df$area2_new[lsoa_flood_long.df$combid == i] <- overalp_area2
  lsoa_flood_long.df$area1_new[lsoa_flood_long.df$combid == i] <- overalp_area1

  # print(progress)
  if(k %% 100==0){
    cat(paste0(k, ","))
  }
  k <- k + 1

}


save(lsoa_flood_long.df, file = "lsoa_flood_long2.RData")

load("lsoa_flood_long2.RData")




#########################################
### Create variables per lsoa and day ###
#########################################

### percentage flooded using corrected data
lsoa_flood_long.df$flood_percentage1 <- NA
lsoa_flood_long.df$flood_percentage2 <- NA
lsoa_flood_long.df$flood_percentage3 <- NA

oo <- which(! lsoa_flood_long.df$combid %in% nids)
lsoa_flood_long.df$flood_percentage1[oo] <- lsoa_flood_long.df$flood_area1[oo]
lsoa_flood_long.df$flood_percentage2[oo] <- lsoa_flood_long.df$flood_area2[oo]
lsoa_flood_long.df$flood_percentage3[oo] <- lsoa_flood_long.df$flood_area3[oo]


### Use newly computed for all overalpping floods (attention: is already sum!)
rm(pi)
s1 <- b1^2 * pi
s2 <- b2^2 * pi
s3 <- b3^2 * pi

oo <- which(lsoa_flood_long.df$combid %in% nids)
lsoa_flood_long.df$flood_percentage1[oo] <- as.numeric(lsoa_flood_long.df$area1_new[oo]) / s1
lsoa_flood_long.df$flood_percentage2[oo] <- as.numeric(lsoa_flood_long.df$area2_new[oo]) / s2
lsoa_flood_long.df$flood_percentage3[oo] <- as.numeric(lsoa_flood_long.df$area3_new[oo]) / s3


summary(lsoa_flood_long.df[, c("flood_percentage1", "flood_percentage2", "flood_percentage3")])


# ### minimal distance to centre
# lsoa_flood_long.df$flood_distcent <- ave(lsoa_flood_long.df$dist, 
#                                      paste(lsoa_flood_long.df$lsoa01, lsoa_flood_long.df$date, sep = "_"),
#                                      FUN = function(x) min(x))


### duration (use average)
lsoa_flood_long.df$flood_duration <- ave(as.numeric(lsoa_flood_long.df$new_duration), 
                                         paste(lsoa_flood_long.df$lsoa01, lsoa_flood_long.df$date, sep = "_"),
                                         FUN = function(x) mean(x))

### indicator for duration approximate
lsoa_flood_long.df$flood_durapprox <- ave(lsoa_flood_long.df$date_approx, 
                                         paste(lsoa_flood_long.df$lsoa01, lsoa_flood_long.df$date, sep = "_"),
                                         FUN = function(x) max(x))




#####################################
### reduce to unique obs per lsoa ###
#####################################

lsoa_flood_day.df <- data.frame(lsoa_flood_long.df[, c("lsoa01", "date", "flood_percentage1", "flood_percentage2", "flood_percentage3", 
                                            "flood_duration", "flood_durapprox")])

lsoa_flood_day.df <- unique(lsoa_flood_day.df)

lsoa_flood_day.df <- lsoa_flood_day.df[order(lsoa_flood_day.df$lsoa01, lsoa_flood_day.df$date),]





#######################################################
### Compute the duration as cumulative days flooded ###
#######################################################

# Instead of the duration of the single flood incidences use the number of consecutive days affected by a flood

### Id for consecutive days
lsoa_flood_day.df <- lsoa_flood_day.df[order(lsoa_flood_day.df$lsoa01, lsoa_flood_day.df$date), ]

lsoa_flood_day.df$tmp_diff <- ave(as.numeric(lsoa_flood_day.df$date),
                                  lsoa_flood_day.df$lsoa01,
                                  FUN = function(x) x - dplyr::lag(x))
table(lsoa_flood_day.df$tmp_diff)

# Count as consecutive if <= 7 days difference
oo <- which(lsoa_flood_day.df$tmp_diff <= 7)
lsoa_flood_day.df$tmp_diff[oo] <- 0
lsoa_flood_day.df$tmp_diff[-oo] <- 1

# Create instance id
lsoa_flood_day.df$fid <- ave(as.numeric(lsoa_flood_day.df$tmp_diff),
                         lsoa_flood_day.df$lsoa01,
                         FUN = function(x) cumsum(x))


### Duration as cumsum per instance
lsoa_flood_day.df$tmp <- 1
lsoa_flood_day.df$flood_duration3 <- ave(lsoa_flood_day.df$tmp,
                                          paste(lsoa_flood_day.df$lsoa01, lsoa_flood_day.df$fid, sep = "_"),
                                          FUN = function(x) cumsum(x))

lsoa_flood_day.df$tmp[lsoa_flood_day.df$flood_percentage2 == 0] <- 0
lsoa_flood_day.df$flood_duration2 <- ave(lsoa_flood_day.df$tmp,
                                         paste(lsoa_flood_day.df$lsoa01, lsoa_flood_day.df$fid, sep = "_"),
                                         FUN = function(x) cumsum(x))
lsoa_flood_day.df$flood_duration2[lsoa_flood_day.df$flood_percentage2 == 0] <- NA

lsoa_flood_day.df$tmp[lsoa_flood_day.df$flood_percentage1 == 0] <- 0
lsoa_flood_day.df$flood_duration1 <- ave(lsoa_flood_day.df$tmp,
                                         paste(lsoa_flood_day.df$lsoa01, lsoa_flood_day.df$fid, sep = "_"),
                                         FUN = function(x) cumsum(x))
lsoa_flood_day.df$flood_duration1[lsoa_flood_day.df$flood_percentage1 == 0] <- NA


summary(lsoa_flood_day.df[, c("flood_duration1", "flood_duration2", "flood_duration3")])

lsoa_flood_day.df$tmp <- lsoa_flood_day.df$tmp_diff <- lsoa_flood_day.df$fid <- NULL

# Drop duration from pure flood incidence
lsoa_flood_day.df$flood_duration <- NULL


### Save
save(lsoa_flood_day.df, file = "Lsoa_flood_daily.RData")










