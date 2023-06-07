
#### Combine Defra Flood Reports and LSOA ####
#### Tobias Ruettenauer ####
#### 2020/ 01 / 31 ####

rm(list=ls())

### Load packages
library(foreign)
library(rgdal)
library(spdep)
library(rgeos)
library(doParallel)
library(foreign)
library(GISTools)
library(cleangeo)
library(fuzzyjoin)


### Working Directory
setwd("C:/work/Forschung/Climate Change_Replication/02_Data")


#################
### Load data ###
#################


load("Lsoa_flood_daily.RData")

bhps.df <- read.dta("all_bhpsukhls.dta")


####################
### Prepare data ###
####################


### Code Missing LSOA
bhps.df$lsoa01[bhps.df$lsoa01 == ""] <- NA


### For UKHLS use derived date
bhps.df$istrtdaty[bhps.df$wave > 18] <- bhps.df$intdaty_dv[bhps.df$wave > 18]
bhps.df$istrtdatm[bhps.df$wave > 18] <- bhps.df$intdatm_dv[bhps.df$wave > 18]
bhps.df$istrtdatd[bhps.df$wave > 18] <- bhps.df$intdatd_dv[bhps.df$wave > 18]


### Interview date (impute for first year)
bhps.df$istrtdaty[is.na(bhps.df$istrtdaty)] <- bhps.df$year[is.na(bhps.df$istrtdaty)]

### Generate interview date
bhps.df$date <- as.Date(paste(bhps.df$istrtdatd, bhps.df$istrtdatm, bhps.df$istrtdaty, sep = "/"),
                              "%d/%m/%Y")

summary(bhps.df$date)


# Still 24 NA dates

### Order columns
vars <- c("pidp", "year", "wave", "lsoa01", "date", "istrtdatd", "istrtdatm", "istrtdaty")
bhps.df <- bhps.df[, c(vars, names(bhps.df)[which(!names(bhps.df) %in% vars)])]




### Distribution of flood percentages
summary(lsoa_flood_day.df$flood_percentage1)
quantile(lsoa_flood_day.df$flood_percentage1, c(seq(0, .9, 0.1), .95, .99))
summary(lsoa_flood_day.df$flood_percentage2)
quantile(lsoa_flood_day.df$flood_percentage2,  c(seq(0, .9, 0.1), .95, .99))
summary(lsoa_flood_day.df$flood_percentage3)
quantile(lsoa_flood_day.df$flood_percentage3, c(seq(0, .9, 0.1), .95, .99))

summary(c(lsoa_flood_day.df$flood_percentage1, lsoa_flood_day.df$flood_percentage2, lsoa_flood_day.df$flood_percentage3))
quantile(c(lsoa_flood_day.df$flood_percentage1, lsoa_flood_day.df$flood_percentage2, lsoa_flood_day.df$flood_percentage3), c(seq(0, .9, 0.1), .95, .99))

par(mfrow = c(1, 3))
hist(lsoa_flood_day.df$flood_percentage3[lsoa_flood_day.df$flood_percentage3 > 0], main = "5km")
hist(lsoa_flood_day.df$flood_percentage2[lsoa_flood_day.df$flood_percentage2 > 0], main = "2km")
hist(lsoa_flood_day.df$flood_percentage1[lsoa_flood_day.df$flood_percentage1 > 0], main = "1km")
par(mfrow = c(1, 1))

summary(lsoa_flood_day.df$flood_percentage2[lsoa_flood_day.df$date >= "2007-01-01"])
quantile(lsoa_flood_day.df$flood_percentage2[lsoa_flood_day.df$date >= "2007-01-01"], c(seq(0, .9, 0.1), .95, .99))

par(mfrow = c(1, 3))
hist(lsoa_flood_day.df$flood_percentage3[lsoa_flood_day.df$flood_percentage3 > 0 & lsoa_flood_day.df$date >= "2007-01-01"], main = "5km")
hist(lsoa_flood_day.df$flood_percentage2[lsoa_flood_day.df$flood_percentage2 > 0 & lsoa_flood_day.df$date >= "2007-01-01"], main = "2km")
hist(lsoa_flood_day.df$flood_percentage1[lsoa_flood_day.df$flood_percentage1 > 0 & lsoa_flood_day.df$date >= "2007-01-01"], main = "1km")
par(mfrow = c(1, 1))


### Distribution of flood areas (total size)
b1 <- 1000
b2 <- 2000
b3 <- 5000
rm(pi)
s1 <- b1^2 * pi
s2 <- b2^2 * pi
s3 <- b3^2 * pi

lsoa_flood_day.df$flood_area3 <- lsoa_flood_day.df$flood_percentage3 * s3
lsoa_flood_day.df$flood_area2 <- lsoa_flood_day.df$flood_percentage2 * s2
lsoa_flood_day.df$flood_area1 <- lsoa_flood_day.df$flood_percentage1 * s1


# Histogram
par(mfrow = c(1, 3))
hist(lsoa_flood_day.df$flood_area3[lsoa_flood_day.df$flood_area3 > 0 & lsoa_flood_day.df$date >= "2007-01-01"], main = "5km")
hist(lsoa_flood_day.df$flood_area2[lsoa_flood_day.df$flood_area2 > 0 & lsoa_flood_day.df$date >= "2007-01-01"], main = "2km")
hist(lsoa_flood_day.df$flood_area1[lsoa_flood_day.df$flood_area1 > 0 & lsoa_flood_day.df$date >= "2007-01-01"], main = "1km")
par(mfrow = c(1, 1))

# Quantiles
summary(lsoa_flood_day.df$flood_area1[lsoa_flood_day.df$date >= "2007-01-01"])
quantile(lsoa_flood_day.df$flood_area1[lsoa_flood_day.df$date >= "2007-01-01"], c(seq(0, .9, 0.1), .95, .99))

summary(lsoa_flood_day.df$flood_area2[lsoa_flood_day.df$date >= "2007-01-01"])
quantile(lsoa_flood_day.df$flood_area2[lsoa_flood_day.df$date >= "2007-01-01"], c(seq(0, .9, 0.1), .95, .99))

summary(lsoa_flood_day.df$flood_area3[lsoa_flood_day.df$date >= "2007-01-01"])
quantile(lsoa_flood_day.df$flood_area3[lsoa_flood_day.df$date >= "2007-01-01"], c(seq(0, .9, 0.1), .95, .99))

# Use one hectar (10,000 sqm) and 50,000 sqm 


### Table per year
by(lsoa_flood_day.df$flood_area3 > 10000, substring(lsoa_flood_day.df$date, 1, 4), FUN = function(x) table(x))
by(lsoa_flood_day.df$flood_area2 > 10000, substring(lsoa_flood_day.df$date, 1, 4), FUN = function(x) table(x))
by(lsoa_flood_day.df$flood_area1 > 10000, substring(lsoa_flood_day.df$date, 1, 4), FUN = function(x) table(x))



#######################################
### Calculate whether in flood area ###
#######################################

# Define cutoffs
areaco <- c(0, 10000, 50000)
names(areaco) <- c("", "3", "10")
intens <- c("", "3", "10")
buf <- c("1", "2", "3")

### Ever flooded
for(b in buf){
  k <- 1
  for(i in intens){
    v1 <- paste0("flood", "_area", b)
    
    n1 <- paste0("flood", i, "_ever", b)
    
    
    floodids <- unique(lsoa_flood_day.df$lsoa01[lsoa_flood_day.df[, v1] > areaco[which(intens == i)]])
    bhps.df[, n1] <- 0 
    bhps.df[which(bhps.df$lsoa01 %in% floodids), n1]  <- 1
    bhps.df[which(is.na(bhps.df$lsoa01)), n1] <- NA 
  }
}




###########################################
### Fuzzy Merge events by lsoa and date ###
###########################################

### --- functions for fuzzy match --- ### 

### Minimal distance match function
mindistmatch_fun <- function(v1, v2, method = "euclidean", dir = c("both", "positive", "negative"),
                             max_dist = 1, distance_col = "distance") {
  if (is.null(dim(v1))) {
    # If the vectors are one-dimensional, turn them into 1-column matrices
    v1 <- t(t(as.numeric(v1)))
    v2 <- t(t(as.numeric(v2)))
  }
  
  l <- dim(v1)[1]
  dir <- match.arg(dir)
  pn <- sign(rowSums(v1 - v2))
  if(dir == "positive"){
    dp <- which(pn >= 0)
    v1 <- v1[dp, , drop = F]; v2 <- v2[dp, , drop = F]
  } else if(dir == "negative"){
    dn <- which(pn < 0)
    v1 <- v1[dn, , drop = F]; v2 <- v2[dn, , drop = F]
  }
  
  if (method == "euclidean") {
    d <- sqrt(rowSums((v1 - v2) ^ 2))
  } else if (method == "manhattan") {
    d <- rowSums(abs(v1 - v2))
  }
  
  if(dir == "positive" && length(dp) < l){
    ds <- d
    d <- rep(NA, l)
    d[dp] <- ds
  } else if(dir == "negative" && length(dn) < l){
    ds <- d
    d <- rep(NA, l)
    d[dn] <- ds
  }
  
  ret <- tibble::tibble(instance = !is.na(d) & d == min(d, na.rm = T) & d <= max_dist)
  if (!is.null(distance_col)) {
    ret[[distance_col]] <- pn * d
  }
  ret
}


### Fuzzy join 2, adopted from https://stackoverflow.com/questions/48008903/combined-fuzzy-and-exact-matching
fuzzy_inner_join2 <- function(x,y,by, match_fun, ...){
  match_fun_equal_lgl <- sapply(match_fun, identical, `==`)
  # columns to use for exact join equivalent
  by_exact = by[match_fun_equal_lgl]
  # columns to use for fuzzy join on relevant subsets of data (for efficiency)
  by_fuzzy = by[!match_fun_equal_lgl]
  # update match_fun
  match_fun <- match_fun[!match_fun_equal_lgl]
  # trim inputs of irrelevant data
  x <- dplyr::semi_join(x,y,by = by_exact)
  y <- dplyr::semi_join(y,x,by = by_exact)
  y <- dplyr::left_join(x, y, by = by_exact) # duplicate y for each x
  # make lists so we have pairs of data frames to fuzzy join together
  x_list <- dplyr::group_split(dplyr::group_by_at(x, c(by_exact, names(by_fuzzy))))
  y_list <- dplyr::group_split(dplyr::group_by_at(y, c(by_exact, names(by_fuzzy))), keep = FALSE)
  # apply fuzzy join on pairs and bind the results
  purrr::map2_dfr(x_list,y_list, fuzzyjoin::fuzzy_join, match_fun = match_fun,
                  by = by_fuzzy, mode = "inner", ...)
}



### --- Perform merging --- ###

### Reduce to relevant flood variables
lsoa_flood_day.tmp <- lsoa_flood_day.df[, c("lsoa01", "date", "flood_percentage1", "flood_percentage2", "flood_percentage3",
                                            "flood_duration1", "flood_duration2", "flood_duration3", "flood_area1", "flood_area2", "flood_area3")]


### Reduce to ids for inner match first, use unique lsoa date combinations
bhps_ids.df <- unique(bhps.df[, c("lsoa01", "date")])
bhps_flood.df <- bhps.df

# Define cutoffs
areaco <- c(0, 10000, 50000)
intens <- c("", "3", "10")
buf <- c("1", "2", "3")


### Fuzzy join past and future in loop
for(b in buf){
  for(i in intens){
    v1 <- paste0("flood_percentage", b)
    v2 <- paste0("flood_duration", b)
    v3 <- paste0("flood_area", b)
    
    p1 <- paste0("past_flood", i, "_percentage", b)
    p2 <- paste0("past_flood", i, "_duration", b)
    p3 <- paste0("past_flood", i, "_area", b)
    p4 <- paste0("past_flood", i, "_date", b)
    p5 <- paste0("past_flood", i, "_tempdist", b)
    
    f1 <- paste0("future_flood", i, "_percentage", b)
    f2 <- paste0("future_flood", i, "_duration", b)
    f3 <- paste0("future_flood", i, "_area", b)
    f4 <- paste0("future_flood", i, "_date", b)
    f5 <- paste0("future_flood", i, "_tempdist", b)
    
    
    
    # reduce to areas with respective flood size
    flood.tmp <- lsoa_flood_day.tmp[which(lsoa_flood_day.tmp[, v3] > areaco[which(intens == i)]), ] 
    flood_past.tmp <- flood.tmp[, c("lsoa01", "date", v1, v2, v3)]
    names(flood_past.tmp)[-c(1)] <- c(p4, p1, p2, p3)
    flood_future.tmp <- flood.tmp[, c("lsoa01", "date", v1, v2, v3)]
    names(flood_future.tmp)[-c(1)] <- c(f4, f1, f2, f3)
    
    # Join past (within 10 years)
    flood_past.df <- fuzzy_inner_join2(tibble::as_tibble(bhps_ids.df),
                                       tibble::as_tibble(flood_past.tmp),
                                       by = c( "lsoa01" = "lsoa01",
                                               "date" = p4),
                                       match_fun = list(`==`, function(x,y) mindistmatch_fun(x, y, max_dist = 365*10, 
                                                                                             dir = "positive",
                                                                                             distance_col = p5)))
    
    
    # Join future (within 10 years)
    flood_future.df <- fuzzy_inner_join2(tibble::as_tibble(bhps_ids.df),
                                         tibble::as_tibble(flood_future.tmp),
                                         by = c( "lsoa01" = "lsoa01",
                                                 "date" = f4),
                                         match_fun = list(`==`, function(x,y) mindistmatch_fun(x,y, max_dist = 365*10, 
                                                                                               dir = "negative",
                                                                                               distance_col = f5)))
    
    ### Left join to full data
    bhps_flood.df <- dplyr::left_join(bhps_flood.df, flood_past.df, by = c("lsoa01", "date"))
    bhps_flood.df <- dplyr::left_join(bhps_flood.df, flood_future.df, by = c("lsoa01", "date"))
    
    print(paste(b, i))
  }
}


### Save
save(bhps_flood.df, file = "bhpsukhls_flood.RData")



