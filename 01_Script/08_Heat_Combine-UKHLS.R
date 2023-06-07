
#### Combine MET office data and LSOA ####
#### Tobias Ruettenauer ####
#### 2020/ 08 / 13 ####

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
library(tibble)
library(fuzzyjoin)


### Working Directory
setwd("C:/work/Forschung/Climate Change_Replication/02_Data")


#################
### Load data ###
#################


load("heatwave_lsoa.RData")

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


### Average heatwave temperature
buf <- c(1, 2, 3)
intens <- ""
for(b in buf){
  k <- 1
  for(i in intens){
    v1 <- paste0("heatwave", b)
    
    id <- paste0("heat_id", b)
    n1 <- paste0("heat_tasmax", b)
    n2 <- paste0("heat_tasmax_diff", b)
    
    ### Gen common id for each heatwave (consecutive days)
    heatwave_lsoa.df <- heatwave_lsoa.df[order(heatwave_lsoa.df$lsoa01, heatwave_lsoa.df$date),]
    heatwave_lsoa.df$tmp <- heatwave_lsoa.df$date
    heatwave_lsoa.df$tmp[heatwave_lsoa.df[, v1 ] == 0] <- NA
    
    heatwave_lsoa.df[, id] <- ave(as.numeric(heatwave_lsoa.df$tmp),
                                    heatwave_lsoa.df$lsoa01,
                                    FUN = function(x) x - dplyr::lag(x, default = 0))
    heatwave_lsoa.df[which(heatwave_lsoa.df[, id] <= 1), id] <- 0 
    heatwave_lsoa.df[which(heatwave_lsoa.df[, id] > 1), id] <- 1 
    heatwave_lsoa.df[, id] <- ave(heatwave_lsoa.df[, id],
                                    heatwave_lsoa.df$lsoa01,
                                    FUN = function(x) cumsum(x))
    
    ### Gen average in tasmax and diff for each heatwave id
    oo <- which(!is.na(heatwave_lsoa.df[, id]))
    heatwave_lsoa.df[, n1] <- NA; heatwave_lsoa.df[, n2] <- NA
    
    heatwave_lsoa.df[oo, n1] <- ave(heatwave_lsoa.df$tasmax[oo],
                                        heatwave_lsoa.df$lsoa01[oo], heatwave_lsoa.df[oo, id],
                                        FUN = function(x) mean(x))
    heatwave_lsoa.df[oo, n2] <- ave(heatwave_lsoa.df$tasmax_diff[oo],
                                             heatwave_lsoa.df$lsoa01[oo], heatwave_lsoa.df[oo, id],
                                             FUN = function(x) mean(x))
    heatwave_lsoa.df$tmp <- NULL
    
  }
}





#######################################
### Calculate whether in heat area ###
#######################################

### Ever heatave affected
buf <- c(1, 2, 3)
intens <- ""
for(b in buf){
  k <- 1
  for(i in intens){
    v1 <- paste0("heatwave", b)
    
    n1 <- paste0("heatwave", "_ever", b)
    
    heatids <- unique(heatwave_lsoa.df$lsoa01[heatwave_lsoa.df[, v1] > 0])
    bhps.df[, n1] <- 0 
    bhps.df[which(bhps.df$lsoa01 %in% heatids), n1]  <- 1
    bhps.df[which(is.na(bhps.df$lsoa01)), n1] <- NA 
    
    
  }
}





#####################################
### Merge events by lsoa and date ###
#####################################

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
# Fuzzy join 2 from https://stackoverflow.com/questions/48008903/combined-fuzzy-and-exact-matching
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


### Fuzzy merge
buf <- c(1, 2, 3)
intens <- ""

lsoa_heat_day.tmp <- heatwave_lsoa.df


# Reduce to ids for inner match first, use unique lsoa date combinations
# XX ATTENTION: Only match after 2000 for now XX
bhps_ids.df <- unique(bhps.df[bhps.df$istrtdaty >= 1990, c("lsoa01", "date")])
bhps_heat.df <- bhps.df


for(b in buf){
  for(i in intens){
    v1 <- paste0("heatwave", b)
    v2 <- paste0("heat_tasmax", b)
    v3 <- paste0("heat_tasmax_diff", b)
    
    p1 <- paste0("past_heatwave", b)
    p2 <- paste0("past_heat", "_tasmax", b)
    p3 <- paste0("past_heat", "_tasmax_diff", b)
    p4 <- paste0("past_heat", "_date", b)
    p5 <- paste0("past_heat", "_tempdist", b)
    
    f1 <- paste0("future_heatwave", b)
    f2 <- paste0("future_heat", "_tasmax", b)
    f3 <- paste0("future_heat", "_tasmax_diff", b)
    f4 <- paste0("future_heat", "_date", b)
    f5 <- paste0("future_heat", "_tempdist", b)
    
    
    
    # reduce to areas with respective heat size
    heat.tmp <- lsoa_heat_day.tmp[which(lsoa_heat_day.tmp[, v1] == 1), ] 
    heat_past.tmp <- heat.tmp[, c("lsoa01", "date", v1, v2, v3)]
    names(heat_past.tmp)[-c(1)] <- c(p4, p1, p2, p3)
    heat_future.tmp <- heat.tmp[, c("lsoa01", "date", v1, v2, v3)]
    names(heat_future.tmp)[-c(1)] <- c(f4, f1, f2, f3)
    
    # Join past (within 3 years)
    heat_past.df <- fuzzy_inner_join2(tibble::as_tibble(bhps_ids.df),
                                       tibble::as_tibble(heat_past.tmp),
                                       by = c( "lsoa01" = "lsoa01",
                                               "date" = p4),
                                       match_fun = list(`==`, function(x,y) mindistmatch_fun(x, y, max_dist = 365*10, 
                                                                                             dir = "positive",
                                                                                             distance_col = p5)))
    
    
    # Join future (within 3 years)
    heat_future.df <- fuzzy_inner_join2(tibble::as_tibble(bhps_ids.df),
                                         tibble::as_tibble(heat_future.tmp),
                                         by = c( "lsoa01" = "lsoa01",
                                                 "date" = f4),
                                         match_fun = list(`==`, function(x,y) mindistmatch_fun(x,y, max_dist = 365*10, 
                                                                                               dir = "negative",
                                                                                               distance_col = f5)))
    
    ### Left join to full data
    bhps_heat.df <- dplyr::left_join(bhps_heat.df, heat_past.df, by = c("lsoa01", "date"))
    bhps_heat.df <- dplyr::left_join(bhps_heat.df, heat_future.df, by = c("lsoa01", "date"))
    
    print(paste(b, i))
  }
}


### Save
save(bhps_heat.df, file = "bhpsukhls_heat.RData")


for(b in buf){
  print(table(!is.na(bhps_heat.df[, paste0("past_heat_tasmax", b)])))
  print(table(!is.na(bhps_heat.df[, paste0("future_heat_tasmax", b)])))
}





### Reduce to relevant heatwave variables
lsoa_heat_day.tmp <- heatwave_lsoa.df[which(heatwave_lsoa.df[, "heatwave"] > 0),
                                      c("lsoa01", "date", "heat_tasmax", "heat_tasmax_diff")]


### Reduce to ids for inner match first, use unique lsoa date combinations
# Only match after 2000  
bhps_ids.df <- unique(bhps.df[bhps.df$istrtdaty >= 2000, c("lsoa01", "date")])

# Fuzzy join past
names(lsoa_heat_day.tmp)[2:ncol(lsoa_heat_day.tmp)] <- paste0("past_", names(lsoa_heat_day.tmp)[2:ncol(lsoa_heat_day.tmp)])
heat_past.df <- fuzzy_inner_join2(tibble::as_tibble(bhps_ids.df),
                          tibble::as_tibble(lsoa_heat_day.tmp),
                          by = c( "lsoa01" = "lsoa01",
                            "date" = "past_date"),
                          match_fun = list(`==`, function(x,y) mindistmatch_fun(x,y, max_dist = 365*3,
                                                                                dir = "positive",
                                                                                distance_col = "past_heat_tempdist")))

# Fuzzy join future
names(lsoa_heat_day.tmp) <- gsub("past", "future", names(lsoa_heat_day.tmp))
heat_future.df <- fuzzy_inner_join2(tibble::as_tibble(bhps_ids.df),
                                  tibble::as_tibble(lsoa_heat_day.tmp),
                                  by = c( "lsoa01" = "lsoa01",
                                          "date" = "future_date"),
                                  match_fun = list(`==`, function(x,y) mindistmatch_fun(x,y, max_dist = 365*3,
                                                                                        dir = "negative",
                                                                                        distance_col = "future_heat_tempdist")))

# Save
save(heat_past.df, heat_future.df, file = "Fuzzy_weather_id.RData")


### Left join to full data
bhps_heat.df <- dplyr::left_join(bhps.df, heat_past.df, by = c("lsoa01", "date"))
bhps_heat.df <- dplyr::left_join(bhps_heat.df, heat_future.df, by = c("lsoa01", "date"))

table(!is.na(bhps_heat.df$past_heat_tasmax))
table(!is.na(bhps_heat.df$future_heat_tasmax))



### Save
save(bhps_heat.df, file = "bhpsukhls_heat.RData")




