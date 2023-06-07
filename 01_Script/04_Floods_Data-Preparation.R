
#### Prepare variables ####
#### Tobias Ruettenauer ####
#### 2020/ 01 / 31 ####

rm(list=ls())

### Load packages
library(rgdal)
library(spdep)
library(rgeos)
library(doParallel)
library(foreign)
library(GISTools)
library(cleangeo)


### Working Directory
setwd("C:/work/Forschung/Climate Change_Replication/02_Data")


#################
### Load data ###
#################


load("bhpsukhls_flood.RData")


# Merge new variables (if necessay)
bhps.df <- read.dta("all_bhpsukhls.dta")
oo <- names(bhps.df)[which(!names(bhps.df) %in% names(bhps_flood.df))]
bhps_flood.df <- merge(bhps_flood.df, bhps.df[, c("pidp", "wave", oo)],
                       by = c("pidp", "wave"), all.x = TRUE)




####################
### Move in date ###
####################

### Repeat last function
repeat_last = function(x, forward = TRUE, maxgap = Inf, na.rm = FALSE) {
  if (!forward) x = rev(x)           # reverse x twice if carrying backward
  ind = which(!is.na(x))             # get positions of nonmissing values
  if (is.na(x[1]) && !na.rm)         # if it begins with NA
    ind = c(1,ind)                 # add first pos
  rep_times = diff(                  # diffing the indices + length yields how often
    c(ind, length(x) + 1) )          # they need to be repeated
  if (maxgap < Inf) {
    exceed = rep_times - 1 > maxgap  # exceeding maxgap
    if (any(exceed)) {               # any exceed?
      ind = sort(c(ind[exceed] + 1, ind))      # add NA in gaps
      rep_times = diff(c(ind, length(x) + 1) ) # diff again
    }
  }
  x = rep(x[ind], times = rep_times) # repeat the values at these indices
  if (!forward) x = rev(x)           # second reversion
  x
}


### Gen person year number
bhps_flood.df <- bhps_flood.df[order(bhps_flood.df$pidp, bhps_flood.df$wave), ]
bhps_flood.df$pynr <- ave(bhps_flood.df$wave,
                                by = bhps_flood.df$pidp,
                                FUN = function(x) 1:length(x))



### Clean plnowy4
bhps_flood.df$plnowy4[which(bhps_flood.df$plnowy4 > 2030)] <- NA

### Clean plnowm
bhps_flood.df$plnowm[which(bhps_flood.df$plnowm > 12)] <- 1


### Impute birth date for households never moved in their life
oo <- which(bhps_flood.df$mvever == 1 & is.na(bhps_flood.df$plnowy4))
bhps_flood.df$plnowy4[oo] <- bhps_flood.df$birthy[oo]
bhps_flood.df$plnowm[oo] <- bhps_flood.df$birthm[oo]

### Replace plnow with mvyr if NA
bhps_flood.df$plnowy4[which(is.na(bhps_flood.df$plnowy4))] <- bhps_flood.df$mvyr[which(is.na(bhps_flood.df$plnowy4))]
bhps_flood.df$plnowm[which(is.na(bhps_flood.df$plnowm))] <- bhps_flood.df$mvmnth[which(is.na(bhps_flood.df$plnowm))]


### Construct Move-in year

# Order date
bhps_flood.df <- bhps_flood.df[order(bhps_flood.df$pidp, bhps_flood.df$year),] 

# Impute last moving date
bhps_flood.df$movein_yr <- ave(bhps_flood.df$plnowy4,
                                       by = bhps_flood.df$pidp,
                                       FUN = function(x) repeat_last(x))


# ### For some obs the first move in year is missing (for instance children)
# # Impute the households pno1 movein value
# bhps_flood.df$tmp <- bhps_flood.df$movein_yr
# bhps_flood.df$tmp[bhps_flood.df$pno != 1] <- NA
# bhps_flood.df$tmp <- ave(bhps_flood.df$tmp,
#                          by = paste0(bhps_flood.df$hidp, "_", bhps_flood.df$wave),
#                          FUN = function(x) max(x, na.rm = TRUE))


# Children make a big part of those not having a move in date
# Could also use parents move in date or birthday
# However, for know use their age 16 survey as move in date.

# For those without moving year, use first observation year at location
bhps_flood.df$combid <- paste0(bhps_flood.df$pidp, "_", bhps_flood.df$lsoa01) 
bhps_flood.df$tmp <- ave(bhps_flood.df$istrtdaty,
                                     by = bhps_flood.df$combid,
                                     FUN = function(x) min(x))

bhps_flood.df$movein_yr[is.na(bhps_flood.df$movein_yr)] <- bhps_flood.df$tmp[is.na(bhps_flood.df$movein_yr)]
bhps_flood.df$tmp <- NULL




### Contruct Move-in month

# Order date
bhps_flood.df <- bhps_flood.df[order(bhps_flood.df$pidp, bhps_flood.df$year),] 

# Impute last moving date
bhps_flood.df$movein_mn <- ave(bhps_flood.df$plnowm,
                                     by = bhps_flood.df$pidp,
                                     FUN = function(x) repeat_last(x))

# For those without moving year, use first observation year at location
bhps_flood.df$combid <- paste0(bhps_flood.df$pidp, bhps_flood.df$lsoa01) 
bhps_flood.df$tmp <- ave(bhps_flood.df$istrtdatm,
                               by = bhps_flood.df$combid,
                               FUN = function(x) min(x))

# If there before 1991 and month is na, use random month
set.seed(1241234)
oo <- which(is.na(bhps_flood.df$movein_mn) & bhps_flood.df$movein_yr < 1991)
bhps_flood.df$movein_mn[oo] <- ave(bhps_flood.df$movein_mn[oo],
                                         by = bhps_flood.df$pidp[oo],
                                         FUN = function(x) rep(sample(c(1:12), 1), length(x)))
  
# Replace remaining with int month
bhps_flood.df$movein_mn[is.na(bhps_flood.df$movein_mn)] <- bhps_flood.df$tmp[is.na(bhps_flood.df$movein_mn)]
bhps_flood.df$tmp <- NULL


### Contruct Move-in Date

bhps_flood.df$movein_date <- as.Date(paste("15", bhps_flood.df$movein_mn, bhps_flood.df$movein_yr, sep = "/"),
                                           "%d/%m/%Y")

summary(bhps_flood.df$movein_date)






############################
### Resctrict the sample ###
############################



### Restrict to England only
bhps_flood.df <- bhps_flood.df[which(bhps_flood.df$gor_dv < 10), ]

### Drop BHPS UKHPS without lsoa or date
bhps_flood.df <- bhps_flood.df[which(!is.na(bhps_flood.df$lsoa01) & !is.na(bhps_flood.df$date)), ]




############################################################
### Generate period (move in and move out) for each lsoa ###
############################################################

bhps_flood.df <- bhps_flood.df[order(bhps_flood.df$pidp, bhps_flood.df$wave), ]

### Indicator for change in lsoa
bhps_flood.df$tmp <- ave(bhps_flood.df$lsoa01,
                         by = bhps_flood.df$pidp,
                         FUN = function(x) dplyr::lag(x))
  
bhps_flood.df$lsoa_change <- 0
oo <- which(bhps_flood.df$lsoa01 != bhps_flood.df$tmp & !is.na(bhps_flood.df$tmp))
bhps_flood.df$lsoa_change[oo] <- 1
bhps_flood.df$tmp <- NULL

### Cumulative number of lsoa spell per person
bhps_flood.df$lsoa_nr <- ave(bhps_flood.df$lsoa_change,
                             by = bhps_flood.df$pidp,
                             FUN = function(x) cumsum(x))

### Move in date for each lsoa spell (use first non-missing)
bhps_flood.df$movein_lsoa_date <- as.Date(ave(as.numeric(bhps_flood.df$movein_date),
                                              by = paste0(bhps_flood.df$pidp, "_", bhps_flood.df$lsoa_nr),
                                              FUN = function(x) x[!is.na(x)][1]), origin = "1970-01-01")




### Clean move in date if lsoa change, but no new move-in date (seems mostly people moving often)
bhps_flood.df$tmp_diff <- ave(as.numeric(bhps_flood.df$movein_lsoa_date),
                              by = bhps_flood.df$pidp,
                              FUN = function(x) x - dplyr::lag(x))

bhps_flood.df$tmp_diff2 <- ave(as.numeric(bhps_flood.df$date),
                               by = bhps_flood.df$pidp,
                               FUN = function(x) x - dplyr::lag(x))

# ids <- bhps_flood.df$pidp[which(bhps_flood.df$tmp_diff <= 65 & bhps_flood.df$lsoa_change == 1)]
# View(bhps_flood.df[bhps_flood.df$pidp %in% ids, c(1:5, 9:11, 36:37, 40:45, 264, 270, 341:ncol(bhps_flood.df))])

# Impute for those having no or negative change in move in
bhps_flood.df$impute_mvin <- 0
oo <- which(bhps_flood.df$tmp_diff <= 0 & bhps_flood.df$lsoa_change == 1)
bhps_flood.df$impute_mvin[oo] <- 1

# Rule: impute random date between the two interviews --------------
set.seed(54654887)
bhps_flood.df$tmp_randdiff <- NA
oo <- which(bhps_flood.df$impute_mvin == 1)
bhps_flood.df$tmp_randdiff[oo] <- as.numeric(lapply(bhps_flood.df$tmp_diff2[oo], 
                                              FUN = function(x) sample(0:x, size = 1)))

oo <- which(bhps_flood.df$impute_mvin == 1)
bhps_flood.df$movein_lsoa_date[oo] <- as.Date(as.numeric(bhps_flood.df$date[oo]) - bhps_flood.df$tmp_randdiff[oo],
                                              origin = "1970-01-01")

bhps_flood.df$tmp_diff <- NULL
bhps_flood.df$tmp_diff2 <- NULL
bhps_flood.df$tmp_randdiff <- NULL



# Replace Move in date for each lsoa spell by the more recent date (fill imputation for following years)
bhps_flood.df$movein_lsoa_date <- as.Date(ave(as.numeric(bhps_flood.df$movein_lsoa_date),
                                              by = paste0(bhps_flood.df$pidp, "_", bhps_flood.df$lsoa_nr),
                                              FUN = function(x) max(x)), origin = "1970-01-01")




### Gen moveout date for each lsoa

# Lead move in date if lsoa change
bhps_flood.df$tmp <- bhps_flood.df$movein_date
bhps_flood.df$tmp[which(bhps_flood.df$lsoa_change == 0)] <- NA

bhps_flood.df$tmp <- as.Date(ave(as.numeric(bhps_flood.df$tmp),
                                 by = bhps_flood.df$pidp,
                                 FUN = function(x) dplyr::lead(x)), origin = "1970-01-01")

# Use last interview data for the final wave of each person
bhps_flood.df$pynr <- ave(bhps_flood.df$pidp,
                          by = bhps_flood.df$pidp,
                          FUN = function(x) 1:length(x))
bhps_flood.df$pN <- ave(bhps_flood.df$pidp,
                        by = bhps_flood.df$pidp,
                        FUN = function(x) length(x))

oo <- which(bhps_flood.df$pN == bhps_flood.df$pynr)
bhps_flood.df$tmp[oo] <- bhps_flood.df$date[oo] 

# For each lsoa spell, use max move out date
bhps_flood.df$moveout_lsoa_date <- as.Date(ave(as.numeric(bhps_flood.df$tmp),
                                              by = paste0(bhps_flood.df$pidp, "_", bhps_flood.df$lsoa_nr),
                                              FUN = function(x) max(x, na.rm = TRUE)), 
                                           origin = "1970-01-01")
bhps_flood.df$tmp <- NULL





# #### Flag inconsistent observations ####
# 
# bhps_flood.df$tmp <- ave(as.numeric(bhps_flood.df$movein_date),
#                          by = bhps_flood.df$pidp,
#                          FUN = function(x) x - dplyr::lag(x))
# bhps_flood.df$movein_date_flag <- 0
# bhps_flood.df$movein_date_flag[bhps_flood.df$tmp < 0] <- 1
# bhps_flood.df$tmp <- NULL
# 
# table(bhps_flood.df$movein_date_flag)
# 
# ids <- bhps_flood.df$pidp[bhps_flood.df$movein_date_flag == 1]
# 
# View(bhps_flood.df[which(bhps_flood.df$pidp %in% ids), c("lsoa01", "pidp", "year", "wave", "date", 
#                                                          "plnowy4", "plnowm", "movein_date", "movein_date_flag", "movest")])
# 
# 
# ### By each combid, replace with first stated date if negative difference
# bhps_flood.df$movein_date_flag <- ave(as.numeric(bhps_flood.df$movein_date_flag),
#                                       by = bhps_flood.df$combid,
#                                       FUN = function(x) cumsum(x))





#########################
### Affected by flood ###
#########################

bhps_flood.df <- bhps_flood.df[order(bhps_flood.df$pidp, bhps_flood.df$wave), ]



### Floods date format
bhps_flood.df$past_flood_date1 <- as.Date(bhps_flood.df$past_flood_date1, origin = "1970-01-01")
bhps_flood.df$future_flood_date1 <- as.Date(bhps_flood.df$future_flood_date1, origin = "1970-01-01")
bhps_flood.df$past_flood10_date1 <- as.Date(bhps_flood.df$past_flood10_date1, origin = "1970-01-01")
bhps_flood.df$future_flood10_date1 <- as.Date(bhps_flood.df$future_flood10_date1, origin = "1970-01-01")
bhps_flood.df$past_flood3_date1 <- as.Date(bhps_flood.df$past_flood3_date1, origin = "1970-01-01")
bhps_flood.df$future_flood3_date1 <- as.Date(bhps_flood.df$future_flood3_date1, origin = "1970-01-01")

bhps_flood.df$past_flood_date2 <- as.Date(bhps_flood.df$past_flood_date2, origin = "1970-01-01")
bhps_flood.df$future_flood_date2 <- as.Date(bhps_flood.df$future_flood_date2, origin = "1970-01-01")
bhps_flood.df$past_flood10_date2 <- as.Date(bhps_flood.df$past_flood10_date2, origin = "1970-01-01")
bhps_flood.df$future_flood10_date2 <- as.Date(bhps_flood.df$future_flood10_date2, origin = "1970-01-01")
bhps_flood.df$past_flood3_date2 <- as.Date(bhps_flood.df$past_flood3_date2, origin = "1970-01-01")
bhps_flood.df$future_flood3_date2 <- as.Date(bhps_flood.df$future_flood3_date2, origin = "1970-01-01")

bhps_flood.df$past_flood_date3 <- as.Date(bhps_flood.df$past_flood_date3, origin = "1970-01-01")
bhps_flood.df$future_flood_date3 <- as.Date(bhps_flood.df$future_flood_date3, origin = "1970-01-01")
bhps_flood.df$past_flood10_date3 <- as.Date(bhps_flood.df$past_flood10_date3, origin = "1970-01-01")
bhps_flood.df$future_flood10_date3 <- as.Date(bhps_flood.df$future_flood10_date3, origin = "1970-01-01")
bhps_flood.df$past_flood3_date3 <- as.Date(bhps_flood.df$past_flood3_date3, origin = "1970-01-01")
bhps_flood.df$future_flood3_date3 <- as.Date(bhps_flood.df$future_flood3_date3, origin = "1970-01-01")



### Calculate whether hit at place of residence ###

# Future date
bhps_flood.df$date_lead <- as.Date(ave(as.numeric(bhps_flood.df$date),
                                       by = bhps_flood.df$pidp,
                                       FUN = function(x) dplyr::lead(x)), 
                                   origin = "1970-01-01")



### Past Floods ###
buf <- c("1", "2", "3") # For each buffer
intens <- c("", "3", "10") # and each intensity (percentage)

for(b in buf){
  for(i in intens){
    newname <- paste0("flood", i, "_affect_past", b)
    pastdate <- paste0("past_flood", i, "_date", b)
    futuredate <- paste0("future_flood", i, "_date", b)
    
    newnamedate <- paste0("flood", i, "_affect_past_date", b)
    newnamedur <- paste0("flood", i, "_affect_past_duration", b)
    newnameper <- paste0("flood", i, "_affect_past_percentage", b)
    
    pastdur <- paste0("past_flood", i, "_duration", b)
    futuredur <- paste0("future_flood", i, "_duration", b)
    pastper <- paste0("past_flood", i, "_percentage", b)
    futureper <- paste0("future_flood", i, "_percentage", b)
    
    bhps_flood.df[, newname] <- 0
    
    # Directly affected affected
    oo <- which(bhps_flood.df[, pastdate] >= bhps_flood.df$movein_lsoa_date)
    bhps_flood.df[, newname][oo] <- 1
    
    # Also save date
    bhps_flood.df$tmp_date <- NA
    bhps_flood.df$tmp_date[oo] <- bhps_flood.df[oo, pastdate]
    
    # Also save duration
    bhps_flood.df$tmp_dur <- NA
    bhps_flood.df$tmp_dur[oo] <- bhps_flood.df[oo, pastdur]
    
    # Also save percentage
    bhps_flood.df$tmp_per <- NA
    bhps_flood.df$tmp_per[oo] <- bhps_flood.df[oo, pastper]
    
    # Affected at old lsoa but moved afterwards
    bhps_flood.df$tmp <- NA 
    oo <- which(bhps_flood.df[, futuredate] > bhps_flood.df$date & 
                  bhps_flood.df[, futuredate] <= bhps_flood.df$moveout_lsoa_date & 
                  bhps_flood.df[, futuredate] < bhps_flood.df$date_lead)
    bhps_flood.df$tmp[oo] <- 1
    
    # Lag and declare as affected
    bhps_flood.df$tmp <- ave(bhps_flood.df$tmp,
                             by = bhps_flood.df$pidp,
                             FUN = function(x) dplyr::lag(x))
    
    bhps_flood.df[, newname][bhps_flood.df$tmp == 1] <- 1
    bhps_flood.df$tmp <- NULL
    
    # Also save date
    bhps_flood.df$tmp_date2 <- NA
    bhps_flood.df$tmp_date2[oo] <- bhps_flood.df[oo, futuredate]
    
    # Lag date 
    bhps_flood.df$tmp_date2 <- ave(bhps_flood.df$tmp_date2,
                                  by = bhps_flood.df$pidp,
                                  FUN = function(x) dplyr::lag(x))
    
    # Also save dur
    bhps_flood.df$tmp_dur2 <- NA
    bhps_flood.df$tmp_dur2[oo] <- bhps_flood.df[oo, futuredur]
    
    # Lag dur 
    bhps_flood.df$tmp_dur2 <- ave(bhps_flood.df$tmp_dur2,
                                   by = bhps_flood.df$pidp,
                                   FUN = function(x) dplyr::lag(x))
    
    # Also save per
    bhps_flood.df$tmp_per2 <- NA
    bhps_flood.df$tmp_per2[oo] <- bhps_flood.df[oo, futureper]
    
    # Lag per 
    bhps_flood.df$tmp_per2 <- ave(bhps_flood.df$tmp_per2,
                                  by = bhps_flood.df$pidp,
                                  FUN = function(x) dplyr::lag(x))
    
    # Use more recent date and duration
    oo <- which(bhps_flood.df$tmp_date2 > bhps_flood.df$tmp_date | is.na(bhps_flood.df$tmp_date))
    bhps_flood.df$tmp_date[oo] <- bhps_flood.df$tmp_date2[oo]
    bhps_flood.df$tmp_dur[oo] <- bhps_flood.df$tmp_dur2[oo]
    bhps_flood.df$tmp_per[oo] <- bhps_flood.df$tmp_per2[oo]
    
    bhps_flood.df$tmp_dur2 <- bhps_flood.df$tmp_date2 <- bhps_flood.df$tmp_per2 <- NA
    
    
    # Fill na
    bhps_flood.df[, newnamedate] <- as.Date(ave(bhps_flood.df$tmp_date,
                                        by = bhps_flood.df$pidp,
                                        FUN = function(x) repeat_last(x, forward = TRUE)),
                                        origin = "1970-01-01")
    bhps_flood.df[, newnamedur] <- ave(bhps_flood.df$tmp_dur,
                                        by = bhps_flood.df$pidp,
                                        FUN = function(x) repeat_last(x, forward = TRUE))
    bhps_flood.df[, newnameper] <- ave(bhps_flood.df$tmp_per,
                                       by = bhps_flood.df$pidp,
                                       FUN = function(x) repeat_last(x, forward = TRUE))
    
    bhps_flood.df$tmp_date <- bhps_flood.df$tmp_dur <- bhps_flood.df$tmp_per <- NA
    
    # Fill NAs after occasion
    bhps_flood.df[, newname] <- ave(bhps_flood.df[, newname],
                                    by = bhps_flood.df$pidp,
                                    FUN = function(x) cummax(x))
    
    
  }
  
}

# View(bhps_flood.df[, c("pidp", "year", "wave", "lsoa01", "date", "movein_lsoa_date", "moveout_lsoa_date", pastdate, futuredate, newname, newnamedate, newnamedur)])


### Future Floods ###
buf <- c("1", "2", "3") # For each buffer
intens <- c("", "3", "10") # and each intensity (percentage)

for(b in buf){
  for(i in intens){
    newname <- paste0("flood", i, "_affect_future", b)
    pastdate <- paste0("past_flood", i, "_date", b)
    futuredate <- paste0("future_flood", i, "_date", b)
    
    newnamedate <- paste0("flood", i, "_affect_future_date", b)
    newnamedur <- paste0("flood", i, "_affect_future_duration", b)
    newnameper <- paste0("flood", i, "_affect_future_percentage", b)
    
    pastdur <- paste0("past_flood", i, "_duration", b)
    futuredur <- paste0("future_flood", i, "_duration", b)
    pastper <- paste0("past_flood", i, "_percentage", b)
    futureper <- paste0("future_flood", i, "_percentage", b)
    
    pastaffect <- paste0("flood", i, "_affect_past", b)
    
    bhps_flood.df[, newname] <- 0
    
    # Directly affected affected
    oo <- which(bhps_flood.df[, futuredate] <= bhps_flood.df$moveout_lsoa_date)
    bhps_flood.df[, newname][oo] <- 1
    
    # Also save date
    bhps_flood.df$tmp_date <- NA
    bhps_flood.df$tmp_date[oo] <- bhps_flood.df[oo, futuredate]
    
    # Also save duration
    bhps_flood.df$tmp_dur <- NA
    bhps_flood.df$tmp_dur[oo] <- bhps_flood.df[oo, futuredur]
    
    # Also save percentage
    bhps_flood.df$tmp_per <- NA
    bhps_flood.df$tmp_per[oo] <- bhps_flood.df[oo, futureper]
    
    # Newly past affected
    bhps_flood.df$tmp <- as.Date(ave(as.numeric(bhps_flood.df$date),
                                     by = paste0(bhps_flood.df$pidp),
                                     FUN = function(x) dplyr::lag(x)), 
                                 origin = "1970-01-01")
    
    bhps_flood.df$tmp2 <- NA
    oo <- which(bhps_flood.df[, pastdate] <= bhps_flood.df$date & 
                  bhps_flood.df[, pastdate] > bhps_flood.df$tmp &
                  bhps_flood.df[, pastdate] >= bhps_flood.df$movein_lsoa_date)
    bhps_flood.df$tmp2[oo] <- 1
    
    
    # Lead and declare as future affected
    bhps_flood.df$tmp2 <- ave(bhps_flood.df$tmp2,
                              by = bhps_flood.df$pidp,
                              FUN = function(x) dplyr::lead(x))
    
    bhps_flood.df[, newname][bhps_flood.df$tmp2 == 1] <- 1
    bhps_flood.df$tmp <- NULL
    bhps_flood.df$tmp2 <- NULL
    
    # Also save date
    bhps_flood.df$tmp_date2 <- NA
    bhps_flood.df$tmp_date2[oo] <- bhps_flood.df[oo, pastdate]
    
    # Lead date 
    bhps_flood.df$tmp_date2 <- ave(bhps_flood.df$tmp_date2,
                                   by = bhps_flood.df$pidp,
                                   FUN = function(x) dplyr::lead(x))
    
    # Also save dur
    bhps_flood.df$tmp_dur2 <- NA
    bhps_flood.df$tmp_dur2[oo] <- bhps_flood.df[oo, pastdur]
    
    # Lead dur 
    bhps_flood.df$tmp_dur2 <- ave(bhps_flood.df$tmp_dur2,
                                  by = bhps_flood.df$pidp,
                                  FUN = function(x) dplyr::lead(x))
    
    # Also save per
    bhps_flood.df$tmp_per2 <- NA
    bhps_flood.df$tmp_per2[oo] <- bhps_flood.df[oo, pastper]
    
    # Lead per 
    bhps_flood.df$tmp_per2 <- ave(bhps_flood.df$tmp_per2,
                                  by = bhps_flood.df$pidp,
                                  FUN = function(x) dplyr::lead(x))
    
    # Use more earlier date and duration
    oo <- which(bhps_flood.df$tmp_date2 < bhps_flood.df$tmp_date | is.na(bhps_flood.df$tmp_date))
    bhps_flood.df$tmp_date[oo] <- bhps_flood.df$tmp_date2[oo]
    bhps_flood.df$tmp_dur[oo] <- bhps_flood.df$tmp_dur2[oo]
    bhps_flood.df$tmp_per[oo] <- bhps_flood.df$tmp_per2[oo]
    
    bhps_flood.df$tmp_dur2 <- bhps_flood.df$tmp_date2  <- bhps_flood.df$tmp_per2 <- NA
    
    # Fill na
    bhps_flood.df[, newnamedate] <- as.Date(ave(bhps_flood.df$tmp_date,
                                                by = bhps_flood.df$pidp,
                                                FUN = function(x) repeat_last(x, forward = FALSE)),
                                            origin = "1970-01-01")
    bhps_flood.df[, newnamedur] <- ave(bhps_flood.df$tmp_dur,
                                       by = bhps_flood.df$pidp,
                                       FUN = function(x) repeat_last(x, forward = FALSE))
    bhps_flood.df[, newnameper] <- ave(bhps_flood.df$tmp_per,
                                       by = bhps_flood.df$pidp,
                                       FUN = function(x) repeat_last(x, forward = FALSE))
    
    bhps_flood.df$tmp_date <- bhps_flood.df$tmp_dur <- bhps_flood.df$tmp_per <- NA
    
    
    # Fill NAs before occasion
    bhps_flood.df[, newname] <- ave(bhps_flood.df[, newname],
                                    by = bhps_flood.df$pidp,
                                    FUN = function(x) rev(cummax(rev(x))))
    
    
    # Correct those occasions with an overlap in migration
    tmp <- ave(bhps_flood.df[, pastaffect],
               by = bhps_flood.df$pidp,
               FUN = function(x) dplyr::lead(x))
    oo <- which(tmp == 1 & bhps_flood.df[, newname] == 0 & bhps_flood.df[, pastaffect] == 0)
    bhps_flood.df[, newname][oo] <- 1
    bhps_flood.df[, newname] <- ave(bhps_flood.df[, newname],
                                    by = bhps_flood.df$pidp,
                                    FUN = function(x) rev(cummax(rev(x))))
    
    # Also correct date and duration and percentage
    bhps_flood.df$tmp_date <- NA
    bhps_flood.df$tmp_date[oo] <- bhps_flood.df[oo + 1, pastdate]
    bhps_flood.df[!is.na(bhps_flood.df$tmp_date), futuredate] <- as.Date(bhps_flood.df$tmp_date[!is.na(bhps_flood.df$tmp_date)],
                                                                         origin = "1970-01-01")
    bhps_flood.df[, futuredate] <- as.Date(ave(as.numeric(bhps_flood.df[, futuredate]),
                                               by = bhps_flood.df$pidp,
                                               FUN = function(x) repeat_last(x, forward = FALSE)), 
                                           origin = "1970-01-01")
    
    bhps_flood.df$tmp_dur <- NA
    bhps_flood.df$tmp_dur[oo] <- bhps_flood.df[oo + 1, pastdur]
    bhps_flood.df[!is.na(bhps_flood.df$tmp_dur), futuredur] <- bhps_flood.df$tmp_dur[!is.na(bhps_flood.df$tmp_dur)]
    bhps_flood.df[, futuredur] <- ave(as.numeric(bhps_flood.df[, futuredur]),
                                               by = bhps_flood.df$pidp,
                                               FUN = function(x) repeat_last(x, forward = FALSE))
    
    bhps_flood.df$tmp_per <- NA
    bhps_flood.df$tmp_per[oo] <- bhps_flood.df[oo + 1, pastper]
    bhps_flood.df[!is.na(bhps_flood.df$tmp_per), futureper] <- bhps_flood.df$tmp_per[!is.na(bhps_flood.df$tmp_per)]
    bhps_flood.df[, futureper] <- ave(as.numeric(bhps_flood.df[, futureper]),
                                      by = bhps_flood.df$pidp,
                                      FUN = function(x) repeat_last(x, forward = FALSE))
    
    
    bhps_flood.df$tmp_date <- bhps_flood.df$tmp_dur <- bhps_flood.df$tmp_per <- NULL
    
  }
}



### Summary Tables 


# 1km
table(bhps_flood.df$flood_affect_past1)
table(bhps_flood.df$flood3_affect_past1)
table(bhps_flood.df$flood10_affect_past1)

table(bhps_flood.df$flood_affect_future1)
table(bhps_flood.df$flood3_affect_future1)
table(bhps_flood.df$flood10_affect_future1)

# 2km
table(bhps_flood.df$flood_affect_past2)
table(bhps_flood.df$flood3_affect_past2)
table(bhps_flood.df$flood10_affect_past2)

table(bhps_flood.df$flood_affect_future2)
table(bhps_flood.df$flood3_affect_future2)
table(bhps_flood.df$flood10_affect_future2)

# 5km
table(bhps_flood.df$flood_affect_past3)
table(bhps_flood.df$flood3_affect_past3)
table(bhps_flood.df$flood10_affect_past3)

table(bhps_flood.df$flood_affect_future3)
table(bhps_flood.df$flood3_affect_future3)
table(bhps_flood.df$flood10_affect_future3)


# Duration
summary(bhps_flood.df[, paste0("flood", rep(intens, length(buf)),
                               "_affect_past_duration", rep(buf, each = length(intens)))])

# Future duration is nor rally meaningful (because mostly fist day is used)

# Percentage
summary(bhps_flood.df[, paste0("flood", rep(intens, length(buf)),
                               "_affect_past_percentage", rep(buf, each = length(intens)))])



################################
### Ever affected by a flood ###
################################

buf <- c("1", "2", "3") # For each buffer
intens <- c("", "3", "10") # and each intensity (percentage)

for(b in buf){
  for(i in intens){
    pname <- paste0("flood", i, "_affect_past", b)
    fname <- paste0("flood", i, "_affect_future", b)
    newname <-  paste0("flood", i, "_affect_ever", b)
    
    bhps_flood.df[, newname] <- ave(bhps_flood.df[, pname] + bhps_flood.df[, fname],
                                    by = bhps_flood.df$pidp,
                                    FUN = function(x) max(x))
    bhps_flood.df[, newname][bhps_flood.df[, newname] > 1] <- 1
  }
}



### Test

nrow(bhps_flood.df[bhps_flood.df$flood_affect_ever1 == 1, ])
nrow(bhps_flood.df[bhps_flood.df$flood_affect_past1 == 1, ])
nrow(bhps_flood.df[bhps_flood.df$flood_affect_future1 == 1 & bhps_flood.df$flood_affect_past1 == 0, ])

nrow(bhps_flood.df[bhps_flood.df$flood3_affect_ever1 == 1, ])
nrow(bhps_flood.df[bhps_flood.df$flood3_affect_past1 == 1, ])
nrow(bhps_flood.df[bhps_flood.df$flood3_affect_future1 == 1 & bhps_flood.df$flood3_affect_past1 == 0, ])

nrow(bhps_flood.df[bhps_flood.df$flood3_affect_ever3 == 1, ])
nrow(bhps_flood.df[bhps_flood.df$flood3_affect_past3 == 1, ])
nrow(bhps_flood.df[bhps_flood.df$flood3_affect_future3 == 1 & bhps_flood.df$flood3_affect_past3 == 0, ])

# ids <- bhps_flood.df$pidp[which(bhps_flood.df$flood_affect_ever == 1 & !(bhps_flood.df$flood_affect_past == 1 | bhps_flood.df$flood_affect_future == 1))]
# View(bhps_flood.df[which(bhps_flood.df$pidp %in% ids), c(1:5, 264, 270, 341:ncol(bhps_flood.df))])




############################################
### Temporal distance to affected floods ###
############################################

bhps_flood.df <- bhps_flood.df[order(bhps_flood.df$pidp, bhps_flood.df$wave), ]


### As above, first identify which floods are actually affecting people
buf <- c("1", "2", "3") # For each buffer
intens <- c("", "3", "10") # and each intensity (percentage)

for(b in buf){
  for(i in intens){
    
    distname <- paste0("tempdist", i, "_affect_past", b)
    distnamef <- paste0("tempdist", i, "_affect_future", b)
    
    pastdate <- paste0("flood", i, "_affect_past_date", b)
    futuredate <- paste0("flood", i, "_affect_future_date", b)
    
    # Temporal distance to interview
    bhps_flood.df[, distname] <- as.numeric(bhps_flood.df$date) - as.numeric(bhps_flood.df[, pastdate])
    bhps_flood.df[, distnamef] <- as.numeric(bhps_flood.df$date) - as.numeric(bhps_flood.df[, futuredate])
  }
}

summary(bhps_flood.df[, paste0("tempdist", rep(intens, length(buf)),
                               "_affect_past", rep(buf, each = length(intens)))])

summary(bhps_flood.df[, paste0("tempdist", rep(intens, length(buf)),
                               "_affect_future", rep(buf, each = length(intens)))])





###################################################################################
### For those living in an lsoa that was ever flooded, feed across entire waves ###
### Indicator of ever lived in an ever flooded lsoa                             ###
###################################################################################

for(b in buf){
  for(i in intens){
    
    n1 <- paste0("flood", i, "_ever", b)
    
    bhps_flood.df[, n1] <- ave(bhps_flood.df[, n1],
                               bhps_flood.df$pidp,
                               FUN = function(x) max(x, na.rm = TRUE))
  }
}



### Save
save(bhps_flood.df, file = "bhpsukhls_flood_final.RData")








####### ---------------------------------------------------------------------- ######
#######         Environmental variables from BHPS US                           ######
####### ---------------------------------------------------------------------- ######


################################
### Attitudes and Intentions ###
################################

### People will be affected by climate change

bhps_flood.df$climatechange30 <- 2 - bhps_flood.df$scopecl30
table(bhps_flood.df$climatechange30)

bhps_flood.df$climatechange200 <- 2 - bhps_flood.df$scopecl200
table(bhps_flood.df$climatechange200)


### Climate change cause heats in UK

bhps_flood.df$climatechange_heats <- 2 - bhps_flood.df$opccb
table(bhps_flood.df$climatechange_heats)


### 1) Climate change influence by behaviour

bhps_flood.df$ccbehavcontributes <- 6 - bhps_flood.df$scenv_bccc
table(bhps_flood.df$ccbehavcontributes)

bhps_flood.df$ccbehavcontributes_w1 <- bhps_flood.df$scenv_ccls - 1 # other way round?!
table(bhps_flood.df$ccbehavcontributes_w1)


### 2) Pay more for env friendly products

bhps_flood.df$ccpaymore <- 6 - bhps_flood.df$scenv_pmep
table(bhps_flood.df$ccpaymore)

bhps_flood.df$ccpaymore_w1 <- 2 - bhps_flood.df$scenv_pmre
table(bhps_flood.df$ccpaymore_w1)


### 3) Soon experience disaster

bhps_flood.df$ccsoondisaster <- 6 - bhps_flood.df$scenv_meds
table(bhps_flood.df$ccsoondisaster)

bhps_flood.df$ccsoondisaster_w1 <- 2 - bhps_flood.df$scenv_dstr
table(bhps_flood.df$ccsoondisaster_w1)


### 4) Climate change been axaggerated

bhps_flood.df$ccexag <- 6 - bhps_flood.df$scenv_crex
table(bhps_flood.df$ccexag)

bhps_flood.df$ccexag_w1 <- 2 - bhps_flood.df$scenv_exag
table(bhps_flood.df$ccexag_w1)


### 5) Climate change beyond control

bhps_flood.df$ccbeyondcontrol <- 6 - bhps_flood.df$scenv_tlat
table(bhps_flood.df$ccbeyondcontrol)

bhps_flood.df$ccbeyondcontrol_w1 <- 2 - bhps_flood.df$scenv_bcon
table(bhps_flood.df$ccbeyondcontrol_w1)


### 6) Climate change to far away to worry

bhps_flood.df$ccnoworry <- 6 - bhps_flood.df$scenv_nowo
table(bhps_flood.df$ccnoworry)

bhps_flood.df$ccnoworry_w1 <- 2 - bhps_flood.df$scenv_futr
table(bhps_flood.df$ccnoworry_w1)


### 7) Changes need to fit lifestyle

bhps_flood.df$ccchangesfitlifestyle <- 6 - bhps_flood.df$scenv_fitl
table(bhps_flood.df$ccchangesfitlifestyle)

bhps_flood.df$ccchangesfitlifestyle_w1 <- 2 - bhps_flood.df$scenv_cfit
table(bhps_flood.df$ccchangesfitlifestyle_w1)


### 8) not worth contributing if others dont

bhps_flood.df$ccnotworthifotherdont <-  6 - bhps_flood.df$scenv_noot
table(bhps_flood.df$ccnotworthifotherdont)

bhps_flood.df$ccnotworthifotherdont_w1 <- 2 - bhps_flood.df$scenv_chwo
table(bhps_flood.df$ccnotworthifotherdont_w1)


### 9) Not worth britain doing anything

bhps_flood.df$ccnotworthbritain <- 6 - bhps_flood.df$scenv_canc
table(bhps_flood.df$ccnotworthbritain)

bhps_flood.df$ccnotworthbritain_w1 <- 2 - bhps_flood.df$scenv_brit
table(bhps_flood.df$ccnotworthbritain_w1)








##################################
### Vote intention green party ###
##################################

table(bhps_flood.df$vote3, useNA = "always")
table(bhps_flood.df$vote4, useNA = "always")

### Combine clossness and vote intention (as variable vote in BHPS)

table(bhps_flood.df$vote3, bhps_flood.df$vote4, useNA = "always")
# Should have receive only one, but not for all cases true

bhps_flood.df$vote_comb <- NA
# Closest to first
oo <- which(bhps_flood.df$vote1 == 1 | bhps_flood.df$vote2 == 1)
bhps_flood.df$vote_comb[oo] <- bhps_flood.df$vote4[oo]
# Vote for tomorrow second
oo <- which(bhps_flood.df$vote1 == 2 & bhps_flood.df$vote2 == 2)
bhps_flood.df$vote_comb[oo] <- bhps_flood.df$vote3[oo]
table(bhps_flood.df$vote_comb, useNA = "always")


# # Set to NA, if different party preferences
# oo <- which(!is.na(bhps_flood.df$vote3) & !is.na(bhps_flood.df$vote4) &
#               bhps_flood.df$vote3 != bhps_flood.df$vote4)
# bhps_flood.df$tmp1 <- bhps_flood.df$vote3
# bhps_flood.df$tmp1[oo] <- NA
# bhps_flood.df$tmp2 <- bhps_flood.df$vote4
# bhps_flood.df$tmp2[oo] <- NA
# 
# # Comine to one variable
# bhps_flood.df$vote_comb <- rowMeans(bhps_flood.df[, c("tmp1", "tmp2")], na.rm = TRUE)
# bhps_flood.df$vote_comb[is.nan(bhps_flood.df$vote_comb)] <- NA
# table(bhps_flood.df$vote_comb, useNA = "always")
# 
# bhps_flood.df[, c("tmp1", "tmp2")] <- NULL



# Recode to green vote preference
bhps_flood.df$vote_green <- bhps_flood.df$vote_comb
bhps_flood.df$vote_green[which(bhps_flood.df$vote_green %in% c(1:5, 7:13))] <- 0
bhps_flood.df$vote_green[which(bhps_flood.df$vote_green %in% c(95, 97))] <- 0 # non-voters and other party
bhps_flood.df$vote_green[which(bhps_flood.df$vote_green %in% c(96))] <- NA # cannot vote to NA
bhps_flood.df$vote_green[which(bhps_flood.df$vote_green %in% c(6))] <- 1 # cannot vote to NA

table(bhps_flood.df$vote_green)


# Recode to not vote
bhps_flood.df$vote_none <- bhps_flood.df$vote_comb
bhps_flood.df$vote_none[which(bhps_flood.df$vote_none %in% c(1:13))] <- 0
bhps_flood.df$vote_none[which(bhps_flood.df$vote_none %in% c(97))] <- 0 # other party
bhps_flood.df$vote_none[which(bhps_flood.df$vote_none %in% c(96))] <- NA # cannot vote to NA
bhps_flood.df$vote_none[which(bhps_flood.df$vote_none %in% c(95))] <- 1 # cannot vote to NA

table(bhps_flood.df$vote_none)




#############################
### Oppinion on politics ###
#############################

table(bhps_flood.df$poleff3)
bhps_flood.df$pol_nocare <- 6 - bhps_flood.df$poleff3
table(bhps_flood.df$pol_nocare)

table(bhps_flood.df$poleff4)
bhps_flood.df$pol_nosay <- 6 - bhps_flood.df$poleff4
table(bhps_flood.df$pol_nocare)

# Index of fustration
bhps_flood.df$pol_frust <- (bhps_flood.df$pol_nocare + bhps_flood.df$pol_nosay)/2
table(bhps_flood.df$pol_frust)




##############################
### Environmental behavior ###
##############################

vars <- paste0("envhabit", c(1:11))

### code 6 as missing

for(i in vars){
  bhps_flood.df[which(bhps_flood.df[, i] == 6), i] <- NA
  print(table(bhps_flood.df[, i]))
}


### Recode that positive is more environ friendly 2 4 5 6 7 8 9 10 11
c <- c(2, 4, 5, 6, 7, 8, 9, 10, 11)

for(i in c){
  var <- paste0("envhabit", i)
  bhps_flood.df[, var] <- 6 - bhps_flood.df[, var]
  print(table(bhps_flood.df[, var]))
}



### Principal components
tmp.df <- na.omit(bhps_flood.df[, vars]) # Only non missing values

pca.env <- princomp(tmp.df, cor=TRUE)
summary(pca.env)

pca.env <- psych::principal(tmp.df, nfactors = 3, rotate = "varimax")
print(pca.env)




### General behavioral index
bhps_flood.df$env_index_home <- (bhps_flood.df$envhabit1 +
                                         bhps_flood.df$envhabit2 +
                                         bhps_flood.df$envhabit3 +
                                         bhps_flood.df$envhabit4 +
                                         bhps_flood.df$envhabit7)/5

cronbach <- psych::alpha(bhps_flood.df[, paste0("envhabit", c(1:4, 7))])
print(cronbach)


bhps_flood.df$env_index_shop <- (bhps_flood.df$envhabit5 +
                                         bhps_flood.df$envhabit6)/2

bhps_flood.df$env_index_trans <- (bhps_flood.df$envhabit8 +
                                         bhps_flood.df$envhabit9)/2

# All available items
bhps_flood.df$env_index_all <- (bhps_flood.df$envhabit1 +
                                        bhps_flood.df$envhabit2 +
                                        bhps_flood.df$envhabit3 +
                                        bhps_flood.df$envhabit4 +
                                        bhps_flood.df$envhabit5 +
                                        bhps_flood.df$envhabit6 +
                                        bhps_flood.df$envhabit7 +
                                        bhps_flood.df$envhabit8 +
                                        bhps_flood.df$envhabit9 +
                                        bhps_flood.df$envhabit10 +
                                        bhps_flood.df$envhabit11)/11

cronbach <- psych::alpha(bhps_flood.df[, paste0("envhabit", c(1:11))])
print(cronbach)

# All items available in BHPS and UKHLS
bhps_flood.df$env_index_all2 <- (bhps_flood.df$envhabit1 +
                                        bhps_flood.df$envhabit2 +
                                        bhps_flood.df$envhabit3 +
                                        bhps_flood.df$envhabit4 +
                                        bhps_flood.df$envhabit5 +
                                        bhps_flood.df$envhabit6 +
                                        bhps_flood.df$envhabit7)/7

cronbach <- psych::alpha(bhps_flood.df[, paste0("envhabit", c(1:7))])
print(cronbach)




### Car use frequency

table(bhps_flood.df$trcarfq)
bhps_flood.df$car_freq <- 9 - bhps_flood.df$trcarfq
table(bhps_flood.df$car_freq)


### Cycling use frequency

table(bhps_flood.df$trbikefq)
bhps_flood.df$bike_freq <- 9 - bhps_flood.df$trbikefq
table(bhps_flood.df$bike_freq)


###########################################
### Prepare socio-demographic variables ###
###########################################

getmode <- function(v) {
  uniqv <- unique(v[!is.na(v)])
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

### Age
table(bhps_flood.df$age_dv)

# Age for bhps
oo <- which(bhps_flood.df$wave <= 18)
bhps_flood.df$age_dv[oo] <- bhps_flood.df$age[oo]

# Cut into bins
cuts <- seq(min(bhps_flood.df$age_dv, na.rm = T), max(bhps_flood.df$age_dv, na.rm = T), 5)
cuts <- c(cuts[-which(cuts > 90)],  max(bhps_flood.df$age_dv, na.rm = T))
bhps_flood.df$age_cat <- cut(bhps_flood.df$age_dv,
                             breaks = cuts)
table(bhps_flood.df$age_cat)

bhps_flood.df$age_cat <- relevel(as.factor(bhps_flood.df$age_cat),
                                 ref = as.numeric(getmode(bhps_flood.df$age_cat)))

### Gender
table(bhps_flood.df$sex)

bhps_flood.df$female <- bhps_flood.df$sex - 1

### Migration to UK
table(bhps_flood.df$bornuk_dv)
bhps_flood.df$migback <- bhps_flood.df$bornuk_dv - 1


### Migration background generation
table(bhps_flood.df$generation)
bhps_flood.df$migback_gen <- bhps_flood.df$generation
bhps_flood.df$migback_gen[bhps_flood.df$migback_gen %in% c(4:6)] <- 0
table(bhps_flood.df$migback_gen)


### Ethnic group
table(bhps_flood.df$ethn_dv)
bhps_flood.df$ethn_dv <- relevel(as.factor(bhps_flood.df$ethn_dv),
                                 ref = as.numeric(getmode(bhps_flood.df$ethn_dv)))


### Ethnic group short
bhps_flood.df$ethn_dv_short <- as.numeric(bhps_flood.df$ethn_dv)
bhps_flood.df$ethn_dv_short[bhps_flood.df$ethn_dv_short %in% c(1:4)] <- 0
bhps_flood.df$ethn_dv_short[bhps_flood.df$ethn_dv_short %in% c(5:8)] <- 1
bhps_flood.df$ethn_dv_short[bhps_flood.df$ethn_dv_short %in% c(9:13)] <- 2
bhps_flood.df$ethn_dv_short[bhps_flood.df$ethn_dv_short %in% c(14:16)] <- 3
bhps_flood.df$ethn_dv_short[bhps_flood.df$ethn_dv_short %in% c(17:97)] <- 4

bhps_flood.df$ethn_dv_short <- as.factor(bhps_flood.df$ethn_dv_short)

table(bhps_flood.df$ethn_dv_short)


### Education
table(bhps_flood.df$hiqual_dv)

bhps_flood.df$hiqual_dv <- relevel(as.factor(bhps_flood.df$hiqual_dv),
                                   ref = as.numeric(getmode(bhps_flood.df$hiqual_dv)))


### Child(ren) in household
table(bhps_flood.df$nkids_dv)
bhps_flood.df$child <- bhps_flood.df$nkids_dv
bhps_flood.df$child[bhps_flood.df$child > 0 ] <- 1
table(bhps_flood.df$child)

### Marital status
table(bhps_flood.df$marstat_dv)

bhps_flood.df$marstat_dv <- relevel(as.factor(bhps_flood.df$marstat_dv),
                                    ref = as.numeric(getmode(bhps_flood.df$marstat_dv)))

### Household total income
summary(bhps_flood.df$fimngrs_dv)
summary(bhps_flood.df$fimnlabgrs_dv)
summary(bhps_flood.df$fihhmngrs_dv)
hist(bhps_flood.df$fihhmngrs_dv)

bhps_flood.df$hhinc <- bhps_flood.df$fihhmngrs_dv/1000

summary(bhps_flood.df$hhinc)

# Drop outlier (below zero or above 100.000 per month )
oo <- which(bhps_flood.df$hhinc < 0 |
              bhps_flood.df$hhinc > 100)
bhps_flood.df$hhinc[oo] <- NA

# Squared income
bhps_flood.df$hhinc_sq <- bhps_flood.df$hhinc^2


### OECD equivalence income (square root method for simplicity)
bhps_flood.df$hhinc_eq <- bhps_flood.df$hhinc / sqrt(bhps_flood.df$hhsize)
bhps_flood.df$hhinc_eq_sq <- bhps_flood.df$hhinc_eq^2

summary(bhps_flood.df$hhinc_eq)


### Income deciles
breaks <- quantile(bhps_flood.df$hhinc[which(bhps_flood.df$wave %in% c(18, 19, 22, 28))], probs = seq(0, 1, by = 0.1), na.rm = TRUE)
bhps_flood.df$hhinc_dec <- cut(bhps_flood.df$hhinc, breaks = breaks, include.lowest = TRUE)



### Home owner
table(bhps_flood.df$tenure_dv, useNA = "always")

bhps_flood.df$owner <- bhps_flood.df$tenure_dv
bhps_flood.df$owner[bhps_flood.df$owner %in% c(1,2)] <- 1
bhps_flood.df$owner[bhps_flood.df$owner %in% c(3:8)] <- 0
table(bhps_flood.df$owner)


### Moving Preference
bhps_flood.df$lkmove <- bhps_flood.df$lkmove - 1



### Political party

# No party support as reference
bhps_flood.df$vote_comb <- relevel(as.factor(bhps_flood.df$vote_comb),
                                   ref = "95")



### Satisfaction with income

# Variable for bhps and US
oo <- which(bhps_flood.df$wave <= 18)
bhps_flood.df$sat_inc <- ifelse(bhps_flood.df$wave <= 18, 
                                bhps_flood.df$lfsat2,
                                bhps_flood.df$sclfsat2)

# Make ordinal                                    
table(bhps_flood.df$sat_inc)
bhps_flood.df$sat_inc[which(bhps_flood.df$sat_inc == 0)] <- NA
bhps_flood.df$sat_inc <- as.factor(bhps_flood.df$sat_inc)


### 4 seasons
bhps_flood.df$season <- NA
bhps_flood.df$season[bhps_flood.df$istrtdatm %in% c(3:5)] <- 1
bhps_flood.df$season[bhps_flood.df$istrtdatm %in% c(6:8)] <- 2
bhps_flood.df$season[bhps_flood.df$istrtdatm %in% c(9:11)] <- 3
bhps_flood.df$season[bhps_flood.df$istrtdatm %in% c(12, 1, 2)] <- 4
table(bhps_flood.df$season)


### Rolling season index

# Year -> previous year if month is Jan or Feb (winter season)
bhps_flood.df$tmp <- bhps_flood.df$istrtdaty
oo <- which(bhps_flood.df$istrtdatm %in% c(1,2))
bhps_flood.df$tmp[oo] <- bhps_flood.df$tmp[oo] - 1

# Gen index based on transformed year and season
bhps_flood.df$season_index <- paste(bhps_flood.df$tmp, bhps_flood.df$season, sep = "_")
table(bhps_flood.df$season_index)
table(bhps_flood.df$season_index[bhps_flood.df$istrtdaty >=2008], bhps_flood.df$istrtdaty[bhps_flood.df$istrtdaty >=2008])

bhps_flood.df$tmp <- NULL




####################################
### Impute vote comp for wave 10 ###
####################################

### Repeat last function
repeat_last = function(x, forward = TRUE, maxgap = Inf, na.rm = FALSE) {
  if (!forward) x = rev(x)           # reverse x twice if carrying backward
  ind = which(!is.na(x))             # get positions of nonmissing values
  if (is.na(x[1]) && !na.rm)         # if it begins with NA
    ind = c(1,ind)                 # add first pos
  rep_times = diff(                  # diffing the indices + length yields how often
    c(ind, length(x) + 1) )          # they need to be repeated
  if (maxgap < Inf) {
    exceed = rep_times - 1 > maxgap  # exceeding maxgap
    if (any(exceed)) {               # any exceed?
      ind = sort(c(ind[exceed] + 1, ind))      # add NA in gaps
      rep_times = diff(c(ind, length(x) + 1) ) # diff again
    }
  }
  x = rep(x[ind], times = rep_times) # repeat the values at these indices
  if (!forward) x = rev(x)           # second reversion
  x
}

# Order
bhps_flood.df <- bhps_flood.df[order(bhps_flood.df$pidp, bhps_flood.df$wave),]

### Use previous wave to impute missing wave 10 vote_comp 
oo <- which(bhps_flood.df$wave == 27 | bhps_flood.df$wave == 28)

# Feed forward
bhps_flood.df$vote_comb_imp <- bhps_flood.df$vote_comb

bhps_flood.df$vote_comb_imp[] <- ave(bhps_flood.df$vote_comb_imp[],
                                     bhps_flood.df$pidp[],
                                     FUN = function(x) repeat_last(x, maxgap = 1))





#########################
### "Personal traits" ###
#########################

### Repeat last function
repeat_last = function(x, forward = TRUE, maxgap = Inf, na.rm = FALSE) {
  if (!forward) x = rev(x)           # reverse x twice if carrying backward
  ind = which(!is.na(x))             # get positions of nonmissing values
  if (is.na(x[1]) && !na.rm)         # if it begins with NA
    ind = c(1,ind)                 # add first pos
  rep_times = diff(                  # diffing the indices + length yields how often
    c(ind, length(x) + 1) )          # they need to be repeated
  if (maxgap < Inf) {
    exceed = rep_times - 1 > maxgap  # exceeding maxgap
    if (any(exceed)) {               # any exceed?
      ind = sort(c(ind[exceed] + 1, ind))      # add NA in gaps
      rep_times = diff(c(ind, length(x) + 1) ) # diff again
    }
  }
  x = rep(x[ind], times = rep_times) # repeat the values at these indices
  if (!forward) x = rev(x)           # second reversion
  x
}



### Partisanship
table(bhps_flood.df$vote4, useNA = "ifany")
bhps_flood.df$partisan <- bhps_flood.df$vote4
# Only >= fairly strong support
oo <- which(bhps_flood.df$vote5 == 3)
bhps_flood.df$partisan[oo] <- 99 # no partisanship
# if answered 2 in vote1 and vote2 code as no partisanship
oo <- which(bhps_flood.df$vote1 == 2 & bhps_flood.df$vote2 == 2)
bhps_flood.df$partisan[oo] <- 99 # no partisanship
table(bhps_flood.df$partisan, useNA = "ifany")
bhps_flood.df$partisan <- factor(bhps_flood.df$partisan)


### Problem solving
bhps_flood.df$problemsolv <- (bhps_flood.df$scghqd + bhps_flood.df$scghqh + bhps_flood.df$scghqf)/3



### Climate change related vars
table(bhps_flood.df$opccb)
bhps_flood.df$ccfloods <- 2 - bhps_flood.df$opccb

table(bhps_flood.df$scenv_bcon)
bhps_flood.df$ccbeyondcontrol <- 6 - bhps_flood.df$scenv_tlat

table(bhps_flood.df$scenv_ccls)
bhps_flood.df$ccothersdontcont <-  6 - bhps_flood.df$scenv_noot

table(bhps_flood.df$scenv_bccc)
bhps_flood.df$ccbehavcontributes <- 6 - bhps_flood.df$scenv_bccc



### Trust vars
table(bhps_flood.df$sctrust, useNA = "ifany")

bhps_flood.df$trust <- bhps_flood.df$sctrust
bhps_flood.df$trust[which(bhps_flood.df$trust == 2 | bhps_flood.df$trust == 3)] <- 0
table(bhps_flood.df$trust)

# Trust measure 2 (attention, in UKHLS 11 points and BHPS 10 points. Use above 5 as trust high)
table(bhps_flood.df$riskb, useNA = "ifany")
table(bhps_flood.df$scriskb, useNA = "ifany")

bhps_flood.df$trust2 <- bhps_flood.df$scriskb
bhps_flood.df$trust2[which(bhps_flood.df$wave == 18)] <- bhps_flood.df$riskb[which(bhps_flood.df$wave == 18)]
bhps_flood.df$trust3 <- bhps_flood.df$trust2
bhps_flood.df$trust2[which(bhps_flood.df$trust2 <= 5)] <- 0
bhps_flood.df$trust2[which(bhps_flood.df$trust2 > 5)] <- 1
table(bhps_flood.df$trust2)



### Cooperation / volunteering
bhps_flood.df$cooperation <- bhps_flood.df$ivcoop
bhps_flood.df$cooperation[which(bhps_flood.df$wave <= 18)] <- bhps_flood.df$iv4[which(bhps_flood.df$wave <= 18)]
table(bhps_flood.df$wave[!is.na(bhps_flood.df$cooperation)])
bhps_flood.df$cooperation <- 6 - bhps_flood.df$cooperation

# Volunteering
bhps_flood.df$voluntary_work <- ifelse(bhps_flood.df$volun == 1, 1, 0)
table(bhps_flood.df$wave[!is.na(bhps_flood.df$voluntary_work)])

bhps_flood.df$voluntary_work_freq <- 10 - bhps_flood.df$volfreq
bhps_flood.df$voluntary_work_freq[bhps_flood.df$voluntary_work_freq %in% c(1:3)] <- 1
bhps_flood.df$voluntary_work_freq[bhps_flood.df$voluntary_work_freq %in% c(4:6)] <- 2
bhps_flood.df$voluntary_work_freq[bhps_flood.df$voluntary_work_freq %in% c(7:9)] <- 3
table(bhps_flood.df$wave[!is.na(bhps_flood.df$voluntary_work_freq)])

oo <- which(bhps_flood.df$voluntary_work == 0 & bhps_flood.df$wave %in% c(20, 22, 24, 26, 28))
bhps_flood.df$voluntary_work_freq[oo] <- 0
table(bhps_flood.df$voluntary_work_freq)

# Charity spending
bhps_flood.df$charity <- ifelse(bhps_flood.df$chargv == 1, 1, 0)
table(bhps_flood.df$charity)
table(bhps_flood.df$wave[!is.na(bhps_flood.df$charity)])

bhps_flood.df$charity_freq <- bhps_flood.df$charfreq
bhps_flood.df$charity_freq[bhps_flood.df$charity == 0] <- 0
table(bhps_flood.df$charity_freq)




######################################
#### Feed forward first wave vars ####
######################################

### For time varying vars, use the first wave information in FE models
vars <- c("climatechange30", "env_index_all2",
          "age_dv", "female", "migback", "ethn_dv_short", "hiqual_dv", "child",  "marstat_dv", "hhinc",
          "vote_comb", "partisan", 
          "ccbehavcontributes_w1", "ccpaymore_w1", "ccsoondisaster_w1", "ccexag_w1", "ccbeyondcontrol_w1", 
          "ccnoworry_w1", "ccchangesfitlifestyle_w1", "ccnotworthifotherdont_w1", "ccnotworthbritain_w1",
          "env_lifestylecur", "env_lifestyleint", "greenisalternative",
          "problemsolv", "trust", "trust2", "trust3", "cooperation")
bhps_flood.df <- bhps_flood.df[order(bhps_flood.df$pidp, bhps_flood.df$wave),]
for(v in vars){
  nv <- paste0("c_", v)
  nv2 <- paste0("a_", v)
  
  oo <- which(bhps_flood.df$wave == 18 | bhps_flood.df$wave == 19)
  bhps_flood.df[oo, nv] <-  bhps_flood.df[oo, v] 
  
  if(is.factor(bhps_flood.df[, v])){
    lev <- levels(droplevels(bhps_flood.df[, nv]))
    bhps_flood.df[, nv] <- as.numeric(bhps_flood.df[, nv])
  }
  
  bhps_flood.df[, nv] <- ave(bhps_flood.df[, nv], bhps_flood.df$pidp, 
                            FUN = function(x) repeat_last(x))
  
  if(is.factor(bhps_flood.df[, v])){
    bhps_flood.df[, nv] <- factor(bhps_flood.df[, nv], labels = lev)
  }
  
  ### Add average
  if(!is.factor(bhps_flood.df[, v])){
    bhps_flood.df[, nv2] <- ave(bhps_flood.df[, v], bhps_flood.df$pidp, 
                               FUN = function(x) mean(x, na.rm = TRUE))
  }

  
}


### For voluntary work and charity use second wave data
vars <- c("voluntary_work", "voluntary_work_freq", "charity", "charity_freq", "pol_nocare", "pol_nosay")
for(v in vars){
  nv <- paste0("c_", v)
  
  oo <- which(bhps_flood.df$wave == 20 | bhps_flood.df$wave == 21)
  bhps_flood.df[oo, nv] <-  bhps_flood.df[oo, v] 
  
  if(is.factor(bhps_flood.df[, v])){
    lev <- levels(droplevels(bhps_flood.df[, nv]))
    bhps_flood.df[, nv] <- as.numeric(bhps_flood.df[, nv])
  }
  
  bhps_flood.df[, nv] <- ave(bhps_flood.df[, nv], bhps_flood.df$pidp, 
                            FUN = function(x) mean(x, na.rm = TRUE))
  
  if(is.factor(bhps_flood.df[, v])){
    bhps_flood.df[, nv] <- factor(bhps_flood.df[, nv], labels = lev)
  }
}




### Self-efficacy (needs backward imputation)

v <- paste0("se", c(1:10))
for(i in v){
  print(table(bhps_flood.df[, i]))
}

bhps_flood.df$selfeff <- rowSums(bhps_flood.df[,v])/10
table(bhps_flood.df$selfeff)

# Feed backward from wave 5
bhps_flood.df$c_selfeff <- ave(bhps_flood.df$selfeff, bhps_flood.df$pidp, 
                               FUN = function(x) repeat_last(x, forward = FALSE))






############
### Save ###
############


save(bhps_flood.df, file = "bhpsukhls_flood_final.RData")





