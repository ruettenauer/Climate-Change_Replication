
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


load("bhpsukhls_heat.RData")



# Merge new variables (if necessay)
bhps.df <- read.dta("all_bhpsukhls.dta")
oo <- names(bhps.df)[which(!names(bhps.df) %in% names(bhps_heat.df))]
bhps_heat.df <- merge(bhps_heat.df, bhps.df[, c("pidp", "wave", oo)],
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
bhps_heat.df <- bhps_heat.df[order(bhps_heat.df$pidp, bhps_heat.df$wave), ]
bhps_heat.df$pynr <- ave(bhps_heat.df$wave,
                                by = bhps_heat.df$pidp,
                                FUN = function(x) 1:length(x))



### Clean plnowy4
bhps_heat.df$plnowy4[which(bhps_heat.df$plnowy4 > 2030)] <- NA

### Clean plnowm
bhps_heat.df$plnowm[which(bhps_heat.df$plnowm > 12)] <- 1


### Impute birth date for households never moved in their life
oo <- which(bhps_heat.df$mvever == 1 & is.na(bhps_heat.df$plnowy4))
bhps_heat.df$plnowy4[oo] <- bhps_heat.df$birthy[oo]
bhps_heat.df$plnowm[oo] <- bhps_heat.df$birthm[oo]

### Replace plnow with mvyr if NA
bhps_heat.df$plnowy4[which(is.na(bhps_heat.df$plnowy4))] <- bhps_heat.df$mvyr[which(is.na(bhps_heat.df$plnowy4))]
bhps_heat.df$plnowm[which(is.na(bhps_heat.df$plnowm))] <- bhps_heat.df$mvmnth[which(is.na(bhps_heat.df$plnowm))]


### Construct Move-in year

# Order date
bhps_heat.df <- bhps_heat.df[order(bhps_heat.df$pidp, bhps_heat.df$year),] 

# Impute last moving date
bhps_heat.df$movein_yr <- ave(bhps_heat.df$plnowy4,
                                       by = bhps_heat.df$pidp,
                                       FUN = function(x) repeat_last(x))


# ### For some obs the first move in year is missing (for instance children)
# # Impute the households pno1 movein value
# bhps_heat.df$tmp <- bhps_heat.df$movein_yr
# bhps_heat.df$tmp[bhps_heat.df$pno != 1] <- NA
# bhps_heat.df$tmp <- ave(bhps_heat.df$tmp,
#                          by = paste0(bhps_heat.df$hidp, "_", bhps_heat.df$wave),
#                          FUN = function(x) max(x, na.rm = TRUE))


# Children make a big part of those not having a move in date
# Could also use parents move in date or birthday
# However, for know use their age 16 survey as move in date.

# For those without moving year, use first observation year at location
bhps_heat.df$combid <- paste0(bhps_heat.df$pidp, "_", bhps_heat.df$lsoa01) 
bhps_heat.df$tmp <- ave(bhps_heat.df$istrtdaty,
                                     by = bhps_heat.df$combid,
                                     FUN = function(x) min(x))

bhps_heat.df$movein_yr[is.na(bhps_heat.df$movein_yr)] <- bhps_heat.df$tmp[is.na(bhps_heat.df$movein_yr)]
bhps_heat.df$tmp <- NULL




### Contruct Move-in month

# Order date
bhps_heat.df <- bhps_heat.df[order(bhps_heat.df$pidp, bhps_heat.df$year),] 

# Impute last moving date
bhps_heat.df$movein_mn <- ave(bhps_heat.df$plnowm,
                                     by = bhps_heat.df$pidp,
                                     FUN = function(x) repeat_last(x))

# For those without moving year, use first observation year at location
bhps_heat.df$combid <- paste0(bhps_heat.df$pidp, bhps_heat.df$lsoa01) 
bhps_heat.df$tmp <- ave(bhps_heat.df$istrtdatm,
                               by = bhps_heat.df$combid,
                               FUN = function(x) min(x))

# If there before 1991 and month is na, use random month
set.seed(1241234)
oo <- which(is.na(bhps_heat.df$movein_mn) & bhps_heat.df$movein_yr < 1991)
bhps_heat.df$movein_mn[oo] <- ave(bhps_heat.df$movein_mn[oo],
                                         by = bhps_heat.df$pidp[oo],
                                         FUN = function(x) rep(sample(c(1:12), 1), length(x)))
  
# Replace remaining with int month
bhps_heat.df$movein_mn[is.na(bhps_heat.df$movein_mn)] <- bhps_heat.df$tmp[is.na(bhps_heat.df$movein_mn)]
bhps_heat.df$tmp <- NULL


### Contruct Move-in Date

bhps_heat.df$movein_date <- as.Date(paste("15", bhps_heat.df$movein_mn, bhps_heat.df$movein_yr, sep = "/"),
                                           "%d/%m/%Y")

summary(bhps_heat.df$movein_date)






############################
### Resctrict the sample ###
############################



# ### Restrict to England only
# bhps_heat.df <- bhps_heat.df[which(bhps_heat.df$gor_dv < 10), ]

### Drop BHPS UKHPS without lsoa or date
bhps_heat.df <- bhps_heat.df[which(!is.na(bhps_heat.df$lsoa01) & !is.na(bhps_heat.df$date)), ]




############################################################
### Generate period (move in and move out) for each lsoa ###
############################################################

bhps_heat.df <- bhps_heat.df[order(bhps_heat.df$pidp, bhps_heat.df$wave), ]

### Indicator for change in lsoa
bhps_heat.df$tmp <- ave(bhps_heat.df$lsoa01,
                         by = bhps_heat.df$pidp,
                         FUN = function(x) dplyr::lag(x))
  
bhps_heat.df$lsoa_change <- 0
oo <- which(bhps_heat.df$lsoa01 != bhps_heat.df$tmp & !is.na(bhps_heat.df$tmp))
bhps_heat.df$lsoa_change[oo] <- 1
bhps_heat.df$tmp <- NULL

### Cumulative number of lsoa spell per person
bhps_heat.df$lsoa_nr <- ave(bhps_heat.df$lsoa_change,
                             by = bhps_heat.df$pidp,
                             FUN = function(x) cumsum(x))

### Move in date for each lsoa spell (use first non-missing)
bhps_heat.df$movein_lsoa_date <- as.Date(ave(as.numeric(bhps_heat.df$movein_date),
                                              by = paste0(bhps_heat.df$pidp, "_", bhps_heat.df$lsoa_nr),
                                              FUN = function(x) x[!is.na(x)][1]), origin = "1970-01-01")




### Clean move in date if lsoa change, but no new move-in date (seems mostly people moving often)
bhps_heat.df$tmp_diff <- ave(as.numeric(bhps_heat.df$movein_lsoa_date),
                              by = bhps_heat.df$pidp,
                              FUN = function(x) x - dplyr::lag(x))

bhps_heat.df$tmp_diff2 <- ave(as.numeric(bhps_heat.df$date),
                               by = bhps_heat.df$pidp,
                               FUN = function(x) x - dplyr::lag(x))

# ids <- bhps_heat.df$pidp[which(bhps_heat.df$tmp_diff <= 65 & bhps_heat.df$lsoa_change == 1)]
# View(bhps_heat.df[bhps_heat.df$pidp %in% ids, c(1:5, 9:11, 36:37, 40:45, 264, 270, 341:ncol(bhps_heat.df))])

# Impute for those having no or negative change in move in
bhps_heat.df$impute_mvin <- 0
oo <- which(bhps_heat.df$tmp_diff <= 0 & bhps_heat.df$lsoa_change == 1)
bhps_heat.df$impute_mvin[oo] <- 1

# Rule: impute random date between the two interviews --------------
set.seed(54654887)
bhps_heat.df$tmp_randdiff <- NA
oo <- which(bhps_heat.df$impute_mvin == 1)
bhps_heat.df$tmp_randdiff[oo] <- as.numeric(lapply(bhps_heat.df$tmp_diff2[oo], 
                                              FUN = function(x) sample(0:x, size = 1)))

oo <- which(bhps_heat.df$impute_mvin == 1)
bhps_heat.df$movein_lsoa_date[oo] <- as.Date(as.numeric(bhps_heat.df$date[oo]) - bhps_heat.df$tmp_randdiff[oo],
                                              origin = "1970-01-01")

bhps_heat.df$tmp_diff <- NULL
bhps_heat.df$tmp_diff2 <- NULL
bhps_heat.df$tmp_randdiff <- NULL



# Replace Move in date for each lsoa spell by the more recent date (fill imputation for following years)
bhps_heat.df$movein_lsoa_date <- as.Date(ave(as.numeric(bhps_heat.df$movein_lsoa_date),
                                              by = paste0(bhps_heat.df$pidp, "_", bhps_heat.df$lsoa_nr),
                                              FUN = function(x) max(x)), origin = "1970-01-01")




### Gen moveout date for each lsoa

# Lead move in date if lsoa change
bhps_heat.df$tmp <- bhps_heat.df$movein_date
bhps_heat.df$tmp[which(bhps_heat.df$lsoa_change == 0)] <- NA

bhps_heat.df$tmp <- as.Date(ave(as.numeric(bhps_heat.df$tmp),
                                 by = bhps_heat.df$pidp,
                                 FUN = function(x) dplyr::lead(x)), origin = "1970-01-01")

# Use last interview data for the final wave of each person
bhps_heat.df$pynr <- ave(bhps_heat.df$pidp,
                          by = bhps_heat.df$pidp,
                          FUN = function(x) 1:length(x))
bhps_heat.df$pN <- ave(bhps_heat.df$pidp,
                        by = bhps_heat.df$pidp,
                        FUN = function(x) length(x))

oo <- which(bhps_heat.df$pN == bhps_heat.df$pynr)
bhps_heat.df$tmp[oo] <- bhps_heat.df$date[oo] 

# For each lsoa spell, use max move out date
bhps_heat.df$moveout_lsoa_date <- as.Date(ave(as.numeric(bhps_heat.df$tmp),
                                              by = paste0(bhps_heat.df$pidp, "_", bhps_heat.df$lsoa_nr),
                                              FUN = function(x) max(x, na.rm = TRUE)), 
                                           origin = "1970-01-01")
bhps_heat.df$tmp <- NULL











#########################
### Affected by heat ###
#########################

bhps_heat.df <- bhps_heat.df[order(bhps_heat.df$pidp, bhps_heat.df$wave), ]




### Calculate whether hit at place of residence ###

# Future date
bhps_heat.df$date_lead <- as.Date(ave(as.numeric(bhps_heat.df$date),
                                       by = bhps_heat.df$pidp,
                                       FUN = function(x) dplyr::lead(x)), 
                                   origin = "1970-01-01")



### Past heats ###
buf <- c(1, 2, 3) # For each buffer
intens <- c("") # and each intensity (percentage)
names(bhps_heat.df)[names(bhps_heat.df) == "past_date"] <- "past_heat_date"
names(bhps_heat.df)[names(bhps_heat.df) == "future_date"] <- "future_heat_date"


for(b in buf){
  for(i in intens){
    newname <- paste0("heat", i, "_affect_past", b)
    pastdate <- paste0("past_heat", i, "_date", b)
    futuredate <- paste0("future_heat", i, "_date", b)
    
    newnamedate <- paste0("heat", i, "_affect_past_date", b)
    # newnamedur <- paste0("heat", i, "_affect_past_duration", b)
    newnameper <- paste0("heat", i, "_affect_past_tasmax", b)
    
    pastdur <- paste0("past_heat", i, "_duration", b)
    futuredur <- paste0("future_heat", i, "_duration", b)
    pastper <- paste0("past_heat", i, "_tasmax", b)
    futureper <- paste0("future_heat", i, "_tasmax", b)
    
    bhps_heat.df[, newname] <- 0
    
    # Directly affected affected
    oo <- which(bhps_heat.df[, pastdate] >= bhps_heat.df$movein_lsoa_date)
    bhps_heat.df[, newname][oo] <- 1
    
    # Also save date
    bhps_heat.df$tmp_date <- NA
    bhps_heat.df$tmp_date[oo] <- bhps_heat.df[oo, pastdate]
    
    # # Also save duration
    # bhps_heat.df$tmp_dur <- NA
    # bhps_heat.df$tmp_dur[oo] <- bhps_heat.df[oo, pastdur]
    
    # Also save percentage
    bhps_heat.df$tmp_per <- NA
    bhps_heat.df$tmp_per[oo] <- bhps_heat.df[oo, pastper]
    
    # Affected at old lsoa but moved afterwards
    bhps_heat.df$tmp <- NA 
    oo <- which(bhps_heat.df[, futuredate] > bhps_heat.df$date & 
                  bhps_heat.df[, futuredate] <= bhps_heat.df$moveout_lsoa_date & 
                  bhps_heat.df[, futuredate] < bhps_heat.df$date_lead)
    bhps_heat.df$tmp[oo] <- 1
    
    # Lag and declare as affected
    bhps_heat.df$tmp <- ave(bhps_heat.df$tmp,
                             by = bhps_heat.df$pidp,
                             FUN = function(x) dplyr::lag(x))
    
    bhps_heat.df[, newname][bhps_heat.df$tmp == 1] <- 1
    bhps_heat.df$tmp <- NULL
    
    # Also save date
    bhps_heat.df$tmp_date2 <- NA
    bhps_heat.df$tmp_date2[oo] <- bhps_heat.df[oo, futuredate]
    
    # Lag date 
    bhps_heat.df$tmp_date2 <- ave(bhps_heat.df$tmp_date2,
                                  by = bhps_heat.df$pidp,
                                  FUN = function(x) dplyr::lag(x))
    
    # # Also save dur
    # bhps_heat.df$tmp_dur2 <- NA
    # bhps_heat.df$tmp_dur2[oo] <- bhps_heat.df[oo, futuredur]
    # 
    # # Lag dur 
    # bhps_heat.df$tmp_dur2 <- ave(bhps_heat.df$tmp_dur2,
    #                                by = bhps_heat.df$pidp,
    #                                FUN = function(x) dplyr::lag(x))
    
    # Also save per
    bhps_heat.df$tmp_per2 <- NA
    bhps_heat.df$tmp_per2[oo] <- bhps_heat.df[oo, futureper]
    
    # Lag per 
    bhps_heat.df$tmp_per2 <- ave(bhps_heat.df$tmp_per2,
                                  by = bhps_heat.df$pidp,
                                  FUN = function(x) dplyr::lag(x))
    
    # Use more recent date and duration
    oo <- which(bhps_heat.df$tmp_date2 > bhps_heat.df$tmp_date | is.na(bhps_heat.df$tmp_date))
    bhps_heat.df$tmp_date[oo] <- bhps_heat.df$tmp_date2[oo]
    # bhps_heat.df$tmp_dur[oo] <- bhps_heat.df$tmp_dur2[oo]
    bhps_heat.df$tmp_per[oo] <- bhps_heat.df$tmp_per2[oo]
    
    bhps_heat.df$tmp_date2 <- bhps_heat.df$tmp_per2 <- NA
    
    
    # Fill na
    bhps_heat.df[, newnamedate] <- as.Date(ave(bhps_heat.df$tmp_date,
                                        by = bhps_heat.df$pidp,
                                        FUN = function(x) repeat_last(x, forward = TRUE)),
                                        origin = "1970-01-01")
    # bhps_heat.df[, newnamedur] <- ave(bhps_heat.df$tmp_dur,
    #                                     by = bhps_heat.df$pidp,
    #                                     FUN = function(x) repeat_last(x, forward = TRUE))
    bhps_heat.df[, newnameper] <- ave(bhps_heat.df$tmp_per,
                                       by = bhps_heat.df$pidp,
                                       FUN = function(x) repeat_last(x, forward = TRUE))
    
    bhps_heat.df$tmp_date <- bhps_heat.df$tmp_per <- NA
    
    # Fill NAs after occasion
    bhps_heat.df[, newname] <- ave(bhps_heat.df[, newname],
                                    by = bhps_heat.df$pidp,
                                    FUN = function(x) cummax(x))
    
    
  }
  
}



### Future heats ###
buf <- c(1, 2, 3) # For each buffer
intens <- c("") # and each intensity (percentage)

for(b in buf){
  for(i in intens){
    newname <- paste0("heat", i, "_affect_future", b)
    pastdate <- paste0("past_heat", i, "_date", b)
    futuredate <- paste0("future_heat", i, "_date", b)
    
    newnamedate <- paste0("heat", i, "_affect_future_date", b)
    # newnamedur <- paste0("heat", i, "_affect_future_duration", b)
    newnameper <- paste0("heat", i, "_affect_future_tasmax", b)
    
    # pastdur <- paste0("past_heat", i, "_duration", b)
    # futuredur <- paste0("future_heat", i, "_duration", b)
    pastper <- paste0("past_heat", i, "_tasmax", b)
    futureper <- paste0("future_heat", i, "_tasmax", b)
    
    pastaffect <- paste0("heat", i, "_affect_past", b)
    
    bhps_heat.df[, newname] <- 0
    
    # Directly affected affected
    oo <- which(bhps_heat.df[, futuredate] <= bhps_heat.df$moveout_lsoa_date)
    bhps_heat.df[, newname][oo] <- 1
    
    # Also save date
    bhps_heat.df$tmp_date <- NA
    bhps_heat.df$tmp_date[oo] <- bhps_heat.df[oo, futuredate]
    
    # # Also save duration
    # bhps_heat.df$tmp_dur <- NA
    # bhps_heat.df$tmp_dur[oo] <- bhps_heat.df[oo, futuredur]
    
    # Also save percentage
    bhps_heat.df$tmp_per <- NA
    bhps_heat.df$tmp_per[oo] <- bhps_heat.df[oo, futureper]
    
    # Newly past affected
    bhps_heat.df$tmp <- as.Date(ave(as.numeric(bhps_heat.df$date),
                                     by = paste0(bhps_heat.df$pidp),
                                     FUN = function(x) dplyr::lag(x)), 
                                 origin = "1970-01-01")
    
    bhps_heat.df$tmp2 <- NA
    oo <- which(bhps_heat.df[, pastdate] <= bhps_heat.df$date & 
                  bhps_heat.df[, pastdate] > bhps_heat.df$tmp &
                  bhps_heat.df[, pastdate] >= bhps_heat.df$movein_lsoa_date)
    bhps_heat.df$tmp2[oo] <- 1
    
    
    # Lead and declare as future affected
    bhps_heat.df$tmp2 <- ave(bhps_heat.df$tmp2,
                              by = bhps_heat.df$pidp,
                              FUN = function(x) dplyr::lead(x))
    
    bhps_heat.df[, newname][bhps_heat.df$tmp2 == 1] <- 1
    bhps_heat.df$tmp <- NULL
    bhps_heat.df$tmp2 <- NULL
    
    # Also save date
    bhps_heat.df$tmp_date2 <- NA
    bhps_heat.df$tmp_date2[oo] <- bhps_heat.df[oo, pastdate]
    
    # Lead date 
    bhps_heat.df$tmp_date2 <- ave(bhps_heat.df$tmp_date2,
                                   by = bhps_heat.df$pidp,
                                   FUN = function(x) dplyr::lead(x))
    
    # # Also save dur
    # bhps_heat.df$tmp_dur2 <- NA
    # bhps_heat.df$tmp_dur2[oo] <- bhps_heat.df[oo, pastdur]
    # 
    # # Lead dur 
    # bhps_heat.df$tmp_dur2 <- ave(bhps_heat.df$tmp_dur2,
    #                               by = bhps_heat.df$pidp,
    #                               FUN = function(x) dplyr::lead(x))
    
    # Also save per
    bhps_heat.df$tmp_per2 <- NA
    bhps_heat.df$tmp_per2[oo] <- bhps_heat.df[oo, pastper]
    
    # Lead per 
    bhps_heat.df$tmp_per2 <- ave(bhps_heat.df$tmp_per2,
                                  by = bhps_heat.df$pidp,
                                  FUN = function(x) dplyr::lead(x))
    
    # Use more earlier date and duration
    oo <- which(bhps_heat.df$tmp_date2 < bhps_heat.df$tmp_date | is.na(bhps_heat.df$tmp_date))
    bhps_heat.df$tmp_date[oo] <- bhps_heat.df$tmp_date2[oo]
    # bhps_heat.df$tmp_dur[oo] <- bhps_heat.df$tmp_dur2[oo]
    bhps_heat.df$tmp_per[oo] <- bhps_heat.df$tmp_per2[oo]
    
    bhps_heat.df$tmp_date2  <- bhps_heat.df$tmp_per2 <- NA
    
    # Fill na
    bhps_heat.df[, newnamedate] <- as.Date(ave(bhps_heat.df$tmp_date,
                                                by = bhps_heat.df$pidp,
                                                FUN = function(x) repeat_last(x, forward = FALSE)),
                                            origin = "1970-01-01")
    # bhps_heat.df[, newnamedur] <- ave(bhps_heat.df$tmp_dur,
    #                                    by = bhps_heat.df$pidp,
    #                                    FUN = function(x) repeat_last(x, forward = FALSE))
    bhps_heat.df[, newnameper] <- ave(bhps_heat.df$tmp_per,
                                       by = bhps_heat.df$pidp,
                                       FUN = function(x) repeat_last(x, forward = FALSE))
    
    bhps_heat.df$tmp_date <- bhps_heat.df$tmp_per <- NA
    
    
    # Fill NAs before occasion
    bhps_heat.df[, newname] <- ave(bhps_heat.df[, newname],
                                    by = bhps_heat.df$pidp,
                                    FUN = function(x) rev(cummax(rev(x))))
    
    
    # Correct those occasions with an overlap in migration
    tmp <- ave(bhps_heat.df[, pastaffect],
               by = bhps_heat.df$pidp,
               FUN = function(x) dplyr::lead(x))
    oo <- which(tmp == 1 & bhps_heat.df[, newname] == 0 & bhps_heat.df[, pastaffect] == 0)
    bhps_heat.df[, newname][oo] <- 1
    bhps_heat.df[, newname] <- ave(bhps_heat.df[, newname],
                                    by = bhps_heat.df$pidp,
                                    FUN = function(x) rev(cummax(rev(x))))
    
    # Also correct date and duration and percentage
    bhps_heat.df$tmp_date <- NA
    bhps_heat.df$tmp_date[oo] <- bhps_heat.df[oo + 1, pastdate]
    bhps_heat.df[!is.na(bhps_heat.df$tmp_date), futuredate] <- as.Date(bhps_heat.df$tmp_date[!is.na(bhps_heat.df$tmp_date)],
                                                                         origin = "1970-01-01")
    bhps_heat.df[, futuredate] <- as.Date(ave(as.numeric(bhps_heat.df[, futuredate]),
                                               by = bhps_heat.df$pidp,
                                               FUN = function(x) repeat_last(x, forward = FALSE)), 
                                           origin = "1970-01-01")
    
    # bhps_heat.df$tmp_dur <- NA
    # bhps_heat.df$tmp_dur[oo] <- bhps_heat.df[oo + 1, pastdur]
    # bhps_heat.df[!is.na(bhps_heat.df$tmp_dur), futuredur] <- bhps_heat.df$tmp_dur[!is.na(bhps_heat.df$tmp_dur)]
    # bhps_heat.df[, futuredur] <- ave(as.numeric(bhps_heat.df[, futuredur]),
    #                                            by = bhps_heat.df$pidp,
    #                                            FUN = function(x) repeat_last(x, forward = FALSE))
    
    bhps_heat.df$tmp_per <- NA
    bhps_heat.df$tmp_per[oo] <- bhps_heat.df[oo + 1, pastper]
    bhps_heat.df[!is.na(bhps_heat.df$tmp_per), futureper] <- bhps_heat.df$tmp_per[!is.na(bhps_heat.df$tmp_per)]
    bhps_heat.df[, futureper] <- ave(as.numeric(bhps_heat.df[, futureper]),
                                      by = bhps_heat.df$pidp,
                                      FUN = function(x) repeat_last(x, forward = FALSE))
    
    
    bhps_heat.df$tmp_date <- bhps_heat.df$tmp_per <- NULL
    
  }
}



### Summary Tables 

table(bhps_heat.df$heat_affect_past)
table(bhps_heat.df$heat_affect_future)





################################
### Ever affected by a heat ###
################################

buf <- c(1, 2, 3) # For each buffer
intens <- c("") # and each intensity (percentage)

for(b in buf){
  for(i in intens){
    pname <- paste0("heat", i, "_affect_past", b)
    fname <- paste0("heat", i, "_affect_future", b)
    newname <-  paste0("heat", i, "_affect_ever", b)
    
    bhps_heat.df[, newname] <- ave(bhps_heat.df[, pname] + bhps_heat.df[, fname],
                                    by = bhps_heat.df$pidp,
                                    FUN = function(x) max(x))
    bhps_heat.df[, newname][bhps_heat.df[, newname] > 1] <- 1
  }
}



### Test

nrow(bhps_heat.df[bhps_heat.df$heat_affect_ever == 1, ])
nrow(bhps_heat.df[bhps_heat.df$heat_affect_past == 1, ])
nrow(bhps_heat.df[bhps_heat.df$heat_affect_future == 1 & bhps_heat.df$heat_affect_past == 0, ])




############################################
### Temporal distance to affected heats ###
############################################

bhps_heat.df <- bhps_heat.df[order(bhps_heat.df$pidp, bhps_heat.df$wave), ]


### As above, first identify which heats are actually affecting people
buf <- c(1, 2, 3) # For each buffer
intens <- c("") # and each intensity (percentage)

for(b in buf){
  for(i in intens){
    
    distname <- paste0("tempdist", i, "_affect_past", b)
    distnamef <- paste0("tempdist", i, "_affect_future", b)
    
    pastdate <- paste0("heat", i, "_affect_past_date", b)
    futuredate <- paste0("heat", i, "_affect_future_date", b)
    
    # Temporal distance to interview
    bhps_heat.df[, distname] <- as.numeric(bhps_heat.df$date) - as.numeric(bhps_heat.df[, pastdate])
    bhps_heat.df[, distnamef] <- as.numeric(bhps_heat.df$date) - as.numeric(bhps_heat.df[, futuredate])
  }
}

summary(bhps_heat.df[, paste0("tempdist", rep(intens, length(buf)),
                               "_affect_past", rep(buf, each = length(intens)))])

summary(bhps_heat.df[, paste0("tempdist", rep(intens, length(buf)),
                               "_affect_future", rep(buf, each = length(intens)))])




###################################################################################
### For those living in an lsoa that was ever heat, feed across entire waves    ###
### Indicator of ever lived in an ever heat lsoa                             ###
###################################################################################
buf <- c(1, 2, 3) # For each buffer
intens <- c("") # and each intensity (percentage)

for(b in buf){
  for(i in intens){
    
    n_old <- paste0("heatwave", i, "_ever", b)
    n1 <- paste0("heat", i, "_ever", b)
    
    # Rename
    names(bhps_heat.df)[which(names(bhps_heat.df) == n_old)] <- n1
    
    # Distribute
    bhps_heat.df[, n1] <- ave(bhps_heat.df[, n1],
                               bhps_heat.df$pidp,
                               FUN = function(x) max(x, na.rm = TRUE))
  }
}







####### ---------------------------------------------------------------------- ######
#######         Environmental variables from BHPS US                           ######
####### ---------------------------------------------------------------------- ######


################################
### Attitudes and Intentions ###
################################

### People will be affected by climate change

bhps_heat.df$climatechange30 <- 2 - bhps_heat.df$scopecl30
table(bhps_heat.df$climatechange30)

bhps_heat.df$climatechange200 <- 2 - bhps_heat.df$scopecl200
table(bhps_heat.df$climatechange200)


### Climate change cause heats in UK

bhps_heat.df$climatechange_heats <- 2 - bhps_heat.df$opccb
table(bhps_heat.df$climatechange_heats)


### 1) Climate change influence by behaviour

bhps_heat.df$ccbehavcontributes <- 6 - bhps_heat.df$scenv_bccc
table(bhps_heat.df$ccbehavcontributes)

bhps_heat.df$ccbehavcontributes_w1 <- bhps_heat.df$scenv_ccls - 1 # other way round?!
table(bhps_heat.df$ccbehavcontributes_w1)


### 2) Pay more for env friendly products

bhps_heat.df$ccpaymore <- 6 - bhps_heat.df$scenv_pmep
table(bhps_heat.df$ccpaymore)

bhps_heat.df$ccpaymore_w1 <- 2 - bhps_heat.df$scenv_pmre
table(bhps_heat.df$ccpaymore_w1)


### 3) Soon experience disaster

bhps_heat.df$ccsoondisaster <- 6 - bhps_heat.df$scenv_meds
table(bhps_heat.df$ccsoondisaster)

bhps_heat.df$ccsoondisaster_w1 <- 2 - bhps_heat.df$scenv_dstr
table(bhps_heat.df$ccsoondisaster_w1)


### 4) Climate change been axaggerated

bhps_heat.df$ccexag <- 6 - bhps_heat.df$scenv_crex
table(bhps_heat.df$ccexag)

bhps_heat.df$ccexag_w1 <- 2 - bhps_heat.df$scenv_exag
table(bhps_heat.df$ccexag_w1)


### 5) Climate change beyond control

bhps_heat.df$ccbeyondcontrol <- 6 - bhps_heat.df$scenv_tlat
table(bhps_heat.df$ccbeyondcontrol)

bhps_heat.df$ccbeyondcontrol_w1 <- 2 - bhps_heat.df$scenv_bcon
table(bhps_heat.df$ccbeyondcontrol_w1)


### 6) Climate change to far away to worry

bhps_heat.df$ccnoworry <- 6 - bhps_heat.df$scenv_nowo
table(bhps_heat.df$ccnoworry)

bhps_heat.df$ccnoworry_w1 <- 2 - bhps_heat.df$scenv_futr
table(bhps_heat.df$ccnoworry_w1)


### 7) Changes need to fit lifestyle

bhps_heat.df$ccchangesfitlifestyle <- 6 - bhps_heat.df$scenv_fitl
table(bhps_heat.df$ccchangesfitlifestyle)

bhps_heat.df$ccchangesfitlifestyle_w1 <- 2 - bhps_heat.df$scenv_cfit
table(bhps_heat.df$ccchangesfitlifestyle_w1)


### 8) not worth contributing if others dont

bhps_heat.df$ccnotworthifotherdont <-  6 - bhps_heat.df$scenv_noot
table(bhps_heat.df$ccnotworthifotherdont)

bhps_heat.df$ccnotworthifotherdont_w1 <- 2 - bhps_heat.df$scenv_chwo
table(bhps_heat.df$ccnotworthifotherdont_w1)


### 9) Not worth britain doing anything

bhps_heat.df$ccnotworthbritain <- 6 - bhps_heat.df$scenv_canc
table(bhps_heat.df$ccnotworthbritain)

bhps_heat.df$ccnotworthbritain_w1 <- 2 - bhps_heat.df$scenv_brit
table(bhps_heat.df$ccnotworthbritain_w1)









##################################
### Vote intention green party ###
##################################

table(bhps_heat.df$vote3, useNA = "always")
table(bhps_heat.df$vote4, useNA = "always")

### Combine clossness and vote intention (as variable vote in BHPS)

table(bhps_heat.df$vote3, bhps_heat.df$vote4, useNA = "always")
# Should have receive only one, but not for all cases true


# Set to NA, if different party preferences
oo <- which(!is.na(bhps_heat.df$vote3) & !is.na(bhps_heat.df$vote4) &
              bhps_heat.df$vote3 != bhps_heat.df$vote4)
bhps_heat.df$tmp1 <- bhps_heat.df$vote3
bhps_heat.df$tmp1[oo] <- NA
bhps_heat.df$tmp2 <- bhps_heat.df$vote4
bhps_heat.df$tmp2[oo] <- NA

# Comine to one variable
bhps_heat.df$vote_comb <- rowMeans(bhps_heat.df[, c("tmp1", "tmp2")], na.rm = TRUE)
bhps_heat.df$vote_comb[is.nan(bhps_heat.df$vote_comb)] <- NA
table(bhps_heat.df$vote_comb, useNA = "always")

bhps_heat.df[, c("tmp1", "tmp2")] <- NULL

# Recode to green vote preference
bhps_heat.df$vote_green <- bhps_heat.df$vote_comb
bhps_heat.df$vote_green[which(bhps_heat.df$vote_green %in% c(1:5, 7:13))] <- 0
bhps_heat.df$vote_green[which(bhps_heat.df$vote_green %in% c(95, 97))] <- 0 # non-voters and other party
bhps_heat.df$vote_green[which(bhps_heat.df$vote_green %in% c(96))] <- NA # cannot vote to NA
bhps_heat.df$vote_green[which(bhps_heat.df$vote_green %in% c(6))] <- 1 # cannot vote to NA

table(bhps_heat.df$vote_green)


# Recode to not vote
bhps_heat.df$vote_none <- bhps_heat.df$vote_comb
bhps_heat.df$vote_none[which(bhps_heat.df$vote_none %in% c(1:13))] <- 0
bhps_heat.df$vote_none[which(bhps_heat.df$vote_none %in% c(97))] <- 0 # other party
bhps_heat.df$vote_none[which(bhps_heat.df$vote_none %in% c(96))] <- NA # cannot vote to NA
bhps_heat.df$vote_none[which(bhps_heat.df$vote_none %in% c(95))] <- 1 # cannot vote to NA

table(bhps_heat.df$vote_none)




###########################
### Opinion on politics ###
###########################

table(bhps_heat.df$poleff3)
bhps_heat.df$pol_nocare <- 6 - bhps_heat.df$poleff3
table(bhps_heat.df$pol_nocare)

table(bhps_heat.df$poleff4)
bhps_heat.df$pol_nosay <- 6 - bhps_heat.df$poleff4
table(bhps_heat.df$pol_nosay)

# Index of fustration
bhps_heat.df$pol_frust <- (bhps_heat.df$pol_nocare + bhps_heat.df$pol_nosay)/2
table(bhps_heat.df$pol_frust)




##############################
### Environmental behavior ###
##############################

vars <- paste0("envhabit", c(1:11))

### code 6 as missing

for(i in vars){
  bhps_heat.df[which(bhps_heat.df[, i] == 6), i] <- NA
  print(table(bhps_heat.df[, i]))
}


### Recode that positive is more environ friendly 2 4 5 6 7 8 9 10 11
c <- c(2, 4, 5, 6, 7, 8, 9, 10, 11)

for(i in c){
  var <- paste0("envhabit", i)
  bhps_heat.df[, var] <- 6 - bhps_heat.df[, var]
  print(table(bhps_heat.df[, var]))
}



### Principal components
tmp.df <- na.omit(bhps_heat.df[, vars]) # Only non missing values

pca.env <- princomp(tmp.df, cor=TRUE)
summary(pca.env)

pca.env <- psych::principal(tmp.df, nfactors = 3, rotate = "varimax")
print(pca.env)




### General behavioral index
bhps_heat.df$env_index_home <- (bhps_heat.df$envhabit1 +
                                         bhps_heat.df$envhabit2 +
                                         bhps_heat.df$envhabit3 +
                                         bhps_heat.df$envhabit4 +
                                         bhps_heat.df$envhabit7)/5

cronbach <- psych::alpha(bhps_heat.df[, paste0("envhabit", c(1:4, 7))])
print(cronbach)


bhps_heat.df$env_index_shop <- (bhps_heat.df$envhabit5 +
                                         bhps_heat.df$envhabit6)/2

bhps_heat.df$env_index_trans <- (bhps_heat.df$envhabit8 +
                                         bhps_heat.df$envhabit9)/2

# All available items
bhps_heat.df$env_index_all <- (bhps_heat.df$envhabit1 +
                                        bhps_heat.df$envhabit2 +
                                        bhps_heat.df$envhabit3 +
                                        bhps_heat.df$envhabit4 +
                                        bhps_heat.df$envhabit5 +
                                        bhps_heat.df$envhabit6 +
                                        bhps_heat.df$envhabit7 +
                                        bhps_heat.df$envhabit8 +
                                        bhps_heat.df$envhabit9 +
                                        bhps_heat.df$envhabit10 +
                                        bhps_heat.df$envhabit11)/11

cronbach <- psych::alpha(bhps_heat.df[, paste0("envhabit", c(1:11))])
print(cronbach)

# All items available in BHPS and UKHLS
bhps_heat.df$env_index_all2 <- (bhps_heat.df$envhabit1 +
                                        bhps_heat.df$envhabit2 +
                                        bhps_heat.df$envhabit3 +
                                        bhps_heat.df$envhabit4 +
                                        bhps_heat.df$envhabit5 +
                                        bhps_heat.df$envhabit6 +
                                        bhps_heat.df$envhabit7)/7

cronbach <- psych::alpha(bhps_heat.df[, paste0("envhabit", c(1:7))])
print(cronbach)



### Car use frequency

table(bhps_heat.df$trcarfq)
bhps_heat.df$car_freq <- 9 - bhps_heat.df$trcarfq
table(bhps_heat.df$car_freq)


### Cycling use frequency

table(bhps_heat.df$trbikefq)
bhps_heat.df$bike_freq <- 9 - bhps_heat.df$trbikefq
table(bhps_heat.df$bike_freq)



###########################################
### Prepare socio-demographic variables ###
###########################################

getmode <- function(v) {
  uniqv <- unique(v[!is.na(v)])
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

### Age
table(bhps_heat.df$age_dv)

# Age for bhps
oo <- which(bhps_heat.df$wave <= 18)
bhps_heat.df$age_dv[oo] <- bhps_heat.df$age[oo]

# Cut into bins
cuts <- seq(min(bhps_heat.df$age_dv, na.rm = T), max(bhps_heat.df$age_dv, na.rm = T), 5)
cuts <- c(cuts[-which(cuts > 90)],  max(bhps_heat.df$age_dv, na.rm = T))
bhps_heat.df$age_cat <- cut(bhps_heat.df$age_dv,
                             breaks = cuts)
table(bhps_heat.df$age_cat)

bhps_heat.df$age_cat <- relevel(as.factor(bhps_heat.df$age_cat),
                                 ref = as.numeric(getmode(bhps_heat.df$age_cat)))

### Gender
table(bhps_heat.df$sex)

bhps_heat.df$female <- bhps_heat.df$sex - 1

### Migration to UK
table(bhps_heat.df$bornuk_dv)
bhps_heat.df$migback <- bhps_heat.df$bornuk_dv - 1


### Migration background generation
table(bhps_heat.df$generation)
bhps_heat.df$migback_gen <- bhps_heat.df$generation
bhps_heat.df$migback_gen[bhps_heat.df$migback_gen %in% c(4:6)] <- 0
table(bhps_heat.df$migback_gen)


### Ethnic group
table(bhps_heat.df$ethn_dv)
bhps_heat.df$ethn_dv <- relevel(as.factor(bhps_heat.df$ethn_dv),
                                 ref = as.numeric(getmode(bhps_heat.df$ethn_dv)))


### Ethnic group short
bhps_heat.df$ethn_dv_short <- as.numeric(bhps_heat.df$ethn_dv)
bhps_heat.df$ethn_dv_short[bhps_heat.df$ethn_dv_short %in% c(1:4)] <- 0
bhps_heat.df$ethn_dv_short[bhps_heat.df$ethn_dv_short %in% c(5:8)] <- 1
bhps_heat.df$ethn_dv_short[bhps_heat.df$ethn_dv_short %in% c(9:13)] <- 2
bhps_heat.df$ethn_dv_short[bhps_heat.df$ethn_dv_short %in% c(14:16)] <- 3
bhps_heat.df$ethn_dv_short[bhps_heat.df$ethn_dv_short %in% c(17:97)] <- 4

bhps_heat.df$ethn_dv_short <- as.factor(bhps_heat.df$ethn_dv_short)

table(bhps_heat.df$ethn_dv_short)


### Education
table(bhps_heat.df$hiqual_dv)

bhps_heat.df$hiqual_dv <- relevel(as.factor(bhps_heat.df$hiqual_dv),
                                   ref = as.numeric(getmode(bhps_heat.df$hiqual_dv)))


### Child(ren) in household
table(bhps_heat.df$nkids_dv)
bhps_heat.df$child <- bhps_heat.df$nkids_dv
bhps_heat.df$child[bhps_heat.df$child > 0 ] <- 1
table(bhps_heat.df$child)

### Marital status
table(bhps_heat.df$marstat_dv)

bhps_heat.df$marstat_dv <- relevel(as.factor(bhps_heat.df$marstat_dv),
                                    ref = as.numeric(getmode(bhps_heat.df$marstat_dv)))

### Household total income
summary(bhps_heat.df$fimngrs_dv)
summary(bhps_heat.df$fimnlabgrs_dv)
summary(bhps_heat.df$fihhmngrs_dv)
hist(bhps_heat.df$fihhmngrs_dv)

bhps_heat.df$hhinc <- bhps_heat.df$fihhmngrs_dv/1000

summary(bhps_heat.df$hhinc)

# Drop outlier (below zero or above 100.000 per month )
oo <- which(bhps_heat.df$hhinc < 0 |
              bhps_heat.df$hhinc > 100)
bhps_heat.df$hhinc[oo] <- NA

# Squared income
bhps_heat.df$hhinc_sq <- bhps_heat.df$hhinc^2


### OECD equivalence income (square root method for simplicity)
bhps_heat.df$hhinc_eq <- bhps_heat.df$hhinc / sqrt(bhps_heat.df$hhsize)
bhps_heat.df$hhinc_eq_sq <- bhps_heat.df$hhinc_eq^2

summary(bhps_heat.df$hhinc_eq)


### Income deciles
breaks <- quantile(bhps_heat.df$hhinc[which(bhps_heat.df$wave %in% c(18, 19, 22, 28))], probs = seq(0, 1, by = 0.1), na.rm = TRUE)
bhps_heat.df$hhinc_dec <- cut(bhps_heat.df$hhinc, breaks = breaks, include.lowest = TRUE)


### Home owner
table(bhps_heat.df$tenure_dv, useNA = "always")

bhps_heat.df$owner <- bhps_heat.df$tenure_dv
bhps_heat.df$owner[bhps_heat.df$owner %in% c(1,2)] <- 1
bhps_heat.df$owner[bhps_heat.df$owner %in% c(3:8)] <- 0
table(bhps_heat.df$owner)


### Moving Preference
bhps_heat.df$lkmove <- bhps_heat.df$lkmove - 1



### Political party

# No party support as reference
bhps_heat.df$vote_comb <- relevel(as.factor(bhps_heat.df$vote_comb),
                                   ref = "95")



### Satisfaction with income

# Variable for bhps and US
oo <- which(bhps_heat.df$wave <= 18)
bhps_heat.df$sat_inc <- ifelse(bhps_heat.df$wave <= 18, 
                                bhps_heat.df$lfsat2,
                                bhps_heat.df$sclfsat2)

# Make ordinal                                    
table(bhps_heat.df$sat_inc)
bhps_heat.df$sat_inc[which(bhps_heat.df$sat_inc == 0)] <- NA
bhps_heat.df$sat_inc <- as.factor(bhps_heat.df$sat_inc)



### 4 seasons
bhps_heat.df$season <- NA
bhps_heat.df$season[bhps_heat.df$istrtdatm %in% c(3:5)] <- 1
bhps_heat.df$season[bhps_heat.df$istrtdatm %in% c(6:8)] <- 2
bhps_heat.df$season[bhps_heat.df$istrtdatm %in% c(9:11)] <- 3
bhps_heat.df$season[bhps_heat.df$istrtdatm %in% c(12, 1, 2)] <- 4
table(bhps_heat.df$season)


### Rolling season index

# Year -> previous year if month is Jan or Feb (winter season)
bhps_heat.df$tmp <- bhps_heat.df$istrtdaty
oo <- which(bhps_heat.df$istrtdatm %in% c(1,2))
bhps_heat.df$tmp[oo] <- bhps_heat.df$tmp[oo] - 1

# Gen index based on transformed year and season
bhps_heat.df$season_index <- paste(bhps_heat.df$tmp, bhps_heat.df$season, sep = "_")
table(bhps_heat.df$season_index)
table(bhps_heat.df$season_index[bhps_heat.df$istrtdaty >=2008], bhps_heat.df$istrtdaty[bhps_heat.df$istrtdaty >=2008])

bhps_heat.df$tmp <- NULL



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
bhps_heat.df <- bhps_heat.df[order(bhps_heat.df$pidp, bhps_heat.df$wave),]

### Use previous wave to impute missing wave 10 vote_comp 
oo <- which(bhps_heat.df$wave == 27 | bhps_heat.df$wave == 28)

# Feed forward
bhps_heat.df$vote_comb_imp <- bhps_heat.df$vote_comb

bhps_heat.df$vote_comb_imp[] <- ave(bhps_heat.df$vote_comb_imp[],
                                    bhps_heat.df$pidp[],
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
table(bhps_heat.df$vote4, useNA = "ifany")
bhps_heat.df$partisan <- bhps_heat.df$vote4
# Only >= fairly strong support
oo <- which(bhps_heat.df$vote5 == 3)
bhps_heat.df$partisan[oo] <- 99 # no partisanship
# if answered 2 in vote1 and vote2 code as no partisanship
oo <- which(bhps_heat.df$vote1 == 2 & bhps_heat.df$vote2 == 2)
bhps_heat.df$partisan[oo] <- 99 # no partisanship
table(bhps_heat.df$partisan, useNA = "ifany")
bhps_heat.df$partisan <- factor(bhps_heat.df$partisan)


### Problem solving
bhps_heat.df$problemsolv <- (bhps_heat.df$scghqd + bhps_heat.df$scghqh + bhps_heat.df$scghqf)/3



### Trust vars
table(bhps_heat.df$sctrust, useNA = "ifany")

bhps_heat.df$trust <- bhps_heat.df$sctrust
bhps_heat.df$trust[which(bhps_heat.df$trust == 2 | bhps_heat.df$trust == 3)] <- 0
table(bhps_heat.df$trust)

# Trust measure 2 (attention, in UKHLS 11 points and BHPS 10 points. Use above 5 as trust)
table(bhps_heat.df$riskb, useNA = "ifany")
table(bhps_heat.df$scriskb, useNA = "ifany")

bhps_heat.df$trust2 <- bhps_heat.df$scriskb
bhps_heat.df$trust2[which(bhps_heat.df$wave == 18)] <- bhps_heat.df$riskb[which(bhps_heat.df$wave == 18)]
bhps_heat.df$trust3 <- bhps_heat.df$trust2
bhps_heat.df$trust2[which(bhps_heat.df$trust2 <= 5)] <- 0
bhps_heat.df$trust2[which(bhps_heat.df$trust2 > 5)] <- 1
table(bhps_heat.df$trust2)


### Cooperation / volunteering
bhps_heat.df$cooperation <- bhps_heat.df$ivcoop
bhps_heat.df$cooperation[which(bhps_heat.df$wave <= 18)] <- bhps_heat.df$iv4[which(bhps_heat.df$wave <= 18)]
table(bhps_heat.df$wave[!is.na(bhps_heat.df$cooperation)])
bhps_heat.df$cooperation <- 6 - bhps_heat.df$cooperation

# Volunteering
bhps_heat.df$voluntary_work <- ifelse(bhps_heat.df$volun == 1, 1, 0)
table(bhps_heat.df$wave[!is.na(bhps_heat.df$voluntary_work)])

bhps_heat.df$voluntary_work_freq <- 10 - bhps_heat.df$volfreq
bhps_heat.df$voluntary_work_freq[bhps_heat.df$voluntary_work_freq %in% c(1:3)] <- 1
bhps_heat.df$voluntary_work_freq[bhps_heat.df$voluntary_work_freq %in% c(4:6)] <- 2
bhps_heat.df$voluntary_work_freq[bhps_heat.df$voluntary_work_freq %in% c(7:9)] <- 3
table(bhps_heat.df$wave[!is.na(bhps_heat.df$voluntary_work_freq)])

oo <- which(bhps_heat.df$voluntary_work == 0 & bhps_heat.df$wave %in% c(20, 22, 24, 26, 28))
bhps_heat.df$voluntary_work_freq[oo] <- 0
table(bhps_heat.df$voluntary_work_freq)

# Charity spending
bhps_heat.df$charity <- ifelse(bhps_heat.df$chargv == 1, 1, 0)
table(bhps_heat.df$charity)
table(bhps_heat.df$wave[!is.na(bhps_heat.df$charity)])

bhps_heat.df$charity_freq <- bhps_heat.df$charfreq
bhps_heat.df$charity_freq[bhps_heat.df$charity == 0] <- 0
table(bhps_heat.df$charity_freq)



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
bhps_heat.df <- bhps_heat.df[order(bhps_heat.df$pidp, bhps_heat.df$wave),]

for(v in vars){
  nv <- paste0("c_", v)
  nv2 <- paste0("a_", v)
  
  
  oo <- which(bhps_heat.df$wave == 18 | bhps_heat.df$wave == 19)
  bhps_heat.df[oo, nv] <-  bhps_heat.df[oo, v] 
  
  if(is.factor(bhps_heat.df[, v])){
    lev <- levels(droplevels(bhps_heat.df[, nv]))
    bhps_heat.df[, nv] <- as.numeric(bhps_heat.df[, nv])
  }
  
  bhps_heat.df[, nv] <- ave(bhps_heat.df[, nv], bhps_heat.df$pidp, 
                            FUN = function(x) repeat_last(x))
  
  if(is.factor(bhps_heat.df[, v])){
    bhps_heat.df[, nv] <- factor(bhps_heat.df[, nv], labels = lev)
  }
  
  ### Add average
  if(!is.factor(bhps_heat.df[, v])){
    bhps_heat.df[, nv2] <- ave(bhps_heat.df[, v], bhps_heat.df$pidp, 
                               FUN = function(x) mean(x, na.rm = TRUE))
  }
  
}


### For voluntary work and charity use second wave data
vars <- c("voluntary_work", "voluntary_work_freq", "charity", "charity_freq", "pol_nocare", "pol_nosay")
for(v in vars){
  nv <- paste0("c_", v)
  
  oo <- which(bhps_heat.df$wave == 20 | bhps_heat.df$wave == 21)
  bhps_heat.df[oo, nv] <-  bhps_heat.df[oo, v] 
  
  if(is.factor(bhps_heat.df[, v])){
    lev <- levels(droplevels(bhps_heat.df[, nv]))
    bhps_heat.df[, nv] <- as.numeric(bhps_heat.df[, nv])
  }
  
  bhps_heat.df[, nv] <- ave(bhps_heat.df[, nv], bhps_heat.df$pidp, 
                            FUN = function(x) mean(x, na.rm = TRUE))
  
  if(is.factor(bhps_heat.df[, v])){
    bhps_heat.df[, nv] <- factor(bhps_heat.df[, nv], labels = lev)
  }
}




### Self-efficacy (needs backward imputation)

v <- paste0("se", c(1:10))
for(i in v){
  print(table(bhps_heat.df[, i]))
}

bhps_heat.df$selfeff <- rowSums(bhps_heat.df[,v])/10
table(bhps_heat.df$selfeff)

# Feed backward from wave 5
bhps_heat.df$c_selfeff <- ave(bhps_heat.df$selfeff, bhps_heat.df$pidp, 
                              FUN = function(x) repeat_last(x, forward = FALSE))







############
### Save ###
############


save(bhps_heat.df, file = "bhpsukhls_heat_final.RData")





