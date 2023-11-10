
#### Analyses ####
#### Tobias Ruettenauer ####
#### 2020/ 01 / 31 ####

rm(list=ls())

### Load packages
#library(rgdal)
#library(spdep)
#library(rgeos)
library(doParallel)
library(foreign)
#library(GISTools)
#library(cleangeo)
library(ggplot2)
library(gridExtra)
library(extrafont)
loadfonts()
library(plm)
library(lfe)
library(sandwich)
library(lmtest)
library(texreg)



#### Extract methods for felm (from github texreg)
extract.felm <- function(model,
                         include.nobs = TRUE,
                         include.rsquared = TRUE,
                         include.adjrs = TRUE,
                         include.fstatistic = FALSE,
                         include.proj.stats = FALSE,
                         include.groups = TRUE,
                         ...) {
  
  s <- summary(model, ...)
  nam <- rownames(s$coefficients)
  co <- s$coefficients[, 1]
  se <- s$coefficients[, 2]
  pval <- s$coefficients[, 4]
  
  gof <- numeric()
  gof.names <- character()
  gof.decimal <- logical()
  if (include.rsquared == TRUE) {
    gof <- c(gof, s$r2)
    gof.decimal <- c(gof.decimal, TRUE)
    if (include.proj.stats == TRUE) {
      gof <- c(gof, s$P.r.squared)
      gof.decimal <- c(gof.decimal, TRUE)
      gof.names <- c(gof.names, "R$^2$ (full model)", "R$^2$ (proj model)")
    } else {
      gof.names <- c(gof.names, "R$^2$")
    }
  }
  if (include.adjrs == TRUE) {
    gof <- c(gof, s$r2adj)
    gof.decimal <- c(gof.decimal, TRUE)
    if (include.proj.stats == TRUE) {
      gof <- c(gof, s$P.adj.r.squared)
      gof.decimal <- c(gof.decimal, TRUE)
      gof.names <- c(gof.names, "Adj. R$^2$ (full model)", "Adj. R$^2$ (proj model)")
    } else {
      gof.names <- c(gof.names, "Adj. R$^2$")
    }
  }
  if (include.fstatistic == TRUE) {
    gof <- c(gof, s$F.fstat[1], s$F.fstat[4])
    gof.decimal <- c(gof.decimal, TRUE, TRUE)
    if (include.proj.stats == TRUE) {
      gof <- c(gof, s$P.fstat[length(s$P.fstat) - 1], s$P.fstat[1])
      gof.decimal <- c(gof.decimal, TRUE, TRUE)
      gof.names <- c(gof.names, "F statistic (full model)",
                     "F (full model): p-value", "F statistic (proj model)",
                     "F (proj model): p-value")
    } else {
      gof.names <- c(gof.names, "F statistic", "F p-value")
    }
  }
  if (include.nobs == TRUE) {
    gof <- c(gof, s$N)
    gof.names <- c(gof.names, "Num. obs.")
    gof.decimal <- c(gof.decimal, FALSE)
  }
  if (include.groups == TRUE && length(s$fe) > 0) {
    grp <- sapply(s$fe, function(x) length(levels(x)))
    grp.names <- paste0("Num. groups: ", names(grp))
    gof <- c(gof, grp)
    gof.names <- c(gof.names, grp.names)
    gof.decimal <- c(gof.decimal, rep(FALSE, length(grp)))
  }
  
  tr <- createTexreg(
    coef.names = nam,
    coef = co,
    se = se,
    pvalues = pval,
    gof.names = gof.names,
    gof = gof,
    gof.decimal = gof.decimal
  )
  return(tr)
}
setMethod("extract", signature = className("felm", "lfe"),
          definition = extract.felm)







### Working Directory
setwd("C:/work/Forschung/Climate Change_Replication/02_Data")


#################
### Load data ###
#################


load("bhpsukhls_heat_final.RData")




### time by region interaction
bhps_heat.df$time_gov <- factor(paste0(bhps_heat.df$istrtdaty, "_", bhps_heat.df$gor_dv))




###---------------------------------###
###       Prepare variables         ###
###---------------------------------###


### heat affected within the last two years

intens <- c("")
buf <- c(1, 2, 3)
for(i in intens){
  for(b in buf){
    v1 <- paste0("heat", i, "_affect_past", b)
    t1 <- paste0("tempdist", i, "_affect_past", b)
    nv <- paste0("heat", 12, "_affect1y_past", b)
    nv1 <- paste0("heat", 6, "_affect1y_past", b)
    nv2 <- paste0("heat", 1, "_affect1y_past", b)
    
    
    
    oo <- which(bhps_heat.df[, t1] <= (365 *  0.25))
    bhps_heat.df[, nv] <- 0
    bhps_heat.df[oo, nv] <- 1
    
    oo <- which(bhps_heat.df[, t1] <= 30)
    bhps_heat.df[, nv1] <- 0
    bhps_heat.df[oo, nv1] <- 1
    
    oo <- which(bhps_heat.df[, t1] <= 14)
    bhps_heat.df[, nv2] <- 0
    bhps_heat.df[oo, nv2] <- 1
    

  }
}






###-----------------------------------###
###       Reduce to final sample      ###
###-----------------------------------###

depvar <- c("climatechange30", "env_index_all2")
vars <- c("heat1_affect1y_past1", "age_cat", "female", "migback", "ethn_dv_short", "hiqual_dv", "child", "marstat_dv", "hhinc_dec",  
                          "istrtdaty", "istrtdatm")

for(d in depvar){
  cc <- complete.cases(bhps_heat.df[, c(d, vars)])
  nc <- ave(cc, bhps_heat.df$pidp, FUN = function(x) length(x[x == TRUE]))
  
  bhps_heat.df[, paste0("sample_", d)] <- 0
  oo <- which(cc == TRUE)
  bhps_heat.df[oo, paste0("sample_", d)] <- 1
  
  bhps_heat.df[, paste0("samplefe_", d)] <- 0
  oo <- which(cc == 1 & nc >= 2)
  bhps_heat.df[oo, paste0("samplefe_", d)] <- 1
}





###-----------------------------------###
###       Linear Models 500m          ###
###-----------------------------------###



### Standard lm + controls ###

depvar <- c("climatechange30", "env_index_all2")
intens <- c("12", "6", "1")
buffer <- c(1, 2, 3)
b <- buffer[1]
lm.l <- vector("list", nrow(expand.grid(intens)) * 2)

cl <- makeCluster(2)
registerDoParallel(cl)

res1 <- foreach(d = depvar) %dopar% {
  library(lfe)
  names(lm.l) <- paste0(d, c(1:length(lm.l)))
  j <- 1
  for(i in intens){
    
    v1 <- paste0("heat", i, "_affect1y_past", b)
    ev1 <- paste0("heat_ever", b)
    
    fm1 <- paste0(d, " ~ ", v1)
    fm2 <- paste0(d, " ~ ", v1, 
                  "+ age_cat + female + migback + ethn_dv_short + hiqual_dv + child + marstat_dv + hhinc_dec")
    
    lm1 <- felm(formula(paste0(fm1, " | season_index | 0 | pidp + lsoa01")),
                data = bhps_heat.df[bhps_heat.df[, paste0("sample_", d)] == 1, ],
                psdef = FALSE)
    
    lm2 <- felm(formula(paste0(fm2, " | season_index | 0 | pidp + lsoa01")),
                data = bhps_heat.df[bhps_heat.df[, paste0("sample_", d)] == 1, ],
                psdef = FALSE)
    
    
    
    lm.l[[j]] <- lm1
    lm.l[[j+1]] <- lm2
    j <- j + 2
  }
  res <- lm.l
  res
}
stopCluster(cl)


# Map
ns <- as.list(paste0("heat", intens, "_affect1y_past", b))
names(ns) <- paste0("heat", intens, "_affect1y_past", b)

# Climate change
screenreg(res1[[1]], custom.coef.map = ns, digits = 3)
# Envir behav
screenreg(res1[[2]], custom.coef.map = ns, digits = 3)






###----------------------------------###
###       Panel Models 500 m         ###
###----------------------------------###


### Create categorical indicator of direction in change

intens <- c("12", "6", "1")
buffer <- c(1, 2, 3)
for(i in intens){
  for(b in buf){
    
    v1 <- paste0("heat", i, "_affect1y_past", b)
    nv <- paste0("pattern", i, "_affect1y_past", b)
    
    oo <- which(!(bhps_heat.df$wave == 18 | bhps_heat.df$wave == 19))
    bhps_heat.df$tmp1 <- bhps_heat.df[, v1]
    bhps_heat.df[oo, "tmp1"] <- NA
    bhps_heat.df[, "tmp1"] <- ave(bhps_heat.df[, "tmp1"],
                                  bhps_heat.df$pidp,
                                  FUN = function(x) mean(x, na.rm = TRUE)) 
    oo <- which(!bhps_heat.df$wave == 22)
    bhps_heat.df$tmp2 <- bhps_heat.df[, v1]
    bhps_heat.df[oo, "tmp2"] <- NA
    bhps_heat.df[, "tmp2"] <- ave(bhps_heat.df[, "tmp2"],
                                  bhps_heat.df$pidp,
                                  FUN = function(x) mean(x, na.rm = TRUE))
    
    bhps_heat.df[, nv] <- paste(bhps_heat.df$tmp1, bhps_heat.df$tmp2, sep = " - ")
    oo <- grep("NaN", bhps_heat.df[, nv])
    bhps_heat.df[oo, nv] <- NA
    bhps_heat.df[, nv] <- as.factor(bhps_heat.df[, nv])
  }
}



### Standard plm + controls ###

depvar <- c("climatechange30", "env_index_all2")
intens <- c("12", "6", "1")
buffer <- c(1, 2, 3)
b <- buffer[1]
plm.l <- vector("list", nrow(expand.grid(intens)) * 2)

cl <- makeCluster(2)
registerDoParallel(cl)

res2 <- foreach(d = depvar) %dopar% {
  library(lfe)
  intens <- c("12", "6", "1")
  names(plm.l) <- paste0(d, c(1:length(plm.l)))
  j <- 1
  for(i in intens){
    
    v1 <- paste0("heat", i, "_affect1y_past", b)
    ev1 <- paste0("heat_ever", b)
    
    fm1 <- paste0(d, " ~ ", v1)
    fm2 <- paste0(d, " ~ ", v1, 
                  "+ age_cat + hiqual_dv + child + marstat_dv + hhinc_dec")
    
    # Estimate via felm with twoway clustered SEs
    bhps_heat.df$interaction <- as.numeric(as.factor(bhps_heat.df$season_index)) * bhps_heat.df[, ev1]
    
    plm1 <- felm(formula(paste0(fm1, " | pidp + season_index + interaction | 0 | pidp + lsoa01")),
                 data = bhps_heat.df[bhps_heat.df[, paste0("samplefe_", d)] == 1, ],
                 psdef = FALSE)
    
    plm2 <- felm(formula(paste0(fm2, " | pidp + season_index + interaction | 0 | pidp + lsoa01")),
                 data = bhps_heat.df[bhps_heat.df[, paste0("samplefe_", d)] == 1, ],
                 psdef = FALSE)
    
    plm.l[[j]] <- plm1
    plm.l[[j+1]] <- plm2
    j <- j + 2
  }
  res <- plm.l
  res
}
stopCluster(cl)


# Map
ns <- as.list(paste0("heat", intens, "_affect1y_past", b))
names(ns) <- paste0("heat", intens, "_affect1y_past", b)

# ns2 <- as.list(paste0("heat", intens, "_affectwo1y_past", b))
# names(ns2) <- paste0("heat", intens, "_affectwo1y_past", b)

# Climate change
screenreg(res2[[1]], custom.coef.map = ns, digits = 3)
# Envir behav
screenreg(res2[[2]], custom.coef.map = ns, digits = 3)






###----------------------------###
###       Export Tables        ###
###----------------------------###

#############
#### OLS ####
#############


### Attitudes

tex_t1 <- texreg(l = res1[[1]],
                 #file="../03_Output/Mod_fe_siting.tex",
                 digits = 3, leading.zero = TRUE,
                 stars = c(0.001, 0.01, 0.05, 0.1),
                 symbol = "\\dagger",
                 label = "tab:lm_belief_heat",
                 caption = "Linear OLS model. Dep. var.: Belief in climate change.",
                 custom.model.names = c("(1)", "(2)", "(3)", "(4)", "(5)", "(6)"),
                 #groups=list("Census cell level"=1:6, "Community level"=7:nb),
                 custom.coef.map = list('heat12_affect1y_past1' =  'Heatwave affected',
                                        'heat6_affect1y_past1' =  'Heatwave affected',
                                        'heat1_affect1y_past1' =  'Heatwave affected'),
                 custom.note = "%stars. Two-sided test. Cluster robust standard errors in parentheses (clustered by person and LSOA).",
                 dcolumn = TRUE, caption.above = TRUE, use.packages = FALSE, include.rmse = FALSE
)


# Customize
tex_t1 <- gsub("D[{].[}][{].[}][{][[:digit:]]+\\.*[[:digit:]]*[}]", "D{.}{.}{2.4}", tex_t1)
tex_t1 <- gsub("\\hline", "\\hline\\\\[-1.2ex]", tex_t1, tex_t1, fixed = TRUE)



head <- c("\\centering\n\\footnotesize\n{\\begin{threeparttable}")
tex_t1 <- gsub("\n\\begin{table}", head, tex_t1, fixed = TRUE)

head2 <- c("\\hline\\\\[-1.2ex] \n  & \\multicolumn{2}{c}{Heatwave 4 months}  & \\multicolumn{2}{c}{Heatwave 1 month}  & \\multicolumn{2}{c}{Heatwave 14 days}  \\\\ 
           \\cmidrule(lr){2-3} \\cmidrule(lr){4-5} \\cmidrule(lr){6-7}") 
tex_t1 <- sub("\\hline\\\\[-1.2ex]", head2, tex_t1, fixed = TRUE)

cont <- paste0("\\hline\\\\[-1.2ex]
  Basic controls & \\multicolumn{1}{c}{Yes} & \\multicolumn{1}{c}{Yes} & \\multicolumn{1}{c}{Yes} & \\multicolumn{1}{c}{Yes} & \\multicolumn{1}{c}{Yes} & \\multicolumn{1}{c}{Yes}   \\\\
  Additional controls & \\multicolumn{1}{c}{No} & \\multicolumn{1}{c}{Yes} & \\multicolumn{1}{c}{No} & \\multicolumn{1}{c}{Yes} & \\multicolumn{1}{c}{No} & \\multicolumn{1}{c}{Yes}  \\\\
 \\hline\\\\[-1.2ex]\n","R$^2$")
tex_t1 <- sub("\\hline\\\\[-1.2ex]\nR$^2$", cont, tex_t1, fixed = TRUE)



bottom <- paste0("\n\\\\hline\\\\\\\\[-1.2ex]\n ",
                 "\\\\end{tabular}\n ",
                 "\\\\begin{tablenotes}\n ",
                 "\\\\item \\\\scriptsize{$^{***}p<0.001$, $^{**}p<0.01$, $^*p<0.05$, $^{\\\\dagger}p<0.1$. Two-sided test. Cluster robust standard errors in parentheses (clustered by person and LSOA). Basic controls: year-season of interview (omitted). Additional controls: age (5 year intervals, omitted), highest education, child(ren), marital status, household income deciles.}\n",
                 "\\\\end{tablenotes}\n",
                 "\\\\label{tab:lm_belief_heat}\n",
                 "\\\\end{center}\n",
                 "\\\\end{threeparttable}\n",
                 "}\n")
l <- gregexpr("\\hline\\\\[-1.2ex]", tex_t1, fixed = TRUE)
l <- l[[1]][length(l[[1]])]
tex_t1 <- substr(tex_t1, 1, (l-2))
tex_t1 <- sub("$", bottom, tex_t1, fixed = FALSE)

# N into middle
l1 <- gregexpr("Num. obs", tex_t1, fixed = TRUE)
l1 <- l1[[1]][length(l1[[1]])]


tmp <- substr(tex_t1, l1, (l-2))
tmp2 <- gsub("([[:digit:]]+)", "\\\\multicolumn{1}{c}{\\1}", tmp)
tex_t1 <- gsub(tmp, tmp2, tex_t1, fixed = TRUE)

# Print
write.table(tex_t1, file = "../03_Output/Table1_heat_attitudes.tex",
            col.names = FALSE, row.names = FALSE, quote = FALSE)



### Behaviour

tex_t2 <- texreg(l = res1[[2]],
                 #file="../03_Output/Mod_fe_siting.tex",
                 digits = 3, leading.zero = TRUE,
                 stars = c(0.001, 0.01, 0.05, 0.1),
                 symbol = "\\dagger",
                 label = "tab:lm_behave_heat",
                 caption = "Linear OLS model. Dep. var.: Environmental behaviour.",
                 custom.model.names = c("(1)", "(2)", "(3)", "(4)", "(5)", "(6)"),
                 #groups=list("Census cell level"=1:6, "Community level"=7:nb),
                 custom.coef.map = list('heat12_affect1y_past1' =  'Heatwave affected',
                                        'heat6_affect1y_past1' =  'Heatwave affected',
                                        'heat1_affect1y_past1' =  'Heatwave affected'),
                 custom.note = "%stars. Two-sided test. Cluster robust standard errors in parentheses (clustered by person and LSOA).",
                 dcolumn = TRUE, caption.above = TRUE, use.packages = FALSE, include.rmse = FALSE
)


# Customize
tex_t2 <- gsub("D[{].[}][{].[}][{][[:digit:]]+\\.*[[:digit:]]*[}]", "D{.}{.}{2.4}", tex_t2)
tex_t2 <- gsub("\\hline", "\\hline\\\\[-1.2ex]", tex_t2, tex_t2, fixed = TRUE)



head <- c("\\centering\n\\footnotesize\n{\\begin{threeparttable}")
tex_t2 <- gsub("\n\\begin{table}", head, tex_t2, fixed = TRUE)

head2 <- c("\\hline\\\\[-1.2ex] \n  & \\multicolumn{2}{c}{Heatwave 4 months}  & \\multicolumn{2}{c}{Heatwave 1 month}  & \\multicolumn{2}{c}{Heatwave 14 days}  \\\\ 
           \\cmidrule(lr){2-3} \\cmidrule(lr){4-5} \\cmidrule(lr){6-7}") 
tex_t2 <- sub("\\hline\\\\[-1.2ex]", head2, tex_t2, fixed = TRUE)

cont <- paste0("\\hline\\\\[-1.2ex]
  Basic controls & \\multicolumn{1}{c}{Yes} & \\multicolumn{1}{c}{Yes} & \\multicolumn{1}{c}{Yes} & \\multicolumn{1}{c}{Yes} & \\multicolumn{1}{c}{Yes} & \\multicolumn{1}{c}{Yes}   \\\\
  Additional controls & \\multicolumn{1}{c}{No} & \\multicolumn{1}{c}{Yes} & \\multicolumn{1}{c}{No} & \\multicolumn{1}{c}{Yes} & \\multicolumn{1}{c}{No} & \\multicolumn{1}{c}{Yes}  \\\\
 \\hline\\\\[-1.2ex]\n","R$^2$")
tex_t2 <- sub("\\hline\\\\[-1.2ex]\nR$^2$", cont, tex_t2, fixed = TRUE)



bottom <- paste0("\n\\\\hline\\\\\\\\[-1.2ex]\n ",
                 "\\\\end{tabular}\n ",
                 "\\\\begin{tablenotes}\n ",
                 "\\\\item \\\\scriptsize{$^{***}p<0.001$, $^{**}p<0.01$, $^*p<0.05$, $^{\\\\dagger}p<0.1$. Two-sided test. Cluster robust standard errors in parentheses (clustered by person and LSOA). Basic controls: year-season of interview (omitted). Additional controls: age (5 year intervals, omitted), highest education, child(ren), marital status, household income deciles.}\n",
                 "\\\\end{tablenotes}\n",
                 "\\\\label{tab:lm_behave_heat}\n",
                 "\\\\end{center}\n",
                 "\\\\end{threeparttable}\n",
                 "}\n")
l <- gregexpr("\\hline\\\\[-1.2ex]", tex_t2, fixed = TRUE)
l <- l[[1]][length(l[[1]])]
tex_t2 <- substr(tex_t2, 1, (l-2))
tex_t2 <- sub("$", bottom, tex_t2, fixed = FALSE)

# N into middle
l1 <- gregexpr("Num. obs", tex_t2, fixed = TRUE)
l1 <- l1[[1]][length(l1[[1]])]


tmp <- substr(tex_t2, l1, (l-2))
tmp2 <- gsub("([[:digit:]]+)", "\\\\multicolumn{1}{c}{\\1}", tmp)
tex_t2 <- gsub(tmp, tmp2, tex_t2, fixed = TRUE)

# Print
write.table(tex_t2, file = "../03_Output/Table2_heat_behav.tex",
            col.names = FALSE, row.names = FALSE, quote = FALSE)




#############
#### FE ####
#############


### Attitudes

tex_t1 <- texreg(l = res2[[1]],
                 #file="../03_Output/Mod_fe_siting.tex",
                 digits = 3, leading.zero = TRUE, fontsize = "scriptsize",
                 stars = c(0.001, 0.01, 0.05, 0.1),
                 symbol = "\\dagger",
                 label = "tab:fe_belief_heat",
                 caption = "Individual fixed effects model. Dep. var.: Belief in climate change.",
                 custom.model.names = c("(1)", "(2)", "(3)", "(4)", "(5)", "(6)"),
                 groups=list("Education (ref: GCSE)" = 2:6, 
                             "Marital (ref: married/civil partner)" = 8:12,
                             "Income (ref: lowest decile)" = 13:21),
                 custom.coef.map = list('heat12_affect1y_past1' =  'Heatwave affected',
                                        'heat6_affect1y_past1' =  'Heatwave affected',
                                        'heat1_affect1y_past1' =  'Heatwave affected',
                                        "hiqual_dv1"            = "Degree",
                                        "hiqual_dv2"            = "Other higher degree",
                                        "hiqual_dv3"            = "A-level etc",
                                        "hiqual_dv5"            = "Other qualification",
                                        "hiqual_dv9"            = "No qualification",
                                        "child"                 = "Child(ren)",
                                        "marstat_dv2"           = "Living as couple",
                                        "marstat_dv3"           = "Widowed",
                                        "marstat_dv4"           = "Divorced",
                                        "marstat_dv5"           = "Separated",
                                        "marstat_dv6"           = "Never married",
                                        "hhinc_dec(1.04,1.5]"   = "Income dec 2",
                                        "hhinc_dec(1.5,1.95]"   = "Income dec 3",
                                        "hhinc_dec(1.95,2.43]"  = "Income dec 4",
                                        "hhinc_dec(2.43,2.96]"  = "Income dec 5",
                                        "hhinc_dec(2.96,3.55]"  = "Income dec 6",
                                        "hhinc_dec(3.55,4.31]"  = "Income dec 7",
                                        "hhinc_dec(4.31,5.27]"  = "Income dec 8",
                                        "hhinc_dec(5.27,6.87]"  = "Income dec 9",
                                        "hhinc_dec(6.87,89.7]"  = "Income dec 10"),
                 custom.note = "%stars. Two-sided test. Cluster robust standard errors in parentheses (clustered by person and LSOA).",
                 dcolumn = TRUE, caption.above = TRUE, use.packages = FALSE, include.proj.stats = FALSE
)


# Customize
tex_t1 <- gsub("D[{].[}][{].[}][{][[:digit:]]+\\.*[[:digit:]]*[}]", "D{.}{.}{2.4}", tex_t1)
tex_t1 <- gsub("\\hline", "\\hline\\\\[-1.2ex]", tex_t1, tex_t1, fixed = TRUE)
tex_t1 <- gsub("pidp", "Person-ID", tex_t1, tex_t1, fixed = TRUE)



head <- c("\\centering\n\\footnotesize\n{\\begin{threeparttable}")
tex_t1 <- gsub("\n\\begin{table}", head, tex_t1, fixed = TRUE)

head2 <- c("\\hline\\\\[-1.2ex] \n  & \\multicolumn{2}{c}{Heatwave 4 months}  & \\multicolumn{2}{c}{Heatwave 1 month}  & \\multicolumn{2}{c}{Heatwave 14 days}  \\\\ 
           \\cmidrule(lr){2-3} \\cmidrule(lr){4-5} \\cmidrule(lr){6-7}") 
tex_t1 <- sub("\\hline\\\\[-1.2ex]", head2, tex_t1, fixed = TRUE)

cont <- paste0("\\hline\\\\[-1.2ex]
  Basic controls & \\multicolumn{1}{c}{Yes} & \\multicolumn{1}{c}{Yes} & \\multicolumn{1}{c}{Yes} & \\multicolumn{1}{c}{Yes} & \\multicolumn{1}{c}{Yes} & \\multicolumn{1}{c}{Yes}   \\\\
  Additional controls & \\multicolumn{1}{c}{No} & \\multicolumn{1}{c}{Yes} & \\multicolumn{1}{c}{No} & \\multicolumn{1}{c}{Yes} & \\multicolumn{1}{c}{No} & \\multicolumn{1}{c}{Yes}  \\\\
 \\hline\\\\[-1.2ex]\n","R$^2$")
tex_t1 <- sub("\\hline\\\\[-1.2ex]\nR$^2$", cont, tex_t1, fixed = TRUE)




bottom <- paste0("\n\\\\hline\\\\\\\\[-1.2ex]\n ",
                 "\\\\end{tabular}\n ",
                 "\\\\end{scriptsize}\n ",
                 "\\\\begin{tablenotes}\n ",
                 "\\\\item \\\\scriptsize{$^{***}p<0.001$, $^{**}p<0.01$, $^*p<0.05$, $^{\\\\dagger}p<0.1$. Two-sided test. Cluster robust standard errors in parentheses (clustered by person and LSOA). Basic controls: year-season of interview (omitted). Additional controls: age (5 year intervals, omitted), highest education, child(ren), marital status, household income deciles.}\n",
                 "\\\\end{tablenotes}\n",
                 "\\\\label{tab:fe_belief_heat}\n",
                 "\\\\end{center}\n",
                 "\\\\end{threeparttable}\n",
                 "}\n")
l <- gregexpr("\\hline\\\\[-1.2ex]", tex_t1, fixed = TRUE)
l <- l[[1]][length(l[[1]])]
tex_t1 <- substr(tex_t1, 1, (l-2))
tex_t1 <- sub("$", bottom, tex_t1, fixed = FALSE)

# N into middle
l1 <- gregexpr("Num. obs", tex_t1, fixed = TRUE)
l1 <- l1[[1]][length(l1[[1]])]


tmp <- substr(tex_t1, l1, (l-2))
tmp2 <- gsub("([[:digit:]]+)", "\\\\multicolumn{1}{c}{\\1}", tmp)
tex_t1 <- gsub(tmp, tmp2, tex_t1, fixed = TRUE)

# Print
write.table(tex_t1, file = "../03_Output/Table1_heat_fe_attitudes.tex",
            col.names = FALSE, row.names = FALSE, quote = FALSE)



### Behaviour

tex_t2 <- texreg(l = res2[[2]],
                 #file="../03_Output/Mod_fe_siting.tex",
                 digits = 3, leading.zero = TRUE, fontsize = "scriptsize",
                 stars = c(0.001, 0.01, 0.05, 0.1),
                 symbol = "\\dagger",
                 label = "tab:fe_behave_heat",
                 caption = "Individual fixed effects model. Dep. var.: Environmental behaviour.",
                 custom.model.names = c("(1)", "(2)", "(3)", "(4)", "(5)", "(6)"),
                 groups=list("Education (ref: GCSE)" = 2:6, 
                             "Marital (ref: married/civil partner)" = 8:12,
                             "Income (ref: lowest decile)" = 13:21),
                 custom.coef.map = list('heat12_affect1y_past1' =  'Heatwave affected',
                                        'heat6_affect1y_past1' =  'Heatwave affected',
                                        'heat1_affect1y_past1' =  'Heatwave affected',
                                        "hiqual_dv1"            = "Degree",
                                        "hiqual_dv2"            = "Other higher degree",
                                        "hiqual_dv3"            = "A-level etc",
                                        "hiqual_dv5"            = "Other qualification",
                                        "hiqual_dv9"            = "No qualification",
                                        "child"                 = "Child(ren)",
                                        "marstat_dv2"           = "Living as couple",
                                        "marstat_dv3"           = "Widowed",
                                        "marstat_dv4"           = "Divorced",
                                        "marstat_dv5"           = "Separated",
                                        "marstat_dv6"           = "Never married",
                                        "hhinc_dec(1.04,1.5]"   = "Income dec 2",
                                        "hhinc_dec(1.5,1.95]"   = "Income dec 3",
                                        "hhinc_dec(1.95,2.43]"  = "Income dec 4",
                                        "hhinc_dec(2.43,2.96]"  = "Income dec 5",
                                        "hhinc_dec(2.96,3.55]"  = "Income dec 6",
                                        "hhinc_dec(3.55,4.31]"  = "Income dec 7",
                                        "hhinc_dec(4.31,5.27]"  = "Income dec 8",
                                        "hhinc_dec(5.27,6.87]"  = "Income dec 9",
                                        "hhinc_dec(6.87,89.7]"  = "Income dec 10"),
                 custom.note = "%stars. Two-sided test. Cluster robust standard errors in parentheses (clustered by person and LSOA).",
                 dcolumn = TRUE, caption.above = TRUE, use.packages = FALSE, include.proj.stats = FALSE
)


# Customize
tex_t2 <- gsub("D[{].[}][{].[}][{][[:digit:]]+\\.*[[:digit:]]*[}]", "D{.}{.}{2.4}", tex_t2)
tex_t2 <- gsub("\\hline", "\\hline\\\\[-1.2ex]", tex_t2, tex_t2, fixed = TRUE)
tex_t2 <- gsub("pidp", "Person-ID", tex_t2, tex_t2, fixed = TRUE)



head <- c("\\centering\n\\footnotesize\n{\\begin{threeparttable}")
tex_t2 <- gsub("\n\\begin{table}", head, tex_t2, fixed = TRUE)

head2 <- c("\\hline\\\\[-1.2ex] \n  & \\multicolumn{2}{c}{Heatwave 4 months}  & \\multicolumn{2}{c}{Heatwave 1 month}  & \\multicolumn{2}{c}{Heatwave 14 days}  \\\\ 
           \\cmidrule(lr){2-3} \\cmidrule(lr){4-5} \\cmidrule(lr){6-7}") 
tex_t2 <- sub("\\hline\\\\[-1.2ex]", head2, tex_t2, fixed = TRUE)

cont <- paste0("\\hline\\\\[-1.2ex]
  Basic controls & \\multicolumn{1}{c}{Yes} & \\multicolumn{1}{c}{Yes} & \\multicolumn{1}{c}{Yes} & \\multicolumn{1}{c}{Yes} & \\multicolumn{1}{c}{Yes} & \\multicolumn{1}{c}{Yes}   \\\\
  Additional controls & \\multicolumn{1}{c}{No} & \\multicolumn{1}{c}{Yes} & \\multicolumn{1}{c}{No} & \\multicolumn{1}{c}{Yes} & \\multicolumn{1}{c}{No} & \\multicolumn{1}{c}{Yes}  \\\\
 \\hline\\\\[-1.2ex]\n","R$^2$")
tex_t2 <- sub("\\hline\\\\[-1.2ex]\nR$^2$", cont, tex_t2, fixed = TRUE)



bottom <- paste0("\n\\\\hline\\\\\\\\[-1.2ex]\n ",
                 "\\\\end{tabular}\n ",
                 "\\\\end{scriptsize}\n ",
                 "\\\\begin{tablenotes}\n ",
                 "\\\\item \\\\scriptsize{$^{***}p<0.001$, $^{**}p<0.01$, $^*p<0.05$, $^{\\\\dagger}p<0.1$. Two-sided test. Cluster robust standard errors in parentheses (clustered by person and LSOA). Basic controls: year-season of interview (omitted). Additional controls: age (5 year intervals, omitted), highest education, child(ren), marital status, household income deciles.}\n",
                 "\\\\end{tablenotes}\n",
                 "\\\\label{tab:fe_behave_heat}\n",
                 "\\\\end{center}\n",
                 "\\\\end{threeparttable}\n",
                 "}\n")
l <- gregexpr("\\hline\\\\[-1.2ex]", tex_t2, fixed = TRUE)
l <- l[[1]][length(l[[1]])]
tex_t2 <- substr(tex_t2, 1, (l-2))
tex_t2 <- sub("$", bottom, tex_t2, fixed = FALSE)

# N into middle
l1 <- gregexpr("Num. obs", tex_t2, fixed = TRUE)
l1 <- l1[[1]][length(l1[[1]])]


tmp <- substr(tex_t2, l1, (l-2))
tmp2 <- gsub("([[:digit:]]+)", "\\\\multicolumn{1}{c}{\\1}", tmp)
tex_t2 <- gsub(tmp, tmp2, tex_t2, fixed = TRUE)

# Print
write.table(tex_t2, file = "../03_Output/Table2_heat_fe_behav.tex",
            col.names = FALSE, row.names = FALSE, quote = FALSE)









###-----------------------------------------###
###                   Plot                  ###
###-----------------------------------------###

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}



##################
### Attitudes  ###
##################


#### Coefficients Plot ####

# Set up df
eff <- data.frame(matrix(NA, ncol = 7, nrow = 12))
colnames(eff) <- c("Variable", "Coefficient", "SE", "t", "p", "N", "modelName")


eff[1:nrow(eff), 1] <- c(rep(c(1:2), times = 3), rep(c(3:4), times = 3))

eff$modelName <- rep(rep(c(1:3), each = 2), 2)

# Plug in values all
for(i in c(1:6)){
  eff[i, 2:6] <- c(summary(res1[[1]][[i]])$coefficients[1, 1:4], length(res1[[1]][[i]]$residuals))
  eff[i+6, 2:6] <- c(summary(res2[[1]][[i]])$coefficients[1, 1:4], length(res2[[1]][[i]]$residuals))
}


### Combine df, make labels

allModelFrame <- data.frame(eff)

allModelFrame$Variable <- factor(allModelFrame$Variable, levels = rev(c(1:4)),
                                 labels = rev(c("POLS",
                                                "POLS\n + additional controls",
                                                "FE",
                                                "FE\n + additional controls")))

allModelFrame$modelName <- factor(allModelFrame$modelName, levels = rev(c(1:3)),
                                  labels = rev(c("Heatwave 4 months", 
                                                 "Heatwave 1 month", 
                                                 "Heatwave 14 days")))



# #  Reorder levels
# allModelFrame$Variable <- factor(allModelFrame$Variable,
#                                levels(allModelFrame$Variable)[rev(c(1,6,2,3,4,5,7))])
# allModelFrame$modelName <- factor(allModelFrame$modelName,
#                                levels(allModelFrame$modelName)[rev(c(1,2))])


# Add to combined df
allModelFrame$eff <- 1
mf_all <- allModelFrame



# Confidence intervals
interval2  <-  -qnorm((1-0.95)/2)  # 95% multiplier

# Colours
# cols <- viridis_pal(begin = 0.2, end = 0.75, direction = 1, option = "magma")(3)
cols <- gg_color_hue(3)


# Plot
zp1 <- ggplot(allModelFrame, aes(colour = modelName, shape = modelName))
zp1 <- zp1 + geom_hline(yintercept = 0, colour = gray(1/2), lty = 2)
zp1 <- zp1 + geom_pointrange(aes(x = Variable, y = Coefficient, ymin = Coefficient - SE*interval2,
                                 ymax = Coefficient + SE*interval2),
                             lwd = 1, position = position_dodge(width = 1/2),
                             fill ="black")
zp1 <- zp1 + coord_flip() + theme_bw()
zp1 <- zp1 + theme(legend.title = element_blank())
zp1 <- zp1 + labs(y = "Coefficients of past heat experience", x = "")
zp1 <- zp1 + scale_shape_manual(values = rev(c(19, 17, 15)), labels = rev(c("Heatwave 4 months", 
                                                                            "Heatwave 1 month",
                                                                            "Heatwave 14 days"))) 
zp1 <- zp1 + scale_color_manual(values = c(cols), labels = rev(c("Heatwave 4 months", 
                                                                 "Heatwave 1 month",
                                                                 "Heatwave 14 days")))
# zp1 <- zp1 + scale_color_hue(direction = -1, h.start=90)
zp1 <- zp1 + theme(text = element_text(family = "CMU Serif", size = 20),
                   legend.position = c(0.95,0.99), legend.justification = c(0.95,0.99),
                   legend.background = element_blank(),
                   legend.box.background = element_rect(colour = "black"),
                   legend.key=element_blank(),
                   axis.text.y=element_text(colour="black"),
                   axis.title.x=element_text(colour="black"),
                   axis.text.x = element_text(colour="black"),
                   #axis.text.x = element_text(size=16),
                   #axis.text.y = element_text(size=16),
                   #axis.title.x = element_text(size=16),
                   #axis.title.y = element_text(size=16)
)
zp1 <- zp1 + ggtitle(element_blank())
zp1 <- zp1 + guides(colour = guide_legend(override.aes = list(linetype = 0), reverse = T),
                    shape = guide_legend(reverse = T))


print(zp1)


cairo_pdf(file="../03_Output/Coefplot1_heat_belief.pdf", width = 9, height = 9, bg = "white", family = "CMU Serif")
par(mar=c(0,0,0,0))
par(mfrow=c(1,1),oma=c(0,0,0,0))
zp1
dev.off()






##################
### Behaviour  ###
##################


#### Coefficients Plot ####

# Set up df
eff <- data.frame(matrix(NA, ncol = 7, nrow = 12))
colnames(eff) <- c("Variable", "Coefficient", "SE", "t", "p", "N", "modelName")


eff[1:nrow(eff), 1] <- c(rep(c(1:2), times = 3), rep(c(3:4), times = 3))

eff$modelName <- rep(rep(c(1:3), each = 2), 2)

# Plug in values all
for(i in c(1:6)){
  eff[i, 2:6] <- c(summary(res1[[2]][[i]])$coefficients[1, 1:4], length(res1[[2]][[i]]$residuals))
  eff[i+6, 2:6] <- c(summary(res2[[2]][[i]])$coefficients[1, 1:4], length(res2[[2]][[i]]$residuals))
}



### Combine df, make labels

allModelFrame <- data.frame(eff)

allModelFrame$Variable <- factor(allModelFrame$Variable, levels = rev(c(1:4)),
                                 labels = rev(c("POLS",
                                                "POLS\n + additional controls",
                                                "FE",
                                                "FE\n + additional controls")))

allModelFrame$modelName <- factor(allModelFrame$modelName, levels = rev(c(1:3)),
                                  labels = rev(c("Heatwave 4 months", 
                                                 "Heatwave 1 month", 
                                                 "Heatwave 14 days")))

# #  Reorder levels
# allModelFrame$Variable <- factor(allModelFrame$Variable,
#                                levels(allModelFrame$Variable)[rev(c(1,6,2,3,4,5,7))])
# allModelFrame$modelName <- factor(allModelFrame$modelName,
#                                levels(allModelFrame$modelName)[rev(c(1,2))])


# Add to combined df
allModelFrame$eff <- 2
mf_all <- rbind(mf_all, allModelFrame)


# Confidence intervals
interval2  <-  -qnorm((1-0.95)/2)  # 95% multiplier

# Colours
# cols <- viridis_pal(begin = 0.2, end = 0.75, direction = 1, option = "magma")(3)
cols <- gg_color_hue(3)



# Plot
zp3 <- ggplot(allModelFrame, aes(colour = modelName, shape = modelName))
zp3 <- zp3 + geom_hline(yintercept = 0, colour = gray(1/2), lty = 2)
zp3 <- zp3 + geom_pointrange(aes(x = Variable, y = Coefficient, ymin = Coefficient - SE*interval2,
                                 ymax = Coefficient + SE*interval2),
                             lwd = 1, position = position_dodge(width = 1/2),
                             fill ="black")
zp3 <- zp3 + coord_flip() + theme_bw()
zp3 <- zp3 + theme(legend.title = element_blank())
zp3 <- zp3 + labs(y = "Coefficients of past heat experience", x = "")
zp3 <- zp3 + scale_shape_manual(values = rev(c(19, 17, 15)), labels = rev(c("Heatwave 4 months", 
                                                                            "Heatwave 1 month",
                                                                            "Heatwave 14 days"))) 
zp3 <- zp3 + scale_color_manual(values = c(cols), labels = rev(c("Heatwave 4 months", 
                                                                 "Heatwave 1 month",
                                                                 "Heatwave 14 days")))
# zp3 <- zp3 + scale_color_hue(direction = -1, h.start=90)
zp3 <- zp3 + theme(text = element_text(family = "CMU Serif", size = 20),
                   legend.position = c(0.95,0.99), legend.justification = c(0.95,0.99),
                   legend.background = element_blank(),
                   legend.box.background = element_rect(colour = "black"),
                   legend.key=element_blank(),
                   axis.text.y=element_text(colour="black"),
                   axis.title.x=element_text(colour="black"),
                   axis.text.x = element_text(colour="black"),
                   #axis.text.x = element_text(size=16),
                   #axis.text.y = element_text(size=16),
                   #axis.title.x = element_text(size=16),
                   #axis.title.y = element_text(size=16)
)
zp3 <- zp3 + ggtitle(element_blank())
zp3 <- zp3 + guides(colour = guide_legend(override.aes = list(linetype = 0), reverse = T),
                    shape = guide_legend(reverse = T))


print(zp3)


cairo_pdf(file="../03_Output/Coefplot3_heat_behav.pdf", width = 9, height = 9, bg = "white", family = "CMU Serif")
par(mar=c(0,0,0,0))
par(mfrow=c(1,1),oma=c(0,0,0,0))
zp3
dev.off()





#####################
### Combine Plots ### 
#####################




#### Combine plots ####

mf_all$eff <- factor(mf_all$eff, levels = c(1, 2),
                     labels = c("Climate change belief\n (0 - 1)",
                                "Pro-environmental behaviour\n (1 - 5)"))

# # Set min and max
# mf_all$y_min <- NA
# mf_all$y_min[which(mf_all$eff %in% c("Direct x1", "Direct x2"))] <- -0.3
# mf_all$y_min[which(mf_all$eff %in% c("Indirect x1", "Indirect x2"))] <- -1.5
# 
# mf_all$y_max <- NA
# mf_all$y_max[which(mf_all$eff %in% c("Direct x1", "Direct x2"))] <- 0.7
# mf_all$y_max[which(mf_all$eff %in% c("Indirect x1", "Indirect x2"))] <- 2.2



# Confidence intervals
interval2 <- -qnorm((1-0.95)/2)  # 95% multiplier

mf_all$lb <- mf_all$Coefficient - mf_all$SE * interval2
mf_all$ub <- mf_all$Coefficient + mf_all$SE * interval2

# Set limits (trim data manually)
# xlim = c(-1.5, 2.2)
# 
# mf_all$ub[which(mf_all$ub > mf_all$y_max)] <- mf_all$y_max[which(mf_all$ub > mf_all$y_max)]
# mf_all$lb[which(mf_all$lb < mf_all$y_min)] <- mf_all$y_min[which(mf_all$lb < mf_all$y_min)]

# Coef Labels
mf_all$lab <- as.character(sprintf("%.3f", round(mf_all$Coefficient, 3)))

mf_all$lab[mf_all$p <= 0.1 & mf_all$p > 0.05] <- paste0(mf_all$lab[mf_all$p <= 0.1 & mf_all$p > 0.05], expression("\u2020"))
mf_all$lab[mf_all$p <= 0.05 & mf_all$p > 0.01] <- paste0(mf_all$lab[mf_all$p <= 0.05 & mf_all$p > 0.01], "*")
mf_all$lab[mf_all$p <= 0.01 & mf_all$p > 0.001] <- paste0(mf_all$lab[mf_all$p <= 0.01 & mf_all$p > 0.001], "**")
mf_all$lab[mf_all$p <= 0.001] <- paste0(mf_all$lab[mf_all$p <= 0.001], "***")

# Number of cases position
mf_all$maxx <- ave(mf_all$ub, 
                   by = mf_all$eff,
                   FUN = function(x) max(x))
mf_all$minx <- ave(mf_all$lb, 
                   by = mf_all$eff,
                   FUN = function(x) min(x))

mf_all$tx <- NA
oo <- which(mf_all$modelName == levels(mf_all$modelName)[3])
mf_all$tx[oo] <- mf_all$maxx[oo] - (mf_all$maxx[oo] - mf_all$minx[oo]) * 0.60
oo <- which(mf_all$modelName == levels(mf_all$modelName)[2])
mf_all$tx[oo] <- mf_all$maxx[oo] - (mf_all$maxx[oo] - mf_all$minx[oo]) * 0.38
oo <- which(mf_all$modelName == levels(mf_all$modelName)[1])
mf_all$tx[oo] <- mf_all$maxx[oo] - (mf_all$maxx[oo] - mf_all$minx[oo]) * 0.15

mf_all$ty <- NA
oo <- which(mf_all$Variable %in% levels(mf_all$Variable)[5:6])
mf_all$ty[oo] <- 4.5 
oo <- which(mf_all$Variable %in% levels(mf_all$Variable)[3:4])
mf_all$ty[oo] <- 2.5 
oo <- which(mf_all$Variable %in% levels(mf_all$Variable)[1:2])
mf_all$ty[oo] <- 0.5 


mf_all$tn <- NA
oo <- which(mf_all$modelName == levels(mf_all$modelName)[3])
mf_all$tn[oo] <- "N ="

oo <- which(duplicated(mf_all[, c("eff", "tn", "tx", "ty")]))
mf_all$tn[oo] <- NA



# Plot
zp_all <- ggplot(mf_all, aes(colour = modelName, shape = modelName))
zp_all <- zp_all + facet_grid(. ~ eff, scales = "free_x")
zp_all <- zp_all + geom_hline(yintercept = 0, colour = scales::alpha("black", 0.3), lty = 2)
zp_all <- zp_all + geom_pointrange(aes(x = Variable, y = Coefficient, ymin = Coefficient - SE*interval2,
                                       ymax = Coefficient + SE*interval2),
                                   lwd = 0.7, position = position_dodge(width = 1/1.8),
                                   fill = "black")
zp_all <- zp_all + geom_text(aes(label = lab,
                                 x = Variable, 
                                 y = Coefficient + ave(Coefficient, eff, FUN = max)* rep(c(0.3, 0.25), each = 12)), 
                             size = 5, show.legend  = FALSE, 
                              vjust = -0.3, position = position_dodge(width = 1/1.8))
zp_all <- zp_all + geom_text(aes(label = N, x = ty, y = tx), 
                             size = 5, show.legend  = FALSE, 
                             hjust = 0, vjust = -0.3)
zp_all <- zp_all + geom_text(data = unique(mf_all[!is.na(mf_all$tn), c("eff", "tn", "tx", "ty", "modelName")]),
                             mapping = aes(label = tn, x = ty, y = tx), size = 5, 
                             show.legend  = FALSE, colour = "black",
                             hjust = 1.3, vjust = -0.3)
zp_all <- zp_all + geom_blank(aes(x = 4.4, y = 0))
zp_all <- zp_all + coord_flip() + theme_bw()
# zp_all <- zp_all + ylim(xlim[1], xlim[2])
# zp_all <- zp_all + geom_blank(aes(y = y_min)) + geom_blank(aes(y = y_max))
zp_all <- zp_all + scale_x_discrete(expand = c(0.05,0.1) )
zp_all <- zp_all + theme(legend.title = element_blank())
zp_all <- zp_all + labs(y = "Coefficients of past heatwave", x = "")
zp_all <- zp_all + scale_shape_manual(values = rev(c(19, 17, 15)), labels = rev(c("Heatwave 4 months", 
                                                                            "Heatwave 1 month",
                                                                            "Heatwave 14 days"))) 
zp_all <- zp_all + scale_color_manual(values = c(cols), labels = rev(c("Heatwave 4 months", 
                                                                 "Heatwave 1 month",
                                                                 "Heatwave 14 days")))
zp_all <- zp_all + theme(text = element_text(family = "CMU Serif", size = 20),
                         # legend.position = c(0.95, 0.99), legend.justification = c(0.95,0.99),
                         legend.position = "bottom",
                         legend.background = element_blank(),
                         legend.box.background = element_rect(colour = "black"),
                         legend.key = element_blank(),
                         axis.text.y = element_text(colour="black", size = 20),
                         axis.title.x = element_text(colour="black"),
                         axis.text.x = element_text(colour="black"),
                         strip.background = element_blank(),
                         strip.text = element_text(size = 20, colour = "black"),
                         #axis.text.x = element_text(size=16),
                         #axis.text.y = element_text(size=16),
                         #axis.title.x = element_text(size=16),
                         #axis.title.y = element_text(size=16)
)
zp_all <- zp_all + ggtitle(element_blank())
zp_all <- zp_all + guides(colour = guide_legend(override.aes = list(linetype = 0), reverse = T),
                          shape = guide_legend(reverse = T))

zp_all <- zp_all + geom_vline(aes(xintercept = 2.5))
zp_all <- zp_all + geom_vline(aes(xintercept = 4.5))


print(zp_all)

#### Combine plots ####


cairo_pdf(file = paste("../03_Output/", "Coefplot_heat_combined_1.pdf", sep=""), width = 13, height = 10, 
          bg = "white", family="CMU Serif")
par(mar = c(0, 0, 0, 0))
par(mfrow = c(1, 1), oma = c(0, 0, 0, 0))
zp_all
dev.off()



### Stepwise plot ###

for(i in c(1:4)){
  
  mf_tmp <- mf_all
  mf_tmp$alpha <- 1
  s <- unlist(list(c(1:6), c(13:18), c(7:12), c(19:24))[1:i])
  mf_tmp$alpha[s] <- 0
  mf_tmp$alpha <- as.factor(mf_tmp$alpha)
  
  zp_all <- ggplot(mf_tmp, aes(colour = modelName, shape = modelName))
  zp_all <- zp_all + facet_grid(. ~ eff, scales = "free_x")
  zp_all <- zp_all + geom_hline(yintercept = 0, colour = scales::alpha("black", 0.3), lty = 2)
  zp_all <- zp_all + geom_pointrange(aes(x = Variable, y = Coefficient, ymin = Coefficient - SE*interval2,
                                         ymax = Coefficient + SE*interval2, alpha = alpha),
                                     lwd = 0.7, position = position_dodge(width = 1/1.8),
                                     fill = "black")
  zp_all <- zp_all + geom_text(aes(label = lab, alpha = alpha,
                                   x = Variable, y = Coefficient), size = 5, show.legend  = FALSE, 
                               hjust = -0.3, vjust = -0.3, position = position_dodge(width = 1/1.8))
  zp_all <- zp_all + geom_text(aes(label = N, x = ty, y = tx), 
                               size = 5, show.legend  = FALSE, 
                               hjust = 0, vjust = -0.3)
  zp_all <- zp_all + geom_text(data = unique(mf_all[!is.na(mf_all$tn), c("eff", "tn", "tx", "ty", "modelName")]),
                               mapping = aes(label = tn, x = ty, y = tx), size = 5, 
                               show.legend  = FALSE, colour = "black", 
                               hjust = 1.3, vjust = -0.3)
  zp_all <- zp_all + coord_flip() + theme_bw()
  # zp_all <- zp_all + ylim(xlim[1], xlim[2])
  # zp_all <- zp_all + geom_blank(aes(y = y_min)) + geom_blank(aes(y = y_max))
  zp_all <- zp_all + scale_x_discrete(expand = c(0.05,0.1) )
  zp_all <- zp_all + theme(legend.title = element_blank())
  zp_all <- zp_all + labs(y = "Coefficients of past flood experience (within 2 years)", x = "")
  zp_all <- zp_all + scale_shape_manual(values = rev(c(19, 17, 15)), labels = rev(c("Heatwave 4 months", 
                                                                                    "Heatwave 1 month",
                                                                                    "Heatwave 14 days"))) 
  zp_all <- zp_all + scale_color_manual(values = c(cols), labels = rev(c("Heatwave 4 months", 
                                                                         "Heatwave 1 month",
                                                                         "Heatwave 14 days")))
  zp_all <- zp_all + scale_alpha_manual(values = c(1, 0), guide = FALSE)
  zp_all <- zp_all + theme(text = element_text(family = "CMU Serif", size = 20),
                           # legend.position = c(0.95, 0.99), legend.justification = c(0.95,0.99),
                           legend.position = "bottom",
                           legend.background = element_blank(),
                           legend.box.background = element_rect(colour = "black"),
                           legend.key = element_blank(),
                           axis.text.y = element_text(colour="black", size = 20),
                           axis.title.x = element_text(colour="black"),
                           axis.text.x = element_text(colour="black"),
                           strip.background = element_blank(),
                           strip.text = element_text(size = 20, colour = "black"),
                           #axis.text.x = element_text(size=16),
                           #axis.text.y = element_text(size=16),
                           #axis.title.x = element_text(size=16),
                           #axis.title.y = element_text(size=16)
  )
  zp_all <- zp_all + ggtitle(element_blank())
  zp_all <- zp_all + guides(colour = guide_legend(override.aes = list(linetype = 0), reverse = T),
                            shape = guide_legend(reverse = T))
  
  zp_all <- zp_all + geom_vline(aes(xintercept = 2.5))
  zp_all <- zp_all + geom_vline(aes(xintercept = 4.5))
  
  cairo_pdf(file = paste("../03_Output/", "Coefplot_heat_combined_1_", i, ".pdf", sep=""), width = 13, height = 10, 
            bg = "white", family="CMU Serif")
  par(mar = c(0, 0, 0, 0))
  par(mfrow = c(1, 1), oma = c(0, 0, 0, 0))
  print(zp_all)
  dev.off()
  
  
}








###--------------------------------------------------------------###
###         COmbined plot with floods and heatwaves              ###
###--------------------------------------------------------------###


### Combine heatwave and flood results
mf_heat <- mf_all
mf_heat$group <- 2

load("Coeftable_combined_flood.RData")
mf_flood <- mf_all
mf_flood$group <- 1



##################
### Plot POLS  ### 
##################

allModelFrame <- rbind(mf_flood, mf_heat)

### Subset to POLS
allModelFrame <- allModelFrame[which(allModelFrame$Variable == "POLS\n + additional controls" | allModelFrame$Variable == "POLS"), ]

### Labels

allModelFrame$Variable <- as.numeric(allModelFrame$Variable)
allModelFrame$Variable <- allModelFrame$Variable - ifelse(allModelFrame$group == 2, 2, 0)
allModelFrame$Variable <- factor(allModelFrame$Variable, labels = c("Heatwaves\n + additional\n controls",
                                                                    "Heatwaves",
                                                                    "Floods\n + additional\n controls",
                                                                    "Floods"))

allModelFrame$group <- factor(allModelFrame$group, labels = c("Floods", "Heatwaves"))

allModelFrame$modelName <- factor(allModelFrame$modelName,levels(allModelFrame$modelName)[c(3, 2, 1, 6, 5, 4)])

allModelFrame$ty <- allModelFrame$ty - ifelse(allModelFrame$group == "Heatwaves", 2, 0)


### Adjust minimum and maximum
allModelFrame$maxx <- ave(allModelFrame$ub, 
                          by = allModelFrame$eff,
                          FUN = function(x) max(x))
allModelFrame$minx <- ave(allModelFrame$lb, 
                          by = allModelFrame$eff,
                          FUN = function(x) min(x))

allModelFrame$tx <- NA
oo <- which(allModelFrame$modelName == levels(allModelFrame$modelName)[1])
allModelFrame$tx[oo] <- allModelFrame$maxx[oo] - (allModelFrame$maxx[oo] - allModelFrame$minx[oo]) * 0.60
oo <- which(allModelFrame$modelName == levels(allModelFrame$modelName)[2])
allModelFrame$tx[oo] <- allModelFrame$maxx[oo] - (allModelFrame$maxx[oo] - allModelFrame$minx[oo]) * 0.38
oo <- which(allModelFrame$modelName == levels(allModelFrame$modelName)[3])
allModelFrame$tx[oo] <- allModelFrame$maxx[oo] - (allModelFrame$maxx[oo] - allModelFrame$minx[oo]) * 0.15
oo <- which(allModelFrame$modelName == levels(allModelFrame$modelName)[4])
allModelFrame$tx[oo] <- allModelFrame$maxx[oo] - (allModelFrame$maxx[oo] - allModelFrame$minx[oo]) * 0.60
oo <- which(allModelFrame$modelName == levels(allModelFrame$modelName)[5])
allModelFrame$tx[oo] <- allModelFrame$maxx[oo] - (allModelFrame$maxx[oo] - allModelFrame$minx[oo]) * 0.38
oo <- which(allModelFrame$modelName == levels(allModelFrame$modelName)[6])
allModelFrame$tx[oo] <- allModelFrame$maxx[oo] - (allModelFrame$maxx[oo] - allModelFrame$minx[oo]) * 0.15

# vjust values
allModelFrame$labvjust <- 0.3 + (0.05 * (nchar(allModelFrame$lab) - max(nchar(allModelFrame$lab))))

### Plot the coefficients


# Plot
zp_all <- ggplot(allModelFrame, aes(colour = modelName, shape = modelName))
zp_all <- zp_all + facet_grid(. ~ eff, scales = "free_x")
zp_all <- zp_all + geom_hline(yintercept = 0, colour = scales::alpha("black", 0.3), lty = 2)
zp_all <- zp_all + geom_pointrange(aes(x = Variable, y = Coefficient, ymin = Coefficient - SE*interval2,
                                       ymax = Coefficient + SE*interval2),
                                   lwd = 0.7, position = position_dodge2(width = 1/1.8, reverse = TRUE),
                                   fill = "black")
zp_all <- zp_all + geom_text(aes(label = lab,
                                 x = Variable, y = Coefficient), size = 6, show.legend  = FALSE, 
                             hjust = - allModelFrame$labvjust, vjust = -0.3, position = position_dodge2(width = 1/1.8, reverse = TRUE))
zp_all <- zp_all + geom_text(aes(label = N, x = ty, y = tx), 
                             size = 6, show.legend  = FALSE, 
                             hjust = 0, vjust = -0.3)
zp_all <- zp_all + geom_text(data = unique(allModelFrame[!is.na(allModelFrame$tn), c("eff", "tn", "tx", "ty", "modelName")]),
                             mapping = aes(label = tn, x = ty, y = tx), size = 6, 
                             show.legend  = FALSE, colour = "black",
                             hjust = 1.3, vjust = -0.3)
zp_all <- zp_all + geom_blank(aes(x = 4.4, y = 0))
zp_all <- zp_all + coord_flip() + theme_bw()
zp_all <- zp_all + scale_x_discrete(expand = c(0.05,0.00) )
zp_all <- zp_all + theme(legend.title = element_blank())
zp_all <- zp_all + labs(y = "POLS coefficients of past experience", x = "")
zp_all <- zp_all + scale_shape_manual(values = rev(c(6, 5, 8, 15, 17, 19))) 
zp_all <- zp_all + scale_color_manual(values = rev(c(gg_color_hue(6)[4:6], gg_color_hue(6)[1:3])))
zp_all <- zp_all + theme(text = element_text(family = "CMU Serif", size = 20),
                         # legend.position = c(0.95, 0.99), legend.justification = c(0.95,0.99),
                         plot.title = element_text(hjust = 0.5),
                         legend.position = "bottom",
                         legend.box="vertical",
                         legend.background = element_blank(),
                         legend.box.background = element_rect(colour = "black"),
                         legend.key = element_blank(),
                         axis.text.y = element_text(colour="black", size = 20),
                         axis.title.x = element_text(colour="black"),
                         axis.text.x = element_text(colour="black"),
                         strip.background = element_blank(),
                         strip.text = element_text(size = 20, colour = "black"),
                         panel.spacing.x = unit(6, "mm"),
                         plot.margin = unit(c(0,6,0,0), "mm")
                         #axis.text.x = element_text(size=16),
                         #axis.text.y = element_text(size=16),
                         #axis.title.x = element_text(size=16),
                         #axis.title.y = element_text(size=16)
)
zp_all <- zp_all + guides(colour = guide_legend(override.aes = list(linetype = 0), reverse = F, nrow=2, byrow=T),
                          shape = guide_legend(reverse = F, nrow=2, byrow=T))

zp_all <- zp_all + geom_vline(aes(xintercept = 2.5))


print(zp_all)

#### Combine plots ####


cairo_pdf(file = paste("../03_Output/", "Coefplot_POLS_floods_heatwaves.pdf", sep=""), width = 13, height = 10, 
          bg = "white", family="CMU Serif")
par(mar = c(0, 0, 0, 0))
par(mfrow = c(1, 1), oma = c(0, 0, 0, 0))
zp_all
dev.off()




################
### Plot FE  ### 
################

allModelFrame <- rbind(mf_flood, mf_heat)

### Subset to FE
allModelFrame <- allModelFrame[which(allModelFrame$Variable == "FE\n + additional controls" | allModelFrame$Variable == "FE"), ]

### Labels

allModelFrame$Variable <- as.numeric(allModelFrame$Variable)
allModelFrame$Variable <- allModelFrame$Variable + ifelse(allModelFrame$group == 1, 2, 0)
allModelFrame$Variable <- factor(allModelFrame$Variable, labels = c("Heatwaves\n + additional\n controls",
                                                                    "Heatwaves",
                                                                    "Floods\n + additional\n controls",
                                                                    "Floods"))

allModelFrame$group <- factor(allModelFrame$group, labels = c("Floods", "Heatwaves"))

allModelFrame$modelName <- factor(allModelFrame$modelName,levels(allModelFrame$modelName)[c(3, 2, 1, 6, 5, 4)])

allModelFrame$ty <- allModelFrame$ty + ifelse(allModelFrame$group == "Floods", 2, 0)


### Adjust minimum and maximum
allModelFrame$maxx <- ave(allModelFrame$ub, 
                    by = allModelFrame$eff,
                    FUN = function(x) max(x))
allModelFrame$minx <- ave(allModelFrame$lb, 
                    by = allModelFrame$eff,
                    FUN = function(x) min(x))

allModelFrame$tx <- NA
oo <- which(allModelFrame$modelName == levels(allModelFrame$modelName)[1])
allModelFrame$tx[oo] <- allModelFrame$maxx[oo] - (allModelFrame$maxx[oo] - allModelFrame$minx[oo]) * 0.60
oo <- which(allModelFrame$modelName == levels(allModelFrame$modelName)[2])
allModelFrame$tx[oo] <- allModelFrame$maxx[oo] - (allModelFrame$maxx[oo] - allModelFrame$minx[oo]) * 0.38
oo <- which(allModelFrame$modelName == levels(allModelFrame$modelName)[3])
allModelFrame$tx[oo] <- allModelFrame$maxx[oo] - (allModelFrame$maxx[oo] - allModelFrame$minx[oo]) * 0.15
oo <- which(allModelFrame$modelName == levels(allModelFrame$modelName)[4])
allModelFrame$tx[oo] <- allModelFrame$maxx[oo] - (allModelFrame$maxx[oo] - allModelFrame$minx[oo]) * 0.60
oo <- which(allModelFrame$modelName == levels(allModelFrame$modelName)[5])
allModelFrame$tx[oo] <- allModelFrame$maxx[oo] - (allModelFrame$maxx[oo] - allModelFrame$minx[oo]) * 0.38
oo <- which(allModelFrame$modelName == levels(allModelFrame$modelName)[6])
allModelFrame$tx[oo] <- allModelFrame$maxx[oo] - (allModelFrame$maxx[oo] - allModelFrame$minx[oo]) * 0.15

### Plot the coefficients


# Plot
zp_all <- ggplot(allModelFrame, aes(colour = modelName, shape = modelName))
zp_all <- zp_all + facet_grid(. ~ eff, scales = "free_x")
zp_all <- zp_all + geom_hline(yintercept = 0, colour = scales::alpha("black", 0.3), lty = 2)
zp_all <- zp_all + geom_pointrange(aes(x = Variable, y = Coefficient, ymin = Coefficient - SE*interval2,
                                       ymax = Coefficient + SE*interval2),
                                   lwd = 0.7, position = position_dodge2(width = 1/1.8, reverse = TRUE),
                                   fill = "black")
zp_all <- zp_all + geom_text(aes(label = lab,
                                 x = Variable, y = Coefficient), size = 6, show.legend  = FALSE, 
                             hjust = -0.3, vjust = -0.3, position = position_dodge2(width = 1/1.8, reverse = TRUE))
zp_all <- zp_all + geom_text(aes(label = N, x = ty, y = tx), 
                             size = 6, show.legend  = FALSE, 
                             hjust = 0, vjust = -0.3)
zp_all <- zp_all + geom_text(data = unique(allModelFrame[!is.na(allModelFrame$tn), c("eff", "tn", "tx", "ty", "modelName")]),
                             mapping = aes(label = tn, x = ty, y = tx), size = 6, 
                             show.legend  = FALSE, colour = "black",
                             hjust = 1.3, vjust = -0.3)
zp_all <- zp_all + geom_blank(aes(x = 4.4, y = 0))
zp_all <- zp_all + coord_flip() + theme_bw()
zp_all <- zp_all + scale_x_discrete(expand = c(0.05,0.00) )
zp_all <- zp_all + theme(legend.title = element_blank())
zp_all <- zp_all + labs(y = "FE coefficients of past experience", x = "")
zp_all <- zp_all + scale_shape_manual(values = rev(c(6, 5, 8, 15, 17, 19))) 
zp_all <- zp_all + scale_color_manual(values = rev(c(gg_color_hue(6)[4:6], gg_color_hue(6)[1:3])))
zp_all <- zp_all + theme(text = element_text(family = "CMU Serif", size = 20),
                         # legend.position = c(0.95, 0.99), legend.justification = c(0.95,0.99),
                         plot.title = element_text(hjust = 0.5),
                         legend.position = "bottom",
                         legend.box="vertical",
                         legend.background = element_blank(),
                         legend.box.background = element_rect(colour = "black"),
                         legend.key = element_blank(),
                         axis.text.y = element_text(colour="black", size = 20),
                         axis.title.x = element_text(colour="black"),
                         axis.text.x = element_text(colour="black"),
                         strip.background = element_blank(),
                         strip.text = element_text(size = 20, colour = "black"),
                         panel.spacing.x = unit(6, "mm"),
                         plot.margin = unit(c(0,6,0,0), "mm")
                         #axis.text.x = element_text(size=16),
                         #axis.text.y = element_text(size=16),
                         #axis.title.x = element_text(size=16),
                         #axis.title.y = element_text(size=16)
)
zp_all <- zp_all + guides(colour = guide_legend(override.aes = list(linetype = 0), reverse = F, nrow=2, byrow=T),
                          shape = guide_legend(reverse = F, nrow=2, byrow=T))

zp_all <- zp_all + geom_vline(aes(xintercept = 2.5))


print(zp_all)






cairo_pdf(file = paste("../03_Output/", "Coefplot_FE_floods_heatwaves.pdf", sep=""), width = 13, height = 10, 
          bg = "white", family="CMU Serif")
par(mar = c(0, 0, 0, 0))
par(mfrow = c(1, 1), oma = c(0, 0, 0, 0))
zp_all
dev.off()





#### Combine plots stepwise ####
zp_all2 <- ggplot(allModelFrame, aes(colour = modelName, shape = modelName))
zp_all2 <- zp_all2 + facet_grid(. ~ eff, scales = "free_x")
zp_all2 <- zp_all2 + geom_hline(yintercept = 0, colour = scales::alpha("black", 0.3), lty = 2)
zp_all2 <- zp_all2 + geom_pointrange(aes(x = Variable, y = Coefficient, ymin = Coefficient - SE*interval2,
                                       ymax = Coefficient + SE*interval2),
                                   lwd = 0.7, position = position_dodge2(width = 1/1.8, reverse = TRUE),
                                   fill = "black", alpha = 0)
zp_all2 <- zp_all2 + geom_text(aes(label = lab,
                                 x = Variable, y = Coefficient), size = 6, show.legend  = FALSE,  alpha = 0,
                             hjust = -0.3, vjust = -0.3, position = position_dodge2(width = 1/1.8, reverse = TRUE))
zp_all2 <- zp_all2 + geom_text(aes(label = N, x = ty, y = tx), 
                             size = 6, show.legend  = FALSE, 
                             hjust = 0, vjust = -0.3)
zp_all2 <- zp_all2 + geom_text(data = unique(allModelFrame[!is.na(allModelFrame$tn), c("eff", "tn", "tx", "ty", "modelName")]),
                             mapping = aes(label = tn, x = ty, y = tx), size = 6, 
                             show.legend  = FALSE, colour = "black",
                             hjust = 1.3, vjust = -0.3)
zp_all2 <- zp_all2 + geom_blank(aes(x = 4.4, y = 0))
zp_all2 <- zp_all2 + coord_flip() + theme_bw()
zp_all2 <- zp_all2 + scale_x_discrete(expand = c(0.05,0.00) )
zp_all2 <- zp_all2 + theme(legend.title = element_blank())
zp_all2 <- zp_all2 + labs(y = "FE coefficients of past experience", x = "")
zp_all2 <- zp_all2 + scale_shape_manual(values = rev(c(6, 5, 8, 15, 17, 19))) 
zp_all2 <- zp_all2 + scale_color_manual(values = rev(c(gg_color_hue(6)[4:6], gg_color_hue(6)[1:3])))
zp_all2 <- zp_all2 + theme(text = element_text(family = "CMU Serif", size = 20),
                         # legend.position = c(0.95, 0.99), legend.justification = c(0.95,0.99),
                         plot.title = element_text(hjust = 0.5),
                         legend.position = "bottom",
                         legend.box="vertical",
                         legend.background = element_blank(),
                         legend.box.background = element_rect(colour = "black"),
                         legend.key = element_blank(),
                         axis.text.y = element_text(colour="black", size = 20),
                         axis.title.x = element_text(colour="black"),
                         axis.text.x = element_text(colour="black"),
                         strip.background = element_blank(),
                         strip.text = element_text(size = 20, colour = "black"),
                         panel.spacing.x = unit(6, "mm"),
                         plot.margin = unit(c(0,6,0,0), "mm")
                         #axis.text.x = element_text(size=16),
                         #axis.text.y = element_text(size=16),
                         #axis.title.x = element_text(size=16),
                         #axis.title.y = element_text(size=16)
)
zp_all2 <- zp_all2 + guides(colour = guide_legend(override.aes = list(linetype = 0, alpha = 1), reverse = F, nrow=2, byrow=T),
                          shape = guide_legend(reverse = F, nrow=2, byrow=T))

zp_all2 <- zp_all2 + geom_vline(aes(xintercept = 2.5))

print(zp_all2)

cairo_pdf(file = paste("../03_Output/", "Coefplot_FE_floods_heatwaves_0.pdf", sep=""), width = 13, height = 10, 
          bg = "white", family="CMU Serif")
par(mar = c(0, 0, 0, 0))
par(mfrow = c(1, 1), oma = c(0, 0, 0, 0))
zp_all2
dev.off()
  
  
  





###-----------------------------------------###
###         Summary statistics              ###
###-----------------------------------------###

bhps_heat.df$index1 <- bhps_heat.df$c_ccexag_w1 + bhps_heat.df$c_ccnoworry_w1 + (1 - bhps_heat.df$c_ccsoondisaster_w1)
bhps_heat.df$index1 <- ifelse(bhps_heat.df$index1 > 0, 1, 0)

bhps_heat.df$index2 <- NA
oo <- which(bhps_heat.df$c_partisan %in% c(2, 3, 4, 5, 6, 8, 9, 11))
bhps_heat.df$index2[oo] <- 1
oo <- which(bhps_heat.df$c_partisan %in% c(95:99))
bhps_heat.df$index2[oo] <- 2
oo <- which(bhps_heat.df$c_partisan %in% c(1, 7, 10, 12, 13, 14))
bhps_heat.df$index2[oo] <- 3

bhps_heat.df$index1 <- as.integer(bhps_heat.df$index1)
bhps_heat.df$index2 <- as.integer(bhps_heat.df$index2)

bhps_heat.df$index2 <- as.factor(bhps_heat.df$index2)




### Sample 1

vars <- c("climatechange30", 
          "heat1_affect1y_past1", "heat6_affect1y_past1", "heat12_affect1y_past1", 
          "age_dv", "female", "migback", "ethn_dv_short", "hiqual_dv", "child", 
          "marstat_dv", "hhinc", "index1", "index2")

sample1 <- bhps_heat.df[which(bhps_heat.df$samplefe_climatechange30 == 1), vars]

# Make model frame
oo <- names(Filter(is.factor, sample1))
sample1[, oo] <- droplevels(sample1[, oo])
for(i in oo){
  tmp <-  sapply(levels(sample1[, i]), function(x) as.integer(x == sample1[, i]))
  colnames(tmp) <- paste0(i, "_", colnames(tmp))
  sample1 <- cbind(sample1[, -which(names(sample1) == i)], tmp)
}


stg1 <- stargazer(sample1[ , ], 
                  #out = "../03_Output/summarystats1.tex", 
                  type = "latex", style = "asr", digits = 3, align = T, 
                  summary.stat  =  c("n", "mean", "sd", "min", "max"),
                  label  =  "tab:desc3",
                  font.size = "scriptsize", 
                  table.placement = "h!",
                  column.sep.width  =  ".2pt" , title  =  "Estimation sample 3 (Heatwaves - Climate change belief)",
                  covariate.labels = c("Climate change belief", 
                                       #"Pro-environmental behaviour", 
                                       "Heatwave affected (14 days)",   
                                       "Heatwave affected (1 month)", 
                                       "Heatwave affected (4 months)", 
                                       "Age",
                                       "Sex (female)", 
                                       "Migration background", 
                                       "Child(ren) in household",
                                       "Household income (in thousand)",  
                                       "High scepticism",
                                       "XX1 \\quad Any White", 
                                       "\\quad Mixed",
                                       "\\quad Asian", 
                                       "\\quad Black", 
                                       "\\quad Other",
                                       "XX2 \\quad GCSE etc",
                                       "\\quad Degree", 
                                       "\\quad Other higher degree", 
                                       "\\quad A-level etc",
                                       "\\quad Other qualification", 
                                       "\\quad No qualification",
                                       "XX3 \\quad Married/Civil partner", 
                                       "\\quad Living as couple", 
                                       "\\quad Widowed/surviving civil partner",
                                       "\\quad Divorced/dissolved civil partner", 
                                       "\\quad Separated (incl. from civil partner)", 
                                       "\\quad Never married",
                                       "XX4 \\quad Left-leaning partisanship", 
                                       "\\quad No partisanship",
                                       "\\quad Right-leaning partisanship"
                                       # "XX4 \\quad None",
                                       # "\\quad Conservative",
                                       # "\\quad Labour",
                                       # "\\quad Lib dem /lib/sdp",
                                       # "\\quad Scot nat",
                                       # "\\quad Plaid cymru",
                                       # "\\quad Green party",
                                       # "\\quad Ulster unionist",
                                       # "\\quad Sdlp",
                                       # "\\quad Alliance party",
                                       # "\\quad Democratic unionist",
                                       # "\\quad Sinn fein",
                                       # "\\quad Cant vote",
                                       # "\\quad Other party"
                                       )
)


### Finalize table
stg1 <- gsub("XX1", "Ethnic background & & & &  &  \\\\\\\\ ",
             stg1,)
stg1 <- gsub("XX2", "Highest education & & & &  &  \\\\\\\\ ",
             stg1)
stg1 <- gsub("XX3", "Marital status & & & &  &  \\\\\\\\ ",
             stg1)
stg1 <- gsub("XX4", "Partisanship & & & &  &  \\\\\\\\ ",
             stg1)

stg1 <- gsub(".000", "",
             stg1)

stg1 <- gsub("D{.}{.}{-3} D{.}{.}{-3} D{.}{.}{-3} D{.}{.}{-3} D{.}{.}{-3}", 
             "D{.}{.}{5.0} D{.}{.}{3.3} D{.}{.}{3.3} D{.}{.}{3.3} D{.}{.}{3.3}",
             stg1, fixed = TRUE)

write.table(stg1, file = "../03_Output/summarystats3.tex", 
            quote = FALSE, col.names = FALSE, row.names = FALSE)




### Sample 2

vars <- c("env_index_all2",
          "heat1_affect1y_past1", "heat6_affect1y_past1", "heat12_affect1y_past1", 
          "age_dv", "female", "migback", "ethn_dv_short", "hiqual_dv", "child", 
          "marstat_dv", "hhinc", "index1", "index2")

sample2 <- bhps_heat.df[which(bhps_heat.df$samplefe_env_index_all2 == 1), vars]

# Make model frame
oo <- names(Filter(is.factor, sample2))
sample2[, oo] <- droplevels(sample2[, oo])
for(i in oo){
  tmp <-  sapply(levels(sample2[, i]), function(x) as.integer(x == sample2[, i]))
  colnames(tmp) <- paste0(i, "_", colnames(tmp))
  sample2 <- cbind(sample2[, -which(names(sample2) == i)], tmp)
}


stg2 <- stargazer(sample2[ , ], 
                  #out = "../03_Output/summarystats1.tex", 
                  type = "latex", style = "asr", digits = 3, align = T, 
                  summary.stat  =  c("n", "mean", "sd", "min", "max"),
                  label  =  "tab:desc4",
                  font.size = "scriptsize", 
                  table.placement = "h!",
                  column.sep.width  =  ".2pt" , title  =  "Estimation sample 4 (Heatwaves - Behaviour)",
                  covariate.labels = c(#"Climate change belief", 
                    "Pro-environmental behaviour", 
                    "Heatwave affected (14 days)",   
                    "Heatwave affected (1 month)", 
                    "Heatwave affected (4 months)", 
                    "Age",
                    "Sex (female)", 
                    "Migration background", 
                    "Child(ren) in household",
                    "Household income (in thousand)",  
                    "High scepticism",
                    "XX1 \\quad Any White", 
                    "\\quad Mixed",
                    "\\quad Asian", 
                    "\\quad Black", 
                    "\\quad Other",
                    "XX2 \\quad GCSE etc",
                    "\\quad Degree", 
                    "\\quad Other higher degree", 
                    "\\quad A-level etc",
                    "\\quad Other qualification", 
                    "\\quad No qualification",
                    "XX3 \\quad Married/Civil partner", 
                    "\\quad Living as couple", 
                    "\\quad Widowed/surviving civil partner",
                    "\\quad Divorced/dissolved civil partner", 
                    "\\quad Separated (incl. from civil partner)", 
                    "\\quad Never married",
                    "XX4 \\quad Left-leaning partisanship", 
                    "\\quad No partisanship",
                    "\\quad Right-leaning partisanship"
                    # "XX4 \\quad None",
                    # "\\quad Conservative",
                    # "\\quad Labour",
                    # "\\quad Lib dem /lib/sdp",
                    # "\\quad Scot nat",
                    # "\\quad Plaid cymru",
                    # "\\quad Green party",
                    # "\\quad Ulster unionist",
                    # "\\quad Sdlp",
                    # "\\quad Alliance party",
                    # "\\quad Democratic unionist",
                    # "\\quad Sinn fein",
                    # "\\quad Cant vote",
                    # "\\quad Other party"
                    )
)


### Finalize table
stg2 <- gsub("XX1", "Ethnic background & & & &  &  \\\\\\\\ ",
             stg2,)
stg2 <- gsub("XX2", "Highest education & & & &  &  \\\\\\\\ ",
             stg2)
stg2 <- gsub("XX3", "Marital status & & & &  &  \\\\\\\\ ",
             stg2)
stg2 <- gsub("XX4", "Partisanship & & & &  &  \\\\\\\\ ",
             stg2)

stg2 <- gsub(".000", "",
             stg2)

stg2 <- gsub("D{.}{.}{-3} D{.}{.}{-3} D{.}{.}{-3} D{.}{.}{-3} D{.}{.}{-3}", 
             "D{.}{.}{5.0} D{.}{.}{3.3} D{.}{.}{3.3} D{.}{.}{3.3} D{.}{.}{3.3}",
             stg2, fixed = TRUE)

write.table(stg2, file = "../03_Output/summarystats4.tex", 
            quote = FALSE, col.names = FALSE, row.names = FALSE)















###-----------------------------------------------------###
###         Effect heterogeneity  - Exaggerated         ###
###-----------------------------------------------------###

bhps_heat.df$index1 <- bhps_heat.df$c_ccexag_w1 + bhps_heat.df$c_ccnoworry_w1 + (1 - bhps_heat.df$c_ccsoondisaster_w1)
bhps_heat.df$index1 <- ifelse(bhps_heat.df$index1 > 0, 1, 0)


### exag ###

depvar <- c("climatechange30", "env_index_all2")
intens <- c("12", "6", "1")
buffer <- c(1, 2, 3)
b <- buffer[1]
plm.l <- vector("list", nrow(expand.grid(intens)) * 1)
plm3.l <- vector("list", nrow(expand.grid(buffer)) * 1) 

res4 <- vector("list", length(depvar)) 

l <- 1
for(d in depvar) {
  library(lfe)
  names(plm.l) <- paste0(d, c(1:length(plm.l)))
  names(plm3.l) <- paste0(d, c(1:length(plm3.l)))
  
  j <- 1
  k <- 1
  for(i in intens){
    
    v1 <- paste0("heat", i, "_affect1y_past", b)
    ev1 <- paste0("heat_ever", b)
    
    fm2 <- paste0(d, " ~ ", v1, 
                  "+ age_cat + hiqual_dv + child + marstat_dv + hhinc_dec")
    
    # Estimate via felm with twoway clustered SEs
    bhps_heat.df$interaction <- as.numeric(as.factor(bhps_heat.df$season_index)) * bhps_heat.df[, ev1]
    
    plm2_1 <- felm(formula(paste0(fm2, " | pidp + season_index + interaction | 0 | pidp + lsoa01")),
                   data = bhps_heat.df[which(bhps_heat.df[, paste0("samplefe_", d)] == 1 & bhps_heat.df[, "index1"] == 0), ],
                   psdef = FALSE)

    plm2_2 <- felm(formula(paste0(fm2, " | pidp + season_index + interaction | 0 | pidp + lsoa01")),
                   data = bhps_heat.df[which(bhps_heat.df[, paste0("samplefe_", d)] == 1 & bhps_heat.df[, "index1"] == 1), ],
                   psdef = FALSE)

    
    
    # Interacted model
    fm3 <- paste0(d, " ~ index1 * (", v1, 
                  "+ age_cat + hiqual_dv + child + marstat_dv + hhinc_dec", 
                  "  + as.factor(season_index)*", ev1, ")")
    
    plm3 <- felm(formula(paste0(fm3, " | pidp | 0 | pidp + lsoa01")),
                 data = bhps_heat.df[which(bhps_heat.df[, paste0("samplefe_", d)] == 1), ], 
                 psdef = FALSE)
    
    
    plm.l[[j]] <- plm2_1
    plm.l[[j+1]] <- plm2_2
    plm3.l[[k]] <- plm3
    j <- j + 2
    k <- k + 1
    
  }
  res4[[l]] <- list(plm.l,  plm3.l)
  l <- l + 1
}



# Map
ns <- as.list(paste0("heat", intens, "_affect1y_past", b))
names(ns) <- paste0("heat", intens, "_affect1y_past", b)


# Climate change
screenreg(res4[[1]][[1]], custom.coef.map = ns, digits = 3)
# Envir behav
screenreg(res4[[2]][[1]], custom.coef.map = ns, digits = 3)


# Interaction
screenreg(res4[[1]][[2]],  digits = 3)

screenreg(res4[[2]][[2]],  digits = 3)





#######################
#### Export tables ####
#######################


### Behaviour

tex_t2 <- texreg(l = res4[[2]][[2]],
                 # override.se = lapply(res4[[2]][[4]], FUN = function(x) x[, 2]),
                 # override.pvalues = lapply(res4[[2]][[4]], FUN = function(x) x[, 4]),
                 #file="../03_Output/Mod_fe_siting.tex",
                 digits = 3, leading.zero = TRUE,
                 stars = c(0.001, 0.01, 0.05, 0.1),
                 symbol = "\\dagger",
                 label = "tab:fe_exag_heat",
                 caption = "Individual fixed effects model. Dep. var.: Pro-environmental behaviour.",
                 custom.model.names = c("(1)", "(2)", "(3)"),
                 #groups=list("Census cell level"=1:6, "Community level"=7:nb),
                 custom.coef.map = list('heat12_affect1y_past1' =  'Heatwave affected',
                                        'heat6_affect1y_past1' =  'Heatwave affected',
                                        'heat1_affect1y_past1' =  'Heatwave affected',
                                        'index1:heat12_affect1y_past1' =  'Heatwave affected $\\times$ sceptic',
                                        'index1:heat6_affect1y_past1' =  'Heatwave affected $\\times$ sceptic',
                                        'index1:heat1_affect1y_past1' =  'Heatwave affected $\\times$ sceptic'),
                 custom.note = "%stars. Two-sided test. Cluster robust standard errors in parentheses (clustered by person and LSOA).",
                 dcolumn = TRUE, caption.above = TRUE, use.packages = FALSE, include.proj.stats = FALSE
)


# Customize
tex_t2 <- gsub("D[{].[}][{].[}][{][[:digit:]]+\\.*[[:digit:]]*[}]", "D{.}{.}{2.4}", tex_t2)
tex_t2 <- gsub("\\hline", "\\hline\\\\[-1.2ex]", tex_t2, tex_t2, fixed = TRUE)
tex_t2 <- gsub("pidp", "Person-ID", tex_t2, tex_t2, fixed = TRUE)



head <- c("\\centering\n\\footnotesize\n{\\begin{threeparttable}")
tex_t2 <- gsub("\n\\begin{table}", head, tex_t2, fixed = TRUE)

head2 <- c("\\hline\\\\[-1.2ex] \n  & \\multicolumn{1}{c}{Heatwave 4 months}  & \\multicolumn{1}{c}{Heatwave 1 month}  & \\multicolumn{1}{c}{Heatwave 14 days}  \\\\ 
           \\cmidrule(lr){2-2} \\cmidrule(lr){3-3} \\cmidrule(lr){4-4}") 
tex_t2 <- sub("\\hline\\\\[-1.2ex]", head2, tex_t2, fixed = TRUE)

cont <- paste0("\\hline\\\\[-1.2ex]
  Basic controls &  \\multicolumn{1}{c}{Yes} &  \\multicolumn{1}{c}{Yes} &  \\multicolumn{1}{c}{Yes}   \\\\
  Additional controls & \\multicolumn{1}{c}{Yes} & \\multicolumn{1}{c}{Yes} & \\multicolumn{1}{c}{Yes}  \\\\
 \\hline\\\\[-1.2ex]\n","R$^2$")
tex_t2 <- sub("\\hline\\\\[-1.2ex]\nR$^2$", cont, tex_t2, fixed = TRUE)



bottom <- paste0("\n\\\\hline\\\\\\\\[-1.2ex]\n ",
                 "\\\\end{tabular}\n ",
                 "\\\\begin{tablenotes}\n ",
                 "\\\\item \\\\scriptsize{$^{***}p<0.001$, $^{**}p<0.01$, $^*p<0.05$, $^{\\\\dagger}p<0.1$. Two-sided test. Cluster robust standard errors in parentheses (clustered by person and LSOA). Basic controls: year-season of interview (omitted). Additional controls: age (5 year intervals, omitted), highest education, child(ren), marital status, household income deciles. All controls are interacted with the initial scepticism.}\n",
                 "\\\\end{tablenotes}\n",
                 "\\\\label{tab:fe_exag_heat}\n",
                 "\\\\end{center}\n",
                 "\\\\end{threeparttable}\n",
                 "}\n")
l <- gregexpr("\\hline\\\\[-1.2ex]", tex_t2, fixed = TRUE)
l <- l[[1]][length(l[[1]])]
tex_t2 <- substr(tex_t2, 1, (l-2))
tex_t2 <- sub("$", bottom, tex_t2, fixed = FALSE)

# N into middle
l1 <- gregexpr("Num. obs", tex_t2, fixed = TRUE)
l1 <- l1[[1]][length(l1[[1]])]


tmp <- substr(tex_t2, l1, (l-2))
tmp2 <- gsub("([[:digit:]]+)", "\\\\multicolumn{1}{c}{\\1}", tmp)
tex_t2 <- gsub(tmp, tmp2, tex_t2, fixed = TRUE)

# Print
write.table(tex_t2, file = "../03_Output/Table2_fe_exag_heat.tex",
            col.names = FALSE, row.names = FALSE, quote = FALSE)





### Belief

tex_t1 <- texreg(l = res4[[1]][[2]],
                 # override.se = lapply(res4[[1]][[4]], FUN = function(x) x[, 2]),
                 # override.pvalues = lapply(res4[[1]][[4]], FUN = function(x) x[, 4]),
                 #file="../03_Output/Mod_fe_siting.tex",
                 digits = 3, leading.zero = TRUE,
                 stars = c(0.001, 0.01, 0.05, 0.1),
                 symbol = "\\dagger",
                 label = "tab:fe_exag_heat2",
                 caption = "Individual fixed effects model. Dep. var.: Belief in climate change.",
                 custom.model.names = c("(1)", "(2)", "(3)"),
                 #groups=list("Census cell level"=1:6, "Community level"=7:nb),
                 custom.coef.map = list('heat12_affect1y_past1' =  'Heatwave affected',
                                        'heat6_affect1y_past1' =  'Heatwave affected',
                                        'heat1_affect1y_past1' =  'Heatwave affected',
                                        'index1:heat12_affect1y_past1' =  'Heatwave affected $\\times$ sceptic',
                                        'index1:heat6_affect1y_past1' =  'Heatwave affected $\\times$ sceptic',
                                        'index1:heat1_affect1y_past1' =  'Heatwave affected $\\times$ sceptic'),
                 custom.note = "%stars. Two-sided test. Cluster robust standard errors in parentheses (clustered by person and LSOA).",
                 dcolumn = TRUE, caption.above = TRUE, use.packages = FALSE, include.proj.stats = FALSE
)


# Customize
tex_t1 <- gsub("D[{].[}][{].[}][{][[:digit:]]+\\.*[[:digit:]]*[}]", "D{.}{.}{2.4}", tex_t1)
tex_t1 <- gsub("\\hline", "\\hline\\\\[-1.2ex]", tex_t1, tex_t1, fixed = TRUE)
tex_t1 <- gsub("pidp", "Person-ID", tex_t1, tex_t1, fixed = TRUE)



head <- c("\\centering\n\\footnotesize\n{\\begin{threeparttable}")
tex_t1 <- gsub("\n\\begin{table}", head, tex_t1, fixed = TRUE)

head2 <- c("\\hline\\\\[-1.2ex] \n  & \\multicolumn{1}{c}{Heatwave 4 months}  & \\multicolumn{1}{c}{Heatwave 1 month}  & \\multicolumn{1}{c}{Heatwave 14 days}  \\\\ 
           \\cmidrule(lr){2-2} \\cmidrule(lr){3-3} \\cmidrule(lr){4-4}") 
tex_t1 <- sub("\\hline\\\\[-1.2ex]", head2, tex_t1, fixed = TRUE)

cont <- paste0("\\hline\\\\[-1.2ex]
  Basic controls &  \\multicolumn{1}{c}{Yes} &  \\multicolumn{1}{c}{Yes} &  \\multicolumn{1}{c}{Yes}   \\\\
  Additional controls & \\multicolumn{1}{c}{Yes} & \\multicolumn{1}{c}{Yes} & \\multicolumn{1}{c}{Yes}  \\\\
 \\hline\\\\[-1.2ex]\n","R$^2$")
tex_t1 <- sub("\\hline\\\\[-1.2ex]\nR$^2$", cont, tex_t1, fixed = TRUE)



bottom <- paste0("\n\\\\hline\\\\\\\\[-1.2ex]\n ",
                 "\\\\end{tabular}\n ",
                 "\\\\begin{tablenotes}\n ",
                 "\\\\item \\\\scriptsize{$^{***}p<0.001$, $^{**}p<0.01$, $^*p<0.05$, $^{\\\\dagger}p<0.1$. Two-sided test. Cluster robust standard errors in parentheses (clustered by person and LSOA). Basic controls: year-season of interview (omitted). Additional controls: age (5 year intervals, omitted), highest education, child(ren), marital status, household income deciles. All controls are interacted with the initial scepticism.}\n",
                 "\\\\end{tablenotes}\n",
                 "\\\\label{tab:fe_exag_heat2}\n",
                 "\\\\end{center}\n",
                 "\\\\end{threeparttable}\n",
                 "}\n")
l <- gregexpr("\\hline\\\\[-1.2ex]", tex_t1, fixed = TRUE)
l <- l[[1]][length(l[[1]])]
tex_t1 <- substr(tex_t1, 1, (l-2))
tex_t1 <- sub("$", bottom, tex_t1, fixed = FALSE)

# N into middle
l1 <- gregexpr("Num. obs", tex_t1, fixed = TRUE)
l1 <- l1[[1]][length(l1[[1]])]


tmp <- substr(tex_t1, l1, (l-2))
tmp2 <- gsub("([[:digit:]]+)", "\\\\multicolumn{1}{c}{\\1}", tmp)
tex_t1 <- gsub(tmp, tmp2, tex_t1, fixed = TRUE)

# Print
write.table(tex_t1, file = "../03_Output/Table2_fe_exag_heat_att.tex",
            col.names = FALSE, row.names = FALSE, quote = FALSE)








#####################
### Plot results  ###
#####################
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}




#### Coefficients Plot Behaviour ####

# Set up df
eff <- data.frame(matrix(NA, ncol = 7, nrow = 6))
colnames(eff) <- c("Variable", "Coefficient", "SE", "t", "p", "N", "modelName")
eff[, 1] <- rep(c(3,2,1), each = 2)

eff$modelName <- rep(c(1:2), times = 3)

# Plug in values all
for(i in c(1:6)){
  eff[i, 2:6] <- c(summary(res4[[2]][[1]][[i]])$coefficients[1, 1:4], length(res4[[2]][[1]][[i]]$residuals))
}


### Combine df, make labels

allModelFrame <- data.frame(eff)

allModelFrame$Variable <- factor(allModelFrame$Variable, levels = c(1:3),
                                 labels = c("Heatwave 14 days", 
                                            "Heatwave 1 month",
                                            "Heatwave 4 months"))

allModelFrame$modelName <- factor(allModelFrame$modelName, levels = rev(c(1:2)),
                                  labels = rev(c("Low scepticism", 
                                                 "High scepticism")))


# Add floods
load("../../Climate-Change_UK/02_Data/Modelframe_exag.RData")

# Combine frames
allModelFrame$eff <- 2
mf_all1 <- rbind(allModelFrame, mf_all)

mf_all1$eff <- factor(mf_all1$eff, labels = c("Floods", "Heatwaves"))

# Exchange Variable and modelName
mn <- mf_all1$modelName
v <- mf_all1$Variable
mf_all1$Variable <- mn
mf_all1$modelName <- v

### Plot the coefficients
# Confidence intervals
interval2 <- -qnorm((1-0.95)/2)  # 95% multiplier

mf_all1$lb <- mf_all1$Coefficient - mf_all1$SE * interval2
mf_all1$ub <- mf_all1$Coefficient + mf_all1$SE * interval2

# Set limits (trim data manually)
# xlim = c(-1.5, 2.2)
# 
# mf_all1$ub[which(mf_all1$ub > mf_all1$y_max)] <- mf_all1$y_max[which(mf_all1$ub > mf_all1$y_max)]
# mf_all1$lb[which(mf_all1$lb < mf_all1$y_min)] <- mf_all1$y_min[which(mf_all1$lb < mf_all1$y_min)]

# Coef Labels
mf_all1$lab <- as.character(sprintf("%.3f", round(mf_all1$Coefficient, 3)))

mf_all1$lab[mf_all1$p <= 0.1 & mf_all1$p > 0.05] <- paste0(mf_all1$lab[mf_all1$p <= 0.1 & mf_all1$p > 0.05], expression("\u2020"))
mf_all1$lab[mf_all1$p <= 0.05 & mf_all1$p > 0.01] <- paste0(mf_all1$lab[mf_all1$p <= 0.05 & mf_all1$p > 0.01], "*")
mf_all1$lab[mf_all1$p <= 0.01 & mf_all1$p > 0.001] <- paste0(mf_all1$lab[mf_all1$p <= 0.01 & mf_all1$p > 0.001], "**")
mf_all1$lab[mf_all1$p <= 0.001] <- paste0(mf_all1$lab[mf_all1$p <= 0.001], "***")

# Number of cases position
mf_all1$maxx <- ave(mf_all1$ub, 
                    by = mf_all1$eff,
                    FUN = function(x) max(x))
mf_all1$minx <- ave(mf_all1$lb, 
                    by = mf_all1$eff,
                    FUN = function(x) min(x))

mf_all1$tx <- NA
oo <- which(mf_all1$modelName == levels(mf_all1$modelName)[3])
mf_all1$tx[oo] <- mf_all1$maxx[oo] - (mf_all1$maxx[oo] - mf_all1$minx[oo]) * 0.60
oo <- which(mf_all1$modelName == levels(mf_all1$modelName)[2])
mf_all1$tx[oo] <- mf_all1$maxx[oo] - (mf_all1$maxx[oo] - mf_all1$minx[oo]) * 0.38
oo <- which(mf_all1$modelName == levels(mf_all1$modelName)[1])
mf_all1$tx[oo] <- mf_all1$maxx[oo] - (mf_all1$maxx[oo] - mf_all1$minx[oo]) * 0.15
oo <- which(mf_all1$modelName == levels(mf_all1$modelName)[6])
mf_all1$tx[oo] <- mf_all1$maxx[oo] - (mf_all1$maxx[oo] - mf_all1$minx[oo]) * 0.60
oo <- which(mf_all1$modelName == levels(mf_all1$modelName)[5])
mf_all1$tx[oo] <- mf_all1$maxx[oo] - (mf_all1$maxx[oo] - mf_all1$minx[oo]) * 0.38
oo <- which(mf_all1$modelName == levels(mf_all1$modelName)[4])
mf_all1$tx[oo] <- mf_all1$maxx[oo] - (mf_all1$maxx[oo] - mf_all1$minx[oo]) * 0.15

mf_all1$ty <- NA
oo <- which(mf_all1$Variable %in% levels(mf_all1$Variable)[2])
mf_all1$ty[oo] <- 1.5 
oo <- which(mf_all1$Variable %in% levels(mf_all1$Variable)[1])
mf_all1$ty[oo] <- 0.5 


mf_all1$tn <- NA
oo <- which(mf_all1$modelName %in% levels(mf_all1$modelName)[c(3, 6)])
mf_all1$tn[oo] <- "N ="

oo <- which(duplicated(mf_all1[, c("modelName", "tn", "tx", "ty")]))
mf_all1$tn[oo] <- NA


# Plot
zp_all_exag <- ggplot(mf_all1, aes(colour = modelName, shape = modelName))
zp_all_exag <- zp_all_exag + facet_grid(. ~ eff, scales = "free_x")
zp_all_exag <- zp_all_exag + geom_hline(yintercept = 0, colour = scales::alpha("black", 0.3), lty = 2)
zp_all_exag <- zp_all_exag + geom_pointrange(aes(x = Variable, y = Coefficient, ymin = Coefficient - SE*interval2,
                                       ymax = Coefficient + SE*interval2),
                                   lwd = 0.7, position = position_dodge(width = 1/1.8),
                                   fill = "black", alpha = 1)
zp_all_exag <- zp_all_exag + geom_text(aes(label = lab,
                                 x = Variable, y = Coefficient), size = 6, show.legend  = FALSE, 
                             hjust = -0.3, vjust = -0.3, position = position_dodge(width = 1/1.8), alpha = 1)
zp_all_exag <- zp_all_exag + geom_text(aes(label = N, x = ty, y = tx), 
                             size = 6, show.legend  = FALSE, 
                             hjust = 0, vjust = -0.3)
zp_all_exag <- zp_all_exag + geom_text(data = unique(mf_all1[!is.na(mf_all1$tn), c("eff", "tn", "tx", "ty", "modelName")]),
                             mapping = aes(label = tn, x = ty, y = tx), size = 6, 
                             show.legend  = FALSE, colour = "black",
                             hjust = 1.3, vjust = -0.3)
zp_all_exag <- zp_all_exag + geom_blank(aes(x = 2.4, y = 0))
zp_all_exag <- zp_all_exag + coord_flip() + theme_bw()
zp_all_exag <- zp_all_exag + scale_x_discrete(expand = c(0.05,0.00) )
zp_all_exag <- zp_all_exag + theme(legend.title = element_blank())
zp_all_exag <- zp_all_exag + labs(y = "FE coefficients of past experience", x = "")
zp_all_exag <- zp_all_exag + scale_shape_manual(values = c(6, 5, 8, 15, 17, 19)) 
zp_all_exag <- zp_all_exag + scale_color_manual(values = c(gg_color_hue(6)[4:6], gg_color_hue(6)[1:3]))
zp_all_exag <- zp_all_exag + theme(text = element_text(family = "CMU Serif", size = 20),
                         # legend.position = c(0.95, 0.99), legend.justification = c(0.95,0.99),
                         plot.title = element_text(hjust = 0.5),
                         legend.position = "bottom",
                         legend.box="vertical",
                         legend.background = element_blank(),
                         legend.box.background = element_rect(colour = "black"),
                         legend.key = element_blank(),
                         axis.text.y = element_text(colour="black", size = 20),
                         axis.title.x = element_text(colour="black"),
                         axis.text.x = element_text(colour="black"),
                         strip.background = element_blank(),
                         strip.text = element_text(size = 20, colour = "black"),
                         panel.spacing.x = unit(6, "mm"),
                         plot.margin = unit(c(0,6,0,0), "mm")
                         #axis.text.x = element_text(size=16),
                         #axis.text.y = element_text(size=16),
                         #axis.title.x = element_text(size=16),
                         #axis.title.y = element_text(size=16)
)
zp_all_exag <- zp_all_exag + ggtitle("Pro-environmental behaviour")
zp_all_exag <- zp_all_exag + guides(colour = guide_legend(override.aes = list(linetype = 0, alpha = 1), reverse = T, nrow=2, byrow=T),
                          shape = guide_legend(reverse = T, nrow=2, byrow=T))

zp_all_exag <- zp_all_exag + geom_vline(aes(xintercept = 1.5))


print(zp_all_exag)

#### Combine plots ####


cairo_pdf(file = paste("../03_Output/", "Coefplot_exag.pdf", sep=""), width = 13, height = 6, 
          bg = "white", family="CMU Serif")
par(mar = c(0, 0, 0, 0))
par(mfrow = c(1, 1), oma = c(0, 0, 0, 0))
zp_all_exag
dev.off()







#### Coefficients Plot Belief ####

# Set up df
eff <- data.frame(matrix(NA, ncol = 7, nrow = 6))
colnames(eff) <- c("Variable", "Coefficient", "SE", "t", "p", "N", "modelName")
eff[, 1] <- rep(c(3,2,1), each = 2)

eff$modelName <- rep(c(1:2), times = 3)

# Plug in values all
for(i in c(1:6)){
  eff[i, 2:6] <- c(summary(res4[[1]][[1]][[i]])$coefficients[1, 1:4], length(res4[[1]][[1]][[i]]$residuals))
}


### Combine df, make labels

allModelFrame <- data.frame(eff)

allModelFrame$Variable <- factor(allModelFrame$Variable, levels = c(1:3),
                                 labels = c("Heatwave 14 days", 
                                            "Heatwave 1 month",
                                            "Heatwave 4 months"))

allModelFrame$modelName <- factor(allModelFrame$modelName, levels = rev(c(1:2)),
                                  labels = rev(c("Low scepticism", 
                                                 "High scepticism")))


# Add floods
load("../../Climate-Change_UK/02_Data/Modelframe_exag_att.RData")

# Combine frames
allModelFrame$eff <- 2
mf_all2 <- rbind(allModelFrame, mf_all)

mf_all2$eff <- factor(mf_all2$eff, labels = c("Floods", "Heatwaves"))

# Exchange Variable and modelName
mn <- mf_all2$modelName
v <- mf_all2$Variable
mf_all2$Variable <- mn
mf_all2$modelName <- v

### Plot the coefficients
# Confidence intervals
interval2 <- -qnorm((1-0.95)/2)  # 95% multiplier

mf_all2$lb <- mf_all2$Coefficient - mf_all2$SE * interval2
mf_all2$ub <- mf_all2$Coefficient + mf_all2$SE * interval2

# Set limits (trim data manually)
# xlim = c(-1.5, 2.2)
# 
# mf_all2$ub[which(mf_all2$ub > mf_all2$y_max)] <- mf_all2$y_max[which(mf_all2$ub > mf_all2$y_max)]
# mf_all2$lb[which(mf_all2$lb < mf_all2$y_min)] <- mf_all2$y_min[which(mf_all2$lb < mf_all2$y_min)]

# Coef Labels
mf_all2$lab <- as.character(sprintf("%.3f", round(mf_all2$Coefficient, 3)))

mf_all2$lab[mf_all2$p <= 0.1 & mf_all2$p > 0.05] <- paste0(mf_all2$lab[mf_all2$p <= 0.1 & mf_all2$p > 0.05], expression("\u2020"))
mf_all2$lab[mf_all2$p <= 0.05 & mf_all2$p > 0.01] <- paste0(mf_all2$lab[mf_all2$p <= 0.05 & mf_all2$p > 0.01], "*")
mf_all2$lab[mf_all2$p <= 0.01 & mf_all2$p > 0.001] <- paste0(mf_all2$lab[mf_all2$p <= 0.01 & mf_all2$p > 0.001], "**")
mf_all2$lab[mf_all2$p <= 0.001] <- paste0(mf_all2$lab[mf_all2$p <= 0.001], "***")

# Number of cases position
mf_all2$maxx <- ave(mf_all2$ub, 
                    by = mf_all2$eff,
                    FUN = function(x) max(x))
mf_all2$minx <- ave(mf_all2$lb, 
                    by = mf_all2$eff,
                    FUN = function(x) min(x))

mf_all2$tx <- NA
oo <- which(mf_all2$modelName == levels(mf_all2$modelName)[3])
mf_all2$tx[oo] <- mf_all2$maxx[oo] - (mf_all2$maxx[oo] - mf_all2$minx[oo]) * 0.60
oo <- which(mf_all2$modelName == levels(mf_all2$modelName)[2])
mf_all2$tx[oo] <- mf_all2$maxx[oo] - (mf_all2$maxx[oo] - mf_all2$minx[oo]) * 0.38
oo <- which(mf_all2$modelName == levels(mf_all2$modelName)[1])
mf_all2$tx[oo] <- mf_all2$maxx[oo] - (mf_all2$maxx[oo] - mf_all2$minx[oo]) * 0.15
oo <- which(mf_all2$modelName == levels(mf_all2$modelName)[6])
mf_all2$tx[oo] <- mf_all2$maxx[oo] - (mf_all2$maxx[oo] - mf_all2$minx[oo]) * 0.60
oo <- which(mf_all2$modelName == levels(mf_all2$modelName)[5])
mf_all2$tx[oo] <- mf_all2$maxx[oo] - (mf_all2$maxx[oo] - mf_all2$minx[oo]) * 0.38
oo <- which(mf_all2$modelName == levels(mf_all2$modelName)[4])
mf_all2$tx[oo] <- mf_all2$maxx[oo] - (mf_all2$maxx[oo] - mf_all2$minx[oo]) * 0.15

mf_all2$ty <- NA
oo <- which(mf_all2$Variable %in% levels(mf_all2$Variable)[2])
mf_all2$ty[oo] <- 1.5 
oo <- which(mf_all2$Variable %in% levels(mf_all2$Variable)[1])
mf_all2$ty[oo] <- 0.5 


mf_all2$tn <- NA
oo <- which(mf_all2$modelName %in% levels(mf_all2$modelName)[c(3, 6)])
mf_all2$tn[oo] <- "N ="

oo <- which(duplicated(mf_all2[, c("modelName", "tn", "tx", "ty")]))
mf_all2$tn[oo] <- NA


# Plot
zp_all_exag2 <- ggplot(mf_all2, aes(colour = modelName, shape = modelName))
zp_all_exag2 <- zp_all_exag2 + facet_grid(. ~ eff, scales = "free_x")
zp_all_exag2 <- zp_all_exag2 + geom_hline(yintercept = 0, colour = scales::alpha("black", 0.3), lty = 2)
zp_all_exag2 <- zp_all_exag2 + geom_pointrange(aes(x = Variable, y = Coefficient, ymin = Coefficient - SE*interval2,
                                         ymax = Coefficient + SE*interval2),
                                     lwd = 0.7, position = position_dodge(width = 1/1.8),
                                     fill = "black", alpha = 1)
zp_all_exag2 <- zp_all_exag2 + geom_text(aes(label = lab,
                                   x = Variable, y = Coefficient), size = 6, show.legend  = FALSE, 
                               hjust = -0.3, vjust = -0.3, position = position_dodge(width = 1/1.8), alpha = 1)
zp_all_exag2 <- zp_all_exag2 + geom_text(aes(label = N, x = ty, y = tx), 
                               size = 6, show.legend  = FALSE, 
                               hjust = 0, vjust = -0.3)
zp_all_exag2 <- zp_all_exag2 + geom_text(data = unique(mf_all2[!is.na(mf_all2$tn), c("eff", "tn", "tx", "ty", "modelName")]),
                               mapping = aes(label = tn, x = ty, y = tx), size = 6, 
                               show.legend  = FALSE, colour = "black",
                               hjust = 1.3, vjust = -0.3)
zp_all_exag2 <- zp_all_exag2 + geom_blank(aes(x = 2.4, y = 0))
zp_all_exag2 <- zp_all_exag2 + coord_flip() + theme_bw()
zp_all_exag2 <- zp_all_exag2 + scale_x_discrete(expand = c(0.05,0.00) )
zp_all_exag2 <- zp_all_exag2 + theme(legend.title = element_blank())
zp_all_exag2 <- zp_all_exag2 + labs(y = "FE coefficients of past experience", x = "")
zp_all_exag2 <- zp_all_exag2 + scale_shape_manual(values = c(6, 5, 8, 15, 17, 19)) 
zp_all_exag2 <- zp_all_exag2 + scale_color_manual(values = c(gg_color_hue(6)[4:6], gg_color_hue(6)[1:3]))
zp_all_exag2 <- zp_all_exag2 + theme(text = element_text(family = "CMU Serif", size = 20),
                           # legend.position = c(0.95, 0.99), legend.justification = c(0.95,0.99),
                           plot.title = element_text(hjust = 0.5),
                           legend.position = "bottom",
                           legend.box="vertical",
                           legend.background = element_blank(),
                           legend.box.background = element_rect(colour = "black"),
                           legend.key = element_blank(),
                           axis.text.y = element_text(colour="black", size = 20),
                           axis.title.x = element_text(colour="black"),
                           axis.text.x = element_text(colour="black"),
                           strip.background = element_blank(),
                           strip.text = element_text(size = 20, colour = "black"),
                           panel.spacing.x = unit(6, "mm"),
                           plot.margin = unit(c(0,6,0,0), "mm")
                           #axis.text.x = element_text(size=16),
                           #axis.text.y = element_text(size=16),
                           #axis.title.x = element_text(size=16),
                           #axis.title.y = element_text(size=16)
)
zp_all_exag2 <- zp_all_exag2 + ggtitle("Climate change belief")
zp_all_exag2 <- zp_all_exag2 + guides(colour = guide_legend(override.aes = list(linetype = 0, alpha = 1), reverse = T, nrow=2, byrow=T),
                            shape = guide_legend(reverse = T, nrow=2, byrow=T))

zp_all_exag2 <- zp_all_exag2 + geom_vline(aes(xintercept = 1.5))


print(zp_all_exag2)

#### Combine plots ####


cairo_pdf(file = paste("../03_Output/", "Coefplot_exag_att.pdf", sep=""), width = 13, height = 6, 
          bg = "white", family="CMU Serif")
par(mar = c(0, 0, 0, 0))
par(mfrow = c(1, 1), oma = c(0, 0, 0, 0))
zp_all_exag2
dev.off()


g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

mylegend <- g_legend(zp_all_exag2)

cairo_pdf(file = paste("../03_Output/", "Coefplot_exag_comb.pdf", sep=""), width = 13, height = 12, 
          bg = "white", family="CMU Serif")
par(mar = c(0, 0, 0, 0))
par(mfrow = c(1, 1), oma = c(0, 0, 0, 0))
grid.arrange(arrangeGrob(zp_all_exag2 + theme(legend.position="none", plot.margin = unit(c(0,6,4,0), "mm")) + labs(y = "", x = ""),
                         zp_all_exag + theme(legend.position="none"),
                         nrow = 2),
             arrangeGrob(ggplot()+geom_blank(aes(1,1)) + cowplot::theme_nothing(), mylegend, ncol = 2, widths = c(1, 10)), 
             nrow = 2, heights = c(10, 1))
dev.off()






### ----- Redo with reduced sample ----- ###


### Reduce to people with zero in first wave

bhps_heat.df <- bhps_heat.df[order(bhps_heat.df$pidp, bhps_heat.df$wave), ]
bhps_heat.df$first_belief <- ave(bhps_heat.df$climatechange30,
                                 bhps_heat.df$pidp,
                                 FUN = function(x) x[1])


### Repeat

depvar <- c("climatechange30")
intens <- c("12", "6", "1")
buffer <- c(1, 2, 3)
b <- buffer[1]
plm.l <- vector("list", nrow(expand.grid(intens)) * 1)
plm3.l <- vector("list", nrow(expand.grid(buffer)) * 1) 

res4 <- vector("list", length(depvar)) 

l <- 1
for(d in depvar) {
  library(lfe)
  names(plm.l) <- paste0(d, c(1:length(plm.l)))
  names(plm3.l) <- paste0(d, c(1:length(plm3.l)))
  
  j <- 1
  k <- 1
  for(i in intens){
    
    v1 <- paste0("heat", i, "_affect1y_past", b)
    ev1 <- paste0("heat_ever", b)
    
    fm2 <- paste0(d, " ~ ", v1, 
                  "+ age_cat + hiqual_dv + child + marstat_dv + hhinc_dec")
    
    # Estimate via felm with twoway clustered SEs
    bhps_heat.df$interaction <- as.numeric(as.factor(bhps_heat.df$season_index)) * bhps_heat.df[, ev1]
    
    plm2_1 <- felm(formula(paste0(fm2, " | pidp + season_index + interaction | 0 | pidp + lsoa01")),
                   data = bhps_heat.df[which(bhps_heat.df[, paste0("samplefe_", d)] == 1 & bhps_heat.df[, "index1"] == 0 & bhps_heat.df$first_belief == 0), ],
                   psdef = FALSE)
    
    plm2_2 <- felm(formula(paste0(fm2, " | pidp + season_index + interaction | 0 | pidp + lsoa01")),
                   data = bhps_heat.df[which(bhps_heat.df[, paste0("samplefe_", d)] == 1 & bhps_heat.df[, "index1"] == 1 & bhps_heat.df$first_belief == 0), ],
                   psdef = FALSE)
    
    
    
    # Interacted model
    fm3 <- paste0(d, " ~ index1 * (", v1, 
                  "+ age_cat + hiqual_dv + child + marstat_dv + hhinc_dec", 
                  "  + as.factor(season_index)*", ev1, ")")
    
    plm3 <- felm(formula(paste0(fm3, " | pidp | 0 | pidp + lsoa01")),
                 data = bhps_heat.df[which(bhps_heat.df[, paste0("samplefe_", d)] == 1 & bhps_heat.df$first_belief == 0), ], 
                 psdef = FALSE)
    
    
    plm.l[[j]] <- plm2_1
    plm.l[[j+1]] <- plm2_2
    plm3.l[[k]] <- plm3
    j <- j + 2
    k <- k + 1
    
  }
  res4[[l]] <- list(plm.l,  plm3.l)
  l <- l + 1
}






#### Coefficients Plot Belief ####

# Set up df
eff <- data.frame(matrix(NA, ncol = 7, nrow = 6))
colnames(eff) <- c("Variable", "Coefficient", "SE", "t", "p", "N", "modelName")
eff[, 1] <- rep(c(3,2,1), each = 2)

eff$modelName <- rep(c(1:2), times = 3)

# Plug in values all
for(i in c(1:6)){
  eff[i, 2:6] <- c(summary(res4[[1]][[1]][[i]])$coefficients[1, 1:4], length(res4[[1]][[1]][[i]]$residuals))
}


### Combine df, make labels

allModelFrame <- data.frame(eff)

allModelFrame$Variable <- factor(allModelFrame$Variable, levels = c(1:3),
                                 labels = c("Heatwave 14 days", 
                                            "Heatwave 1 month",
                                            "Heatwave 4 months"))

allModelFrame$modelName <- factor(allModelFrame$modelName, levels = rev(c(1:2)),
                                  labels = rev(c("Low scepticism", 
                                                 "High scepticism")))


# Add floods
load("../../Climate-Change_UK/02_Data/Subset0_Modelframe_exag_att.RData")

# Combine frames
allModelFrame$eff <- 2
mf_all2 <- rbind(allModelFrame, mf_all)

mf_all2$eff <- factor(mf_all2$eff, labels = c("Floods", "Heatwaves"))

# Exchange Variable and modelName
mn <- mf_all2$modelName
v <- mf_all2$Variable
mf_all2$Variable <- mn
mf_all2$modelName <- v

### Plot the coefficients
# Confidence intervals
interval2 <- -qnorm((1-0.95)/2)  # 95% multiplier

mf_all2$lb <- mf_all2$Coefficient - mf_all2$SE * interval2
mf_all2$ub <- mf_all2$Coefficient + mf_all2$SE * interval2

# Set limits (trim data manually)
# xlim = c(-1.5, 2.2)
# 
# mf_all2$ub[which(mf_all2$ub > mf_all2$y_max)] <- mf_all2$y_max[which(mf_all2$ub > mf_all2$y_max)]
# mf_all2$lb[which(mf_all2$lb < mf_all2$y_min)] <- mf_all2$y_min[which(mf_all2$lb < mf_all2$y_min)]

# Coef Labels
mf_all2$lab <- as.character(sprintf("%.3f", round(mf_all2$Coefficient, 3)))

mf_all2$lab[mf_all2$p <= 0.1 & mf_all2$p > 0.05] <- paste0(mf_all2$lab[mf_all2$p <= 0.1 & mf_all2$p > 0.05], expression("\u2020"))
mf_all2$lab[mf_all2$p <= 0.05 & mf_all2$p > 0.01] <- paste0(mf_all2$lab[mf_all2$p <= 0.05 & mf_all2$p > 0.01], "*")
mf_all2$lab[mf_all2$p <= 0.01 & mf_all2$p > 0.001] <- paste0(mf_all2$lab[mf_all2$p <= 0.01 & mf_all2$p > 0.001], "**")
mf_all2$lab[mf_all2$p <= 0.001] <- paste0(mf_all2$lab[mf_all2$p <= 0.001], "***")

# Number of cases position
mf_all2$maxx <- ave(mf_all2$ub, 
                    by = mf_all2$eff,
                    FUN = function(x) max(x))
mf_all2$minx <- ave(mf_all2$lb, 
                    by = mf_all2$eff,
                    FUN = function(x) min(x))

mf_all2$tx <- NA
oo <- which(mf_all2$modelName == levels(mf_all2$modelName)[3])
mf_all2$tx[oo] <- mf_all2$maxx[oo] - (mf_all2$maxx[oo] - mf_all2$minx[oo]) * 0.60
oo <- which(mf_all2$modelName == levels(mf_all2$modelName)[2])
mf_all2$tx[oo] <- mf_all2$maxx[oo] - (mf_all2$maxx[oo] - mf_all2$minx[oo]) * 0.38
oo <- which(mf_all2$modelName == levels(mf_all2$modelName)[1])
mf_all2$tx[oo] <- mf_all2$maxx[oo] - (mf_all2$maxx[oo] - mf_all2$minx[oo]) * 0.15
oo <- which(mf_all2$modelName == levels(mf_all2$modelName)[6])
mf_all2$tx[oo] <- mf_all2$maxx[oo] - (mf_all2$maxx[oo] - mf_all2$minx[oo]) * 0.60
oo <- which(mf_all2$modelName == levels(mf_all2$modelName)[5])
mf_all2$tx[oo] <- mf_all2$maxx[oo] - (mf_all2$maxx[oo] - mf_all2$minx[oo]) * 0.38
oo <- which(mf_all2$modelName == levels(mf_all2$modelName)[4])
mf_all2$tx[oo] <- mf_all2$maxx[oo] - (mf_all2$maxx[oo] - mf_all2$minx[oo]) * 0.15

mf_all2$ty <- NA
oo <- which(mf_all2$Variable %in% levels(mf_all2$Variable)[2])
mf_all2$ty[oo] <- 1.5 
oo <- which(mf_all2$Variable %in% levels(mf_all2$Variable)[1])
mf_all2$ty[oo] <- 0.5 


mf_all2$tn <- NA
oo <- which(mf_all2$modelName %in% levels(mf_all2$modelName)[c(3, 6)])
mf_all2$tn[oo] <- "N ="

oo <- which(duplicated(mf_all2[, c("modelName", "tn", "tx", "ty")]))
mf_all2$tn[oo] <- NA


# Plot
zp_sub0_exag2 <- ggplot(mf_all2, aes(colour = modelName, shape = modelName))
zp_sub0_exag2 <- zp_sub0_exag2 + facet_grid(. ~ eff, scales = "free_x")
zp_sub0_exag2 <- zp_sub0_exag2 + geom_hline(yintercept = 0, colour = scales::alpha("black", 0.3), lty = 2)
zp_sub0_exag2 <- zp_sub0_exag2 + geom_pointrange(aes(x = Variable, y = Coefficient, ymin = Coefficient - SE*interval2,
                                                     ymax = Coefficient + SE*interval2),
                                                 lwd = 0.7, position = position_dodge(width = 1/1.8),
                                                 fill = "black", alpha = 1)
zp_sub0_exag2 <- zp_sub0_exag2 + geom_text(aes(label = lab,
                                               x = Variable, y = Coefficient), size = 6, show.legend  = FALSE, 
                                           hjust = -0.3, vjust = -0.3, position = position_dodge(width = 1/1.8), alpha = 1)
zp_sub0_exag2 <- zp_sub0_exag2 + geom_text(aes(label = N, x = ty, y = tx), 
                                           size = 6, show.legend  = FALSE, 
                                           hjust = 0, vjust = -0.3)
zp_sub0_exag2 <- zp_sub0_exag2 + geom_text(data = unique(mf_all2[!is.na(mf_all2$tn), c("eff", "tn", "tx", "ty", "modelName")]),
                                           mapping = aes(label = tn, x = ty, y = tx), size = 6, 
                                           show.legend  = FALSE, colour = "black",
                                           hjust = 1.3, vjust = -0.3)
zp_sub0_exag2 <- zp_sub0_exag2 + geom_blank(aes(x = 2.4, y = 0))
zp_sub0_exag2 <- zp_sub0_exag2 + coord_flip() + theme_bw()
zp_sub0_exag2 <- zp_sub0_exag2 + scale_x_discrete(expand = c(0.05,0.00) )
zp_sub0_exag2 <- zp_sub0_exag2 + theme(legend.title = element_blank())
zp_sub0_exag2 <- zp_sub0_exag2 + labs(y = "FE coefficients of past experience", x = "")
zp_sub0_exag2 <- zp_sub0_exag2 + scale_shape_manual(values = c(6, 5, 8, 15, 17, 19)) 
zp_sub0_exag2 <- zp_sub0_exag2 + scale_color_manual(values = c(gg_color_hue(6)[4:6], gg_color_hue(6)[1:3]))
zp_sub0_exag2 <- zp_sub0_exag2 + theme(text = element_text(family = "CMU Serif", size = 20),
                                       # legend.position = c(0.95, 0.99), legend.justification = c(0.95,0.99),
                                       plot.title = element_text(hjust = 0.5),
                                       legend.position = "bottom",
                                       legend.box="vertical",
                                       legend.background = element_blank(),
                                       legend.box.background = element_rect(colour = "black"),
                                       legend.key = element_blank(),
                                       axis.text.y = element_text(colour="black", size = 20),
                                       axis.title.x = element_text(colour="black"),
                                       axis.text.x = element_text(colour="black"),
                                       strip.background = element_blank(),
                                       strip.text = element_text(size = 20, colour = "black"),
                                       panel.spacing.x = unit(6, "mm"),
                                       plot.margin = unit(c(0,6,0,0), "mm")
                                       #axis.text.x = element_text(size=16),
                                       #axis.text.y = element_text(size=16),
                                       #axis.title.x = element_text(size=16),
                                       #axis.title.y = element_text(size=16)
)
zp_sub0_exag2 <- zp_sub0_exag2 + ggtitle("Climate change concern")
zp_sub0_exag2 <- zp_sub0_exag2 + guides(colour = guide_legend(override.aes = list(linetype = 0, alpha = 1), reverse = T, nrow=2, byrow=T),
                                        shape = guide_legend(reverse = T, nrow=2, byrow=T))

zp_sub0_exag2 <- zp_sub0_exag2 + geom_vline(aes(xintercept = 1.5))


print(zp_sub0_exag2)

#### Combine plots ####


cairo_pdf(file = paste("../03_Output/", "Subset0_Coefplot_exag_att.pdf", sep=""), width = 13, height = 6, 
          bg = "white", family="CMU Serif")
par(mar = c(0, 0, 0, 0))
par(mfrow = c(1, 1), oma = c(0, 0, 0, 0))
zp_sub0_exag2
dev.off()


g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

mylegend <- g_legend(zp_sub0_exag2)

# cairo_pdf(file = paste("../03_Output/", "Subset0_Coefplot_exag_comb.pdf", sep=""), width = 13, height = 12, 
#           bg = "white", family="CMU Serif")
# par(mar = c(0, 0, 0, 0))
# par(mfrow = c(1, 1), oma = c(0, 0, 0, 0))
# grid.arrange(arrangeGrob(zp_sub0_exag2 + theme(legend.position="none", plot.margin = unit(c(0,6,4,0), "mm")) + labs(y = "", x = ""),
#                          zp_all_exag + theme(legend.position="none"),
#                          nrow = 2),
#              arrangeGrob(ggplot()+geom_blank(aes(1,1)) + cowplot::theme_nothing(), mylegend, ncol = 2, widths = c(1, 10)), 
#              nrow = 2, heights = c(10, 1))
# dev.off()





















###------------------------------------------------------###
###         Effect heterogeneity  - partisanship         ###
###------------------------------------------------------###

bhps_heat.df$index2 <- ifelse(bhps_heat.df$c_partisan %in% c(1, 7, 10, 12, 13, 14), 1, 0)


bhps_heat.df$index2 <- NA
oo <- which(bhps_heat.df$c_partisan %in% c(2, 3, 4, 5, 6, 8, 9, 11))
bhps_heat.df$index2[oo] <- 1
oo <- which(bhps_heat.df$c_partisan %in% c(95:99))
bhps_heat.df$index2[oo] <- 2
oo <- which(bhps_heat.df$c_partisan %in% c(1, 7, 10, 12, 13, 14))
bhps_heat.df$index2[oo] <- 3

bhps_heat.df$index2 <- as.factor(bhps_heat.df$index2)


### parti ###

depvar <- c("climatechange30", "env_index_all2")
intens <- c("12", "6", "1")
buffer <- c(1, 2, 3)
b <- buffer[1]
plm.l <- vector("list", nrow(expand.grid(intens)) * 1)
plm3.l <- vector("list", nrow(expand.grid(buffer)) * 1) 

res4 <- vector("list", length(depvar)) 

l <- 1
for(d in depvar) {
  library(lfe)
  names(plm.l) <- paste0(d, c(1:length(plm.l)))
  names(plm3.l) <- paste0(d, c(1:length(plm3.l)))
  
  j <- 1
  k <- 1
  for(i in intens){
    
    v1 <- paste0("heat", i, "_affect1y_past", b)
    ev1 <- paste0("heat_ever", b)
    
    fm2 <- paste0(d, " ~ ", v1, 
                  "+ age_cat + hiqual_dv + child + marstat_dv + hhinc_dec")
    
    # Estimate via felm with twoway clustered SEs
    bhps_heat.df$interaction <- as.numeric(as.factor(bhps_heat.df$season_index)) * bhps_heat.df[, ev1]
    
    plm2_1 <- felm(formula(paste0(fm2, " | pidp + season_index + interaction | 0 | pidp + lsoa01")),
                   data = bhps_heat.df[which(bhps_heat.df[, paste0("samplefe_", d)] == 1 &  bhps_heat.df[, "index2"] == 1), ],
                   psdef = FALSE)

    plm2_2 <- felm(formula(paste0(fm2, " | pidp + season_index + interaction | 0 | pidp + lsoa01")),
                   data = bhps_heat.df[which(bhps_heat.df[, paste0("samplefe_", d)] == 1 &  bhps_heat.df[, "index2"] == 2), ],
                   psdef = FALSE)
    
    plm2_3 <- felm(formula(paste0(fm2, " | pidp + season_index + interaction | 0 | pidp + lsoa01")),
                   data = bhps_heat.df[which(bhps_heat.df[, paste0("samplefe_", d)] == 1 &  bhps_heat.df[, "index2"] == 3), ],
                   psdef = FALSE)
    
    
    
    # Interacted model
    fm3 <- paste0(d, " ~ index2 * (", v1, 
                  "+ age_cat + hiqual_dv + child + marstat_dv + hhinc_dec", 
                  "  + as.factor(season_index)*", ev1, ")")
    
    plm3 <- felm(formula(paste0(fm3, " | pidp | 0 | pidp + lsoa01")),
                 data = bhps_heat.df[which(bhps_heat.df[, paste0("samplefe_", d)] == 1  ), ], 
                 psdef = FALSE)
    
    
    plm.l[[j]] <- plm2_1
    plm.l[[j+1]] <- plm2_2
    plm.l[[j+2]] <- plm2_3
    plm3.l[[k]] <- plm3
    j <- j + 3
    k <- k + 1
    
  }
  res4[[l]] <- list(plm.l,  plm3.l)
  l <- l + 1
}



# Map
ns <- as.list(paste0("heat", intens, "_affect1y_past", b))
names(ns) <- paste0("heat", intens, "_affect1y_past", b)


# Climate change
screenreg(res4[[1]][[1]], custom.coef.map = ns, digits = 3)
# Envir behav
screenreg(res4[[2]][[1]], custom.coef.map = ns, digits = 3)


# Interaction
screenreg(res4[[1]][[2]],  digits = 3)

screenreg(res4[[2]][[2]],  digits = 3)





#######################
#### Export tables ####
#######################


### Behaviour

tex_t2 <- texreg(l = res4[[2]][[2]],
                 # override.se = lapply(res4[[2]][[4]], FUN = function(x) x[, 2]),
                 # override.pvalues = lapply(res4[[2]][[4]], FUN = function(x) x[, 4]),
                 #file="../03_Output/Mod_fe_siting.tex",
                 digits = 3, leading.zero = TRUE,
                 stars = c(0.001, 0.01, 0.05, 0.1),
                 symbol = "\\dagger",
                 label = "tab:fe_parti_heat",
                 caption = "Individual fixed effects model. Dep. var.: Pro-environmental behaviour.",
                 custom.model.names = c("(1)", "(2)", "(3)"),
                 #groups=list("Census cell level"=1:6, "Community level"=7:nb),
                 custom.coef.map = list('heat12_affect1y_past1' =  'Heatwave affected',
                                        'heat6_affect1y_past1' =  'Heatwave affected',
                                        'heat1_affect1y_past1' =  'Heatwave affected',
                                        'index22:heat12_affect1y_past1' =  'Heatwave affected $\\times$ no partisan',
                                        'index22:heat6_affect1y_past1' =  'Heatwave affected $\\times$ no partisan',
                                        'index22:heat1_affect1y_past1' =  'Heatwave affected $\\times$ no partisan',
                                        'index23:heat12_affect1y_past1' =  'Heatwave affected $\\times$ right partisan',
                                        'index23:heat6_affect1y_past1' =  'Heatwave affected $\\times$ right partisan',
                                        'index23:heat1_affect1y_past1' =  'Heatwave affected $\\times$ right partisan'),
                 custom.note = "%stars. Two-sided test. Cluster robust standard errors in parentheses (clustered by person and LSOA).",
                 dcolumn = TRUE, caption.above = TRUE, use.packages = FALSE, include.proj.stats = FALSE
)


# Customize
tex_t2 <- gsub("D[{].[}][{].[}][{][[:digit:]]+\\.*[[:digit:]]*[}]", "D{.}{.}{2.4}", tex_t2)
tex_t2 <- gsub("\\hline", "\\hline\\\\[-1.2ex]", tex_t2, tex_t2, fixed = TRUE)
tex_t2 <- gsub("pidp", "Person-ID", tex_t2, tex_t2, fixed = TRUE)



head <- c("\\centering\n\\footnotesize\n{\\begin{threeparttable}")
tex_t2 <- gsub("\n\\begin{table}", head, tex_t2, fixed = TRUE)

head2 <- c("\\hline\\\\[-1.2ex] \n  & \\multicolumn{1}{c}{Heatwave 4 months}  & \\multicolumn{1}{c}{Heatwave 1 month}  & \\multicolumn{1}{c}{Heatwave 14 days}  \\\\ 
           \\cmidrule(lr){2-2} \\cmidrule(lr){3-3} \\cmidrule(lr){4-4}") 
tex_t2 <- sub("\\hline\\\\[-1.2ex]", head2, tex_t2, fixed = TRUE)

cont <- paste0("\\hline\\\\[-1.2ex]
  Basic controls &  \\multicolumn{1}{c}{Yes} &  \\multicolumn{1}{c}{Yes} &  \\multicolumn{1}{c}{Yes}   \\\\
  Additional controls & \\multicolumn{1}{c}{Yes} & \\multicolumn{1}{c}{Yes} & \\multicolumn{1}{c}{Yes}  \\\\
 \\hline\\\\[-1.2ex]\n","R$^2$")
tex_t2 <- sub("\\hline\\\\[-1.2ex]\nR$^2$", cont, tex_t2, fixed = TRUE)



bottom <- paste0("\n\\\\hline\\\\\\\\[-1.2ex]\n ",
                 "\\\\end{tabular}\n ",
                 "\\\\begin{tablenotes}\n ",
                 "\\\\item \\\\scriptsize{$^{***}p<0.001$, $^{**}p<0.01$, $^*p<0.05$, $^{\\\\dagger}p<0.1$. Two-sided test. Cluster robust standard errors in parentheses (clustered by person and LSOA). Basic controls: year-season of interview (omitted). Additional controls: age (5 year intervals, omitted), highest education, child(ren), marital status, household income deciles. All controls are interacted with the initial partisanship.}\n",
                 "\\\\end{tablenotes}\n",
                 "\\\\label{tab:fe_parti_heat}\n",
                 "\\\\end{center}\n",
                 "\\\\end{threeparttable}\n",
                 "}\n")
l <- gregexpr("\\hline\\\\[-1.2ex]", tex_t2, fixed = TRUE)
l <- l[[1]][length(l[[1]])]
tex_t2 <- substr(tex_t2, 1, (l-2))
tex_t2 <- sub("$", bottom, tex_t2, fixed = FALSE)

# N into middle
l1 <- gregexpr("Num. obs", tex_t2, fixed = TRUE)
l1 <- l1[[1]][length(l1[[1]])]


tmp <- substr(tex_t2, l1, (l-2))
tmp2 <- gsub("([[:digit:]]+)", "\\\\multicolumn{1}{c}{\\1}", tmp)
tex_t2 <- gsub(tmp, tmp2, tex_t2, fixed = TRUE)

# Print
write.table(tex_t2, file = "../03_Output/Table2_fe_parti_heat.tex",
            col.names = FALSE, row.names = FALSE, quote = FALSE)





### Belief

tex_t1 <- texreg(l = res4[[1]][[2]],
                 # override.se = lapply(res4[[1]][[4]], FUN = function(x) x[, 2]),
                 # override.pvalues = lapply(res4[[1]][[4]], FUN = function(x) x[, 4]),
                 #file="../03_Output/Mod_fe_siting.tex",
                 digits = 3, leading.zero = TRUE,
                 stars = c(0.001, 0.01, 0.05, 0.1),
                 symbol = "\\dagger",
                 label = "tab:fe_parti_heat2",
                 caption = "Individual fixed effects model. Dep. var.: Belief in climate change.",
                 custom.model.names = c("(1)", "(2)", "(3)"),
                 #groups=list("Census cell level"=1:6, "Community level"=7:nb),
                 custom.coef.map = list('heat12_affect1y_past1' =  'Heatwave affected',
                                        'heat6_affect1y_past1' =  'Heatwave affected',
                                        'heat1_affect1y_past1' =  'Heatwave affected',
                                        'index22:heat12_affect1y_past1' =  'Heatwave affected $\\times$ no partisan',
                                        'index22:heat6_affect1y_past1' =  'Heatwave affected $\\times$ no partisan',
                                        'index22:heat1_affect1y_past1' =  'Heatwave affected $\\times$ no partisan',
                                        'index23:heat12_affect1y_past1' =  'Heatwave affected $\\times$ right partisan',
                                        'index23:heat6_affect1y_past1' =  'Heatwave affected $\\times$ right partisan',
                                        'index23:heat1_affect1y_past1' =  'Heatwave affected $\\times$ right partisan'),
                 custom.note = "%stars. Two-sided test. Cluster robust standard errors in parentheses (clustered by person and LSOA).",
                 dcolumn = TRUE, caption.above = TRUE, use.packages = FALSE, include.proj.stats = FALSE
)


# Customize
tex_t1 <- gsub("D[{].[}][{].[}][{][[:digit:]]+\\.*[[:digit:]]*[}]", "D{.}{.}{2.4}", tex_t1)
tex_t1 <- gsub("\\hline", "\\hline\\\\[-1.2ex]", tex_t1, tex_t1, fixed = TRUE)
tex_t1 <- gsub("pidp", "Person-ID", tex_t1, tex_t1, fixed = TRUE)



head <- c("\\centering\n\\footnotesize\n{\\begin{threeparttable}")
tex_t1 <- gsub("\n\\begin{table}", head, tex_t1, fixed = TRUE)

head2 <- c("\\hline\\\\[-1.2ex] \n  & \\multicolumn{1}{c}{Heatwave 4 months}  & \\multicolumn{1}{c}{Heatwave 1 month}  & \\multicolumn{1}{c}{Heatwave 14 days}  \\\\ 
           \\cmidrule(lr){2-2} \\cmidrule(lr){3-3} \\cmidrule(lr){4-4}") 
tex_t1 <- sub("\\hline\\\\[-1.2ex]", head2, tex_t1, fixed = TRUE)

cont <- paste0("\\hline\\\\[-1.2ex]
  Basic controls &  \\multicolumn{1}{c}{Yes} &  \\multicolumn{1}{c}{Yes} &  \\multicolumn{1}{c}{Yes}   \\\\
  Additional controls & \\multicolumn{1}{c}{Yes} & \\multicolumn{1}{c}{Yes} & \\multicolumn{1}{c}{Yes}  \\\\
 \\hline\\\\[-1.2ex]\n","R$^2$")
tex_t1 <- sub("\\hline\\\\[-1.2ex]\nR$^2$", cont, tex_t1, fixed = TRUE)



bottom <- paste0("\n\\\\hline\\\\\\\\[-1.2ex]\n ",
                 "\\\\end{tabular}\n ",
                 "\\\\begin{tablenotes}\n ",
                 "\\\\item \\\\scriptsize{$^{***}p<0.001$, $^{**}p<0.01$, $^*p<0.05$, $^{\\\\dagger}p<0.1$. Two-sided test. Cluster robust standard errors in parentheses (clustered by person and LSOA). Basic controls: year-season of interview (omitted). Additional controls: age (5 year intervals, omitted), highest education, child(ren), marital status, household income deciles. All controls are interacted with the initial partisanship.}\n",
                 "\\\\end{tablenotes}\n",
                 "\\\\label{tab:fe_parti_heat2}\n",
                 "\\\\end{center}\n",
                 "\\\\end{threeparttable}\n",
                 "}\n")
l <- gregexpr("\\hline\\\\[-1.2ex]", tex_t1, fixed = TRUE)
l <- l[[1]][length(l[[1]])]
tex_t1 <- substr(tex_t1, 1, (l-2))
tex_t1 <- sub("$", bottom, tex_t1, fixed = FALSE)

# N into middle
l1 <- gregexpr("Num. obs", tex_t1, fixed = TRUE)
l1 <- l1[[1]][length(l1[[1]])]


tmp <- substr(tex_t1, l1, (l-2))
tmp2 <- gsub("([[:digit:]]+)", "\\\\multicolumn{1}{c}{\\1}", tmp)
tex_t1 <- gsub(tmp, tmp2, tex_t1, fixed = TRUE)

# Print
write.table(tex_t1, file = "../03_Output/Table2_fe_parti_heat_att.tex",
            col.names = FALSE, row.names = FALSE, quote = FALSE)








#####################
### Plot results  ###
#####################
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}




#### Coefficients Plot Behaviour ####

# Set up df
eff <- data.frame(matrix(NA, ncol = 7, nrow = 9))
colnames(eff) <- c("Variable", "Coefficient", "SE", "t", "p", "N", "modelName")
eff[, 1] <- rep(c(3,2,1), each = 3)

eff$modelName <- rep(c(1:3), times = 3)

# Plug in values all
for(i in c(1:9)){
  eff[i, 2:6] <- c(summary(res4[[2]][[1]][[i]])$coefficients[1, 1:4], length(res4[[2]][[1]][[i]]$residuals))
}


### Combine df, make labels

allModelFrame <- data.frame(eff)

allModelFrame$Variable <- factor(allModelFrame$Variable, levels = c(1:3),
                                 labels = c("Heatwave 14 days", 
                                            "Heatwave 1 month",
                                            "Heatwave 4 months"))

allModelFrame$modelName <- factor(allModelFrame$modelName, levels = rev(c(1:3)),
                                  labels = rev(c("Left-leaning \n partisanship",
                                                 "No partisanship",
                                                 "Right-leaning \n partisanship")))


# Add floods
load("../../Climate-Change_UK/02_Data/Modelframe_parti.RData")

# Combine frames
allModelFrame$eff <- 2
mf_all1 <- rbind(allModelFrame, mf_all)

mf_all1$eff <- factor(mf_all1$eff, labels = c("Floods", "Heatwaves"))

# Exchange Variable and modelName
mn <- mf_all1$modelName
v <- mf_all1$Variable
mf_all1$Variable <- mn
mf_all1$modelName <- v

### Plot the coefficients
# Confidence intervals
interval2 <- -qnorm((1-0.95)/2)  # 95% multiplier

mf_all1$lb <- mf_all1$Coefficient - mf_all1$SE * interval2
mf_all1$ub <- mf_all1$Coefficient + mf_all1$SE * interval2

# Set limits (trim data manually)
# xlim = c(-1.5, 2.2)
# 
# mf_all1$ub[which(mf_all1$ub > mf_all1$y_max)] <- mf_all1$y_max[which(mf_all1$ub > mf_all1$y_max)]
# mf_all1$lb[which(mf_all1$lb < mf_all1$y_min)] <- mf_all1$y_min[which(mf_all1$lb < mf_all1$y_min)]

# Coef Labels
mf_all1$lab <- as.character(sprintf("%.3f", round(mf_all1$Coefficient, 3)))

mf_all1$lab[mf_all1$p <= 0.1 & mf_all1$p > 0.05] <- paste0(mf_all1$lab[mf_all1$p <= 0.1 & mf_all1$p > 0.05], expression("\u2020"))
mf_all1$lab[mf_all1$p <= 0.05 & mf_all1$p > 0.01] <- paste0(mf_all1$lab[mf_all1$p <= 0.05 & mf_all1$p > 0.01], "*")
mf_all1$lab[mf_all1$p <= 0.01 & mf_all1$p > 0.001] <- paste0(mf_all1$lab[mf_all1$p <= 0.01 & mf_all1$p > 0.001], "**")
mf_all1$lab[mf_all1$p <= 0.001] <- paste0(mf_all1$lab[mf_all1$p <= 0.001], "***")

# Number of cases position
mf_all1$maxx <- ave(mf_all1$ub, 
                    by = mf_all1$eff,
                    FUN = function(x) max(x))
mf_all1$minx <- ave(mf_all1$lb, 
                    by = mf_all1$eff,
                    FUN = function(x) min(x))

mf_all1$tx <- NA
oo <- which(mf_all1$modelName == levels(mf_all1$modelName)[3])
mf_all1$tx[oo] <- mf_all1$maxx[oo] - (mf_all1$maxx[oo] - mf_all1$minx[oo]) * 0.60
oo <- which(mf_all1$modelName == levels(mf_all1$modelName)[2])
mf_all1$tx[oo] <- mf_all1$maxx[oo] - (mf_all1$maxx[oo] - mf_all1$minx[oo]) * 0.38
oo <- which(mf_all1$modelName == levels(mf_all1$modelName)[1])
mf_all1$tx[oo] <- mf_all1$maxx[oo] - (mf_all1$maxx[oo] - mf_all1$minx[oo]) * 0.15
oo <- which(mf_all1$modelName == levels(mf_all1$modelName)[6])
mf_all1$tx[oo] <- mf_all1$maxx[oo] - (mf_all1$maxx[oo] - mf_all1$minx[oo]) * 0.60
oo <- which(mf_all1$modelName == levels(mf_all1$modelName)[5])
mf_all1$tx[oo] <- mf_all1$maxx[oo] - (mf_all1$maxx[oo] - mf_all1$minx[oo]) * 0.38
oo <- which(mf_all1$modelName == levels(mf_all1$modelName)[4])
mf_all1$tx[oo] <- mf_all1$maxx[oo] - (mf_all1$maxx[oo] - mf_all1$minx[oo]) * 0.15

mf_all1$ty <- NA
oo <- which(mf_all1$Variable %in% levels(mf_all1$Variable)[3])
mf_all1$ty[oo] <- 2.5 
oo <- which(mf_all1$Variable %in% levels(mf_all1$Variable)[2])
mf_all1$ty[oo] <- 1.5 
oo <- which(mf_all1$Variable %in% levels(mf_all1$Variable)[1])
mf_all1$ty[oo] <- 0.5 


mf_all1$tn <- NA
oo <- which(mf_all1$modelName %in% levels(mf_all1$modelName)[c(3, 6)])
mf_all1$tn[oo] <- "N ="

oo <- which(duplicated(mf_all1[, c("modelName", "tn", "tx", "ty")]))
mf_all1$tn[oo] <- NA


# Plot
zp_all_parti <- ggplot(mf_all1, aes(colour = modelName, shape = modelName))
zp_all_parti <- zp_all_parti + facet_grid(. ~ eff, scales = "free_x")
zp_all_parti <- zp_all_parti + geom_hline(yintercept = 0, colour = scales::alpha("black", 0.3), lty = 2)
zp_all_parti <- zp_all_parti + geom_pointrange(aes(x = Variable, y = Coefficient, ymin = Coefficient - SE*interval2,
                                       ymax = Coefficient + SE*interval2),
                                   lwd = 0.7, position = position_dodge(width = 1/1.8),
                                   fill = "black", alpha = 1)
zp_all_parti <- zp_all_parti + geom_text(aes(label = lab,
                                 x = Variable, y = Coefficient), size = 6, show.legend  = FALSE, 
                             hjust = -0.3, vjust = -0.3, position = position_dodge(width = 1/1.8), alpha = 1)
zp_all_parti <- zp_all_parti + geom_text(aes(label = N, x = ty, y = tx), 
                             size = 6, show.legend  = FALSE, 
                             hjust = 0, vjust = -0.3)
zp_all_parti <- zp_all_parti + geom_text(data = unique(mf_all1[!is.na(mf_all1$tn), c("eff", "tn", "tx", "ty", "modelName")]),
                             mapping = aes(label = tn, x = ty, y = tx), size = 6, 
                             show.legend  = FALSE, colour = "black",
                             hjust = 1.3, vjust = -0.3)
zp_all_parti <- zp_all_parti + geom_blank(aes(x = 3.4, y = 0))
zp_all_parti <- zp_all_parti + coord_flip() + theme_bw()
zp_all_parti <- zp_all_parti + scale_x_discrete(expand = c(0.05,0.00) )
zp_all_parti <- zp_all_parti + theme(legend.title = element_blank())
zp_all_parti <- zp_all_parti + labs(y = "FE coefficients of past experience", x = "")
zp_all_parti <- zp_all_parti + scale_shape_manual(values = c(6, 5, 8, 15, 17, 19)) 
zp_all_parti <- zp_all_parti + scale_color_manual(values = c(gg_color_hue(6)[4:6], gg_color_hue(6)[1:3]))
zp_all_parti <- zp_all_parti + theme(text = element_text(family = "CMU Serif", size = 20),
                         # legend.position = c(0.95, 0.99), legend.justification = c(0.95,0.99),
                         plot.title = element_text(hjust = 0.5),
                         legend.position = "bottom",
                         legend.box="vertical",
                         legend.background = element_blank(),
                         legend.box.background = element_rect(colour = "black"),
                         legend.key = element_blank(),
                         axis.text.y = element_text(colour="black", size = 20),
                         axis.title.x = element_text(colour="black"),
                         axis.text.x = element_text(colour="black"),
                         strip.background = element_blank(),
                         strip.text = element_text(size = 20, colour = "black"),
                         panel.spacing.x = unit(6, "mm"),
                         plot.margin = unit(c(0,6,0,0), "mm")
                         #axis.text.x = element_text(size=16),
                         #axis.text.y = element_text(size=16),
                         #axis.title.x = element_text(size=16),
                         #axis.title.y = element_text(size=16)
)
zp_all_parti <- zp_all_parti + ggtitle("Pro-environmental behaviour")
zp_all_parti <- zp_all_parti + guides(colour = guide_legend(override.aes = list(linetype = 0, alpha = 1), reverse = T, nrow=2, byrow=T),
                          shape = guide_legend(reverse = T, nrow=2, byrow=T))

zp_all_parti <- zp_all_parti + geom_vline(aes(xintercept = 1.5))
zp_all_parti <- zp_all_parti + geom_vline(aes(xintercept = 2.5))

print(zp_all_parti)

#### Combine plots ####


cairo_pdf(file = paste("../03_Output/", "Coefplot_parti.pdf", sep=""), width = 13, height = 9, 
          bg = "white", family="CMU Serif")
par(mar = c(0, 0, 0, 0))
par(mfrow = c(1, 1), oma = c(0, 0, 0, 0))
zp_all_parti
dev.off()







#### Coefficients Plot Belief ####

# Set up df
eff <- data.frame(matrix(NA, ncol = 7, nrow = 9))
colnames(eff) <- c("Variable", "Coefficient", "SE", "t", "p", "N", "modelName")
eff[, 1] <- rep(c(3,2,1), each = 3)

eff$modelName <- rep(c(1:3), times = 3)

# Plug in values all
for(i in c(1:9)){
  eff[i, 2:6] <- c(summary(res4[[1]][[1]][[i]])$coefficients[1, 1:4], length(res4[[1]][[1]][[i]]$residuals))
}


### Combine df, make labels

allModelFrame <- data.frame(eff)

allModelFrame$Variable <- factor(allModelFrame$Variable, levels = c(1:3),
                                 labels = c("Heatwave 14 days", 
                                            "Heatwave 1 month",
                                            "Heatwave 4 months"))

allModelFrame$modelName <- factor(allModelFrame$modelName, levels = rev(c(1:3)),
                                  labels = rev(c("Left-leaning \n partisanship",
                                                 "No partisanship",
                                                 "Right-leaning \n partisanship")))


# Add floods
load("../../Climate-Change_UK/02_Data/Modelframe_parti_att.RData")

# Combine frames
allModelFrame$eff <- 2
mf_all2 <- rbind(allModelFrame, mf_all)

mf_all2$eff <- factor(mf_all2$eff, labels = c("Floods", "Heatwaves"))

# Exchange Variable and modelName
mn <- mf_all2$modelName
v <- mf_all2$Variable
mf_all2$Variable <- mn
mf_all2$modelName <- v

### Plot the coefficients
# Confidence intervals
interval2 <- -qnorm((1-0.95)/2)  # 95% multiplier

mf_all2$lb <- mf_all2$Coefficient - mf_all2$SE * interval2
mf_all2$ub <- mf_all2$Coefficient + mf_all2$SE * interval2

# Set limits (trim data manually)
# xlim = c(-1.5, 2.2)
# 
# mf_all2$ub[which(mf_all2$ub > mf_all2$y_max)] <- mf_all2$y_max[which(mf_all2$ub > mf_all2$y_max)]
# mf_all2$lb[which(mf_all2$lb < mf_all2$y_min)] <- mf_all2$y_min[which(mf_all2$lb < mf_all2$y_min)]

# Coef Labels
mf_all2$lab <- as.character(sprintf("%.3f", round(mf_all2$Coefficient, 3)))

mf_all2$lab[mf_all2$p <= 0.1 & mf_all2$p > 0.05] <- paste0(mf_all2$lab[mf_all2$p <= 0.1 & mf_all2$p > 0.05], expression("\u2020"))
mf_all2$lab[mf_all2$p <= 0.05 & mf_all2$p > 0.01] <- paste0(mf_all2$lab[mf_all2$p <= 0.05 & mf_all2$p > 0.01], "*")
mf_all2$lab[mf_all2$p <= 0.01 & mf_all2$p > 0.001] <- paste0(mf_all2$lab[mf_all2$p <= 0.01 & mf_all2$p > 0.001], "**")
mf_all2$lab[mf_all2$p <= 0.001] <- paste0(mf_all2$lab[mf_all2$p <= 0.001], "***")

# Number of cases position
mf_all2$maxx <- ave(mf_all2$ub, 
                    by = mf_all2$eff,
                    FUN = function(x) max(x))
mf_all2$minx <- ave(mf_all2$lb, 
                    by = mf_all2$eff,
                    FUN = function(x) min(x))

mf_all2$tx <- NA
oo <- which(mf_all2$modelName == levels(mf_all2$modelName)[3])
mf_all2$tx[oo] <- mf_all2$maxx[oo] - (mf_all2$maxx[oo] - mf_all2$minx[oo]) * 0.60
oo <- which(mf_all2$modelName == levels(mf_all2$modelName)[2])
mf_all2$tx[oo] <- mf_all2$maxx[oo] - (mf_all2$maxx[oo] - mf_all2$minx[oo]) * 0.38
oo <- which(mf_all2$modelName == levels(mf_all2$modelName)[1])
mf_all2$tx[oo] <- mf_all2$maxx[oo] - (mf_all2$maxx[oo] - mf_all2$minx[oo]) * 0.15
oo <- which(mf_all2$modelName == levels(mf_all2$modelName)[6])
mf_all2$tx[oo] <- mf_all2$maxx[oo] - (mf_all2$maxx[oo] - mf_all2$minx[oo]) * 0.60
oo <- which(mf_all2$modelName == levels(mf_all2$modelName)[5])
mf_all2$tx[oo] <- mf_all2$maxx[oo] - (mf_all2$maxx[oo] - mf_all2$minx[oo]) * 0.38
oo <- which(mf_all2$modelName == levels(mf_all2$modelName)[4])
mf_all2$tx[oo] <- mf_all2$maxx[oo] - (mf_all2$maxx[oo] - mf_all2$minx[oo]) * 0.15

mf_all2$ty <- NA
oo <- which(mf_all2$Variable %in% levels(mf_all2$Variable)[3])
mf_all2$ty[oo] <- 2.5 
oo <- which(mf_all2$Variable %in% levels(mf_all2$Variable)[2])
mf_all2$ty[oo] <- 1.5 
oo <- which(mf_all2$Variable %in% levels(mf_all2$Variable)[1])
mf_all2$ty[oo] <- 0.5 


mf_all2$tn <- NA
oo <- which(mf_all2$modelName %in% levels(mf_all2$modelName)[c(3, 6)])
mf_all2$tn[oo] <- "N ="

oo <- which(duplicated(mf_all2[, c("modelName", "tn", "tx", "ty")]))
mf_all2$tn[oo] <- NA


# Plot
zp_all_parti2 <- ggplot(mf_all2, aes(colour = modelName, shape = modelName))
zp_all_parti2 <- zp_all_parti2 + facet_grid(. ~ eff, scales = "free_x")
zp_all_parti2 <- zp_all_parti2 + geom_hline(yintercept = 0, colour = scales::alpha("black", 0.3), lty = 2)
zp_all_parti2 <- zp_all_parti2 + geom_pointrange(aes(x = Variable, y = Coefficient, ymin = Coefficient - SE*interval2,
                                         ymax = Coefficient + SE*interval2),
                                     lwd = 0.7, position = position_dodge(width = 1/1.8),
                                     fill = "black", alpha = 1)
zp_all_parti2 <- zp_all_parti2 + geom_text(aes(label = lab,
                                   x = Variable, y = Coefficient), size = 6, show.legend  = FALSE, 
                               hjust = -0.3, vjust = -0.3, position = position_dodge(width = 1/1.8), alpha = 1)
zp_all_parti2 <- zp_all_parti2 + geom_text(aes(label = N, x = ty, y = tx), 
                               size = 6, show.legend  = FALSE, 
                               hjust = 0, vjust = -0.3)
zp_all_parti2 <- zp_all_parti2 + geom_text(data = unique(mf_all2[!is.na(mf_all2$tn), c("eff", "tn", "tx", "ty", "modelName")]),
                               mapping = aes(label = tn, x = ty, y = tx), size = 6, 
                               show.legend  = FALSE, colour = "black",
                               hjust = 1.3, vjust = -0.3)
zp_all_parti2 <- zp_all_parti2 + geom_blank(aes(x = 3.4, y = 0))
zp_all_parti2 <- zp_all_parti2 + coord_flip() + theme_bw()
zp_all_parti2 <- zp_all_parti2 + scale_x_discrete(expand = c(0.05,0.00) )
zp_all_parti2 <- zp_all_parti2 + theme(legend.title = element_blank())
zp_all_parti2 <- zp_all_parti2 + labs(y = "FE coefficients of past experience", x = "")
zp_all_parti2 <- zp_all_parti2 + scale_shape_manual(values = c(6, 5, 8, 15, 17, 19)) 
zp_all_parti2 <- zp_all_parti2 + scale_color_manual(values = c(gg_color_hue(6)[4:6], gg_color_hue(6)[1:3]))
zp_all_parti2 <- zp_all_parti2 + theme(text = element_text(family = "CMU Serif", size = 20),
                           # legend.position = c(0.95, 0.99), legend.justification = c(0.95,0.99),
                           plot.title = element_text(hjust = 0.5),
                           legend.position = "bottom",
                           legend.box="vertical",
                           legend.background = element_blank(),
                           legend.box.background = element_rect(colour = "black"),
                           legend.key = element_blank(),
                           axis.text.y = element_text(colour="black", size = 20),
                           axis.title.x = element_text(colour="black"),
                           axis.text.x = element_text(colour="black"),
                           strip.background = element_blank(),
                           strip.text = element_text(size = 20, colour = "black"),
                           panel.spacing.x = unit(6, "mm"),
                           plot.margin = unit(c(0,6,0,0), "mm")
                           #axis.text.x = element_text(size=16),
                           #axis.text.y = element_text(size=16),
                           #axis.title.x = element_text(size=16),
                           #axis.title.y = element_text(size=16)
)
zp_all_parti2 <- zp_all_parti2 + ggtitle("Climate change belief")
zp_all_parti2 <- zp_all_parti2 + guides(colour = guide_legend(override.aes = list(linetype = 0, alpha = 1), reverse = T, nrow=2, byrow=T),
                            shape = guide_legend(reverse = T, nrow=2, byrow=T))

zp_all_parti2 <- zp_all_parti2 + geom_vline(aes(xintercept = 2.5))
zp_all_parti2 <- zp_all_parti2 + geom_vline(aes(xintercept = 1.5))


print(zp_all_parti2)

#### Combine plots ####


cairo_pdf(file = paste("../03_Output/", "Coefplot_parti_att.pdf", sep=""), width = 13, height = 6, 
          bg = "white", family="CMU Serif")
par(mar = c(0, 0, 0, 0))
par(mfrow = c(1, 1), oma = c(0, 0, 0, 0))
zp_all_parti2
dev.off()


g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

mylegend <- g_legend(zp_all_parti2)

cairo_pdf(file = paste("../03_Output/", "Coefplot_parti_comb.pdf", sep=""), width = 13, height = 12, 
          bg = "white", family="CMU Serif")
par(mar = c(0, 0, 0, 0))
par(mfrow = c(1, 1), oma = c(0, 0, 0, 0))
grid.arrange(arrangeGrob(zp_all_parti2 + theme(legend.position="none", plot.margin = unit(c(0,6,4,0), "mm")) + labs(y = "", x = ""),
                         zp_all_parti + theme(legend.position="none"),
                         nrow = 2),
             arrangeGrob(ggplot()+geom_blank(aes(1,1)) + cowplot::theme_nothing(), mylegend, ncol = 2, widths = c(1, 10)), 
             nrow = 2, heights = c(10, 1))
dev.off()







###----------------------------------###
###         Re-arrange plots         ###
###----------------------------------###


### Beliefs

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

mylegend <- g_legend(zp_all_parti2)

cairo_pdf(file = paste("../03_Output/", "Coefplot_interaction_att.pdf", sep=""), width = 13, height = 14, 
          bg = "white", family="CMU Serif")
par(mar = c(0, 0, 0, 0))
par(mfrow = c(1, 1), oma = c(0, 0, 0, 0))
grid.arrange(arrangeGrob(zp_all_parti2 + theme(legend.position="none", plot.margin = unit(c(0,6,4,0), "mm")) + labs(y = "", x = ""),
                         zp_all_exag2 + theme(legend.position="none") + ggtitle(element_blank()),
                         nrow = 2, heights = c(6, 4)),
             arrangeGrob(ggplot()+geom_blank(aes(1,1)) + cowplot::theme_nothing(), mylegend, ncol = 2, widths = c(1, 10)), 
             nrow = 2, heights = c(10, 1))
dev.off()


### Behaviour

mylegend <- g_legend(zp_all_parti)

cairo_pdf(file = paste("../03_Output/", "Coefplot_interaction_behav.pdf", sep=""), width = 13, height = 14, 
          bg = "white", family="CMU Serif")
par(mar = c(0, 0, 0, 0))
par(mfrow = c(1, 1), oma = c(0, 0, 0, 0))
grid.arrange(arrangeGrob(zp_all_parti + theme(legend.position="none", plot.margin = unit(c(0,6,4,0), "mm")) + labs(y = "", x = ""),
                         zp_all_exag + theme(legend.position="none") + ggtitle(element_blank()),
                         nrow = 2, heights = c(6, 4)),
             arrangeGrob(ggplot()+geom_blank(aes(1,1)) + cowplot::theme_nothing(), mylegend, ncol = 2, widths = c(1, 10)), 
             nrow = 2, heights = c(10, 1))
dev.off()








### ----- Redo with reduced sample ----- ###


### Reduce to people with zero in first wave

bhps_heat.df <- bhps_heat.df[order(bhps_heat.df$pidp, bhps_heat.df$wave), ]
bhps_heat.df$first_belief <- ave(bhps_heat.df$climatechange30,
                                 bhps_heat.df$pidp,
                                 FUN = function(x) x[1])


### Repeat

depvar <- c("climatechange30")
intens <- c("12", "6", "1")
buffer <- c(1, 2, 3)
b <- buffer[1]
plm.l <- vector("list", nrow(expand.grid(intens)) * 1)
plm3.l <- vector("list", nrow(expand.grid(buffer)) * 1) 

res4 <- vector("list", length(depvar)) 

l <- 1
for(d in depvar) {
  library(lfe)
  names(plm.l) <- paste0(d, c(1:length(plm.l)))
  names(plm3.l) <- paste0(d, c(1:length(plm3.l)))
  
  j <- 1
  k <- 1
  for(i in intens){
    
    v1 <- paste0("heat", i, "_affect1y_past", b)
    ev1 <- paste0("heat_ever", b)
    
    fm2 <- paste0(d, " ~ ", v1, 
                  "+ age_cat + hiqual_dv + child + marstat_dv + hhinc_dec")
    
    # Estimate via felm with twoway clustered SEs
    bhps_heat.df$interaction <- as.numeric(as.factor(bhps_heat.df$season_index)) * bhps_heat.df[, ev1]
    
    plm2_1 <- felm(formula(paste0(fm2, " | pidp + season_index + interaction | 0 | pidp + lsoa01")),
                   data = bhps_heat.df[which(bhps_heat.df[, paste0("samplefe_", d)] == 1 &  bhps_heat.df[, "index2"] == 1 & bhps_heat.df$first_belief == 0), ],
                   psdef = FALSE)
    
    plm2_2 <- felm(formula(paste0(fm2, " | pidp + season_index + interaction | 0 | pidp + lsoa01")),
                   data = bhps_heat.df[which(bhps_heat.df[, paste0("samplefe_", d)] == 1 &  bhps_heat.df[, "index2"] == 2 & bhps_heat.df$first_belief == 0), ],
                   psdef = FALSE)
    
    plm2_3 <- felm(formula(paste0(fm2, " | pidp + season_index + interaction | 0 | pidp + lsoa01")),
                   data = bhps_heat.df[which(bhps_heat.df[, paste0("samplefe_", d)] == 1 &  bhps_heat.df[, "index2"] == 3 & bhps_heat.df$first_belief == 0), ],
                   psdef = FALSE)
    
    
    
    # Interacted model
    fm3 <- paste0(d, " ~ index2 * (", v1, 
                  "+ age_cat + hiqual_dv + child + marstat_dv + hhinc_dec", 
                  "  + as.factor(season_index)*", ev1, ")")
    
    plm3 <- felm(formula(paste0(fm3, " | pidp | 0 | pidp + lsoa01")),
                 data = bhps_heat.df[which(bhps_heat.df[, paste0("samplefe_", d)] == 1   & bhps_heat.df$first_belief == 0), ], 
                 psdef = FALSE)
    
    
    plm.l[[j]] <- plm2_1
    plm.l[[j+1]] <- plm2_2
    plm.l[[j+2]] <- plm2_3
    plm3.l[[k]] <- plm3
    j <- j + 3
    k <- k + 1
    
  }
  res4[[l]] <- list(plm.l,  plm3.l)
  l <- l + 1
}



#### Coefficients Plot Belief ####

# Set up df
eff <- data.frame(matrix(NA, ncol = 7, nrow = 9))
colnames(eff) <- c("Variable", "Coefficient", "SE", "t", "p", "N", "modelName")
eff[, 1] <- rep(c(3,2,1), each = 3)

eff$modelName <- rep(c(1:3), times = 3)

# Plug in values all
for(i in c(1:9)){
  eff[i, 2:6] <- c(summary(res4[[1]][[1]][[i]])$coefficients[1, 1:4], length(res4[[1]][[1]][[i]]$residuals))
}


### Combine df, make labels

allModelFrame <- data.frame(eff)

allModelFrame$Variable <- factor(allModelFrame$Variable, levels = c(1:3),
                                 labels = c("Heatwave 14 days", 
                                            "Heatwave 1 month",
                                            "Heatwave 4 months"))

allModelFrame$modelName <- factor(allModelFrame$modelName, levels = rev(c(1:3)),
                                  labels = rev(c("Left-leaning \n partisanship",
                                                 "No partisanship",
                                                 "Right-leaning \n partisanship")))


# Add floods
load("../../Climate-Change_UK/02_Data/Subset0_Modelframe_parti_att.RData")

# Combine frames
allModelFrame$eff <- 2
mf_all2 <- rbind(allModelFrame, mf_all)

mf_all2$eff <- factor(mf_all2$eff, labels = c("Floods", "Heatwaves"))

# Exchange Variable and modelName
mn <- mf_all2$modelName
v <- mf_all2$Variable
mf_all2$Variable <- mn
mf_all2$modelName <- v

### Plot the coefficients
# Confidence intervals
interval2 <- -qnorm((1-0.95)/2)  # 95% multiplier

mf_all2$lb <- mf_all2$Coefficient - mf_all2$SE * interval2
mf_all2$ub <- mf_all2$Coefficient + mf_all2$SE * interval2

# Set limits (trim data manually)
# xlim = c(-1.5, 2.2)
# 
# mf_all2$ub[which(mf_all2$ub > mf_all2$y_max)] <- mf_all2$y_max[which(mf_all2$ub > mf_all2$y_max)]
# mf_all2$lb[which(mf_all2$lb < mf_all2$y_min)] <- mf_all2$y_min[which(mf_all2$lb < mf_all2$y_min)]

# Coef Labels
mf_all2$lab <- as.character(sprintf("%.3f", round(mf_all2$Coefficient, 3)))

mf_all2$lab[mf_all2$p <= 0.1 & mf_all2$p > 0.05] <- paste0(mf_all2$lab[mf_all2$p <= 0.1 & mf_all2$p > 0.05], expression("\u2020"))
mf_all2$lab[mf_all2$p <= 0.05 & mf_all2$p > 0.01] <- paste0(mf_all2$lab[mf_all2$p <= 0.05 & mf_all2$p > 0.01], "*")
mf_all2$lab[mf_all2$p <= 0.01 & mf_all2$p > 0.001] <- paste0(mf_all2$lab[mf_all2$p <= 0.01 & mf_all2$p > 0.001], "**")
mf_all2$lab[mf_all2$p <= 0.001] <- paste0(mf_all2$lab[mf_all2$p <= 0.001], "***")

# Number of cases position
mf_all2$maxx <- ave(mf_all2$ub, 
                    by = mf_all2$eff,
                    FUN = function(x) max(x))
mf_all2$minx <- ave(mf_all2$lb, 
                    by = mf_all2$eff,
                    FUN = function(x) min(x))

mf_all2$tx <- NA
oo <- which(mf_all2$modelName == levels(mf_all2$modelName)[3])
mf_all2$tx[oo] <- mf_all2$maxx[oo] - (mf_all2$maxx[oo] - mf_all2$minx[oo]) * 0.60
oo <- which(mf_all2$modelName == levels(mf_all2$modelName)[2])
mf_all2$tx[oo] <- mf_all2$maxx[oo] - (mf_all2$maxx[oo] - mf_all2$minx[oo]) * 0.38
oo <- which(mf_all2$modelName == levels(mf_all2$modelName)[1])
mf_all2$tx[oo] <- mf_all2$maxx[oo] - (mf_all2$maxx[oo] - mf_all2$minx[oo]) * 0.15
oo <- which(mf_all2$modelName == levels(mf_all2$modelName)[6])
mf_all2$tx[oo] <- mf_all2$maxx[oo] - (mf_all2$maxx[oo] - mf_all2$minx[oo]) * 0.60
oo <- which(mf_all2$modelName == levels(mf_all2$modelName)[5])
mf_all2$tx[oo] <- mf_all2$maxx[oo] - (mf_all2$maxx[oo] - mf_all2$minx[oo]) * 0.38
oo <- which(mf_all2$modelName == levels(mf_all2$modelName)[4])
mf_all2$tx[oo] <- mf_all2$maxx[oo] - (mf_all2$maxx[oo] - mf_all2$minx[oo]) * 0.15

mf_all2$ty <- NA
oo <- which(mf_all2$Variable %in% levels(mf_all2$Variable)[3])
mf_all2$ty[oo] <- 2.5 
oo <- which(mf_all2$Variable %in% levels(mf_all2$Variable)[2])
mf_all2$ty[oo] <- 1.5 
oo <- which(mf_all2$Variable %in% levels(mf_all2$Variable)[1])
mf_all2$ty[oo] <- 0.5 


mf_all2$tn <- NA
oo <- which(mf_all2$modelName %in% levels(mf_all2$modelName)[c(3, 6)])
mf_all2$tn[oo] <- "N ="

oo <- which(duplicated(mf_all2[, c("modelName", "tn", "tx", "ty")]))
mf_all2$tn[oo] <- NA


# Plot
zp_sub0_parti2 <- ggplot(mf_all2, aes(colour = modelName, shape = modelName))
zp_sub0_parti2 <- zp_sub0_parti2 + facet_grid(. ~ eff, scales = "free_x")
zp_sub0_parti2 <- zp_sub0_parti2 + geom_hline(yintercept = 0, colour = scales::alpha("black", 0.3), lty = 2)
zp_sub0_parti2 <- zp_sub0_parti2 + geom_pointrange(aes(x = Variable, y = Coefficient, ymin = Coefficient - SE*interval2,
                                                       ymax = Coefficient + SE*interval2),
                                                   lwd = 0.7, position = position_dodge(width = 1/1.8),
                                                   fill = "black", alpha = 1)
zp_sub0_parti2 <- zp_sub0_parti2 + geom_text(aes(label = lab,
                                                 x = Variable, y = Coefficient), size = 6, show.legend  = FALSE, 
                                             hjust = -0.3, vjust = -0.3, position = position_dodge(width = 1/1.8), alpha = 1)
zp_sub0_parti2 <- zp_sub0_parti2 + geom_text(aes(label = N, x = ty, y = tx), 
                                             size = 6, show.legend  = FALSE, 
                                             hjust = 0, vjust = -0.3)
zp_sub0_parti2 <- zp_sub0_parti2 + geom_text(data = unique(mf_all2[!is.na(mf_all2$tn), c("eff", "tn", "tx", "ty", "modelName")]),
                                             mapping = aes(label = tn, x = ty, y = tx), size = 6, 
                                             show.legend  = FALSE, colour = "black",
                                             hjust = 1.3, vjust = -0.3)
zp_sub0_parti2 <- zp_sub0_parti2 + geom_blank(aes(x = 3.4, y = 0))
zp_sub0_parti2 <- zp_sub0_parti2 + coord_flip() + theme_bw()
zp_sub0_parti2 <- zp_sub0_parti2 + scale_x_discrete(expand = c(0.05,0.00) )
zp_sub0_parti2 <- zp_sub0_parti2 + theme(legend.title = element_blank())
zp_sub0_parti2 <- zp_sub0_parti2 + labs(y = "FE coefficients of past experience", x = "")
zp_sub0_parti2 <- zp_sub0_parti2 + scale_shape_manual(values = c(6, 5, 8, 15, 17, 19)) 
zp_sub0_parti2 <- zp_sub0_parti2 + scale_color_manual(values = c(gg_color_hue(6)[4:6], gg_color_hue(6)[1:3]))
zp_sub0_parti2 <- zp_sub0_parti2 + theme(text = element_text(family = "CMU Serif", size = 20),
                                         # legend.position = c(0.95, 0.99), legend.justification = c(0.95,0.99),
                                         plot.title = element_text(hjust = 0.5),
                                         legend.position = "bottom",
                                         legend.box="vertical",
                                         legend.background = element_blank(),
                                         legend.box.background = element_rect(colour = "black"),
                                         legend.key = element_blank(),
                                         axis.text.y = element_text(colour="black", size = 20),
                                         axis.title.x = element_text(colour="black"),
                                         axis.text.x = element_text(colour="black"),
                                         strip.background = element_blank(),
                                         strip.text = element_text(size = 20, colour = "black"),
                                         panel.spacing.x = unit(6, "mm"),
                                         plot.margin = unit(c(0,6,0,0), "mm")
                                         #axis.text.x = element_text(size=16),
                                         #axis.text.y = element_text(size=16),
                                         #axis.title.x = element_text(size=16),
                                         #axis.title.y = element_text(size=16)
)
zp_sub0_parti2 <- zp_sub0_parti2 + ggtitle("Climate change belief")
zp_sub0_parti2 <- zp_sub0_parti2 + guides(colour = guide_legend(override.aes = list(linetype = 0, alpha = 1), reverse = T, nrow=2, byrow=T),
                                          shape = guide_legend(reverse = T, nrow=2, byrow=T))

zp_sub0_parti2 <- zp_sub0_parti2 + geom_vline(aes(xintercept = 2.5))
zp_sub0_parti2 <- zp_sub0_parti2 + geom_vline(aes(xintercept = 1.5))


print(zp_sub0_parti2)

#### Combine plots ####


cairo_pdf(file = paste("../03_Output/", "Subset0_Coefplot_parti_att.pdf", sep=""), width = 13, height = 6, 
          bg = "white", family="CMU Serif")
par(mar = c(0, 0, 0, 0))
par(mfrow = c(1, 1), oma = c(0, 0, 0, 0))
zp_sub0_parti2
dev.off()


g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

mylegend <- g_legend(zp_sub0_parti2)

# cairo_pdf(file = paste("../03_Output/", "Coefplot_parti_comb.pdf", sep=""), width = 13, height = 12, 
#           bg = "white", family="CMU Serif")
# par(mar = c(0, 0, 0, 0))
# par(mfrow = c(1, 1), oma = c(0, 0, 0, 0))
# grid.arrange(arrangeGrob(zp_all_parti2 + theme(legend.position="none", plot.margin = unit(c(0,6,4,0), "mm")) + labs(y = "", x = ""),
#                          zp_all_parti + theme(legend.position="none"),
#                          nrow = 2),
#              arrangeGrob(ggplot()+geom_blank(aes(1,1)) + cowplot::theme_nothing(), mylegend, ncol = 2, widths = c(1, 10)), 
#              nrow = 2, heights = c(10, 1))
# dev.off()



###----------------------------------###
###         Re-arrange plots         ###
###----------------------------------###


### Beliefs

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

mylegend <- g_legend(zp_sub0_parti2)

cairo_pdf(file = paste("../03_Output/", "Subset0_Coefplot_interaction_att.pdf", sep=""), width = 13, height = 14, 
          bg = "white", family="CMU Serif")
par(mar = c(0, 0, 0, 0))
par(mfrow = c(1, 1), oma = c(0, 0, 0, 0))
grid.arrange(arrangeGrob(zp_sub0_parti2 + theme(legend.position="none", plot.margin = unit(c(0,6,4,0), "mm")) + labs(y = "", x = ""),
                         zp_sub0_exag2 + theme(legend.position="none") + ggtitle(element_blank()),
                         nrow = 2, heights = c(6, 4)),
             arrangeGrob(ggplot()+geom_blank(aes(1,1)) + cowplot::theme_nothing(), mylegend, ncol = 2, widths = c(1, 10)), 
             nrow = 2, heights = c(10, 1))
dev.off()












