
#### Descriptives ####
#### Tobias Ruettenauer ####
#### 2020/ 01 / 31 ####

rm(list=ls())

### Load packages
library(rgdal)
library(spdep)
library(rgeos)
library(doParallel)
library(future.apply)
plan(multisession)
library(foreign)
library(GISTools)
library(cleangeo)
library(ggplot2)
library(extrafont)
loadfonts()

library(MatchIt)
library(Matching)
library(CBPS)


### Working Directory
setwd("C:/work/Forschung/Climate Change_Replication/02_Data")


#################
### Load data ###
#################


load("bhpsukhls_heat_final.RData")





####################################################
### Number of cases that are affected every year ###
####################################################

table(bhps_heat.df$istrtdaty[!is.na(bhps_heat.df$climatechange30)])

table(bhps_heat.df$istrtdaty[!is.na(bhps_heat.df$orga3)])


table(bhps_heat.df$istrtdaty[bhps_heat.df$heat_affect_past == 1])



#####################
### Define macros ###
#####################

depvar <- c("climatechange30", "env_index_all2")
intens <- c("")
buffer <- c(1, 2, 3)

i <- intens[1]
b <- buffer[1]
# dv <- depvar[1]

v1 <- paste0("heat", i, "_affect_past", b)
t1 <- paste0("tempdist", i, "_affect_past", b)
f1 <- paste0("tempdist", i, "_affect_future", b)


depvarchars <- list(climatechange30 = list(lab = "Belief in climate change",
                                           leg = "% belief in climate change",
                                           multipl = 100,
                                           file = "Timepath1_belief"),
                    orgm3 = list(lab = "Environmental organisation",
                                 leg = "% in environmental organisation",
                                 multipl = 100,
                                 file = "Timepath2_orga"),
                    env_index_all2 = list(lab = "Pro-environmental behaviour",
                                          leg = "Pro-environmental behaviour (1-5)",
                                          multipl = 1,
                                          file = "Timepath3_behav"),
                    pol_frust = list(lab = "Political frustration",
                                     leg = "Political frustration (1-5)",
                                     multipl = 1,
                                     file = "Timepath4_pol"))


### Matching
bhps_heat.df$date_num <- as.numeric(bhps_heat.df$date)
fm_exact <- treatment ~ istrtdaty + istrtdatm 
fm_nb <- treatment ~ istrtdatd 
fm_maha <- treatment ~ date_num
ratio <- 3
caliper <- 5

basis <- 10
  
set.seed(456789123)  
  

#### ------- Start looping over depvars ------- ####

sample.list <- vector(mode = "list", length = length(depvar))
names(sample.list) <- depvar

for(dv in depvar){
  
  chars <- depvarchars[dv][[1]]
  
  #######################################
  ### Temporal distance in categories ###
  #######################################
  
  
  ### Rule: Tempdist to past event, if over 3 years, tempdist to future event (dont care beyond 1 years)
  bhps_heat.df$tempdist <- bhps_heat.df[, t1]
  bhps_heat.df$tempdist[bhps_heat.df$tempdist > (365*0.5)] <- NA
  
  oo <- which(is.na(bhps_heat.df$tempdist))
  bhps_heat.df$tempdist[oo] <- bhps_heat.df[, f1][oo]
  bhps_heat.df$tempdist[bhps_heat.df$tempdist < (-365*0.5)] <- NA
  
  summary(bhps_heat.df$tempdist)
  
  
  ### Make 3monthly bins
  min <- min(bhps_heat.df$tempdist, na.rm = TRUE)
  max <- max(bhps_heat.df$tempdist, na.rm = TRUE)
  cat1 <- split(c(min:0), ceiling(seq_along(c(min:0))/(365*0.5/12)))
  cat2 <- split(c(1:max), ceiling(seq_along(c(1:max))/(365*0.5/12)))
  names(cat1) <- as.character(as.numeric(names(cat1)) - 13)
  cat <- c(cat1, cat2)
  
  cuts <- c(min -0.5, unlist(lapply(cat, FUN = function(x) max(x) + 0.5 )))
  
  bhps_heat.df$tempdist_cat <- cut(bhps_heat.df$tempdist,
                                    breaks = cuts,
                                    labels = paste0(c(-12:12)))
  
  table(bhps_heat.df$tempdist_cat)
  
  
  
  
  ###########################
  ### Match control group ###
  ###########################
  
  bhps_heat.df$treatment <- NA
  bhps_heat.df$treatment[!is.na(bhps_heat.df$tempdist_cat)] <- 1
  bhps_heat.df$treatment[is.na(bhps_heat.df$tempdist_cat)] <- 0
  table(bhps_heat.df$treatment)
  
  ### Use nearest nb matching with exact match on year and month
  vars <- unique(c("pidp", "year", "wave",  
                   all.vars(terms(fm_nb)), all.vars(terms(fm_exact)), all.vars(terms(fm_maha)), dv))
  data <- bhps_heat.df[complete.cases(bhps_heat.df[, vars]), c(vars, "tempdist")]
  
  
  # cl <- makeCluster(3)
  m.out <-  matchit(fm_maha, data = data, method = "nearest", ratio = ratio,
                    distance = "mahalanobis", m.order = "random",
                    caliper = c("date_num" = caliper), std.caliper = FALSE,
                    verbose = TRUE)
  m.nb.mat <- data.frame(m.out$match.matrix)
  weights <- data.frame(id = names(m.out$weights), w = m.out$weights)
  
  ## Merge data
  id <- c("pidp", "year", "wave")
  names <- colnames(m.nb.mat)
  m.nb.mat$id <- rownames(m.nb.mat)
  data$id <- rownames(data)
  m.nb.mat <- merge(m.nb.mat, data[, c("id", id, "tempdist")], 
                    by = "id", all.x = TRUE)
  
  ## Reshape to long based on match ids
  data.matched <- reshape(m.nb.mat, direction = "long",
                          varying = c("id", names), 
                          timevar = "X", v.names = "id")
  data.matched <- data.matched[which(!is.na(data.matched$id)),]
  data.matched <- merge(data.matched, data[, c("id", "treatment", "istrtdatd", "istrtdaty", "istrtdatm", dv)], 
                        by = "id", all.x = TRUE)
  data.matched <- merge(data.matched, weights, by = "id", all.x = TRUE)

  
  ### Bins
  data.matched$tempdist_cat <- cut(data.matched$tempdist,
                                   breaks = cuts,
                                   labels = paste0(c(-12:12)))
  
  table(data.matched$treatment, data.matched$tempdist_cat)
  
  
  
  
  ##############################
  ### Aggregate data to bins ###
  ##############################
  
  ### Treatment group
  df <- data.matched[data.matched$treatment == 1, ]
  treat.n <- aggregate(list(N = df[, dv]), list(time = factor(df$tempdist_cat)),
                       FUN = function(x) length(x))
  treat.df <- aggregate(list(CC30 = (df[, dv] * chars$multipl)), list(time = factor(df$tempdist_cat)),
                        FUN = function(x) mean(x, na.rm = TRUE))
  treat.sd <- aggregate(list(sd = (df[, dv] * chars$multipl)), list(time = factor(df$tempdist_cat)),
                        FUN = function(x) sd(x, na.rm = TRUE))
  treat.sd$sem <- treat.sd$sd/sqrt(treat.n$N)
  
  # Combine
  treat.df <- merge(treat.df, treat.sd, by = "time")
  
  # upper and lower bound
  treat.df$lb <- treat.df$CC30 - treat.df$sem
  treat.df$ub <- treat.df$CC30 + treat.df$sem
  
  treat.df$treatment <- 1
  
  
  ### Control
  df <- data.matched[data.matched$treatment == 0, ]
  control.n <- aggregate(list(N = df[, dv]), list(time = factor(df$tempdist_cat)),
                         FUN = function(x) length(x))
  control.df <- aggregate(list(CC30 = (df[, dv] * chars$multipl)), list(time = factor(df$tempdist_cat)),
                          FUN = function(x) mean(x, na.rm = TRUE))
  control.sd <- aggregate(list(sd = (df[, dv] * chars$multipl)), list(time = factor(df$tempdist_cat)),
                          FUN = function(x) sd(x, na.rm = TRUE))
  control.sd$sem <- control.sd$sd/sqrt(control.n$N)
  
  # Combine
  control.df <- merge(control.df, control.sd, by = "time")
  
  # upper and lower bound
  control.df$lb <- control.df$CC30 - control.df$sem
  control.df$ub <- control.df$CC30 + control.df$sem
  
  control.df$treatment <- 0
  
  
  ### Combine
  combined.df <- rbind(treat.df, control.df)
  combined.df <- combined.df[order(combined.df$treatment, combined.df$time), ]
  
  
  
  
  ############
  ### Plot ###
  ############
  
  ### Tempdist in month for scale x
  data.matched$time <- data.matched$tempdist / (365*0.5/12)
  data.matched$CC30 <- data.matched[, dv]
  
  data.matched$before <- 0
  data.matched$before[as.numeric(data.matched$time) < 0] <- 1
  combined.df$before <- 0
  combined.df$before[as.numeric(combined.df$tempdist_cat) <= 12] <- 1
  combined.df$time <- as.numeric(as.character(combined.df$time))
  
  ### Adjust scales Factors
  combined.df$treatment <- as.factor(combined.df$treatment)
  data.matched$treatment <- as.factor(data.matched$treatment)
  combined.df$treatment <- factor(combined.df$treatment, level = c(0, 1),
                                  labels = c("Control group", "Treatment group"))
  data.matched$treatment <- factor(data.matched$treatment, level = c(0, 1),
                                   labels = c("Control group", "Treatment group"))
  
  data.matched$time <- data.matched$time + 12
  combined.df$time <- combined.df$time + 12
  
  
  
  labels <- c("(-6)",
              "(-5.5)",
              "(-5)",
              "(-4.5)",
              "(-4)",
              "(-3.5)",
              "(-3)",
              "(-2.5)",
              "(-2)",
              "(-1.5)",
              "(-1)",
              "(-0.5)",
              "0",
              "(0.5)",
              "(1)",
              "(1.5)",
              "(2)",
              "(2.5)",
              "(3)",
              "(3.5)",
              "(4)",
              "(4.5)",
              "(5)",
              "(5.5)",
              "(6)")
  
  
  ### X axes labels to bins
  my_breaks <- c(0:24)
  my_labels <- labels
  my_labels[c(seq(2, 12, 2), seq(14, 25, 2))] <- ""
  
  
  # basis <- 6
  # knots <- (cuts + abs(min(cuts))) / (365*0.5/12)
  # fm_ga <- formula(paste0(dv, " ~ s(time, bs = 'cr', k =", basis, ", fx = TRUE)"))
  # fit1 <- mgcv::gam(fm_ga, 
  #             data = data.matched[which(data.matched$treatment == "Synthetic control group" & data.matched$before == 1), ], 
  #             method = "REML")
  
  
  ### Plot
  combined.df$before <- ifelse(combined.df$time <= 12, 1, 0)
  zp1 <- ggplot(combined.df, aes(time, CC30)) +
    geom_point( aes(x = time, y = CC30, 
                    group = treatment, colour = treatment, shape = treatment, fill = treatment, alpha = as.factor(paste0(treatment, before))), 
                size = 4, stroke = 1.3) +
    stat_smooth(data = data.matched, mapping = aes(x = time, y = CC30 * chars$multipl, group = paste0(treatment, before), 
                                                   colour = treatment, fill = treatment),
                method = "gam", formula = y ~ s(x, bs = 'ts', k = basis, fx = TRUE), se = TRUE, level = 0.95, fullrange = FALSE,
                size = 1, show.legend = FALSE) +
    stat_smooth(data = data.matched, mapping = aes(x = time, y = CC30 * chars$multipl, group = paste0(treatment, before), 
                                                   colour = treatment, fill = treatment),
                method = "gam", formula = y ~ s(x, bs = 'ts', k = basis, fx = TRUE), se = FALSE,
                size = 1, show.legend = FALSE) +
    scale_shape_manual(values = c(22, 24), name = "treatment") +
    scale_alpha_manual(values = c(1,1,1,1), guide = "none") +
    scale_color_manual(values = colorspace::darken(c("#d8b365", "#5ab4ac"), .3), name = "treatment") +
    scale_fill_manual(values = alpha(c("#d8b365", "#5ab4ac"), .8), name = "treatment") +
    theme_bw() +
    geom_vline(xintercept = 11.5, linetype = "dashed") +
    geom_vline(xintercept = 12.5, linetype = "dashed") +
    annotate("rect", xmin = 11.5, xmax = 12.5, ymin = -Inf, ymax = Inf,
             alpha = .2) +
    labs(title = chars$lab) +
    xlab(expression(paste("before heatwave   ", symbol('\254'), "    ", 'Month intervals', "    ", symbol('\256'), "   after heatwave"))) +
    ylab(chars$leg) +
    geom_blank(aes(x = -0.3, y = mean(CC30))) +
    # ylim(NA, max(c(combined.df$ub_sm, combined.df$CC30)) ) + 
    scale_x_continuous(breaks = my_breaks, labels = my_labels, minor_breaks = my_breaks, expand = c(0.03, 0.02)) + 
    theme(legend.key = element_blank(), legend.title = element_blank(),
          legend.position = c(0.95,0.95), legend.justification = c(0.95,0.95),
          legend.background = element_blank(),
          legend.box.background = element_rect(colour = "black"),
          plot.title = element_text(hjust = 0.5),
          text = element_text(family = "CM Roman", size = 24),
          axis.text.x = element_text(size = 20, colour = "black", margin = margin(t = 5, r = 0 , b = 0, l = 0)),
          axis.text.y = element_text(size = 20, colour = "black", margin = margin(t = 0, r = 5 , b = 0, l = 0)),
          axis.title.x = element_text(size = 24, colour = "black", margin = margin(t = 10, r = 0 , b = 0, l = 0)),
          axis.title.y = element_text(size = 24, colour = "black", margin = margin(t = 0, r = 10 , b = 0, l = 0))) +
    guides(colour = guide_legend(override.aes = list(linetype = 0), reverse = T),
           shape = guide_legend(reverse = T), fill = guide_legend(reverse = T))  
  
  # Coordinates of window
  x1 <- layer_scales(zp1)$x$range$range[2]
  x2 <- layer_scales(zp1)$x$range$range[1]
  y1 <- layer_scales(zp1)$y$range$range[2]
  y2 <- layer_scales(zp1)$y$range$range[1]
  r <- y1 - y2
  
  ### Add N Obs
  zp1 <- zp1 +   geom_text(data = treat.n, show.legend = FALSE,
                           aes(label = N, x =  c(0:24), y = y2 + 0.06*r),
                           family = "CM Roman", size = 6, colour = colorspace::darken(c("#d8b365", "#5ab4ac"), .3)[2]) +
    annotate("text", x =  -0.8, y = y2 + 0.06*r, label = "N[T]", parse = TRUE,
             family = "CM Roman", size = 6, colour = colorspace::darken(c("#d8b365", "#5ab4ac"), .3)[2])
  
  control.n$N <- round(control.n$N, 0)
  zp1 <- zp1 +   geom_text(data = control.n, show.legend = FALSE,
                           aes(label = N, x =  c(0:24), y = y2 + 0.0*r),
                           family = "CM Roman", size = 6, colour = colorspace::darken(c("#d8b365", "#5ab4ac"), .3)[1]) +
    annotate("text", x =  -0.8, y = y2 + 0.0*r, label = "N[C]", parse = TRUE,
             family = "CM Roman", size = 6, colour = colorspace::darken(c("#d8b365", "#5ab4ac"), .3)[1]) +
    geom_hline(yintercept = y2 + 0.10*r)
  
  
  print(zp1)
  
  
  ### Export
  cairo_pdf(file = paste("../03_Output/", chars$file, ".pdf", sep = ""), width = 15, height = 9, 
            bg = "white", family = "CM Roman")
  par(mar = c(0,0,0,0))
  par(mfrow = c(1,1), oma= c(0,0,0,0))
  print(zp1)
  dev.off()
  # svg(file = paste("../03_Output/", chars$file, ".svg", sep = ""), width = 15, height = 9, 
  #     bg = "white", family = "CM Roman")
  # par(mar = c(0,0,0,0))
  # par(mfrow = c(1,1), oma= c(0,0,0,0))
  # print(zp1)
  # dev.off()
  
  
  ### Save the samples
  sample.list[[dv]] <- data.matched
  
  ### Empty graph
  if(dv == depvar[2]){
    combined.df$before <- ifelse(combined.df$time <= 12, 1, 0)
    fc <- alpha(c("#d8b365", "#5ab4ac"), .8)
    zp1 <- ggplot(combined.df, aes(time, CC30)) +
      geom_point( aes(x = time, y = CC30, 
                      group = treatment, colour = treatment, shape = treatment, fill = treatment, alpha = as.factor(paste0(treatment, before))), 
                  size = 4, stroke = 1.3) +
      geom_smooth(data = data.matched, mapping = aes(x = time, y = CC30 * chars$multipl, group = paste0(treatment, before), 
                                                     colour = treatment, fill = treatment, alpha = as.factor(paste0("X", treatment, before))),
                  method = "gam", formula = y ~ s(x, bs = 'ts', k = basis, fx = TRUE), se = TRUE, level = 0.95, fullrange = FALSE,
                  size = NA, show.legend = FALSE) +
      stat_smooth(data = data.matched, mapping = aes(x = time, y = CC30 * chars$multipl, group = paste0(treatment, before), 
                                                     colour = treatment, fill = treatment, alpha = as.factor(paste0(treatment, before))),
                  method = "gam", formula = y ~ s(x, bs = 'ts', k = basis, fx = TRUE), se = FALSE,
                  size = 1, show.legend = FALSE, geom="line") +
      scale_shape_manual(values = c(22, 24), name = "treatment") +
      scale_color_manual(values = colorspace::darken(c("#d8b365", "#5ab4ac"), .3), name = "treatment") +
      scale_alpha_manual(values = c(1,0,1,0, 0.4,0,0.4,0), guide = "none") +
      scale_fill_manual(values = fc) +
      theme_bw() +
      geom_vline(xintercept = 11.5, linetype = "dashed") +
      geom_vline(xintercept = 12.5, linetype = "dashed") +
      annotate("rect", xmin = 11.5, xmax = 12.5, ymin = -Inf, ymax = Inf,
               alpha = .2) +
      labs(title = chars$lab) +
      xlab(expression(paste("before heatwave   ", symbol('\254'), "    ", 'Month intervals', "    ", symbol('\256'), "   after heatwave"))) +
      ylab(chars$leg) +
      geom_blank(aes(x = -0.3, y = mean(CC30))) +
      # ylim(NA, max(c(combined.df$ub_sm, combined.df$CC30)) ) + 
      scale_x_continuous(breaks = my_breaks, labels = my_labels, minor_breaks = my_breaks, expand = c(0.03, 0.02)) + 
      theme(legend.key = element_blank(), legend.title = element_blank(),
            legend.position = c(0.95,0.95), legend.justification = c(0.95,0.95),
            legend.background = element_blank(),
            legend.box.background = element_rect(colour = "black"),
            plot.title = element_text(hjust = 0.5),
            text = element_text(family = "CM Roman", size = 24),
            axis.text.x = element_text(size = 20, colour = "black", margin = margin(t = 5, r = 0 , b = 0, l = 0)),
            axis.text.y = element_text(size = 20, colour = "black", margin = margin(t = 0, r = 5 , b = 0, l = 0)),
            axis.title.x = element_text(size = 24, colour = "black", margin = margin(t = 10, r = 0 , b = 0, l = 0)),
            axis.title.y = element_text(size = 24, colour = "black", margin = margin(t = 0, r = 10 , b = 0, l = 0))) +
      guides(colour = guide_legend(override.aes = list(linetype = 0), reverse = T),
             shape = guide_legend(reverse = T), fill = guide_legend(reverse = T))  
    
    # Coordinates of window
    x1 <- layer_scales(zp1)$x$range$range[2]
    x2 <- layer_scales(zp1)$x$range$range[1]
    y1 <- layer_scales(zp1)$y$range$range[2]
    y2 <- layer_scales(zp1)$y$range$range[1]
    r <- y1 - y2
    
    ### Add N Obs
    zp1 <- zp1 +   geom_text(data = treat.n, show.legend = FALSE,
                             aes(label = N, x =  c(0:24), y = y2 + 0.06*r),
                             family = "CM Roman", size = 6, colour = colorspace::darken(c("#d8b365", "#5ab4ac"), .3)[2]) +
      annotate("text", x =  -0.8, y = y2 + 0.06*r, label = "N[T]", parse = TRUE,
               family = "CM Roman", size = 6, colour = colorspace::darken(c("#d8b365", "#5ab4ac"), .3)[2])
    
    control.n$N <- round(control.n$N, 0)
    zp1 <- zp1 +   geom_text(data = control.n, show.legend = FALSE,
                             aes(label = N, x =  c(0:24), y = y2 + 0.0*r),
                             family = "CM Roman", size = 6, colour = colorspace::darken(c("#d8b365", "#5ab4ac"), .3)[1]) +
      annotate("text", x =  -0.8, y = y2 + 0.0*r, label = "N[C]", parse = TRUE,
               family = "CM Roman", size = 6, colour = colorspace::darken(c("#d8b365", "#5ab4ac"), .3)[1]) +
      geom_hline(yintercept = y2 + 0.10*r)
    
    
    print(zp1)
    
    
    ### Export
    cairo_pdf(file = paste("../03_Output/", chars$file, "_0.pdf", sep = ""), width = 15, height = 9, 
              bg = "white", family = "CM Roman")
    par(mar = c(0,0,0,0))
    par(mfrow = c(1,1), oma= c(0,0,0,0))
    print(zp1)
    dev.off()
  }
  
}



##################################
### Simple OLS on matched data ###
##################################


for(dv in depvar){
  df <- sample.list[[dv]]
  df$after <- 1 - df$before
  df  <- merge(df, bhps_heat.df[, -which(names(bhps_heat.df) %in% c(dv, "tempdist", "tempdist_cat", "treatment"))], 
                             by = c("pidp", "year", "wave"),
                             all.x = TRUE, all.y = FALSE)
  fm <- paste0(dv, " ~ treatment +  treatment")
  lm <- lm(fm, data = df[, ])
  cat("\n", dv, "\n")
  print(summary(lm))
}




