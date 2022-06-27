## LDtandem analysis
## output.R
## N. Mizumoto

# Note
#-------------------------------------------------------------------------#
# This script performs statistical analsys + visualization.
#-------------------------------------------------------------------------#

rm(list = ls())

{
  ## Packages
  {
    require(coxme)
    require(lme4)
    require(car)
    require(multcomp)
    require(ggplot2)
    require("survminer")
    require(viridis)
    require(PupillometryR)
    require(Rmisc)
    library(extrafont)
    font_import(pattern="PT")
    loadfonts()
  }
  
  outputAll()
}

#-------------------------------------------------------------------------#
outputAll <- function(){
  outputTandemDuration()
  outputSepDuration()
  outputAcceleration()
  outputSpeedChange()
  outputMovementComparison()
}
#-------------------------------------------------------------------------#

#-------------------------------------------------------------------------#
#-------------------------------------------------------------------------#
data_summrize <- function(d, variable, cname){
  se  <-  function(x){
    y  <-  x[!is.na(x)]  #  remove  the  missing  values
    sqrt(var(as.vector(y))/length(y))
  }
  l <- length(cname)
  if(l>=2){
    df <- data.frame(as.vector(tapply(d[,cname][,1], d[,cname], unique)))
    for(i in 2:l){
      df <- data.frame(df,
                       as.vector(tapply(d[,cname][,i], d[,cname], unique))
      )
    }
  } else {
    df <- data.frame(as.vector(tapply(d[,cname], d[,cname], unique)))
  }
  df <- data.frame(df,
                   as.vector(tapply(d[,variable], d[,cname], mean)),
                   as.vector(tapply(d[,variable], d[,cname], sd)),
                   as.vector(tapply(d[,variable], d[,cname], se)))
  colnames(df) <- c(cname, paste0(variable, ".", c("mean", "sd", "se")))
  return(df)
}
#-------------------------------------------------------------------------#

#-------------------------------------------------------------------------#
# Tandem duration
#-------------------------------------------------------------------------#
outputTandemDuration <- function(){
  idir = "data/rda/"
  odir = "output/"
  load(paste0(idir,"res.tan.time.rda"))

  fit1 <- survfit(Surv(Tan.time, Cens) ~ Treat, 
          type = "kaplan-meier", 
          data = res.tan.time[res.tan.time$Species=="Rs",])
  ## Rs
  ggsurvplot(
    fit = fit1,
    conf.int = T,
    xlab = "Duration (sec)", 
    ylab = "Tandem Prob",
    xlim = c(0,1800),
    palette = viridis(10,option = "A")[c(1,8)],
    legend = c(0.8,0.8),
    ggtheme = theme_bw(),
    data=res.tan.time[res.tan.time$Species=="Rs",]
  )
  ggsave(paste0(odir, "Rs-Tandem-duration.pdf"), width=4, height = 3)
  
  f.name = "TandemDuration.txt"
  sink(paste0(odir, f.name))
  cat("R. spe.----------------------------------------\n")
  cat("coxme(Surv(Tan.time, Cens) ~ Treat + (1|Source/Video))\n")
  m <- coxme(Surv(Tan.time, Cens) ~ Treat + (1|Source/Video), data = res.tan.time[res.tan.time$Species=="Rs",])
  print(summary(m))
  print(Anova(m))
  
  cat("\n-----------------------------------------------\n")
  
  ## Cf
  ggsurvplot(
    fit = survfit(Surv(Tan.time, Cens) ~ Treat, 
                  type = "kaplan-meier", 
                  data = res.tan.time[res.tan.time$Species=="Cf",]),
    conf.int = T,
    xlab = "Duration (sec)", 
    ylab = "Tandem Prob",
    xlim = c(0,1800),
    palette = viridis(10,option = "A")[c(1,8)],
    legend = c(0.8,0.8),
    ggtheme = theme_bw(),
    data = res.tan.time[res.tan.time$Species=="Cf",]
  )
  ggsave(paste0(odir, "Cf-Tandem-duration.pdf"), width=4, height = 3)
  
  cat("C. for.----------------------------------------\n")
  cat("coxme(Surv(Tan.time, Cens) ~ Treat + (1|Source/Video))\n")
  m <- coxme(Surv(Tan.time, Cens) ~ Treat + (1|Source/Video), data = res.tan.time[res.tan.time$Species=="Cf",])
  print(summary(m))
  print(Anova(m))
  sink()
}
#-------------------------------------------------------------------------#
  
#-------------------------------------------------------------------------#
# Separation duration
#-------------------------------------------------------------------------#
outputSepDuration <- function(){
  idir = "data/rda/"
  odir = "output/"
  load(paste0(idir,"res.sep.time.rda"))
  
  ## Rs
  ggsurvplot(
    fit = survfit(Surv(Sep.time, Cens) ~ Treat, 
                  type = "kaplan-meier", 
                  data = res.sep.time[res.sep.time$Species=="Rs",]),
    conf.int = T,
    xlab = "Duration (sec)", 
    ylab = "Tandem Prob",
    xlim = c(0,100),
    palette = viridis(10,option = "A")[c(1,8)],
    legend = c(0.8,0.8),
    ggtheme = theme_bw(),
    res.sep.time[res.sep.time$Species=="Rs",]
    
  )
  ggsave(paste0(odir, "Rs-Sep-duration.pdf"), width=4, height = 2.6, family="PT Sans")
  
  f.name = "SepDuration.txt"
  sink(paste0(odir, f.name))
  cat("R. spe.----------------------------------------\n")
  cat("coxme(Surv(Sep.time, Cens) ~ Treat + (1|Source/Video))\n")
  
  m <- coxme(Surv(Sep.time, Cens) ~ Treat + (1|Source/Video), data = res.sep.time[res.sep.time$Species=="Rs",])
  print(summary(m))
  print(Anova(m))
  
  ## Cf
  ggsurvplot(
    fit = survfit(Surv(Sep.time, Cens) ~ Treat, 
                  type = "kaplan-meier", 
                  data = res.sep.time[res.sep.time$Species=="Cf",]),
    conf.int = T,
    xlab = "Duration (sec)", 
    ylab = "Tandem Prob",
    xlim = c(0,60),
    palette = viridis(10,option = "A")[c(1,8)],
    legend = c(0.8,0.8),
    ggtheme = theme_bw(),
    data = res.sep.time[res.sep.time$Species=="Cf",]
  )
  ggsave(paste0(odir, "Cf-Sep-duration.pdf"), width=4, height = 2.6, family="PT Sans")
  
  cat("C. for.----------------------------------------\n")
  cat("coxme(Surv(Tan.time, Cens) ~ Treat + (1|Source/Video))\n")
  m <- coxme(Surv(Sep.time, Cens) ~ Treat + (1|Source/Video), data = res.sep.time[res.sep.time$Species=="Cf",])
  print(summary(m))
  print(Anova(m))
  sink()
}
#-------------------------------------------------------------------------#

#-------------------------------------------------------------------------#
# Acceleration during tandem
#-------------------------------------------------------------------------#
outputAcceleration <- function(){
  idir = "data/rda/"
  odir = "output/"
  load(paste0(idir,"Res.tandem.rda"))
  
  # statistics 
  {
    Res.tandem.acc <- Res.tandem[Res.tandem$scheme == "t",]
    df.Cf          <- Res.tandem.acc[Res.tandem.acc$species == "Cf",]
    df.Cf.D        <- df.Cf[df.Cf$treat == "Dark",]
    df.Cf.L        <- df.Cf[df.Cf$treat == "Light",]
    df.Cf.D.F      <- df.Cf.D[df.Cf.D$sex == "F",]
    df.Cf.D.M      <- df.Cf.D[df.Cf.D$sex == "M",]
    df.Cf.L.F      <- df.Cf.L[df.Cf.L$sex == "F",]
    df.Cf.L.M      <- df.Cf.L[df.Cf.L$sex == "M",]
  
    f.name <- "Acceleration.txt"
    sink(paste0(odir, f.name))
    cat("Acceleration correaltion\n")
    cat("-------------------------------------\n")
    cat("C. for.\n")
    cat("-Dark Female\n")
    r <- lmer(acc~pair.dist + (1|source/name), data=df.Cf.D.F)
    print(summary(r))
    cat("-Dark Male\n")
    r <- lmer(acc~pair.dist + (1|source/name), data=df.Cf.D.M)
    print(summary(r))
    cat("-Light Female\n")
    r <- lmer(acc~pair.dist + (1|source/name), data=df.Cf.L.F)
    print(summary(r))
    cat("-Light Male\n")
    r <- lmer(acc~pair.dist + (1|source/name), data=df.Cf.L.M)
    print(summary(r))
  
    df.Rs     <- Res.tandem.acc[Res.tandem.acc$species == "Rs",]
    df.Rs.D   <- df.Rs[df.Rs$treat == "Dark",]
    df.Rs.L   <- df.Rs[df.Rs$treat == "Light",]
    df.Rs.D.F <- df.Rs.D[df.Rs.D$sex == "F",]
    df.Rs.D.M <- df.Rs.D[df.Rs.D$sex == "M",]
    df.Rs.L.F <- df.Rs.L[df.Rs.L$sex == "F",]
    df.Rs.L.M <- df.Rs.L[df.Rs.L$sex == "M",]
    
    cat("-------------------------------------\n")
    cat("R. spe.\n")
    cat("-Dark Female\n")
    r <- lmer(acc~pair.dist + (1|source/name), data=df.Rs.D.F)
    print(summary(r))
    cat("-Dark Male\n")
    r <- lmer(acc~pair.dist + (1|source/name), data=df.Rs.D.M)
    print(summary(r))
    cat("-Light Female\n")
    r <- lmer(acc~pair.dist + (1|source/name), data=df.Rs.L.F)
    print(summary(r))
    cat("-Light Male\n")
    r <- lmer(acc~pair.dist + (1|source/name), data=df.Rs.L.M)
    print(summary(r))
    
    sink()
  }
  
  # plots
  {
    df <- Res.tandem.acc[Res.tandem.acc$species == "Cf" & Res.tandem.acc$treat == "Dark",]
    set.seed(101); df <- df[sample(1:length(df[,1]), 5000),]
    ggplot(data = df, aes(x=pair.dist, y=acc*5, group=sex, col=sex)) +
      #geom_point(alpha=0.3, size=0.25) + 
      stat_smooth(method = "lm") +
      xlim(c(6, 12)) + ylim(c(-5,5)*5) +
      theme_bw() +
      theme(aspect.ratio = 0.75, legend.position = 'none') +
      scale_color_viridis(discrete=T, option="D", end=0.5) +
      xlab("Distance between a pair (mm)") +
      ylab("Acceleration (mm/sec2)")
    ggsave(paste0(odir, "Cf-Dark-Acceleration.pdf"), width=4, height = 3)
    
    ggplot(data = df, aes(x=pair.dist, y=acc*5, group=sex, col=sex)) +
      geom_point(alpha=0.3, size=0.25) + 
      xlim(c(6, 12)) + ylim(c(-5,5)*5) +
      theme_bw() +
      theme(aspect.ratio = 0.75, legend.position = 'none',
            panel.background = element_blank(),
            panel.grid = element_blank(),
            rect = element_rect(fill = "transparent")) +
      scale_color_viridis(discrete=T, option="D", end=0.5) +
      xlab("Distance between a pair (mm)") +
      ylab("Acceleration (mm/sec2)")
    ggsave(paste0(odir, "Cf-Dark-Acceleration.png"), width=4, height = 3,  bg = "transparent")
    
    
    df <- Res.tandem.acc[Res.tandem.acc$species == "Cf" & Res.tandem.acc$treat == "Light",]
    set.seed(101); df <- df[sample(1:length(df[,1]), 5000),]
    ggplot(data = df, aes(x=pair.dist, y=acc*5, group=sex, col=sex)) +
      #geom_point(alpha=0.3, size=0.25) + 
      stat_smooth(method = "lm") +
      xlim(c(6, 12)) + ylim(c(-5,5)*5) +
      theme_bw() +
      theme(aspect.ratio = 0.75, legend.position = 'none') +
      scale_color_viridis(discrete=T, option="D", end=0.5) +
      xlab("Distance between a pair (mm)") +
      ylab("Acceleration (mm/sec2)")
    ggsave(paste0(odir, "Cf-Light-Acceleration.pdf"), width=4, height = 3)
    
    ggplot(data = df, aes(x=pair.dist, y=acc*5, group=sex, col=sex)) +
      geom_point(alpha=0.3, size=0.25) + 
      xlim(c(6, 12)) + ylim(c(-5,5)*5) +
      theme_bw() +
      theme(aspect.ratio = 0.75, legend.position = 'none',
            panel.background = element_blank(),
            panel.grid = element_blank(),
            rect = element_rect(fill = "transparent")) +
      scale_color_viridis(discrete=T, option="D", end=0.5) +
      xlab("Distance between a pair (mm)") +
      ylab("Acceleration (mm/sec2)")
    ggsave(paste0(odir, "Cf-Light-Acceleration.png"), width=4, height = 3,  bg = "transparent")
    
    
    df <- Res.tandem.acc[Res.tandem.acc$species == "Rs" & Res.tandem.acc$treat == "Dark",]
    set.seed(101); df <- df[sample(1:length(df[,1]), 5000),]
    ggplot(data = df, aes(x=pair.dist, y=acc*5, group=sex, col=sex)) +
      #geom_point(alpha=0.3, size=0.25) + 
      stat_smooth(method = "lm") +
      xlim(c(3, 10)) + ylim(c(-4,4)*5) +
      theme_bw() +
      theme(aspect.ratio = 0.75, legend.position = 'none') +
      scale_color_viridis(discrete=T, option="D", end=0.5) +
      xlab("Distance between a pair (mm)") +
      ylab("Acceleration (mm/sec2)")
    ggsave(paste0(odir, "Rs-Dark-Acceleration.pdf"), width=4, height = 3)
    
    ggplot(data = df, aes(x=pair.dist, y=acc*5, group=sex, col=sex)) +
      geom_point(alpha=0.3, size=0.25) + 
      xlim(c(3, 10)) + ylim(c(-4,4)*5) +
      theme_bw() +
      theme(aspect.ratio = 0.75, legend.position = 'none',
            panel.background = element_blank(),
            panel.grid = element_blank(),
            rect = element_rect(fill = "transparent")) +
      scale_color_viridis(discrete=T, option="D", end=0.5) +
      xlab("Distance between a pair (mm)") +
      ylab("Acceleration (mm/sec2)")
    ggsave(paste0(odir, "Rs-Dark-Acceleration.png"), width=4, height = 3,  bg = "transparent")
    
    
    df <- Res.tandem.acc[Res.tandem.acc$species == "Rs" & Res.tandem.acc$treat == "Light",]
    set.seed(101); df <- df[sample(1:length(df[,1]), 5000),]
    ggplot(data = df, aes(x=pair.dist, y=acc*5, group=sex, col=sex)) +
      #geom_point(alpha=0.3, size=0.25) + 
      stat_smooth(method = "lm") +
      xlim(c(3, 10)) + ylim(c(-4,4)*5) +
      theme_bw() +
      theme(aspect.ratio = 0.75, legend.position = 'none') +
      scale_color_viridis(discrete=T, option="D", end=0.5) +
      xlab("Distance between a pair (mm)") +
      ylab("Acceleration (mm/sec2)")
    ggsave(paste0(odir, "Rs-Light-Acceleration.pdf"), width=4, height = 3)
    
    ggplot(data = df, aes(x=pair.dist, y=acc*5, group=sex, col=sex)) +
      geom_point(alpha=0.3, size=0.25) + 
      xlim(c(3, 10)) + ylim(c(-4,4)*5) +
      theme_bw() +
      theme(aspect.ratio = 0.75, legend.position = 'none',
            panel.background = element_blank(),
            panel.grid = element_blank(),
            rect = element_rect(fill = "transparent")) +
      scale_color_viridis(discrete=T, option="D", end=0.5) +
      xlab("Distance between a pair (mm)") +
      ylab("Acceleration (mm/sec2)")
    ggsave(paste0(odir, "Rs-Light-Acceleration.png"), width=4, height = 3,  bg = "transparent")
  }
}
#-------------------------------------------------------------------------#

#-------------------------------------------------------------------------#
# Speed comparison between before and after separation
#-------------------------------------------------------------------------#
outputSpeedChange <- function(){
  idir = "data/rda/"
  odir = "output/"
  load(paste0(idir,"res.tan.time.rda"))
  load(paste0(idir,"res.sep.time.rda"))
  load(paste0(idir,"Res.tandem.rda"))
  load(paste0(idir,"tandem-ind-speed.rda"))
  
  f.name <- "SpeedChange.txt"
  sink(paste0(odir, f.name))
  
  {
    # plot
    {
      # Data prep
      {
        df <- rbind(
          data.frame(
            res.tan.time[,1:4],
            sex="F",
            scheme = "1t",
            mean.speed = res.tan.time$f.speed*5
          ),
          data.frame(
            res.tan.time[,1:4],
            sex="M",
            scheme = "1t",
            mean.speed = res.tan.time$m.speed*5
          ),
          data.frame(
            res.sep.time[,1:4],
            sex="F",
            scheme = "2s",
            mean.speed = res.sep.time$f.speed*5
          ),
          data.frame(
            res.sep.time[,1:4],
            sex="M",
            scheme = "2s",
            mean.speed = res.sep.time$m.speed*5
          )
        )
        df.Cf <- df[df$Species == "Cf",]
        df.Rs <- df[df$Species == "Rs",]
      }
      
      # Plots
      {
        # Cf Dark
        df.Cf.Dark <- df.Cf[df.Cf$Treat=="Dark",]
        sumrepdat <- summarySE(df.Cf.Dark, measurevar = "mean.speed",
                               groupvars=c("sex", "scheme"))
        ggplot(df.Cf.Dark, aes(x = scheme, y = mean.speed, fill = sex, colour = sex)) +
          geom_flat_violin(aes(fill = sex),position = position_nudge(x = .05, y = 0), 
                           adjust = 1.5, trim = FALSE, alpha = .5, colour = NA)+
          geom_point(aes(x = as.numeric(as.factor(scheme))-.05, y = mean.speed),
                     position = position_jitter(width = .05), size = .25, shape = 20, alpha=0.5)+
          geom_line(data = sumrepdat, aes(x = as.numeric(as.factor(scheme))+.02, y = mean.speed),
                    linetype = 3)+
          geom_point(data = sumrepdat, aes(x = as.numeric(as.factor(scheme))+.02, y = mean.speed), shape = 21) +
          geom_errorbar(data = sumrepdat, aes(x = as.numeric(as.factor(scheme))+.02,
                                        y = mean.speed, ymin = mean.speed-ci, ymax = mean.speed+ci), width = .05)+
          scale_color_viridis(discrete=T, option="D", end=0.5)+
          scale_fill_viridis(discrete=T, option="D", end=0.5)+
          theme_bw() + theme(aspect.ratio = 4/3, legend.position="none") +
          ylim(c(-10,50))
        ggsave(paste0(odir, "Cf-Dark-SeparationSpeed.pdf"), width=3, height = 3, family="PT Sans")
        
        # Cf Light
        df.Cf.Light <- df.Cf[df.Cf$Treat=="Light",]
        sumrepdat <- summarySE(df.Cf.Light, measurevar = "mean.speed",
                               groupvars=c("sex", "scheme"))
        ggplot(df.Cf.Light, aes(x = scheme, y = mean.speed, fill = sex, colour = sex)) +
          geom_flat_violin(aes(fill = sex),position = position_nudge(x = .05, y = 0), 
                           adjust = 1.5, trim = FALSE, alpha = .5, colour = NA)+
          geom_point(aes(x = as.numeric(as.factor(scheme))-.05, y = mean.speed),
                     position = position_jitter(width = .05), size = .25, shape = 20, alpha=0.5)+
          geom_line(data = sumrepdat, aes(x = as.numeric(as.factor(scheme))+.02, y = mean.speed),
                    linetype = 3)+
          geom_point(data = sumrepdat, aes(x = as.numeric(as.factor(scheme))+.02, y = mean.speed), shape = 21) +
          geom_errorbar(data = sumrepdat, aes(x = as.numeric(as.factor(scheme))+.02,
                                              y = mean.speed, ymin = mean.speed-ci, ymax = mean.speed+ci), width = .05)+
          scale_color_viridis(discrete=T, option="D", end=0.5)+
          scale_fill_viridis(discrete=T, option="D", end=0.5)+
          theme_bw() + theme(aspect.ratio = 4/3, legend.position="none") +
          ylim(c(-10,50))
        ggsave(paste0(odir, "Cf-Light-SeparationSpeed.pdf"), width=3, height = 3, family="PT Sans")
        
        # Rs Dark
        df.Rs.Dark <- df.Rs[df.Rs$Treat=="Dark",]
        sumrepdat <- summarySE(df.Rs.Dark, measurevar = "mean.speed",
                               groupvars=c("sex", "scheme"))
        ggplot(df.Rs.Dark, aes(x = scheme, y = mean.speed, fill = sex, colour = sex)) +
          geom_flat_violin(aes(fill = sex),position = position_nudge(x = .05, y = 0), 
                           adjust = 1.5, trim = FALSE, alpha = .5, colour = NA)+
          geom_point(aes(x = as.numeric(as.factor(scheme))-.05, y = mean.speed),
                     position = position_jitter(width = .05), size = .25, shape = 20, alpha=0.5)+
          geom_line(data = sumrepdat, aes(x = as.numeric(as.factor(scheme))+.02, y = mean.speed),
                    linetype = 3)+
          geom_point(data = sumrepdat, aes(x = as.numeric(as.factor(scheme))+.02, y = mean.speed), shape = 21) +
          geom_errorbar(data = sumrepdat, aes(x = as.numeric(as.factor(scheme))+.02,
                                              y = mean.speed, ymin = mean.speed-ci, ymax = mean.speed+ci), width = .05)+
          scale_color_viridis(discrete=T, option="D", end=0.5)+
          scale_fill_viridis(discrete=T, option="D", end=0.5)+
          theme_bw() + theme(aspect.ratio = 4/3, legend.position="none") +
          ylim(c(-5,35))
        ggsave(paste0(odir, "Rs-Dark-SeparationSpeed.pdf"), width=3, height = 3, family="PT Sans")
        
        # Rs Light
        df.Rs.Light <- df.Rs[df.Rs$Treat=="Light",]
        sumrepdat <- summarySE(df.Rs.Light, measurevar = "mean.speed",
                               groupvars=c("sex", "scheme"))
        ggplot(df.Rs.Light, aes(x = scheme, y = mean.speed, fill = sex, colour = sex)) +
          geom_flat_violin(aes(fill = sex),position = position_nudge(x = .05, y = 0), 
                           adjust = 1.5, trim = FALSE, alpha = .5, colour = NA)+
          geom_point(aes(x = as.numeric(as.factor(scheme))-.05, y = mean.speed),
                     position = position_jitter(width = .05), size = .25, shape = 20, alpha=0.5)+
          geom_line(data = sumrepdat, aes(x = as.numeric(as.factor(scheme))+.02, y = mean.speed),
                    linetype = 3)+
          geom_point(data = sumrepdat, aes(x = as.numeric(as.factor(scheme))+.02, y = mean.speed), shape = 21) +
          geom_errorbar(data = sumrepdat, aes(x = as.numeric(as.factor(scheme))+.02,
                                              y = mean.speed, ymin = mean.speed-ci, ymax = mean.speed+ci), width = .05)+
          scale_color_viridis(discrete=T, option="D", end=0.5)+
          scale_fill_viridis(discrete=T, option="D", end=0.5)+
          theme_bw() + theme(aspect.ratio = 4/3, legend.position="none") +
          ylim(c(-5,35))
        ggsave(paste0(odir, "Rs-Light-SeparationSpeed.pdf"), width=3, height = 3, family="PT Sans")
      }
    }
    
    # statistics 
    {
      cat("-----Comparison between before and after separation-----\n")
      tandem <- Res.tandem[Res.tandem$scheme == "t",]
      sep <- Res.tandem[Res.tandem$scheme == "s",]
      df <- rbind(
        data.frame(
          res.tandem.speed.ind[res.tandem.speed.ind$sex=="T-M",][,1:4],
          scheme = "1t",
          speed.diff = abs((apply(tapply(tandem$dis, tandem[,4:5], mean), 1, diff)))
        ),
        data.frame(
          res.tandem.speed.ind[res.tandem.speed.ind$sex=="T-M",][,1:4],
          scheme = "1s",
          speed.diff = abs((apply(tapply(sep$dis, sep[,4:5], mean), 1, diff)))
        )
      )
      
      cat("Compare speed difference between sexes\n")
      cat("equation: lmer(speed.diff ~ scheme + (1|source/name))\n")
      cat("-- C. for., Light\n")
      r <- lmer(speed.diff ~ scheme + (1|source/name), data=df[df$species=="Cf" & df$treat=="Light",])
      print(Anova(r))
      
      cat("-- C. for., Dark\n")
      r <- lmer(speed.diff ~ scheme + (1|source/name), data=df[df$species=="Cf" & df$treat=="Dark",])
      print(Anova(r))
      
      cat("-- R. spe., Light\n")
      r <- lmer(speed.diff ~ scheme + (1|source/name), data=df[df$species=="Rs" & df$treat=="Light",])
      print(Anova(r))
      
      cat("-- R. spe., Dark\n")
      r <- lmer(speed.diff ~ scheme + (1|source/name), data=df[df$species=="Rs" & df$treat=="Dark",])
      print(Anova(r))
    }
  }
  
  sink()
}
#-------------------------------------------------------------------------#

#-------------------------------------------------------------------------#
# Solo-Tandem movement comparison
#-------------------------------------------------------------------------#
outputMovementComparison <- function(){
  idir = "data/rda/"
  odir = "output/"
  load(paste0(idir,"tandem-ind-speed.rda"))
  load(paste0(idir,"solo-ind-speed.rda"))
  f.name = "MovementComparison.txt"
  sink(paste0(odir, f.name))
  {
    {
      df <- rbind(
        res.tandem.speed.ind[,1:11],
        res.solo.speed.ind[,c(1:7,9:12)]
      )
      df[,6:9] <- df[,6:9]*5
      df <- df[df$sex!="T-M",]
      df.Cf <- df[df$species=="Cf",]
      df.Rs <- df[df$species=="Rs",]
    }
    
    ## Speed 
    {
      cat("---------- Speed Comparison ----------\n")
      df.Cf.sum <- data_summrize(df.Cf, "move.mean", c("sex", "treat"))
      df.Rs.sum <- data_summrize(df.Rs, "move.mean", c("sex", "treat"))
      # plot
      {
        ggplot(df.Cf, aes(x=sex, y=move.mean, fill=treat))+
          geom_violin(alpha=0.4)+
          geom_dotplot(binaxis = 'y', stackdir = 'center', binwidth = 1, dotsize = .75, alpha = 1, 
                       position=position_dodge(.9))+
          scale_fill_viridis(discrete=T, option="B") +
          ylim(c(0,50)) + 
          ylab("move.mean speed (mm/sec)") +
          theme_bw()+
          theme(aspect.ratio = 3/5, legend.position = "none") +
          geom_errorbar(data=df.Cf.sum, 
                        aes(x=sex, y=move.mean.mean, ymin=move.mean.mean-move.mean.se,
                            ymax=move.mean.mean+move.mean.se), width=0.1, 
                        position=position_dodge(0.9), col=2) +
          geom_point(data=df.Cf.sum, aes(x=sex, y=move.mean.mean),
                     col=2,  position=position_dodge(0.9), size=2)+
          geom_line(data=df.Cf.sum, aes(x=as.numeric(as.factor(treat))/2.2-0.7+as.numeric(as.factor(sex)),
                                        y=move.mean.mean, group=sex),
                     col=2,  position=position_dodge(0.9), size=2)
        ggsave(paste0(odir, "Cf-speed.pdf"), width=4, height = 4, family="PT Sans")
            
        ggplot(df.Rs, aes(x=sex, y=move.mean, fill=treat))+
          geom_violin(alpha=0.4)+
          geom_dotplot(binaxis = 'y', stackdir = 'center', binwidth = 1, dotsize = .45, alpha = 1, 
                       position=position_dodge(.9))+
          scale_fill_viridis(discrete=T, option="B") +
          ylim(c(0,30)) + 
          ylab("move.mean speed (mm/sec)") +
          theme_bw()+
          theme(aspect.ratio = 3/5, legend.position = "none") +
          geom_errorbar(data=df.Rs.sum, 
                        aes(x=sex, y=move.mean.mean, ymin=move.mean.mean-move.mean.se,
                            ymax=move.mean.mean+move.mean.se), width=0.1, 
                        position=position_dodge(0.9), col=2) +
          geom_point(data=df.Rs.sum, aes(x=sex, y=move.mean.mean),
                     col=2,  position=position_dodge(0.9), size=2)+
          geom_line(data=df.Rs.sum, aes(x=as.numeric(as.factor(treat))/2.2-0.7+as.numeric(as.factor(sex)),
                                        y=move.mean.mean, group=sex),
                    col=2,  position=position_dodge(0.9), size=2)
        ggsave(paste0(odir, "Rs-speed.pdf"), width=4, height = 4, family="PT Sans")
      }
      # statistic
      {
        cat("equation: lmer(move.mean ~ treat + sex + (1|source))\n")
        cat("----- C. for. -----\n")
        r <- lmer(move.mean ~ treat + sex + (1|source), data=df.Cf)
        print(Anova(r))
        multicomparison<-glht(r,linfct=mcp(sex="Tukey"))
        print(summary(multicomparison))
        
        cat("\n")
        cat("----- R. spe. -----\n")
        r <- lmer(move.mean ~ treat + sex + (1|source), data=df.Rs)
        print(Anova(r))
      }
    }
    
    ## Pausing duration comparison
    {
      cat("\n\n---------- Pause duration Comparison ----------\n")
      
      df.Cf.sum <- data_summrize(df.Cf, "pause.duration", c("sex", "treat"))
      df.Rs.sum <- data_summrize(df.Rs, "pause.duration", c("sex", "treat"))
      
      # plot
      {
        ggplot(df.Cf, aes(x=sex, y=pause.duration, fill=treat))+
          geom_violin(alpha=0.4)+
          geom_dotplot(binaxis = 'y', stackdir = 'center', binwidth = 0.01, dotsize = .75, alpha = 1, 
                       position=position_dodge(.9))+
          scale_fill_viridis(discrete=T, option="B") +
          ylim(c(0,1)) + 
          ylab("prop of pause") +
          theme_bw()+
          theme(aspect.ratio = 3/5, legend.position = "none") +
          geom_errorbar(data=df.Cf.sum, 
                        aes(x=sex, y=pause.duration.mean, ymin=pause.duration.mean-pause.duration.se,
                            ymax=pause.duration.mean+pause.duration.se), width=0.1, 
                        position=position_dodge(0.9), col=2) +
          geom_point(data=df.Cf.sum, aes(x=sex, y=pause.duration.mean),
                     col=2,  position=position_dodge(0.9), size=2)+
          geom_line(data=df.Cf.sum, aes(x=as.numeric(as.factor(treat))/2.2-0.7+as.numeric(as.factor(sex)),
                                        y=pause.duration.mean, group=sex),
                    col=2,  position=position_dodge(0.9), size=2)
        ggsave(paste0(odir, "Cf-pause.pdf"), width=4, height = 4, family="PT Sans")
        
        ggplot(df.Rs, aes(x=sex, y=pause.duration, fill=treat))+
          geom_violin(alpha=0.4)+
          geom_dotplot(binaxis = 'y', stackdir = 'center', binwidth = 0.01, dotsize = .75, alpha = 1, 
                       position=position_dodge(.9))+
          scale_fill_viridis(discrete=T, option="B") +
          ylim(c(0,1)) + 
          ylab("pause.duration speed (mm/sec)") +
          theme_bw()+
          theme(aspect.ratio = 3/5, legend.position = "none") +
          geom_errorbar(data=df.Rs.sum, 
                        aes(x=sex, y=pause.duration.mean, ymin=pause.duration.mean-pause.duration.se,
                            ymax=pause.duration.mean+pause.duration.se), width=0.1, 
                        position=position_dodge(0.9), col=2) +
          geom_point(data=df.Rs.sum, aes(x=sex, y=pause.duration.mean),
                     col=2,  position=position_dodge(0.9), size=2)+
          geom_line(data=df.Rs.sum, aes(x=as.numeric(as.factor(treat))/2.2-0.7+as.numeric(as.factor(sex)),
                                        y=pause.duration.mean, group=sex),
                    col=2,  position=position_dodge(0.9), size=2)
        ggsave(paste0(odir, "Rs-pause.pdf"), width=4, height = 4, family="PT Sans")
      }
      
      # statistics
      {
        cat("equation: lmer(pause.duration ~ treat + sex + (1|source))\n")
        cat("----- C. for. -----\n")
        r <- lmer(pause.duration ~ treat + sex + (1|source), data=df.Cf)
        print(Anova(r))
        
        cat("\n")
        cat("----- R. spe. -----\n")
        r <- lmer(pause.duration ~ treat + sex + (1|source), data=df.Rs)
        print(Anova(r))
      }
    }
    
    ## Sinuosity comparison
    {
      cat("\n\n---------- Sinuosity Comparison ----------\n")
      
      df.Cf.sum <- data_summrize(df.Cf, "angle.var", c("sex", "treat"))
      df.Rs.sum <- data_summrize(df.Rs, "angle.var", c("sex", "treat"))
      
      # plot
      {
        ggplot(df.Cf, aes(x=sex, y=angle.var, fill=treat))+
          geom_violin(alpha=0.4)+
          geom_dotplot(binaxis = 'y', stackdir = 'center', binwidth = 0.01, dotsize = .75, alpha = 1, 
                       position=position_dodge(.9))+
          scale_fill_viridis(discrete=T, option="B") +
          ylim(c(3,3.5)) + 
          ylab("prop of pause") +
          theme_bw()+
          theme(aspect.ratio = 3/5, legend.position = "none") +
          geom_errorbar(data=df.Cf.sum, 
                        aes(x=sex, y=angle.var.mean, ymin=angle.var.mean-angle.var.se,
                            ymax=angle.var.mean+angle.var.se), width=0.1, 
                        position=position_dodge(0.9), col=2) +
          geom_point(data=df.Cf.sum, aes(x=sex, y=angle.var.mean),
                     col=2,  position=position_dodge(0.9), size=2)+
          geom_line(data=df.Cf.sum, aes(x=as.numeric(as.factor(treat))/2.2-0.7+as.numeric(as.factor(sex)),
                                        y=angle.var.mean, group=sex),
                    col=2,  position=position_dodge(0.9), size=2)
        ggsave(paste0(odir, "Cf-sinuosity.pdf"), width=4, height = 4, family="PT Sans")
        
        
        ggplot(df.Rs, aes(x=sex, y=angle.var, fill=treat))+
          geom_violin(alpha=0.4)+
          geom_dotplot(binaxis = 'y', stackdir = 'center', binwidth = 0.01, dotsize = .75, alpha = 1, 
                       position=position_dodge(.9))+
          scale_fill_viridis(discrete=T, option="B") +
          ylim(c(3,3.6)) + 
          ylab("angle.var speed (mm/sec)") +
          theme_bw()+
          theme(aspect.ratio = 3/5, legend.position = "none") +
          geom_errorbar(data=df.Rs.sum, 
                        aes(x=sex, y=angle.var.mean, ymin=angle.var.mean-angle.var.se,
                            ymax=angle.var.mean+angle.var.se), width=0.1, 
                        position=position_dodge(0.9), col=2) +
          geom_point(data=df.Rs.sum, aes(x=sex, y=angle.var.mean),
                     col=2,  position=position_dodge(0.9), size=2)+
          geom_line(data=df.Rs.sum, aes(x=as.numeric(as.factor(treat))/2.2-0.7+as.numeric(as.factor(sex)),
                                        y=angle.var.mean, group=sex),
                    col=2,  position=position_dodge(0.9), size=2)
        ggsave(paste0(odir, "Rs-sinuosity.pdf"), width=4, height = 4, family="PT Sans")
      }
      
      # statistics
      {
        cat("equation: lmer(angle.var ~ treat + sex + (1|source))\n")
        cat("----- C. for. -----\n")
        r <- lmer(angle.var ~ treat + sex + (1|source), data=df.Cf)
        print(Anova(r))
        
        cat("\n")
        cat("----- R. spe. -----\n")
        r <- lmer(angle.var ~ treat + sex + (1|source), data=df.Rs[df.Rs$angle.var>3,])
        multicomparison<-glht(r,linfct=mcp(sex="Tukey"))
        print(Anova(r))
        print(summary(multicomparison))
      }
      
    }
  }
  sink()
}
#-------------------------------------------------------------------------#