## LDtandem analysis
## dataProcessing.R
## N. Mizumoto

# Note
#-------------------------------------------------------------------------#
# This script processes all the raw data for output results.
# Raw data area at data/raw/
# and processed data will be storaged data/rda/
# The processed data will be used in outputs.R for displaying results.
#-------------------------------------------------------------------------#

rm(list = ls())

{
  ## Packages
  {
    library(data.table)
    library(stringr)
    require(silvermantest) # for h.crit
  }
  
  ## constants
  {
    fps                   <- 5
    arena.size            <- 90 # mm
    move_threshold        <- c(0.94, 0.7)
    names(move_threshold) <- c("Cf", "Rs")
    
    inter_thresh          <- c(10,7)        # threshold for interaction
    names(inter_thresh)   <- c("Cf", "Rs")

    minsec                <- 10 # minimum time for tandem runs
    minsepsec             <- 3
    thresh_dis            <- 50 # minumum moved distance for tandem
  }
  
  processAll()
  
}
#-------------------------------------------------------------------------#

processAll <- function(){
  soloPreprocess()
  tandemPreprocess()
  soloProcessForOutput()
  tandemProcessForOutput()
}

#-------------------------------------------------------------------------#
# This function convert .csv to .rda for solo experiments
# Also, convert pixel to mm. 
# Note 1
#   The scale factor can be variable among cameras (video = v1-7)
#   Assume that termites walk thoroughly on the arena, and obtain the maximum
#   range of the movement for each camera, which can be equivalent to arena.size
# Note 2
#   Video length and realtime do not match in the following videos
#   (may be due to raspberry pi's spec?)
#     Rs A L 1-5; Rs B L 1-7; Rs C L 1-5; Rs D L 1
#   I estimated the position for realtime frame assuming that termites move
#   straightly at constant speed between two successive frames
#-------------------------------------------------------------------------#
soloPreprocess <- function(){
  print("soloPreprocess")
  
  idir <- "data/rawdata/solo/"
  odir <- "data/rda/"
  
  ## data
  rawdata  <- list.files(idir, full.names = TRUE)
  dataname <- list.files(idir, full.names = FALSE)
  
  ## get scale *Note 1
  {
    get.scale <- NULL
    for(v in 1:length(rawdata)){
      dn <- dataname[v]
      # file info
      d <- data.frame(fread(rawdata[v], header=T))
      species = substr(dn, 1, 2)
      if(species == "Cf"){
        source  = paste0("d", substr(dn, 4, 7))
        
        LDtreat = substr(dn, 9, 13)
        if(LDtreat == "Dark-"){LDtreat = "Dark"}
        
        id = (substr(dn, regexpr("solo", dn)[1]+5, regexpr("solo", dn)[1]+6))
        if(str_detect(dn, "v5")){ 
          video = "v5"
        } else if(str_detect(dn, "Light")){ 
          video = "v2"
        } else if(str_detect(dn, "Dark")){ 
          video = "v4"
        }
        
      } else {
        LDtreat = substr(dn, 6, 10)
        if(LDtreat == "Dark-"){LDtreat = "Dark"}
        source = substr(dn, 4, 4)
        id = (substr(dn, regexpr("solo", dn)[1]+5, regexpr("solo", dn)[1]+5))
        video = (substr(dn, regexpr(".csv", dn)[1]-2, regexpr(".csv", dn)[1]-1))
      }
      
      print(paste0(v, " ", species, " ", LDtreat, " ", source, " ", id, " video-", video))
      
      
      # scaling
      d$position <- d$position/25
      d          <- d[seq(1,length(d[,1]),25/fps),]
      d          <- d[d$position < 60*30,]
      m.xmin     <- min(d[,4], na.rm=T)
      m.xmax     <- max(d[,4], na.rm=T)
      m.ymin     <- min(d[,5], na.rm=T)
      m.ymax     <- max(d[,5], na.rm=T)
      f.xmin     <- min(d[,2], na.rm=T)
      f.xmax     <- max(d[,2], na.rm=T)
      f.ymin     <- min(d[,3], na.rm=T)
      f.ymax     <- max(d[,3], na.rm=T)
      
      df <- data.frame(species, LDtreat, source, id, video, 
                       fx = f.xmax-f.xmin, fy = f.ymax-f.ymin,
                       mx = m.xmax-m.xmin, my = m.ymax-m.ymin)
      get.scale <- rbind(get.scale, df)
    }
    
    m.xscale <- tapply(get.scale$mx, get.scale[,c("species","video")], max, na.rm=T)
    m.yscale <- tapply(get.scale$my, get.scale[,c("species","video")], max)
    f.xscale <- tapply(get.scale$fx, get.scale[,c("species","video")], max)
    f.yscale <- tapply(get.scale$fy, get.scale[,c("species","video")], max)
  }
  
  ## scaling
  {
    Species = LD = Source = ID = NULL
    d.loc = NULL
    for(v in 1:length(rawdata)){
      dn <- dataname[v]
      
      # file info
      {
        d <- data.frame(fread(rawdata[v], header=T))
        
        species = substr(dn, 1, 2)
        if(species == "Cf"){
          source = paste0("d", substr(dn, 4, 7))
          LDtreat = substr(dn, 9, 13)
          if(LDtreat == "Dark-"){LDtreat = "Dark"}
          id = as.numeric(substr(dn, regexpr("solo", dn)[1]+5, regexpr("solo", dn)[1]+6))
          if(str_detect(dn, "v5")){ 
            video = "v5"
          } else if(str_detect(dn, "Light")){
            video = "v2"
          } else if(str_detect(dn, "Dark")){ 
            video = "v4"
          }
        } else {
          LDtreat = substr(dn, 6, 10)
          if(LDtreat == "Dark-"){LDtreat = "Dark"}
          source = substr(dn, 4, 4)
          id = as.numeric(substr(dn, regexpr("solo", dn)[1]+5, regexpr("solo", dn)[1]+5))
          video = (substr(dn, regexpr(".csv", dn)[1]-2, regexpr(".csv", dn)[1]-1))
        }
        print(paste0(v, " ", species, " ", LDtreat, " ", source, id, " video-", video))
      }
      
      # scaling
      {
        if(LDtreat=="Light"&&
           species=="Rs"&&
           (source=="A"||source=="B"||source=="C"||(source=="D"&&id==1))){
          # Note 2
          
          d$position <- d$position/25 * 1.354
          
          Lpoint     <- floor(max(d$position) * fps)
          data       <- matrix(0, ncol=5, nrow=Lpoint)
          point      <- seq(0, max(d$position), 1/fps)
          
          count <- 2
          
          data[1,] <- as.numeric(d[1,])
          for(i in 2:length(d[,1])){
            if(d$position[i] > point[count]){
              data[count,1] <- point[count]
              data[count,2] <- d[i-1,2] + (d[i,2]-d[i-1,2])*((point[count]-d$position[i-1])/(d$position[i]-d$position[i-1]))
              data[count,3] <- d[i-1,3] + (d[i,3]-d[i-1,3])*((point[count]-d$position[i-1])/(d$position[i]-d$position[i-1]))
              data[count,4] <- d[i-1,4] + (d[i,4]-d[i-1,4])*((point[count]-d$position[i-1])/(d$position[i]-d$position[i-1]))
              data[count,5] <- d[i-1,5] + (d[i,5]-d[i-1,5])*((point[count]-d$position[i-1])/(d$position[i]-d$position[i-1]))
              count <- count + 1
              if(count==length(point)){break}
            }
          }
          
          d <- data.frame(data)
          names(d) <- c("position","x0","y0","x1","y1")
          d <- d[d$position < 60*30,]
          
        } else {
          d$position <- d$position/25
          d <- d[seq(1,length(d[,1]),25/fps),]
          d <- d[d$position < 60*30,]
          
        }
        f.xmin <- (min(d[,2]) + max(d[,2]))/2
        f.ymin <- (min(d[,3]) + max(d[,3]))/2
        m.xmin <- (min(d[,4]) + max(d[,4]))/2
        m.ymin <- (min(d[,5]) + max(d[,5]))/2
        d[,2] <- (d[,2] - f.xmin) * arena.size / f.xscale[species,video]
        d[,4] <- (d[,4] - m.xmin) * arena.size / m.xscale[species,video]
        d[,3] <- (d[,3] - f.ymin) * arena.size / f.yscale[species,video]
        d[,5] <- (d[,5] - m.ymin) * arena.size / m.yscale[species,video]
      }
      
      # make an array
      d.video = array(0, dim=c(dim(d)[1],3,2))
      dimnames(d.video) = list(NULL, c("time", "x", "y"), NULL)
      
      for(i in 1:2){
        d.video[,1,i] = d[,1]
        d.video[,2,i] = d[,i*2]
        d.video[,3,i] = d[,i*2+1]
      }
      
      print(dim(d.video))
     
      Species = c(Species, species)
      LD = c(LD, LDtreat)
      Source = c(Source, source)
      ID = c(ID, id)
      d.loc <- c(d.loc, list(d.video))
    }
  }
  
  ## output
  
  df.all             <- (list(list(Species, LD, Source, ID, 1:length(rawdata)), d.loc))
  names(df.all)      <- c("meta", "coord")
  names(df.all$meta) <- c("Species", "LD", "Source", "ID", "Unique")
  
  save(df.all, file = paste0(odir,"solo.all.rda"))
}
#-------------------------------------------------------------------------#

#-------------------------------------------------------------------------#
# This function convert .csv to .rda for tandem experiments
# Also, convert pixel to mm. 
# Note 1
#   The scale factor can be variable among cameras (video = v1-7)
#   Assume that termites walk thoroughly on the arena, and obtain the maximum
#   range of the movement for each camera, which can be equivalent to arena.size
# Note 2
#   Video length and realtime do not match in the following videos
#   (may be due to raspberry pi's spec?)
#     Rs A D 2-
#   I estimated the position for realtime frame assuming that termites move
#   straightly at constant speed between two successive frames
#-------------------------------------------------------------------------#
tandemPreprocess <- function(){
  print("tandemPreprocess")
  
  idir <- "data/rawdata/tandem/"
  odir <- "data/rda/"
  
  ## data
  rawdata  <- list.files(idir, full.names = TRUE)
  dataname <- list.files(idir, full.names = FALSE)
  
  ## get scale *Note 1
  {
    get.scale <- NULL
    for(v in 1:length(rawdata)){
      dn <- dataname[v]
      
      # file info
      {
        d <- data.frame(fread(rawdata[v], header=T))
        species = substr(dn, 1, 2)
        if(species == "Cf"){
          source = paste0("d", substr(dn, 4, 7))
          LDtreat = substr(dn, 9, 13)
          if(LDtreat == "Dark-"){LDtreat = "Dark"}
          id = (substr(dn, regexpr("tandem", dn)[1]+7, regexpr("tandem", dn)[1]+8))
          if(str_detect(dn, "v5")){ video = "v5"} else if(str_detect(dn, "Light")){ video = "v2"} else if(str_detect(dn, "Dark")){ video = "v4"}
        } else {
          LDtreat = substr(dn, 6, 10)
          if(LDtreat == "Dark-"){LDtreat = "Dark"}
          source = substr(dn, 4, 4)
          id = (substr(dn, regexpr("tandem", dn)[1]+7, regexpr("tandem", dn)[1]+7))
          video = (substr(dn, regexpr(".csv", dn)[1]-2, regexpr(".csv", dn)[1]-1))
        }
        
        print(paste0(v, " ", species, " ", LDtreat, " ", source, " ", id, " video-", video))
      }
      
      # scaling
      {
        d$position <- d$position/25
        d          <- d[seq(1,length(d[,1]),25/fps),]
        d          <- d[d$position < 60*30,]
        xmin       <- min(c(d[,2],d[,4]))
        xmax       <- max(c(d[,2],d[,4]))
        ymin       <- min(c(d[,3],d[,5]))
        ymax       <- max(c(d[,3],d[,5]))
        
        df         <- data.frame(species, LDtreat, source, id, video, x = xmax-xmin, y = ymax-ymin)
        get.scale  <- rbind(get.scale, df)
      }
    }
    
    xscale <- tapply(get.scale$x, get.scale$video, max)
    yscale <- tapply(get.scale$y, get.scale$video, max)
  }
  
  ## scaling
  {
    Species = LD = Source = ID = NULL
    d.loc = NULL
    for(v in 1:length(rawdata)){
      dn <- dataname[v]
      # file info
      {
        d <- data.frame(fread(rawdata[v], header=T))
        species = substr(dn, 1, 2)
        if(species == "Cf"){
          source = paste0("d", substr(dn, 4, 7))
          LDtreat = substr(dn, 9, 13)
          if(LDtreat == "Dark-"){LDtreat = "Dark"}
          id = (substr(dn, regexpr("tandem", dn)[1]+7, regexpr("tandem", dn)[1]+8))
          if(str_detect(dn, "v5")){ video = "v5"} else if(str_detect(dn, "Light")){ video = "v2"} else if(str_detect(dn, "Dark")){ video = "v4"}
        } else {
          LDtreat = substr(dn, 6, 10)
          if(LDtreat == "Dark-"){LDtreat = "Dark"}
          source = substr(dn, 4, 4)
          id = (substr(dn, regexpr("tandem", dn)[1]+7, regexpr("tandem", dn)[1]+7))
          video = (substr(dn, regexpr(".csv", dn)[1]-2, regexpr(".csv", dn)[1]-1))
        }
        
        print(paste0(v, " ", species, " ", LDtreat, " ", source, " ", id, " video-", video))
      }
      
      # scaling
      {
        if(species=="Rs" && LDtreat=="L" && source == "D" && id>1){ #*Note2
          d$position <- d$position/25 * 1.354
          
          Lpoint <- floor(max(d$position) * fps)
          data <- matrix(0, ncol=5, nrow=Lpoint)
          point <- seq(0, max(d$position), 1/fps)
          
          count <- 2
          
          data[1,] <- as.numeric(d[1,])
          for(i in 2:length(d[,1])){
            if(d$position[i] > point[count]){
              data[count,1] <- point[count]
              data[count,2] <- d[i-1,2] + (d[i,2]-d[i-1,2])*((point[count]-d$position[i-1])/(d$position[i]-d$position[i-1]))
              data[count,3] <- d[i-1,3] + (d[i,3]-d[i-1,3])*((point[count]-d$position[i-1])/(d$position[i]-d$position[i-1]))
              data[count,4] <- d[i-1,4] + (d[i,4]-d[i-1,4])*((point[count]-d$position[i-1])/(d$position[i]-d$position[i-1]))
              data[count,5] <- d[i-1,5] + (d[i,5]-d[i-1,5])*((point[count]-d$position[i-1])/(d$position[i]-d$position[i-1]))
              count <- count + 1
              if(count==length(point)){break}
            }
          }
          d <- data.frame(data)
          names(d) <- c("position","x0","y0","x1","y1")
          d <- d[d$position < 60*30,]
        } else {
          d$position <- d$position/25
          d <- d[seq(1,length(d[,1]),25/fps),]
          d <- d[d$position < 60*30,]
        }
      
      
        xmin  <- (min(c(d[,2],d[,4])) + max(c(d[,2],d[,4])))/2
        ymin  <- (min(c(d[,3],d[,5])) + max(c(d[,3],d[,5])))/2
        d[,2] <- (d[,2] - xmin) * arena.size / xscale[video]
        d[,4] <- (d[,4] - xmin) * arena.size / xscale[video]
        d[,3] <- (d[,3] - ymin) * arena.size / yscale[video]
        d[,5] <- (d[,5] - ymin) * arena.size / yscale[video]
      }
    
      # make an array
      {
        d.video = array(0, dim=c(dim(d)[1],3,2))
        dimnames(d.video) = list(NULL, c("time", "x", "y"), NULL)
        
        for(i in 1:2){
          d.video[,1,i] = d[,1]
          d.video[,2,i] = d[,i*2]
          d.video[,3,i] = d[,i*2+1]
        }
        
        Species = c(Species, species)
        LD      = c(LD, LDtreat)
        Source  = c(Source, source)
        ID      = c(ID, id)
        d.loc   = c(d.loc, list(d.video))
      }
    }
  }
  
  ## output
  df.all <- (list(list(Species, LD, Source, ID, 1:length(rawdata)), d.loc))
  names(df.all) <- c("meta", "coord")
  names(df.all$meta) <- c("Species", "LD", "Source", "ID", "Unique")
  
  save(df.all, file = paste0(odir,"tandem.all.rda"))
}
#-------------------------------------------------------------------------#

#-------------------------------------------------------------------------#
#-------------------------------------------------------------------------#
angle_cal2 <- function(X, Y, Length){
  Ax <- (X[3:Length-1] - X[3:Length-2])
  Ay <- (Y[3:Length-1] - Y[3:Length-2])
  return(atan2(Ay,Ax)-pi/2)
}
#-------------------------------------------------------------------------#

#-------------------------------------------------------------------------#
# function to obtain movement parameters for solo experiments
#-------------------------------------------------------------------------#
soloProcessForOutput <- function(){
  print("soloProcessForOutput")
  iodir = "data/rda/"
  load(paste0(iodir, "solo.all.rda"))
  
  n <- length(df.all$meta$Unique)
  Res = res.speed.ind <- NULL
  for(v in 1:n){
    ## read data
    d <- df.all[[2]][[v]]
    dims = dim(d)
  
    ## get metadata
    species <-df.all$meta[[1]][v]
    treat <-  df.all$meta[[2]][v]
    source <- df.all$meta[[3]][v]
    id <- paste0(source, "-", df.all$meta[[4]][v])
    name = paste(species, treat, id, sep="-")
    print(name)
    
    ## Speed
    moved.dis.each <- sqrt(diff(d[,2,1:2])^2+diff(d[,3,1:2])^2)
    pause <- moved.dis.each < move_threshold[species]
    Pause_duration <- apply(pause, 2, sum)
    
    ### Get individual-level data
    ## 1. representative of speed
    ## 2. pause duration
    if(T){
      for(i in 1:2){
        if(sum(is.na(moved.dis.each[,i]))==0){
          vec <- moved.dis.each[,i]
          
          # get largest peak in a range of 1~12mm/sec
          h0 <- h.crit(vec, 2)
          y <- density(vec)$y
          x <- density(vec)$x
          d1 <- diff(y)
          signs <- diff(d1 / abs(d1))
          ind <- which(signs == -2)
          #truehist(vec, xlim=c(0,15), main=name)
          #points(x,y, xlim=c(0,7), type="l")
          #points(x[ind], y[ind], col=2, pch=19)
          #text(x[ind]+0.5, y[ind], round(x[ind],2))
          #text(12,0.1, round(mod.speed,2))
          mod.speed <- x[ind][ y[ind] == max(y[ind][x[ind]>1 & x[ind] < 12]) ]
          
          # get moving value
          move.speed <- (vec[!pause[,i]])
          
          df <- data.frame(species, treat, source, name, sex = c("F","M")[i], 
                           mean = mean(vec), median = median(vec), mod = mod.speed,
                           move.mean = mean(move.speed), move.median = median(move.speed),
                           angle.var = var(angle_cal2(d[,2,i],d[,3,i],dim(d)[1])),
                           pause.duration = Pause_duration[i]/length(pause))
          res.speed.ind <- rbind(res.speed.ind, df)
        }
      }
    }
    
    ### Output
    sex <- c("F", "M")
    for(i in 1:2){
      if(sum(is.na(moved.dis.each[,i]))==0){
        df <- data.frame(species, treat, source, name, sex = sex[i], 
                         time=1:dim(moved.dis.each)[1]+1, 
                         dis = moved.dis.each[,i],
                         acc = c(diff(moved.dis.each[,i]), NA),
                         pause = pause[,i])
        Res <- rbind(Res, df)
      }
    }
  }
  
  res.solo.speed.ind <- res.speed.ind
  save(res.solo.speed.ind, file = paste0(iodir,"solo-ind-speed.rda"))
  
  Res.solo <- Res
  save(Res.solo, file = paste0(iodir,"Res.solo.rda"))
}
#-------------------------------------------------------------------------#

#-------------------------------------------------------------------------#
#-------------------------------------------------------------------------#
tandem.smoothing  <- function(vec, min.sec){
  if(sum(vec)>0){
    timing <- which(vec)[c(T, diff(which(vec))>1)]
    end <- which(vec)[c(diff(which(vec))>1,T)]
    for(fi in 1:length(timing)){
      if(length( vec[timing[fi]:end[fi]]) < min.sec ){
        vec[timing[fi]:end[fi]] <- F
      }
    }
  }
  return(vec)
}
tandem.smoothing2 <- function(vec, xvec, min.sec){
  if(sum(vec)>0){
    timing <- which(vec)[c(T, diff(which(vec))>1)]
    end <- which(vec)[c(diff(which(vec))>1,T)]
    for(fi in 1:length(timing)){
      if(any(xvec[timing[fi]:end[fi]])){
        vec[timing[fi]:end[fi]] <- T
      } else if(length( vec[timing[fi]:end[fi]]) < min.sec ){
        vec[timing[fi]:end[fi]] <- F
      }
    }
  }
  return(vec)
}
#-------------------------------------------------------------------------#

#-------------------------------------------------------------------------#
# function to obtain tandem/separation duration and movement parameters
# for tandem experiments
#-------------------------------------------------------------------------#
tandemProcessForOutput <- function(){
  print("tandemProcessForOutput")
  iodir = "data/rda/"
  load(paste0(iodir, "tandem.all.rda"))
  
  n <- length(df.all$meta$Unique)
  
  res.speed.ind = res.tan.time = res.sep.time = Res = NULL
  
  for(v in 1:n){
    
    d <- df.all[[2]][[v]]
    dims = dim(d)
    d = d[1:dims[1],,]
    
    species <-df.all$meta[[1]][v]
    treat <-  df.all$meta[[2]][v]
    source <- df.all$meta[[3]][v]
    id <- paste0(source, "-", df.all$meta[[4]][v])
    name = paste(species, treat, id, sep="-")
    print(name)
    
    ## 1. identify tandem
    {
      dist <- sqrt((d[,2,1] - d[,2,2])^2 + (d[,3,1] - d[,3,2])^2)
      
      # interaction: distance between pair < inter_thresh mm
      interaction = dist < inter_thresh[species] & dist > 0
      tandem = interaction
      
      ## not non interaction& less than minsepsec seconds = tandem
      ## non interaction: >= 2*threshold mm
      noninteraction = dist >= inter_thresh[species]*2 & dist != 0
      if(sum(!tandem)>0){
        tandem = !tandem.smoothing2(!tandem, noninteraction, minsepsec*fps)
      }
      
      # smoothing (tandem needs to last 10 seconds)
      if(sum(tandem)>0){
        tandem = tandem.smoothing(tandem, minsec*fps)
      }
      
      ## tandem: both need to move more than thresh_dis or thresh_speed during interaction
      moved.dis.each <- sqrt(diff(d[,2,1:2])^2+diff(d[,3,1:2])^2)
      moved.dis.each2 <- rbind(0, moved.dis.each)
      
      for(i in 1:2) {
        if (sum(tandem) > 0) {
          tan.timing <- which(tandem)[c(T, diff(which(tandem)) > 1)]
          tan.end <- which(tandem)[c(diff(which(tandem)) > 1, T)]
          for (j in 1:length(tan.timing)) {
            if (sum(moved.dis.each2[tan.timing[j]:tan.end[j], i]) < thresh_dis) {
              tandem[tan.timing[j]:tan.end[j]] <- F
            }
            if (sum(diag(moved.dis.each2[tan.timing[j]:tan.end[j], 3 - i])) < thresh_dis) {
              tandem[tan.timing[j]:tan.end[j]] <- F
            }
          }
        }
      }
    }
    
    ## 2. identify separation 
    {
      scheme = duration = tandemduration <- tandem
      scheme[interaction] <- "i"
      scheme[!interaction] <- "r"
      scheme[tandem] <- "t"
      
      if (sum(tandem) > 0) {
        tan.end <- which(tandem)[c(diff(which(tandem)) > 1, T)]
        if (tan.end[1] == length(tandem)) {
          next
          
        }
        sep.begin <- tan.end[tan.end < length(tandem)] + 1
        
        ## sep until next interaction
        tan.timing <- which(interaction)[c(T, diff(which(interaction)) > 1)]
        if (length(tan.timing) == 1) {
          sep.end <- NULL
        } else {
          sep.end <- tan.timing[2:(length(tan.timing))] - 1
        }
        if (length(sep.begin) - length(sep.end) == 1) {
          sep.end <- c(sep.end, length(tandem))
        }
        for (j in 1:length(sep.begin)) {
          for (k in 1:length(sep.end)) {
            if (sep.end[k] > sep.begin[j]) {
              scheme[sep.begin[j]:sep.end[k]] <- "s"
              break
              
            }
          }
        }
      }
    }

    ## 3. Output
    
    sex <- rep(c("f","m"), each=1)
    ## tandem survival time
    {
      Video = Species = Treat = Source = Tan.time = Cens = f.speed = m.speed <- NULL
      if(sum(scheme=="t")>0){
        tan.timing <- which(tandem)[c(T, diff(which(tandem))>1)]
        tan.end <- which(tandem)[c(diff(which(tandem))>1,T)]
        for(j in 1:length(tan.timing)){
          Tan.time = c(Tan.time, tan.end[j] - tan.timing[j] + 1)
          f.speed = c(f.speed, mean(moved.dis.each2[(tan.timing[j]+1):tan.end[j],1]))
          m.speed = c(m.speed, mean(moved.dis.each2[(tan.timing[j]+1):tan.end[j],2]))
          
          if(tan.timing[j]==1 || tan.end[j]==length(tandem)){
            Cens = c(Cens, 0)
          } else {
            Cens = c(Cens, 1)
          }
          
          Video = c(Video, name)
          Species = c(Species, species)
          Treat = c(Treat, treat)
          Source = c(Source, source)
        }
      }
      df <- data.frame(Video, Species, Treat, Source, Tan.time= Tan.time/5, Cens, 
                       f.speed, m.speed)
      res.tan.time <- rbind(res.tan.time, df)
    }
    
    ## sep search survival time
    {
      if(sum(scheme=="s")>0){
        Video =  Species = Treat = Source = Sep.time = Cens <- NULL
        f.speed =  m.speed = sep.max.dis <- NULL
        if(sum(scheme=="s")>0){
          sep.timing <- which(scheme=="s")[c(T, diff(which(scheme=="s"))>1)]
          sep.end <- which(scheme=="s")[c(diff(which(scheme=="s"))>1,T)]
          for(j in 1:length(sep.timing)){
            Sep.time = c(Sep.time, sep.end[j] - sep.timing[j] + 1)
            if(sep.timing[j]==1 || sep.end[j]==length(tandem)){
              Cens = c(Cens, 0)
            } else {
              Cens = c(Cens, 1)
            }
            f.speed <- c(f.speed, mean(moved.dis.each2[sep.timing[j]:sep.end[j],1]))
            m.speed <- c(m.speed, mean(moved.dis.each2[sep.timing[j]:sep.end[j],2]))
            sep.max.dis <- c(sep.max.dis, max(dist[sep.timing[j]:sep.end[j]]))
            Video = c(Video, name)
            Species = c(Species, species)
            Treat = c(Treat, treat)
            Source = c(Source, source)
          }
        }
        df <- data.frame(Video, Species, 
                         Treat, Source, Sep.time=Sep.time/5, Cens,
                         f.speed, m.speed, sep.max.dis)
        res.sep.time <- rbind(res.sep.time, df)
      }
    }
    
    ## speed
    for(i in 1:2){
      if(sum(is.na(moved.dis.each[,i]))==0){
        df <- data.frame(species, treat, source, name, sex = c("F","M")[i], 
                         dis = moved.dis.each[,i], 
                         acc = c(NA, diff(moved.dis.each[,i])),
                         pair.dist = dist[2:dims[1]-1],
                         scheme = scheme[2:dims[1]])
        Res <- rbind(Res, df)
      }
    }
    
    ## individual speed
    ## 1. representative of speed
    ## 2. pause duration
    ## 3. distance between female and male during tandem
    ## 4. variance of turning angle
    {
      #par(mfrow=c(1,2))
      # use female speed as tandem speed
      f.vec <- na.omit(moved.dis.each[,1][scheme=="t"])
      m.vec <- na.omit(moved.dis.each[,2][scheme=="t"])
      if(length(f.vec>0)){
        
        # get moving value
        f.pause <- f.vec < move_threshold[species]
        f.move.speed <- (f.vec[!f.pause])
        m.pause <- m.vec < move_threshold[species]
        m.move.speed <- (m.vec[!m.pause])
        
        df <- data.frame(species, treat, source, name,  sex = "T-F",
                         mean = mean(f.vec), median = median(f.vec), 
                         move.mean = mean(f.move.speed), 
                         move.median = median(f.move.speed),
                         pause.duration = sum(f.pause)/sum(tandem),
                         angle.var = var(angle_cal2(d[,2,1],d[,3,1],dim(d)[1])),
                         tandem.duration = sum(tandem)/5,
                         tandem.distance = mean(dist[tandem]))
        res.speed.ind <- rbind(res.speed.ind, df)
        
        df <- data.frame(species, treat, source, name,  sex = "T-M",
                         mean = mean(m.vec), median = median(m.vec), 
                         move.mean = mean(m.move.speed), 
                         move.median = median(m.move.speed),
                         pause.duration = sum(m.pause)/sum(tandem),
                         angle.var = var(angle_cal2(d[,2,1],d[,3,1],dim(d)[1])),
                         tandem.duration = sum(tandem)/5,
                         tandem.distance = mean(dist[tandem]))
        res.speed.ind <- rbind(res.speed.ind, df)
        
        
      } else {
        print("error"); break;
        df <- data.frame(species, treat, source, name,  sex = "T",
                         #mean = mean(f.vec), median = median(f.vec), mod = mod.speed,
                         #move.mean = mean(move.speed), 
                         f.move.median = NA,
                         m.move.median = NA,
                         f.pause.duration = NA,
                         m.pause.duration = NA,
                         f.angle.var <- NA,
                         m.angle.var <- NA,
                         tandem.duration = sum(tandem)/5,
                         tandem.distance = NA)
        res.speed.ind <- rbind(res.speed.ind, df)
      }
    }
    
    ## save
    {
      save(res.tan.time,          file = paste0(iodir,"res.tan.time.rda"))
      save(res.sep.time,          file = paste0(iodir,"res.sep.time.rda"))
      
      res.tandem.speed.ind  <- res.speed.ind
      save(res.tandem.speed.ind,  file = paste0(iodir,"tandem-ind-speed.rda"))
      
      Res.tandem <- Res
      save(Res.tandem,            file = paste0(iodir,"Res.tandem.rda"))
    }
  }
}
#-------------------------------------------------------------------------#