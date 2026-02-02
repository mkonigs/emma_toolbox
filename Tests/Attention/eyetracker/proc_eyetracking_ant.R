#+eval=FALSE

# to fix the source problem when installing packages
options(repos=list(CRAN="http://cran.rstudio.com/"))

####
#Title   : fancy_parallel_Bart_3.1.R
#Topic   : Prepare gaze data for multiple subjects in parallel
#Version : 3.1
#Author  : Bart Relou 
#Date    : 2018-10-17
#Email   : 
#Links   : -
#Remarks : 

# Set up work environment --------
# Clear workspace
#rm(list = ls())

# Install packeges as required
#if (!require(installr)){
#  install.packages("installr")
#  require(installr)
#}

if (!require(devtools)){
  install.packages("devtools")
  require(devtools)
}

library(devtools)

#if (find_rtools() != TRUE){
#  
#  install.Rtools()
#}

if (!require(scales)){
install.packages("scales")
}

if (!require(lazyeval)){
install.packages("lazyeval")
}

#if (!require(nbs)){
#install.packages("nbs")
#}

if (!require(imputeTS)){
  install.packages("imputeTS")
  require(imputeTS)
}

#if (!require(Rtools)){
#  install.packages("Rtools")
#  require(Rtools)
#}

#if (!require(nbs)){
#  devtools::install_github("mikeniemant/nbs")
#}

if (!require(readr)){
  install.packages("readr")
  require(readr)
}

if (!require(dplyr)){
  install.packages("dplyr")
  require(dplyr)
}

if (!require(zoo)){
  install.packages("zoo")
  require(zoo)
}

if (!require(ggplot2)){
  install.packages("ggplot2")
  require(ggplot2)
}

if (!require(doParallel)){
  install.packages("doParallel")
  require(doParallel)
}

if (!require(gazepath)){
  install.packages("gazepath")
  require(gazepath)
}

if (!require(tidyverse)){
  install.packages("tidyverse")
  require(tidyverse)
}

if (!require(lubridate)){
  install.packages("lubridate")
  require(lubridate)
}

if (!require(mime)){
  install.packages("mime")
  require(mime)
}

#load packages
library(devtools)
library(imputeTS)
#library(Rtools)
#library(nbs)
library(readr)
library(dplyr)
library(zoo)
library(ggplot2)
library(doParallel)
library(gazepath)
library(tidyverse)
library(lubridate)
library(mime)

# Create PATHS object
OFFLINE <- T
PATH1 <- c("../../../Data/Raw/Attention/")
PATH2 <- c("../../../Data/Raw/Attention/eyetracker/")


if(OFFLINE) {
  main.path <- PATH1
} else { 
  main.path <- PATH2
}

PATHS <- main.path


# Define global variables
TEST.BOOL <- T
SAVE.BOOL <- T

# Conversion from px to cm
width_cm <- 34.52
width_px <- 1550
height_cm <- 19.43
height_px <- 900

toptargetpos_px <-  (0.5*height_px) + 50 #target appears 50 pixels above or below center of the screen
bottomtargetpos_px <- (0.5*height_px) - 50
width_pxcm <- width_cm / width_px
height_pxcm <- height_cm / height_px  

# Define target box coordinates
minx <- (width_px/2)-10
maxx <- (width_px/2)+10

mintopy <- toptargetpos_px-10 
maxtopy <- toptargetpos_px+10

minboty <- bottomtargetpos_px+10
maxboty <- bottomtargetpos_px-10

# Define vector of all subject IDs
#ant.subject.v <- list.files(PATHS$data)[grep(pattern = "attention-subject-*", 
#                                           x = list.files(PATHS$data))]
#gaze.subject.v <- list.files(PATHS$data)[grep(pattern = "gazedata_", 
#                                              x = list.files(PATHS$data))]


ant.subject.v <- paste(PATH1, "attention-subject-",subj,".csv", sep="")
gaze.subject.v <- paste(PATH2, "gazedata_",subj, sep="")

# Extract string names so user is left with subject id numbers only                                  
#subject.ids.ant <- sort(as.numeric(gsub("[A-z \\.\\ -]","", ant.subject.v)))
#subject.ids.gaze <- sort(as.numeric(gsub("[A-z \\.\\ -]","", gaze.subject.v)))

# only analyse the data files where both GAZE & ANT are available
#subject.ids.v <- subject.ids.gaze[which(subject.ids.gaze %in% subject.ids.ant == T)]
#subject.ids.v <- subject.ids.v[which(subject.ids.v > 100000 & subject.ids.v < 1000000)]

# Complete list of subject.ids.v based on multiple factors. Use the above to analyse only subject ids within range.
#subject <- read.table(file = "/Users/bart/Documents/VU Amsterdam/Master Research Project/data analyse/Bestanden voor Bart/subject.ids.v.txt",sep = "", header = T)
#subject.ids.v <- subject$x

#if(TEST.BOOL) {
#  subject.ids.v <- subject.ids.v[1:2]  
#}

# Check if all files of the defined subjects exist
#for(subject.id in subject.ids.v) {
ANT.DATA.FILENAME <- paste0(PATH1, "attention-subject-", subj, ".csv")
GAZE.DATA.FILENAME <- paste0(PATH2, "gazedata_", subj, ".txt")
  
  # Sanity check: check if both files exist
  if(file.exists(ANT.DATA.FILENAME)) {
    print(paste0("ANT data found for subject: ", subj))
    ant.exists <- 1
  } else{
    ant.exists <-0
  }

if(file.exists(GAZE.DATA.FILENAME)) {
  print(paste0("Gaze data found for subject: ", subj))
  gaze.exists <- 1
} else{
  gaze.exists <-0
}

all.exists <- ant.exists * gaze.exists

if(all.exists ==1 ) {
  print(paste0("All data found. Starting analysis for subject: ", subj))
} 


# Define all functions that are required for the looping ----
prepareGazeData <- function(GAZE.DATA.FILENAME){
  if(TEST.BOOL) {
    cat("Prepare gaze data")
  }

  print(getwd())
  # Determine when validation is reported
  report.val.L <- grep("deviationXLeft", readLines(GAZE.DATA.FILENAME))
  report.val.R <- grep("deviationXRight", readLines(GAZE.DATA.FILENAME))
  
  # Import validation data
  
  if (length(report.val.L)== 0){
  print("No validation data")
    deviation.x.L <- NA
    deviation.x.R <- NA
    deviation.y.L <- NA
    deviation.y.R <- NA
    
  }else {
    print("Validation data present")
    #Import left validation
    dev.L <- suppressWarnings(suppressMessages(read_table2(file = GAZE.DATA.FILENAME, 
                                                         col_names = F,
                                                         progress = F )))
    dev.L <- dev.L[report.val.L,]
    
    deviation.x.L <- dev.L[1,2]
    deviation.y.L <- dev.L[1,4]
    
    #Import right validation
    dev.R <- suppressWarnings(suppressMessages(read_table2(file = GAZE.DATA.FILENAME, 
                                                           col_names = F,
                                                           progress = F )))
    dev.R <- dev.R[report.val.R,]
    
    deviation.x.R <- dev.R[1,2]
    deviation.y.R <- dev.R[1,4]
    
  }
  
  deviation <- data.frame(deviation.x.L, deviation.y.L, deviation.x.R, deviation.y.R)
  
  #assign colnames to deviation
  colnames(deviation) <- c("deviation.x.L", "deviation.y.L", "deviation.x.R", "deviation.y.R")
  
  
  # Determine row when relevant GAZE data starts
  first.index <- grep("Gaze Data", readLines(GAZE.DATA.FILENAME))[1]


  
  # Import gaze data
  G.t <- suppressWarnings(suppressMessages(read_table2(file = GAZE.DATA.FILENAME, 
                                                       col_names = F,
                                                       skip = first.index-1, 
                                                       progress = F )))
  
  
  
  
  # Add space to OpenSesameTimeStamp so ncols is 19 instead of 18
  test.me <- ncol(G.t) <= 18
  if(test.me) {
    cat("Adding space infront of OST if necessary")
    
    tx <- readLines(GAZE.DATA.FILENAME)
    tx2 <- gsub(pattern = "OpenSesameTimeStamp", replacement = " OpenSesameTimeStamp", x = tx)
    writeLines(tx2, con = paste0(PATH2, "gazedata_", subject.id, ".txt"))
    
    # Determine row when relevant GAZE data starts
    first.index <- grep("IviewTimestamp", readLines(GAZE.DATA.FILENAME))[1]
    
    # Import gaze data
    G.t <- suppressWarnings(suppressMessages(read_table2(file = GAZE.DATA.FILENAME, 
                                                         col_names = F,
                                                         skip = first.index-1, 
                                                         progress = F)))
  }
  
  
  # Subset relevant columns for gaze data and transform object to data frame
  G.df <- as.data.frame(G.t[, c(5, 7, 9, 10, 12, 13, 15, 17, 19)])
  
  colnames(G.df) <- c("Iview","OST", "GLX", "GLY", "GRX", 
                      "GRY", "DisR", "PupL", "PupR")
  
  
  # Change NA's to 0 in Iview column so data isn't deleted 
  G.df$Iview <- as.numeric(G.df$Iview)
  
  for (i in 1:nrow(G.df)){
    if (is.na(G.df$Iview[i]) == TRUE){
      if (i < 3){
        
        G.df$Iview[i] <- G.df$Iview[i-1] + 4000
        
        } else {
      
        G.df$Iview[i] <- G.df$Iview[i-1] + (G.df$Iview[i-1] - G.df$Iview[i-2])
        
        }
    }
}
  #G.df$Iview[is.na(G.df$Iview)] <- 0
  #G.df$Iview <- na_interpolation(G.df$Iview)
  
  # Create timeline by extracting first index from every timepoint and mutate to seconds
  # Lines with na dont influence this, since Iview in that case was also NA. This is not the case 
  # with all the OST NA's.
  G.df$timeline <- (G.df$Iview - G.df$Iview[1])/1000000
  
  # Calculate Hz
  Hz.v <- round(1/diff(G.df$timeline))
  G.df$Hz <- c(Hz.v[1], Hz.v)
  
  G.df[is.na(G.df)] <- 0
  
  # Replace missing values of one eye with the coordinates of the other, so mean GX and GY isnt influenced
  l_0 = G.df$GLX == 0 & G.df$GLY == 0
  r_0 = G.df$GRX == 0 & G.df$GRY == 0
  
  #also for pupils
  pupl_0 = G.df$PupL == 0 
  pupr_0 = G.df$PupR == 0 
  
  G.df[l_0 & ! r_0, c("GLX", "GLY")] = G.df[l_0 & ! r_0, c("GRX", "GRY")]
  G.df[r_0 & ! l_0, c("GRX", "GRY")] = G.df[r_0 & ! l_0, c("GLX", "GLY")]
  
  G.df$PupL[G.df$PupL == 0] <- G.df$PupR[G.df$PupL == 0]
  G.df$PupR[G.df$PupR == 0] <- G.df$PupL[G.df$PupR == 0]
  
  
  # interpolate missing values
  G.df[G.df==0] <- NA
  
  missing_data <<- round((sum(is.na(G.df$GLX))/length(G.df$GLX))*100,1)
  
  G.df$PupL <- na_interpolation(G.df$PupL)
  G.df$PupR <- na_interpolation(G.df$PupR)
  G.df$GLX <- na_interpolation(G.df$GLX)
  G.df$GLY <- na_interpolation(G.df$GLY)
  G.df$GRX <- na_interpolation(G.df$GRX)
  G.df$GRY <- na_interpolation(G.df$GRY)
  G.df$DisR <- na_interpolation(G.df$DisR)
  
  G.df <- na.omit(G.df)
  
  # calculate mean X and Y coordinate 
  G.df <- G.df %>%
    mutate(GX = rowMeans(.[, c("GLX", "GRX")]), 
           GY = rowMeans(.[, c("GLY","GRY")]))
  
  #also for pupils
  G.df <- G.df %>%
    mutate(Pup = rowMeans(.[, c("PupL", "PupR")]))
  
  # Mutate distance from mm to cm
  G.df$DisR <- G.df$DisR / 10
  
  # interpolate DisR when DisR == 0
  DisR <- G.df[, "DisR"]
  DisR[DisR == 0] <- NA
  DisR <- na.approx(DisR, na.rm = F)
  G.df[, "DisR"] <- DisR
  
  # Correlation between eyes
  G.df$cor.x <- cor(G.df$GLX, G.df$GRX, use = "complete.obs")
  G.df$cor.y <- cor(G.df$GLY, G.df$GRY, use = "complete.obs")
  
  # Remove rows with NA's caused by strings "Start test trial #'
  #G.df <- na.omit(G.df$Iview)

  # Add deviation to G.df
  G.df <- data.frame(deviation, G.df)  

  
  # Return result
  return(G.df)
}


 prepareAntData <- function(ANT.DATA.FILENAME) {
  if(TEST.BOOL) {
    cat("Prepare ant data")
  }
  
  # Import ANT data as A.df
  A.t <- read.csv(file = ANT.DATA.FILENAME, 
                   header = T, 
                   sep = ",")
  
  # Subset relevant columns
  A.df <- A.t[, c("time_fixation_2", "time_cue_2", "time_fixation_3", 
                  "time_target_2", "time_new_logger", "cue_condition", 
                  "target_condition", "target_pos", "response_time_new_srbox_2", 
                  "correct_new_srbox_2")]
  
  # Factorize the target position factor
  A.df$target_pos <- factor(A.df$target_pos, labels = c("down","up"))

  # Factorize the cue condition
  A.df$cue_condition <- factor(A.df$cue_condition, labels = c("central", "no cue", "spatial"))
  
  # Change column names
  names(A.df) <- c("time.fixation.1", "time.cue", "time.fixation.2",
                   "time.target", "time.new.logger",
                   "cue.condition", "target.condition", "target.position", 
                   "response.time", "correct.response")
                   
  
  # Delete rows with NA's (test trials)
  A.df <- na.omit(A.df)
  
  # Return result
  return(A.df)
}

combineGazeAnt <- function(G.df, A.df) {
  if(TEST.BOOL) {
    cat("Combine ant + gaze")
  }
  
  G.df$block <- "isi"
  
  # Iterate through all rows of A.df
  for(trial.i in seq(nrow(A.df))) {
    
    # Define start time point
    start.trial.tp <- A.df[trial.i, "time.cue"]-250 #fixation of a trial is 250 ms, so cue - 250 is start fixation period
    cue.tp <- A.df[trial.i, "time.cue"]
    fix.2.tp <- A.df[trial.i, "time.fixation.2"] #fixation.2
    target.tp <- A.df[trial.i, "time.target"]
    resp.tp <- target.tp + A.df[trial.i, "response.time"] #Target period ands after response(time)
    end.trial.tp <- target.tp + 1700
    #end.trial.tp <- target.tp + A.df[trial.i, "response.time"] #Target period ands after response(time)
    
    # Find the row index with the closest time point
    start.trial.gaze.ri <- which(abs(G.df$OST - start.trial.tp) %in%
                                   min(abs(G.df$OST - start.trial.tp)))[1]
    
    cue.gaze.ri <- which(abs(G.df$OST - cue.tp) %in% 
                           min(abs(G.df$OST - cue.tp)))[1]
    
    fix.2.gaze.ri <- which(abs(G.df$OST - fix.2.tp) %in% 
                             min(abs(G.df$OST - fix.2.tp)))[1]
    
    target.gaze.ri <- which(abs(G.df$OST - target.tp) %in% 
                              min(abs(G.df$OST - target.tp)))[1]
    
    resp.trial.gaze.ri <- which(abs(G.df$OST - resp.tp) %in% 
                                 min(abs(G.df$OST - resp.tp)))[1] 
    
    end.trial.gaze.ri <- which(abs(G.df$OST - end.trial.tp) %in% 
                                 min(abs(G.df$OST - end.trial.tp)))[1] 
    
    # Define number, trial, target condition and target position for this trial
    G.df[start.trial.gaze.ri:end.trial.gaze.ri, "trial.number"] <- trial.i
    G.df[start.trial.gaze.ri:end.trial.gaze.ri, "trial.block"] <- ceiling(trial.i/72)
    G.df[start.trial.gaze.ri:end.trial.gaze.ri, "trial.block.number"] <- ((trial.i-1)%%72)+1
    G.df[start.trial.gaze.ri:end.trial.gaze.ri, "cue.condition"] <- A.df[trial.i, "cue.condition"]
    G.df[start.trial.gaze.ri:end.trial.gaze.ri, "target.condition"] <- A.df[trial.i, "target.condition"]
    G.df[start.trial.gaze.ri:end.trial.gaze.ri, "target.position"] <- A.df[trial.i, "target.position"]
    G.df[start.trial.gaze.ri:end.trial.gaze.ri, "response.time"] <- A.df[trial.i, "response.time"]
    G.df[start.trial.gaze.ri:end.trial.gaze.ri, "correct.response"] <- A.df[trial.i, "correct.response"]
    
    # Define block, by first filling all entries with fixation
    

    G.df[start.trial.gaze.ri:end.trial.gaze.ri, "block"] <- "fixation.1"
    
    # Get cue and target time points
    start.cue.tp <- A.df[trial.i, "time.cue"]
    start.fix.2.tp <- A.df[trial.i, "time.fixation.2"]
    start.target.tp <- A.df[trial.i, "time.target"]
    
    # Find row index of time points that lies closests to the cue and target tps
    start.cue.gaze.ri <- which(abs(G.df$OST - start.cue.tp) %in% 
                                 min(abs(G.df$OST - start.cue.tp)))[1]
    start.fix.2.gaze.ri <- which(abs(G.df$OST - start.fix.2.tp) %in% 
                                   min(abs(G.df$OST - start.fix.2.tp)))[1]
    start.target.gaze.ri <- which(abs(G.df$OST - start.target.tp) %in% 
                                    min(abs(G.df$OST - start.target.tp)))[1]
    
    # Fill missing entries
    G.df[start.cue.gaze.ri:start.fix.2.gaze.ri, "block"] <- "cue"
    G.df[start.fix.2.gaze.ri:start.target.gaze.ri, "block"] <- "fixation.2"
    G.df[start.target.gaze.ri:resp.trial.gaze.ri, "block"] <- "target"
    G.df[resp.trial.gaze.ri:end.trial.gaze.ri, "block"] <- "isi"
    
    # Sanity check
    # View(G.df[(start.trial.gaze.ri-1):(end.trial.gaze.ri+1), ])
  }
  
  # Remove rows with NA's 
  #G.df <- na.omit(G.df)
  
  return(G.df)
}

calculateVelocity <- function(G.df) {
  if(TEST.BOOL) {
    cat("Calculate velocity")
  }
  # raw calculation of velocity; includes high frequency noise 
  # Calculate dx dy for every measurement in px and transform into cm
  xdiff <- diff(G.df$GX) * width_pxcm
  ydiff <- diff(G.df$GY) * height_pxcm
  
  # Calculate length of displacement 
  disp.length <- sqrt(xdiff^2 + ydiff^2) 
  
  # make short variable name for distance to the screen
  disscreen <- G.df$DisR
  
  #calculate velocity per second 
  velocity <- atan(disp.length / disscreen[1:length(disp.length)]) * (180/pi) * 1/diff(G.df$timeline)
  
  # Add the velocity to the G.df 
  G.df$velocity <- c(velocity[1], velocity)
  
  # Compute acceleration
  timediff <- diff(G.df$timeline)
  timediff <- c(timediff[1],timediff)
  velocity.diff <- diff(G.df$velocity) 
  velocity.diff <- c(velocity.diff[1], velocity.diff)
  
  # Compute acceleration
  G.df$acceleration <- velocity.diff / timediff
  
  # Calculate amplitude 
  amplitude <- atan(disp.length / disscreen[1:length(disp.length)]) * (180/pi)
  G.df$amplitude <- c(amplitude[1], amplitude)
  # Return G.df
  return(G.df)
}

 calculateVelocity2 <- function(G.df) {
  if(TEST.BOOL) {
    cat("Calculate low pass filtered velocity (central differentiation)")
  }
  # determine amount of lag. 2 is i-1, i+1
  lag = 2
  
  # Calculate dx dy for every measurement in px and transform into cm
  xdiff <- base::diff(G.df$GX, lag = lag) * width_pxcm
  ydiff <- base::diff(G.df$GY, lag = lag) * height_pxcm
  
  # Take the distance to the screen
  disscreen <- G.df$DisR
  
  # Calculate length of displacement on screen
  disp.length <- sqrt(xdiff^2 + ydiff^2) 
  
  # take the difference between timepoints (delta t) between i-1 and i+1
  difftimeline <- base::diff(G.df$timeline, lag = lag)
  
  # Calculate the angle in radians the eye has made, convert to degrees and multiple by Hz 
  #to get the velocity in degrees/sec 
  velocity <- atan(disp.length / disscreen[1:length(disp.length)]) * (180/pi) * (1/difftimeline)
  velocity <- c(rep(velocity[3], lag), velocity)
  
  # Add velocity to the G.df data.frame 
  G.df$velocity2 <- velocity
  
  # Compute acceleration
  velocity.diff <- diff(G.df$velocity2) 
  velocity.diff <- c(velocity.diff[1], velocity.diff)
  
  # Compute acceleration
  G.df$acceleration2 <- velocity.diff / c(rep(difftimeline[3], lag), difftimeline)
  
  # Calculate amplitude 
  amplitude <- atan(disp.length / disscreen[1:length(disp.length)]) * (180/pi)
  G.df$amplitude <- c(rep(amplitude[1], lag), amplitude)
  
  # Return G.df
  return(G.df)
}
butterfilter <- function(G.df) {
  #WATCH OUT! signal packages interferes with dplyr package. Unload signal package when done
  bf <- butter(4, 0.8, "low", plane = "z") #4th order low-pass butterworth filter with freq of 80 Hz 
  
  # replace velocity by filtered velocity 
  G.df$velocity <- filtfilt(bf, G.df$velocity) 
  
  # Return result
  return(G.df)
  
  
  
  
}
prepareFinalData <- function(G.df) {
  if(TEST.BOOL) {
    cat("Finalize data")
  }
  
  # Remove 5 data frames before and after blinks
  zeros <- G.df 
  zeros <- rowSums(G.df[, c("GLX", "GLY", "GRX", "GRY")]) == 0
  zeros <- c(T,T,T,T,T,T, zeros, T,T,T,T,T,T)
  rolling <- diff(cumsum(zeros), 13) 
  G.df <- G.df[rolling == 0, ]
  
  # Remove blinks and flickers from the eye 
  G.df <- G.df %>% 
    filter(GLX > 0 & GRX > 0 & GLY > 0 & GRY > 0) %>% 
    filter(velocity2 <= 700 & acceleration2 < 15000)
  
  # Remove the incomplete rows
  G.df <-  G.df[(rowSums(is.na(G.df[,5:ncol(G.df)]))==0),]
  
  # Now remove the rows with Iview == 0
  G.df <- G.df[which(G.df$Iview != 0), ]
  
  # Return result
  return(G.df)
}

# detecting of fixation and saccade 
fixationANDsaccade <-  function(speed, thres_vel, thres_dur, Hz){
    thres_dur <- thres_dur * (Hz / 1000)
    fixsac <- ifelse(speed > thres_vel, 's', 'f')
    ## Set a minimum saccade duration of 10 ms
    rle <- rle(fixsac)
    fixsac[cumsum(rle$lengths)[which(rle$lengths < 10 * (Hz/1000) & rle$values == 's')]] <- 'f'
    
    while(length(which(rle$lengths < 10 * (Hz/1000) & rle$values == 's')) != 0){
      rle <- rle(fixsac)
      fixsac[cumsum(rle$lengths)[which(rle$lengths < 10 * (Hz/1000) & rle$values == 's')]] <- 'f'
    }
    classify <- numeric()
    for(i in 1:length(rle$values)){
      if(is.na(rle$values[i])){
        classify <- c(classify, rep(NA, rle$lengths[i]))
      } else{
        if(rle$values[i] == 'f' & rle$lengths[i] >= thres_dur){
          classify <- c(classify, rep('f', rle$lengths[i]))
        }
        if(rle$values[i] == 'f' & rle$lengths[i] < thres_dur){
          classify <- c(classify, rep('u', rle$lengths[i]))
        }
        if(rle$values[i] == 's'){
          classify <- c(classify, rep('s', rle$lengths[i]))
        }
      }
    }
    return(classify)
  }

# Calculate threshold estimation & determine position of saccades and length 
SaccadeDetection <- function(G.df) {
  if(TEST.BOOL) {
    cat("Threshhold estimation and saccade identifcation")
  }
  
  # Set up empty vector for thresholds
  saccthreshold.v <- NA
  
  # Set up empty column in G.df 
  G.df$sacc <- NA
  
  for(trial.i in unique(G.df$trial.number)) {
    speed <- G.df$velocity2[G.df$trial.number == trial.i] 
    Hz <- round(mean(G.df$Hz[G.df$trial.number == trial.i], trim=.05))
    saccthreshold.v <- suppressWarnings(suppressMessages(gazepath:::Mould_vel(speed, Hz, plot = F)))
    saccthreshold.v <- ifelse(saccthreshold.v <= 10 | is.na(saccthreshold.v), 10, saccthreshold.v) #set minimum of 10 to filter out microsaccades 
    G.df$sacc[G.df$trial.number == trial.i] <- gazepath:::fixationANDsaccade(speed, 
                                                                  thres_vel = saccthreshold.v,
                                                                  thres_dur = 10,  
                                                                  Hz)
  }
  
  #detach("package:plyr", unload=TRUE) 
  
  # Calculate the duration of the saccade and fixation periods 
  run.length <- rle(G.df$sacc == "s")
  G.df$indexer <- rep(1:length(run.length$lengths), run.length$lengths)
  G.df <- G.df %>%
    group_by(indexer) %>%
    mutate(duration = ((Iview[n()] - Iview[1])/1000)+(1000/Hz)) #Added 4 ms, because n() looks at start time last saccade, not end time(n+1)
  
  G.df <- ungroup(G.df)
  
  return(G.df)
}

# Calculate all variables --
calculateAll <- function(G.df, subject.id) {
  if(TEST.BOOL) {
    cat(paste0("calculating all variables for: ", subject.id, " "))
  }
  
  # Correlation between left and right eye
  correlation.df <- data.frame(cor.x = G.df$cor.x[1], cor.y = G.df$cor.y[1])
  
  # number of saccades made in total
  nr_sac <- G.df %>% 
    filter(sacc == "s", block == "target", correct.response == 1,
           cue.condition %in% c("no cue", "central")) %>%
    distinct(indexer, .keep_all = T) %>% 
    summarise(nr_sac_total = length(unique(indexer)))
  
  # calculate nr of saccades per target.condition 
  nr_sac_condition <- G.df %>% 
    filter(sacc == "s", correct.response == 1, cue.condition %in% c("central", "no cue"), block == "target") %>% 
    distinct(indexer, .keep_all = T) %>% 
    summarise(sac.con = sum(target.condition == "congruent"), sac.incon = sum(target.condition == "incongruent"),
              sac.neut = sum(target.condition == "neutral"))
  
  # count number of trials per trial condition 
  nr_trials_condition <- G.df %>% 
    filter(sacc == "s", correct.response == 1, cue.condition %in% c("central", "no cue"), block == "target") %>% 
    select(trial.number, target.condition) %>% 
    group_by(trial.number) %>% 
    distinct(trial.number, .keep_all = T) %>% 
    summarise(con.trials = sum(target.condition == "congruent"),
              incon.trials = sum(target.condition == "incongruent"),
              neut.trials = sum(target.condition == "neutral"))
    
  nr_trials_condition <- summarise(nr_trials_condition, congruent.trials.n = sum(con.trials), incongruent.trials.n = sum(con.trials),
              neutral.trials.n = sum(neut.trials))
  
  # Ratio between nr of saccades per trial.condition
  new_name <- nr_sac_condition
  colnames(new_name) <- c("sac.con.ratio", "sac.incon.ratio", "sac.neut.ratio")
  sac_trial_ratio <- new_name/nr_trials_condition #had to change name, otherwise same variables as nr_sac_condition

  # Stability of the GY eye position during the second fixation block with no cue or central cue conditions
  # Outcome [1]
  stability <- G.df %>%
    filter(cue.condition %in% c("no cue", "central"), block == "fixation.2", sacc == "f") %>%
    group_by(trial.number) %>%
    summarise(stabilityGY = sd(GY), stabilityGX = sd(GX)) %>%
    summarise(n_trials.1 = length(unique(trial.number)), sdfixationGY.1 = mean(stabilityGY), 
              sdfixationGX.1 = mean(stabilityGX))
  
  # Alternative stability looking at stability of X&Y during the whole trial
  # Outcome [2]
  meanstability <- G.df %>% 
    group_by(trial.number) %>% 
    summarise(sdGY = sd(GY), sdGX = sd(GX)) %>% 
    summarise(overallGY.2 = mean(sdGY), overallGX.2 = mean(sdGX))
  # mean per target condition
  
  stabilitytarget <- G.df %>%
    group_by(trial.number, target.condition) %>%
    summarise(stabilityGY = sd(GY), stabilityGX = sd(GX)) %>%
    group_by(target.condition) %>% 
    summarise(sdGY.2 = mean(stabilityGY), sdGX.2 = mean(stabilityGX)) %>% 
    gather(position, sddev, -target.condition) %>% 
    unite(target.position, c("target.condition", "position"), sep = ".", remove = TRUE) %>% 
    spread(target.position, sddev) 
  

  
  # ggplot(stability) +
  #   geom_col(aes(cue.condition, stabilityGY)) +
  #   geom_col(aes(cue.condition, stabilityGX, colours = "red")) +
  #   theme_minimal()
  
  # Calculate the avarage velocity during a saccade 
  ## only look at real saccades towards targets instead of saccades towards cue(which dissapears before eye can reach)
  # Outcome [3]
  meanSacVel <- G.df %>%
    filter(sacc == "s", block == "target", correct.response == 1,
           cue.condition %in% c("no cue", "central")) %>%
    group_by(trial.number, indexer,target.condition) %>%
    summarise(SacVel = mean(velocity2, na.rm = T))
    
  meanSacVel <- ungroup(meanSacVel)
    
  
  meanSacVel <- summarise(meanSacVel, n_trials.3 = length(unique(trial.number)), n.saccades.3 = length(trial.number),
              meanSacVel.3 = mean(SacVel))
  
  
  # ggplot(meanSacVel, mapping = aes(x = target.condition, y = meanSacVel, colour = cue.condition)) +
  #   geom_point(size = 3) +
  #   theme_bw(base_size = 13) 
  
  # Calculate average peak velocity of the saccades
  ## only look at real saccades towards targets, so not in the fixation.2 block
  # Outcome [4]
  meanPeakSacVel <- G.df %>% 
    filter(correct.response == 1, sacc == "s", block == "target", cue.condition %in% c("no cue", "central")) %>% 
    group_by(cue.condition, trial.number, target.condition) %>% 
    summarise(max_v = max(velocity2))

  meanPeakSacVel <- ungroup(meanPeakSacVel)  
         
  meanPeakSacVel <- summarise(meanPeakSacVel, n_trials.4 = length(unique(trial.number)), avgPeakVelSac.4 = mean(max_v))
  
  # Calculate average peak acceleration of the saccades
  # Outcome [5]
  meanPeakAcc <- G.df %>% 
    filter(correct.response == 1, sacc == "s", block == "target", cue.condition %in% c("no cue", "central")) %>% 
    group_by(trial.number) %>% 
    summarise(max_a = max(acceleration2))
  
  meanPeakAcc <- ungroup(meanPeakAcc)
    
  meanPeakAcc <- summarise(meanPeakAcc, n_trials.5 = length(unique(trial.number)), mean_peak_acc.5 = mean(max_a))
  
  #Calculate mean saccade duration and nr of saccades per target condition
  # Outcome [6]
  mean_sac_duration <- G.df %>% 
    filter(sacc == "s", correct.response == 1, cue.condition %in% c("central", "no cue"), block == "target") %>% 
    group_by(trial.number) %>% 
    summarise(meandur = mean(duration)) %>% 
    summarise(n_trials.6 = length(unique(trial.number)), mean_sac_duration.6 = mean(meandur))
  
  
  # Calculate reaction times when there is no cue or a central cue (both fixations stay in the middle)
  eye_rt <- G.df %>%
    filter(cue.condition %in% c("central", "no cue"), correct.response == 1) %>%
    group_by(trial.number, target.condition, cue.condition) %>%
    summarise(Rt.eye = (Iview[block == "target" & sacc == "s"][1] - Iview[block == "target"][1])/1000) %>%
    filter(Rt.eye > 100) %>% 
    ungroup()
  
  # extract physical reaction time of pressing button
  button_rt <- G.df %>% group_by(trial.number) %>% 
    group_by(trial.number, target.condition) %>%
    summarise(Rt.but = response.time[1])
  
  # Mean button response 
  button_rt_mean.7 <- mean(button_rt$Rt.but)
  
  # Mean button response time per target
  button_rt_target.7 <- button_rt %>% 
    group_by(target.condition) %>% 
    summarise(mean.Rt.but = mean(Rt.but)) %>% 
    gather(processtime, ms, -target.condition) %>% 
    unite(target, c("target.condition", "processtime")) %>% 
    spread(target, ms)
  
  # Calculate the processing time 
  # Outcome [7]
  # mean processtime 
  
  mean_processtime <- inner_join(button_rt, eye_rt, by = "trial.number") %>%
    group_by(trial.number) %>% 
    summarise(ptime = Rt.but - Rt.eye) %>% 
    summarise(mean_processtime = mean(ptime))
  
  # Processtime per target condition
  target_processtime <- inner_join(button_rt, eye_rt, by = "trial.number") %>%
    group_by(trial.number, target.condition.x, cue.condition) %>%
    summarise(ptime = Rt.but - Rt.eye)%>%
    group_by(target.condition.x) %>%
    summarise(meanptime.7 = mean(ptime)) %>% 
    gather(processtime, ms, -target.condition.x) %>% 
    unite(target, c("target.condition.x", "processtime")) %>% 
    spread(target, ms)
  
  # Calculate mean eye_rt
  eye_rt_mean.7 <- mean(eye_rt$Rt.eye)
  
  # Calculate mean eye_rt per target condition 
  eye_rt_target.7 <- eye_rt %>% 
    group_by(target.condition) %>% 
    summarise(eye_rt_target.7 = mean(Rt.eye)) %>% 
    gather(responsetime, ms, -target.condition) %>% 
    unite(target, c("target.condition", "responsetime")) %>% 
    spread(target, ms)
    
  
  # ggplot(processtime, aes(x = factor(trial.block), y = meanptime, colour = target.condition)) + 
  #   geom_point(size = 3) +
  #   theme_bw(base_size = 13) +
  #   labs(x = "trial block", y = "mean processing time [ms]", title = "Processing time of the brain")
  
  # Stability of the fixation period on target 
  # Outcome [8]
  fixation_stab <- G.df %>%
    filter(correct.response == 1, block == "target") %>% 
    group_by(trial.number) %>%
    summarise(sd_fix_GY = sd(GY[sacc == "f" & indexer == max(indexer)]), sd_fix_GX = sd(GX[sacc == "f" & indexer == max(indexer)])) %>% 
    summarise(mean_sd_GY.8 = mean(sd_fix_GY, na.rm = T), mean_sd_GX.8 = mean(sd_fix_GX, na.rm = T))
  
  # Calculate time before eye is within 10 procent of the target {only 144 trials, since spatial is deleted}
  eye_on_target <- G.df %>%
    filter(correct.response == 1, block == "target", cue.condition %in% c("no cue", "central")) %>% 
    mutate(ontarget = (target.position == "up" & GY >= toptargetpos_px-((toptargetpos_px-(0.5*height_px))*0.10)| 
                         (target.position == "down" & GY <= bottomtargetpos_px+(((0.5*height_px)-bottomtargetpos_px)*0.10)))) %>% # -5 and +5 is 10% of 50px.
    group_by(trial.number, cue.condition, target.condition) %>%
    summarise(time_on_target = Iview[ontarget == T][1]-Iview[block == "target"][1]) 
  
  # Calculate how many times target isnt reached
  didnt_reach_target = sum(is.na(eye_on_target$time_on_target))
  
  # Calculate summary of eye_on_target   
  eye_on_target <- eye_on_target %>% 
    group_by(cue.condition, target.condition) %>%
    summarise(meanTOT = mean(time_on_target, na.rm = T))

# Plot relation between amplitude and velocity
  # ggplot(total_amplitude, mapping = aes(x = totalamplitude, y = meanvelocity)) +
  #   geom_point() +
  #   geom_smooth(method = 'lm') +
  #   theme(text = element_text(size = 20)) +
  #   labs(x = "amplitude [deg]", y = "mean velocity [deg/s]", title = 'Saccade velocity-amplitude relationship plot')
  
  # Calculate total number of saccades 
  #sac_before_target <- G.df %>%
  #  filter(block == "target" & cue.condition %in% c("central", "no cue")) %>%
  #  distinct(indexer, .keep_all = T) %>%
  #  mutate(ontarget = (target.position == "up" & GY >= toptargetpos_px-((50/900)*height_px*.10)| 
  #                       (target.position == "down" & GY <= bottomtargetpos_px+((50/900)*height_px*.10))))
  
  # Calculate ratio between final eye_position of first saccade and target position
  # Outcome [9]
  eye_target_ratio_first <- G.df %>%
    filter(correct.response == 1, sacc == "s", block == "target", cue.condition %in% c("no cue", "central")) %>% 
    group_by(trial.number) %>% 
    summarise(pxpos = last(GY[sacc == "s" & indexer == min(indexer)]),
              ratio = if(first(target.position) == "up") ((pxpos-(0.5*height_px))/(toptargetpos_px-(0.5*height_px))) else ((pxpos-(0.5*height_px))/(bottomtargetpos_px-(0.5*height_px)))) %>%
    summarise(eye_ratio_first.9 = mean(ratio))
  
  # Calculate ratio between final eye_position of last saccade and target position
  # Outcome [10]
  eye_target_ratio_final <- G.df %>%
    filter(correct.response == 1, sacc == "s", block == "target", cue.condition %in% c("no cue", "central")) %>% 
    group_by(trial.number) %>% 
    summarise(pxpos = last(GY[sacc == "s"]),
              ratio = if(first(target.position) == "up") ((pxpos-(0.5*height_px))/(toptargetpos_px-(0.5*height_px))) else ((pxpos-(0.5*height_px))/(bottomtargetpos_px-(0.5*height_px)))) %>%
    summarise(eye_ratio_last.10 = mean(ratio))
  
  # Calculate total amplitude of each saccade
  # Outcome [11]
  total_amplitude <- G.df %>%
    filter(correct.response == 1, sacc == "s", block == "target", cue.condition %in% c("no cue", "central")) %>% 
    group_by(indexer, trial.number) %>%
    summarise(totalamplitude = sum(amplitude), totalduration = duration[1], meanvelocity = mean(velocity2)) 
  
  # Calculate linear model of the relation between amplitude and velocity
  # Outcome [12]
  modelfit <- lm(meanvelocity ~ totalamplitude, data = total_amplitude)
  x <- summary(modelfit)$coefficients[2,1]
  intercept <- summary(modelfit)$coefficients[1,1]
  linearmodel <- data.frame(coefficient.11 = x, intercept.11 = intercept)
  
  deviation <- data.frame(G.df$deviation.x.L[1], G.df$deviation.y.L[1], G.df$deviation.x.R[1], G.df$deviation.y.R[1])
  colnames(deviation) <- c("deviation.x.L", "deviation.y.L", "deviation.x.R", "deviation.y.R")
  
  output.df <- data.frame(subject.id, deviation, missing_data, correlation.df, nr_sac, nr_sac_condition, nr_trials_condition, 
                     sac_trial_ratio, stability, meanstability, stabilitytarget, meanSacVel, meanPeakSacVel, meanPeakAcc, 
                     mean_sac_duration, eye_rt_mean.7, eye_rt_target.7, button_rt_mean.7, button_rt_target.7, 
                     mean_processtime, target_processtime, fixation_stab, eye_target_ratio_first, eye_target_ratio_final,
                     linearmodel)
  
  
  # Is the subject.id a patient or not?
  #output.df <- output.df %>% 
  #  mutate(patient = if_else(subject.id < 100000, 0, 1))
  
  #outputall.df <- output.df %>% 
  #  mutate(sac_trial_ratio = nr_sac_total/n_trials.1)
  
  return(output.df)
}


# Run parallel code ----------------------------------------------------
#outputall.df <- foreach(i = seq(length(unique(subject.ids.v))), 
#                  .combine = rbind, 
#                  .packages = c("readr", "dplyr", "zoo", "tidyverse"), 
#                  .inorder = T,
#                  .multicombine = FALSE, #combine data frames into one list
#                  .verbose = T) %do%
#  {
    
    # Define subject id
    subject.id <- subj
    
    # Define file names
    ANT.DATA.FILENAME <- paste0(PATH1, "attention-subject-", subject.id, ".csv")
    GAZE.DATA.FILENAME <- paste0(PATH2, "gazedata_", subject.id, ".txt")
    
    # Import and prepare GAZE data
    G.df <- prepareGazeData(GAZE.DATA.FILENAME)
    
    # # Compute total test time
    # start <- first(G.df$timeline[G.df$trial.number == min(G.df$trial.number)])
    # end <- last(G.df$timeline[G.df$trial.number == max(G.df$trial.number)])
    # DA.df$total.test.duration <- (end-start)/60
    
    # Import and prepare ANT data
    A.df <- prepareAntData(ANT.DATA.FILENAME)
    
    # # Calculate single value for Hz
    # DA.df$avg.HZ <- round(mean(G.df$Hz, na.rm = T))
    
    # Combine both data frames
    G.df <- combineGazeAnt(G.df, A.df)
    
    # Compute velocity / accerlation
  #  G.df <- calculateVelocity(G.df)
    
    # Calculate velocity 2 (two directional bilateral filter)
    G.df <- calculateVelocity2(G.df)
    
    # Filter data with butterworth filter 
    ## pick velocity2 or butterworth, not both 
    ### G.df <- butterfilter(G.df)
    
    # Finalize data frame
    G.df <- prepareFinalData(G.df)
    
    # Detect saccades and fixations
    G.df <- SaccadeDetection(G.df)
    
    # Calculate all variables
    output.df <- calculateAll(G.df, subject.id)

 # }
 
