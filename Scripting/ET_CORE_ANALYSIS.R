
# settings
options(repos = c(CRAN = "http://cran.rstudio.com/"))
analysis <- "all" 

#set working directory
source(paste0(path.expand("~/"), "emma_toolbox/Scripting/setwd.R"))

# load packages
source("Scripting/ipak.R")
package_list <- c("rmarkdown", "knitr", "lubridate", "readr", "plyr")
ipak(package_list)

# get subj nr
args = commandArgs(trailingOnly=TRUE)
subjects <- args

#subjects <- "110001"

for (subj in subjects){

  print("Starting Emma Toolbox Analysis for subject:")
  print(subj)
  
  # create directory
  source(paste0(path.expand("~/"), "emma_toolbox/Scripting/createdir.R"))
  
  if(analysis!="domains"){
    
  ### ATTENTION ###
  testname <- "attention"
  source("Scripting/process_test.R")
  testname <- "attention"
  source("Scripting/process_database.R")
    
  ### VISUAL MEMORY ###
  testname <- "visualmem"
  source("Scripting/process_test.R")
  testname <- "visualmem"
  source("Scripting/process_database.R")

  ### VERBAL MEMORY ###
  testname <- "verbalmem"
  source("Scripting/process_test.R")
  testname <- "verbalmem"
  source("Scripting/process_database.R")
    
  ### VISUAL WM ###
  testname <- "visualwm"
  source("Scripting/process_test.R")
  testname <- "visualwm"
  source("Scripting/process_database.R")
  
  ### VERBAL WM ###
  testname <- "verbalwm"
  source("Scripting/process_test.R")
  testname <- "verbalwm"
  source("Scripting/process_database.R")
  
  ### VISUOMOTOR ###
  testname <- "visuomotor"
  source("Scripting/process_test.R")
  testname <- "visuomotor"
  source("Scripting/process_database.R")}
  
  ### DOMAINS ###
  # get center
  center <- read_file("../emma_toolbox_data/Settings/center.txt")
  center <- gsub("[\r\n]", "", center)
  
  # get pop
  pop <- read_file("../emma_toolbox_data/Settings/pop.txt")
  pop <- gsub("[\r\n]", "", pop)
  
  # start analysis
  source("Scripting/subscripts/domains.R")
  
  ### CHANGE ###
  #source("Scripting/change.R")
  
}

# exit
print("Analysis finished. Do not forget to back-up the data")
print("Press ENTER to exit")
drive <- readLines(file("stdin"),1)

