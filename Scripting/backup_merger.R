# set backup environment
bu <- TRUE

#set working directory
source(paste0(path.expand("~/"), "emma_toolbox/Scripting/setwd.R"))
base <- getwd()

# backup drive
print("Which drive is the backup connected to?, For example, G: ")
drive <- readLines(file("stdin"),1)
print(drive)

shell(paste0('robocopy ../emma_toolbox_data ', drive, '/emma_toolbox_data /E /XO', sep=""))
shell(paste0('robocopy ../emma_toolbox/Scripting ', drive, '/emma_toolbox/Scripting /E /XO', sep=""))

# Move to backup environment
setwd(paste0(drive,"/emma_toolbox/Scripting"))

# merge
source("subscripts/merger.R")    
    
# compile database 
source("subscripts/compile_database.R")

# get center
center <- read_file("../emma_toolbox_data/Settings/center.txt")
center <- gsub("[\r\n]", "", center)

# get pop
pop <- read_file("../emma_toolbox_data/Settings/pop.txt")
pop <- gsub("[\r\n]", "", pop)

# repeated measures
print("Insert base number for neurocognitive monitoring: ")
subj <- readLines(file("stdin"),1)
print(paste0("Starting neurocognitive monitoring for:", subj))

# SGH or CM?
print("Neurocognitive monitoring for SGH (press 1) or CM (press 2)? ")
mode <- readLines(file("stdin"),1)


if(mode == "1"){
  print(paste0("Starting neurocognitive monitoring in SGH for: ", subj))  
  
  source("Scripting/subscripts/repeated_measures.R")

}

if(mode == "2"){
  
  print(paste0("Starting neurocognitive monitoring in CM for:", subj))  
  
  source("Scripting/subscripts/repeated_measures_CM.R")
  
}


shell(paste0('robocopy ../emma_toolbox_data ', drive, '/emma_toolbox_data /E /XO', sep=""))
shell(paste0('robocopy ../emma_toolbox/Scripting ', drive, '/emma_toolbox/Scripting /E /XO', sep=""))


# close
print("Backup & Analysis finished")
print("Press ENTER to exit")
drive <- readLines(file("stdin"),1)