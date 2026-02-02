# load packages

source("ipak.R")
package_list <- c("rmarkdown", "knitr", "lubridate", "readr", "plyr")
ipak(package_list)

# get center
center <- read_file("../../emma_toolbox_data/Settings/center.txt")
center <- gsub("[\r\n]", "", center)

# go to folder 
setwd ("../../emma_toolbox_data/Databases")

attention_dataset <-read.table("Attention/ant_database.txt", header=TRUE, sep=" ")

if (center == "KNVB"){
  attention_dataset_cm <- attention_dataset
  attention_dataset_cm <- attention_dataset_cm[attention_dataset_cm$subj < 110000,]
  
  # separate CM observations
  attention_dataset <- attention_dataset[attention_dataset$subj > 110000,]
}

# remove data with subject = 0
attention_dataset <- attention_dataset[attention_dataset$subj > 0,]

if (file.exists("VisualMemory/visualmt_database.txt")){
  visualmt_dataset <-read.table("VisualMemory/visualmt_database.txt", header=TRUE, sep=" ")
  ET_DATABASE <- merge(unique(attention_dataset), unique(visualmt_dataset), by.x=c("subj", "test"), by.y=c("subj", "test"), all.x = T)
}

if (file.exists("VerbalMemory/verbalmt_database.txt")){
  verbalmt_dataset <-read.table("VerbalMemory/verbalmt_database.txt", header=TRUE, sep=" ")
  ET_DATABASE <- merge(unique(ET_DATABASE), unique(verbalmt_dataset), by.x=c("subj", "test"), by.y=c("subj", "test"), all.x = T)
}

if(file.exists("VisualWM/visualwm_database.txt")){
  visualWM_dataset <-read.table("VisualWM/visualwm_database.txt", header=TRUE, sep=" ")
  ET_DATABASE <- merge(unique(ET_DATABASE), unique(visualWM_dataset), by.x=c("subj", "test"), by.y=c("subj", "test"), all.x = T)
}

if(file.exists("VerbalWM/verbalwm_database.txt")){
  verbalWM_dataset <-read.table("VerbalWM/verbalwm_database.txt", header=TRUE, sep=" ")
  ET_DATABASE <- merge(unique(ET_DATABASE), unique(verbalWM_dataset), by.x=c("subj", "test"), by.y=c("subj", "test"), all.x = T)
}

if(file.exists("VisuoMotor/visuomotor_database.txt")){
  visuomotor_dataset <-read.table("VisuoMotor/visuomotor_database.txt", header=TRUE, sep=" ")
  ET_DATABASE <- merge(unique(ET_DATABASE), unique(visuomotor_dataset), by.x=c("subj", "test"), by.y=c("subj", "test"),all.x = T)
}

# write results
write.table(ET_DATABASE, "ET_DATABASE.txt", row.names = F)

if (center == "KNVB"){
  write.table(attention_dataset_cm, "ET_DATABASE_CM.txt", row.names = F)
}

dups <- ET_DATABASE[ET_DATABASE$subj %in% ET_DATABASE$subj[duplicated(ET_DATABASE$subj)],]

print("ET DATABASE compiled")
setwd("../../emma_toolbox/Scripting")

