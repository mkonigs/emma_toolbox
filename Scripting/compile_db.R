setwd ("Data/Databases")

attention_dataset <-read.table("Attention/ant_database.txt", header=TRUE, sep=" ")

if (center == "KNVB"){
  attention_dataset <- attention_dataset[attention_dataset$subj > 110000,]
}

if (file.exists("Attention/eyetracker/general/eyetracker_database.txt")){
  eyetracker_dataset <- read.table("Attention/eyetracker/general/eyetracker_database.txt", header=TRUE, sep=" ")
  eyetracker_dataset <- dplyr::rename(eyetracker_dataset, subj = subject.id)
}

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

write.table(ET_DATABASE, "ET_DATABASE.txt", row.names = F)

setwd("../../")