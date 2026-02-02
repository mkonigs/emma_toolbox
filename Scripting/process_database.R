
source("Scripting/get_emmadate.R")
source("Scripting/names.R")

write.table(eval(parse(text=paste0(db_name, "_data"))), paste0("../emma_toolbox_data/Files/", subj, "/", emmadate, "/Databases/", db_name, "_data_", subj, "_",emmadate,".txt"), row.names = F)
write.table(eval(parse(text=paste0(db_name, "_data"))), paste0("../emma_toolbox_data/Databases/", foldername, "/", db_name, "_data_", subj, "_",emmadate,".txt"), row.names = F)

if(dir.exists(paste0("../emma_toolbox_data/Databases/", foldername))){
  
  setwd (paste0("../emma_toolbox_data/Databases/", foldername))
  
  if (exists(paste0(db_name, "_database.txt"))){
    file.remove(paste0(db_name, "_database.txt"))
  }
  
  file_list <- list.files(pattern=paste0(db_name, "_data"))
  
  for (file in file_list){
    
    # if the merged dataset doesn't exist, create it
    if (!exists("dataset")){
      dataset <- read.table(file, header=TRUE, sep=" ")
    } else{
    
    # if the merged dataset does exist, append to it
      temp_dataset <-read.table(file, header=TRUE, sep=" ")
      dataset<-rbind.fill(dataset, temp_dataset)
      rm(temp_dataset)
    }
    
  }
  
  testname1 <- as.Date(dataset$test)
  testname2 <- mdy(dataset$test)
  testname1[is.na(testname1)] <- testname2[is.na(testname1)]
  dataset$testname <- testname1
  
  write.table(dataset[2:(length(file_list)+1),], paste0(db_name, "_database.txt"), row.names = F)
  rm(dataset)
  
  setwd("../../../emma_toolbox")
}
