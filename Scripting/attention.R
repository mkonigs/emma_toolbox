### ATTENTION ###

if(file.exists(paste("../emma_toolbox_data/Raw/Attention/attention-subject-", subj, ".csv", sep=""))){
  print("Attention file present")
  
  source("Scripting/get_emmadate.R")
  
  file.copy(paste("../emma_toolbox_data/Raw/Attention/attention-subject-", subj, ".csv", sep=""), paste0("../emma_toolbox_data/Files/", subj, "/", emmadate, "/Raw"))
  file.copy(paste("../emma_toolbox_data/Raw/Attention/attention-subject-", subj, ".csv", sep=""), paste("../emma_toolbox_data/Raw/Attention/attention-subject-", subj,"_",emmadate, ".csv", sep=""))
  
  rmarkdown::render("Scripting/subscripts/process_attention.Rmd",
                    output_file = paste0("../../../emma_toolbox_data/Files/", subj, "/", emmadate, "/Reports/Report - Attention.html", sep=""))
  
  source("Scripting/get_emmadate.R")
  
  write.table(ant_data, paste0("../emma_toolbox_data/Files/", subj, "/", emmadate, "/Databases/ant_data_", subj, "_",emmadate,".txt", sep=""), row.names = F)
  write.table(ant_data, paste0("../emma_toolbox_data/Databases/Attention/ant_data_", subj, "_",emmadate,".txt", sep=""), row.names = F)
  
  setwd("../../")
  
}

if(dir.exists("../emma_toolbox_data/Databases/Attention")){
  
  setwd ("../emma_toolbox_data/Databases/Attention")
  
  if (!exists("dataset")){
    file.remove("ant_database.txt")
  }
  
  file_list <- list.files(pattern="ant_data")
  
  for (file in file_list){
    
    # if the merged dataset doesn't exist, create it
    if (!exists("dataset")){
      dataset <- read.table(file, header=TRUE, sep=" ")
    }
    
    
    # if the merged dataset does exist, append to it
    if (exists("dataset")){
      temp_dataset <-read.table(file, header=TRUE, sep=" ")
      dataset<-rbind.fill(dataset, temp_dataset)
      rm(temp_dataset)
    }
    
  }
  
  test1 <- as.Date(dataset$test)
  test2 <- mdy(dataset$test)
  test1[is.na(test1)] <- test2[is.na(test1)]
  dataset$test <- test1
  
  write.table(dataset[2:(length(file_list)+1),], "ant_database.txt", row.names = F)
  rm(dataset)
  
  setwd("../../")
}
