### VERBAL WM ###

analysis_type <- "Emma Toolbox"


if(file.exists(paste("Data/Raw/VerbalWM/verbalwm-subject-", subj, ".csv", sep=""))){
  print("VerbalWM file present")
  
  file.copy(paste("Data/Raw/VerbalWM/verbalwm-subject-", subj, ".csv", sep=""), paste0("Data/", subj, "/", emmadate, "/Raw"))
  file.copy(paste("Data/Raw/VerbalWM/verbalwm-subject-", subj, ".csv", sep=""), paste("Data/Raw/VerbalWM/verbalwm-subject-", subj,"_",emmadate, ".csv", sep=""))
  
  
  setwd ("Tests/VerbalWM")
  
  emmadate <- read_file("../../Data/emmadate.txt")
  emmadate <- gsub("\n", "", emmadate)
  emmadate <- gsub("\r", "", emmadate)
  
  rmarkdown::render("VerbalWM.Rmd",
                    output_file = paste0("../../Data/", subj, "/", emmadate, "/Reports/Report - Verbal Working Memory.html", sep=""))
  
  emmadate <- read_file("../../Data/emmadate.txt")
  emmadate <- gsub("\n", "", emmadate)
  emmadate <- gsub("\r", "", emmadate)
  
  
  write.table(verbalwm_data, paste0("../../Data/", subj, "/", emmadate, "/Databases/verbalwm_data_", subj, "_",emmadate,".txt", sep=""), row.names = F)
  write.table(verbalwm_data, paste0("../../Data/Databases/VerbalWM/verbalwm_data_", subj, "_",emmadate,".txt", sep=""), row.names = F)
  
  setwd ("../../")
}

if(dir.exists("Data/Databases/VerbalWM")){
  
  setwd ("Data/Databases/VerbalWM")
  
  if (!exists("dataset")){
    file.remove("verbalwm_database.txt")
  }
  
  file_list <- list.files(pattern="verbalwm_data_")
  
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
  
  write.table(dataset[2:(length(file_list)+1),], "verbalwm_database.txt", row.names = F)
  rm(dataset)
  
  
  setwd ("../../../")
  
}
