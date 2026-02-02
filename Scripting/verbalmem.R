### VERBAL MEMORY ###

analysis_type <- "Emma Toolbox"

if(file.exists(paste("Data/Raw/VerbalMemory/verbalmem-subject-", subj, ".csv", sep=""))){
  print("Verbal Memory file present")
  
  file.copy(paste("Data/Raw/VerbalMemory/verbalmem-subject-", subj, ".csv", sep=""), paste0("Data/", subj, "/", emmadate, "/Raw"))
  file.copy(paste("Data/Raw/VerbalMemory/verbalmem-subject-", subj, ".csv", sep=""), paste("Data/Raw/VerbalMemory/verbalmem-subject-", subj,"_",emmadate, ".csv", sep=""))
  
  
  setwd ("Tests/VerbalMemory")
  
  emmadate <- read_file("../../Data/emmadate.txt")
  emmadate <- gsub("\n", "", emmadate)
  emmadate <- gsub("\r", "", emmadate)
  
  rmarkdown::render("Verbal_Memory.Rmd",
                    output_file = paste0("../../Data/", subj, "/", emmadate, "/Reports/Report - Verbal Memory.html", sep=""))
  
  emmadate <- read_file("../../Data/emmadate.txt")
  emmadate <- gsub("\n", "", emmadate)
  emmadate <- gsub("\r", "", emmadate)
  
  write.table(verbalmt_data, paste0("../../Data/", subj, "/", emmadate, "/Databases/verbalmt_data_", subj, "_",emmadate,".txt", sep=""), row.names = F)
  write.table(verbalmt_data, paste0("../../Data/Databases/VerbalMemory/verbalmt_data_", subj, "_",emmadate,".txt", sep=""), row.names = F)
  
  setwd ("../../")
  
}

if (dir.exists("Data/Databases/VerbalMemory")){
  
  setwd ("Data/Databases/VerbalMemory")
  
  if (!exists("dataset")){
    file.remove("verbalmt_database.txt")
  }
  
  file_list <- list.files(pattern="verbalmt_data_")
  
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
  
  write.table(dataset[2:(length(file_list)+1),], "verbalmt_database.txt", row.names = F)
  rm(dataset)
  
  
  setwd ("../../../")
}
