### VISUOMOTOR ###

analysis_type <- "Emma Toolbox"

if(file.exists(paste("../emma_toolbox_data/Raw/VisuoMotor/visuomotor-subject-", subj, ".csv", sep=""))){
  print("Visuomotor file present")
  
  file.copy(paste("../emma_toolbox_data/Raw/VisuoMotor/visuomotor-subject-", subj, ".csv", sep=""), paste0("Data/", subj, "/", emmadate, "/Raw"))
  file.copy(paste("../emma_toolbox_data/Raw/VisuoMotor/visuomotor-subject-", subj, ".csv", sep=""), paste("../emma_toolbox_data/Raw/VisuoMotor/visuomotor-subject-", subj,"_",emmadate, ".csv", sep=""))
  
  
  setwd ("Tests/VisuoMotor")
  
  emmadate <- read_file("../../../emma_toolbox_data/emmadate.txt")
  emmadate <- gsub("\n", "", emmadate)
  emmadate <- gsub("\r", "", emmadate)
  
  rmarkdown::render("Visuomotor.Rmd",
                    output_file = paste0("../../../emma_toolbox_data/", subj, "/", emmadate, "/Reports/Report - Visuomotor.html", sep=""))
  
  emmadate <- read_file("../../../emma_toolbox_data/emmadate.txt")
  emmadate <- gsub("\n", "", emmadate)
  emmadate <- gsub("\r", "", emmadate)
  
  write.table(visuomotor_data, paste0("../../../emma_toolbox_data/", subj, "/", emmadate, "/Databases/visuomotor_data_", subj, "_",emmadate,".txt", sep=""), row.names = F)
  write.table(visuomotor_data, paste0("../../../emma_toolbox_data/Databases/VisuoMotor/visuomotor_data_", subj, "_",emmadate,".txt", sep=""), row.names = F)
  
  setwd ("../../")
}

if(dir.exists("Data/Databases/Visuomotor")){
  
  setwd ("Data/Databases/VisuoMotor")
  
  if (!exists("dataset")){
    file.remove("visuomotor_database.txt")
  }
  
  file_list <- list.files(pattern="visuomotor_data_")
  
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
  
  write.table(dataset[2:(length(file_list)+1),], "visuomotor_database.txt", row.names = F)
  rm(dataset)
  
  setwd("../../../")
  
}