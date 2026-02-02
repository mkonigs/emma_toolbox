
# load packages
package_list <- c("lubridate", "plyr", "readr", "knitr", "rmarkdown")
ipak(package_list)

setwd ("../../emma_toolbox_data/Databases/Attention")

if (!exists("dataset")){
  file.remove("ant_database.txt")
}

file_list <- list.files(pattern="ant_data")
file_list <- file_list[sapply(file_list, file.size) > 0]

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

if (dir.exists("../VisualMemory")){
  setwd ("../VisualMemory")
  
  if (!exists("dataset")){
    file.remove("visualmt_database.txt")
  }
  
  file_list <- list.files(pattern="visualmt_data_")
  file_list <- file_list[sapply(file_list, file.size) > 0]
  
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
  
  write.table(dataset[2:(length(file_list)+1),], "visualmt_database.txt", row.names = F)
  rm(dataset)
}

if (dir.exists("../VerbalMemory")){
  setwd ("../VerbalMemory")
  
  if (!exists("dataset")){
    file.remove("verbalmt_database.txt")
  }
  
  file_list <- list.files(pattern="verbalmt_data_")
  file_list <- file_list[sapply(file_list, file.size) > 0]
  
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
}

if (dir.exists("../VisualWM")){
  setwd ("../VisualWM")
  
  if (!exists("dataset")){
    file.remove("visualwm_database.txt")
  }
  
  file_list <- list.files(pattern="visualwm_data_")
  file_list <- file_list[sapply(file_list, file.size) > 0]
  
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
  
  write.table(dataset[2:(length(file_list)+1),], "visualwm_database.txt", row.names = F)
  rm(dataset)
}

if (dir.exists("../verbalWM")){
  
  setwd ("../VerbalWM")
  
  if (!exists("dataset")){
    file.remove("verbalwm_database.txt")
  }
  
  file_list <- list.files(pattern="verbalwm_data_")
  file_list <- file_list[sapply(file_list, file.size) > 0]
  
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
}

if (dir.exists("../VisuoMotor")){
  setwd ("../VisuoMotor")
  
  if (!exists("dataset")){
    file.remove("visuomotor_database.txt")
  }
  
  file_list <- list.files(pattern="visuomotor_data_")
  file_list <- file_list[sapply(file_list, file.size) > 0]
  
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
}

setwd ("../../../emma_toolbox")
