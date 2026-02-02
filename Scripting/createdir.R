# convert character into integer
print(paste("Starting Emma Toolbox Analysis for:", subj))
print(getwd())

if(file.exists(paste("../emma_toolbox_data/Raw/Attention/attention-subject-", subj, ".csv", sep=""))){
  setwd("Tests/Attention")
  all_data <- read_csv(paste("../../../emma_toolbox_data/Raw/Attention/attention-subject-",subj,".csv",sep=""), col_types = cols())
  time <- unique(all_data$datetime)
  emmadate <- c(strsplit(time," "))
  emmadate <- unlist(emmadate)
  emmadate <- emmadate[1]
  emmadate <- gsub("/", "", emmadate)
  
  setwd("../")
}

dir.create(file.path(paste0("../../emma_toolbox_data/Files/", subj)))
dir.create(file.path(paste0("../../emma_toolbox_data/Files/", subj, "/", emmadate)))
dir.create(file.path(paste0("../../emma_toolbox_data/Files/", subj, "/", emmadate, "/Reports")))
dir.create(file.path(paste0("../../emma_toolbox_data/Files/", subj, "/", emmadate, "/Databases")))
dir.create(file.path(paste0("../../emma_toolbox_data/Files/", subj, "/", emmadate, "/Raw")))
write(emmadate, "../../emma_toolbox_data/emmadate.txt")

setwd("../")