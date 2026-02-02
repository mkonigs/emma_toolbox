### PROCESS ###

source("Scripting/names.R")

if(file.exists(paste0("../emma_toolbox_data/Raw/", foldername, "/", testname, "-subject-", subj, ".csv"))){
  print(paste0(foldername, " file present"))
  
  source("Scripting/get_emmadate.R")
  
  file.copy(paste0("../emma_toolbox_data/Raw/", foldername, "/", testname, "-subject-", subj, ".csv"), paste0("../emma_toolbox_data/Files/", subj, "/", emmadate, "/Raw"))
  file.copy(paste0("../emma_toolbox_data/Raw/", foldername, "/", testname, "-subject-", subj, ".csv"), paste0("../emma_toolbox_data/Raw/", foldername, "/", testname, "-subject-", subj,"_",emmadate, ".csv"))
  
  rmarkdown::render(paste0("Scripting/subscripts/process_", testname, ".Rmd"),
                    output_file = paste0("../../../emma_toolbox_data/Files/", subj, "/", emmadate, "/Reports/Report - ", foldername, ".html"))
  
}



