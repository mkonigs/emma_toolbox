
### CHANGE ###

analysis_type <- "Emma Toolbox"


if(dir.exists("Data/Databases/Change")){
  
  setwd ("Tests/Change")
  
  emmadate <- read_file("../../Data/emmadate.txt")
  emmadate <- gsub("\n", "", emmadate)
  emmadate <- gsub("\r", "", emmadate)
  
  rmarkdown::render("ET_Change.Rmd",
                    output_file = paste0("../../Data/", subj, "/", emmadate, "/Reports/Emma Toolbox - Change.html", sep=""))
  
  emmadate <- read_file("../../Data/emmadate.txt")
  emmadate <- gsub("\n", "", emmadate)
  emmadate <- gsub("\r", "", emmadate)
  
  #write.table(et_domains, paste0("../../Data/", subj, "/", emmadate, "/Databases/domains_data_", subj, ".txt", sep=""), row.names = F)
  #write.table(et_domains, paste0("../../Data/Databases/Domains/domains_data_", subj, "_",emmadate, ".txt", sep=""), row.names = F)
  
  setwd ("../../")
}