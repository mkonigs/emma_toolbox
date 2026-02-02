options(repos = c(CRAN = "http://cran.rstudio.com/"))

# load packages
source("Scripting/ipak.R")
package_list <- c("mice", "reshape2", "readr", "psych", "ggplot2", "lubridate", "kableExtra", "dplyr")
ipak(package_list)


all_data <- read_csv(paste("../emma_toolbox_data/Raw/Attention/attention-subject-", subj,".csv",sep=""), col_types = cols())

Test <- c()
Details <- c()

if (exists("datetime", all_data) == TRUE){
  time <- unique(all_data$datetime)
  
  emmadate <- c(strsplit(time," "))
  emmadate <- unlist(emmadate)
  emmadate <- emmadate[1]
  emmadate <- gsub("/", "", emmadate)
  
  Test <- c(Test, "Date and time")
  Details <- c(Details, time[1])
}

date <- c(strsplit(time," "))
date <- unlist(date)
date <- date[1]
test <- mdy(date)

# get patient data
et_exp <- read.table("../emma_toolbox_data/Databases/ET_DATABASE.txt", header = TRUE, sep=" ")

# format subject numbers so that repeated measures have identical subject number, depending on in backup env or not
if (exists("bu")){
  et_exp$subj <- as.numeric(substr(et_exp$subj, 1, 6))
  subj_base <- as.numeric(substr(subj, 1, 6))
} else{
  subj_base <- subj
}

# select the (repeated) measurement(s) associated with the subject
et_exp <- et_exp[which(et_exp$subj %in% subj_base),]
et_exp$group <- 1

# get control data
et_ctrl <- read.table(paste0("../emma_toolbox_data/Norms/", pop, ".txt"), header = TRUE)
et_ctrl$group <- 0

# combine exp with ctrl data
et_ctrl <- et_ctrl[,which(colnames(et_ctrl) %in% colnames(et_exp))]
et_exp <- et_exp[,which(colnames(et_exp) %in% colnames(et_ctrl))]

et_raw <- rbind(et_exp, et_ctrl)

#PRESELECT
et_raw <- unique(et_raw)
et_raw$subj <- paste0(et_raw$subj, "_", et_raw$test)

# PREPROCESS
rmarkdown::render("Scripting/subscripts/subscripts/preprocess.Rmd")

# IMPUTE
rmarkdown::render("Scripting/subscripts/subscripts/impute.Rmd")

# APPLY PCA
rmarkdown::render("Scripting/subscripts/subscripts//apply_pca.Rmd")

# APPLY NORM REGRESSION
ctrl_demo <- read.table(paste0("../emma_toolbox_data/Norms/", pop, "_age.txt"), header = TRUE)
et_domains_ <- et_domains
rmarkdown::render("Scripting/subscripts/subscripts/norm_scoring.Rmd")

# GENERATE REPORT
if (exists("bu")){
  
  if(bu==FALSE){
    
    rmarkdown::render("Scripting/subscripts/subscripts/report_generator.Rmd", output_file = paste0("../../../../emma_toolbox_data/Files/", subj, "/", emmadate, "/Reports/Emma Toolbox - Report.html", sep=""))}
  
  if(bu==TRUE){
    
    rmarkdown::render("Scripting/subscripts/subscripts/report_generator.Rmd", output_file = paste0(base, "../emma_toolbox_data/Files/", subj, "/Emma Toolbox - Report.html"))
    }
  
}else{
  
  rmarkdown::render("Scripting/subscripts/subscripts/report_generator.Rmd", output_file = paste0("../../../../emma_toolbox_data/Files/", subj, "/", emmadate, "/Reports/Emma Toolbox - Report.html", sep=""))
  
}
