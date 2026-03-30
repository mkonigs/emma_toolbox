
### CHANGE ###

vars <- c("mrt_neutral", "sd_neutral", "acc_neutral", "alerting_mrt", "orienting_mrt", "executive_mrt", "executive_acc", "tau", "encoding_sum.x",   "encoding_sum.y", "vooruit_corr.x", "achteruit_corr.x", "vooruit_corr.y", "achteruit_corr.y",  "overall_random_m", "motor_response")

dirs <- c("neg", "neg", "pos", "neg", "neg", "neg", "pos", "neg", "neg", "pos", "pos", "pos", "pos", "pos", "neg", "neg")

des <- c("Processing speed", "Processing variability", "Processing accuracy", "Alerting attention", "Orienting attention", "Interference control speed", "Interference control accuracy", "Attention consistency", "Visual encoding", "Verbal encoding", "Visual short-term memory", "Visual working memory", "Verbal short-term memory", "Verbal working memory", "Dynamic visuomotor integration", "Motor response inhibition")

trrs <- c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)

scales <-c("Mean of RT (ms)", "SD of RT (ms)", "Accuracy (proportion)", "Change in RT (ms)", "Change in RT (ms)", "Change in accuracy (proportion)", "Change in RT (ms)", "Tau (ms)", "Total displacement (distance from target)", "Direct recall (correct)", "Correct", "Correct", "Correct", "Correct", "Mean Deviance (pixels)", "Latency (ms)") 

summary <- data.frame(vars, dirs, des, trrs, scales)


#Prepare norms for controls
CTRL <- read.table(paste0("../emma_toolbox_data/Norms/", pop, ".txt"), sep=" ", header = TRUE)

for (a in 1:ncol(CTRL)){
  CTRL[,a][is.infinite(CTRL[,a])] <- NA 
}

# import data from local database
ET_DATABASE <- read.table("../emma_toolbox_data/Databases/ET_DATABASE.txt", sep=" ", header = TRUE)

# select baseline number
subj_base <- substr(subj, start=1, stop=6)
ET_subj <- ET_DATABASE[grep(subj_base, ET_DATABASE$subj),]

# GENERATE REPORT
if (exists("bu")){
  
  if(bu==FALSE){
    
    rmarkdown::render("Scripting/subscripts/subscripts/report_generator_change.Rmd", output_file = paste0("../../../../emma_toolbox_data/Files/", subj, "/", emmadate, "/Reports/Emma Toolbox - Change.html", sep=""))}
  
  if(bu==TRUE){
    
    rmarkdown::render("Scripting/subscripts/subscripts/report_generator_change.Rmd", output_file = paste0(base, "/../emma_toolbox_data/Files/", subj, "/Emma Toolbox - Change.html"))
  }
  
}else{
  
  rmarkdown::render("Scripting/subscripts/subscripts/report_generator_change.Rmd", output_file = paste0("../../../../emma_toolbox_data/Files/", subj, "/", emmadate, "/Reports/Emma Toolbox - Change.html", sep=""))
  
}
