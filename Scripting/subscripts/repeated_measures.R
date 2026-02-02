
# open subject folder in Explorer
#shell(paste0("explorer ..\..\emma_toolbox_data\Files\\", subj), intern=FALSE)

### DOMAINS ###
# get center
center <- read_file("../../emma_toolbox_data/Settings/center.txt")
center <- gsub("[\r\n]", "", center)

# get pop
pop <- read_file("../../emma_toolbox_data/Settings/pop.txt")
pop <- gsub("[\r\n]", "", pop)

# start analysis
source("Scripting/subscripts/domains.R")

print("Domain analysis finished")
