#This code is used to create a heatmap to get the TRBJ distrubution of the patient samples under stimulated and un-stimulated conditions
#We upload the excel data and then clean the data for the relevant information. Then we find the expanded clones (defined in terms of TCRusage and CDR3a and CDR3b sequences)
#For the expanded clones, we pool together the data for all the patients at stimulatd and un-stimulated conditions
#We calculate the fraction of different TRBJ and plot them as a heatmap with x-axis giving the different patients at stimulated and un-stimulated conditions and the y-axis giving the fractions of TRAJ expression 

library(readxl)
my_data <- read_xlsx("HIV data_Pooled_Final.xlsx", col_types = "text", col_names = TRUE)#The data is read
new_data1<-my_data[,-c(1,3,4,6,8,11,12,16,17,33,34,35)]#The relevant columns are selected

library(dplyr)
new_data1<-new_data1%>%filter(TRBV != "NA", TRBJ != "NA", TRAV != "NA", TRAJ !="NA", CDR3a !="NA",  CDR3b !="NA", Stimulation != "BSV18")#Get rid of columns with NA values for TCRusage and BSV18 stimulation condition
new_data1<-new_data1%>%filter(TRAV !="TRAV1-1")#Get rid of data with TRAV1-1 because MAIT cells are defined as TRAV1-2
new_data1$TCRusage <-do.call(paste0, new_data1[c("TRAV","TRAJ","TRBV","TRBJ")])#Create a new column called TCRusage which put together TRAV,TRAJ, TRBV, TRBJ

library(dplyr)
test<- count(new_data1, Subject, Stimulation, TRBJ, TCRusage,CDR3a,CDR3b)#Count the number of cells with a given TCR usage and CDR3a and CDR3b sequences


#For expanded clones

final200<-test %>% filter(n != 1)# We are interested in the expanded clones so we get rid of cells that appear once with a given TCR usage and CDR3a and CDR3b sequences
a=grep("HIV", final200$Subject)#In the data set for expanded clones, any expression that has "HIV" in it is chosen
final200$Subject[a]<-"HIV"#Any expression containing "HIV" is renamed as "HIV"


b=grep("EC", final200$Subject)#In the data set for expanded clones, any expression that has "EC" in it is chosen
final200$Subject[b]<-"EC"#Any expression containing "EC" is renamed as "EC"

c=grep("HD", final200$Subject) #In the data set for expanded clones, any expression that has "HD" in it is chosen
final200$Subject[c]<-"HD"#Any expression containing "HD" is renamed as "HD"


final500<-acast(final200, TRBJ~Subject+Stimulation, value.var = "n", fun.aggregate = sum)# Calculate the number of times one kind of TRBJ occurs for different patients at stimulated and un-stimulated conditions


final600<- apply(final500,2,function(x){x/sum(x)})#Calculate the fractions of different TRBJ sequence for different patients at stimulated and un-stimulated conditions
final700<-final600[,c(4,3,2,1,6,5)]#Re-arrange the data 
colnames(final700)<-c("HD no stim","HD stim","EC no stim","EC stim","PR no stim","PR stim")

#Make the heatmap for the processed data
library(pheatmap)
result <-pheatmap(final700, color = colorRampPalette (c("white","navy","firebrick3"))(100), cluster_rows = TRUE, cluster_cols = FALSE,clustering_distance_rows = "euclidean", clustering_method = "complete",gaps_col =c(2,4))




