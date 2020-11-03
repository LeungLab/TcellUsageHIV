#This code is used to create a heatmap to get the TRBV distrubution of the patient samples under stimulated and un-stimulated conditions
#We upload the excel data and then clean the data for the relevant information. Then we find the expanded clones (defined in terms of TCRusage and CDR3a and CDR3b sequences)
#For the expanded clones, we calculate the fraction of different TRBV and plot them as a heatmap with x-axis giving the different patients at stimulated and un-stimulated conditions and the y-axis giving the fractions of TRAJ expression 

library(readxl)
my_data <- read_xlsx("HIV data_Pooled_Final.xlsx", col_types = "text", col_names = TRUE)#The data is read
new_data1<-my_data[,-c(1,3,4,6,8,11,12,16,17,33,34,35)]#The relevant columns are selected

library(dplyr)
new_data1<-new_data1%>%filter(TRBV != "NA", TRBJ != "NA", TRAV != "NA", TRAJ !="NA", CDR3a !="NA",  CDR3b !="NA", Stimulation != "BSV18")#Get rid of columns with NA values for TCRusage and BSV18 stimulation condition
new_data1<-new_data1%>%filter(TRAV !="TRAV1-1")#Get rid of data with TRAV1-1 because MAIT cells are defined as TRAV1-2
new_data1$TCRusage <-do.call(paste0, new_data1[c("TRAV","TRAJ","TRBV","TRBJ")])#Create a new column called TCRusage which put together TRAV,TRAJ, TRBV, TRBJ

library(dplyr)
test<- count(new_data1, Subject, Stimulation, TRBV, TCRusage,CDR3a,CDR3b)#Count the number of cells with a given TCR usage and CDR3a and CDR3b sequences



#For expanded clones
library(reshape2)
final200<-test %>% filter(n != 1)# We are interested in the expanded clones so we get rid of cells that appear once with a given TCR usage and CDR3a and CDR3b sequences

final500<-acast(final200, TRBV~Subject+Stimulation, value.var = "n", fun.aggregate = sum)# Calculate the number of times one kind of TRBV occurs for different patients at stimulated and un-stimulated conditions

final600<- apply(final500,2,function(x){x/sum(x)})#Calculate the fractions of different TRBV sequence for different patients at stimulated and un-stimulated conditions

final700<-final600[,c(8,7,10,9,11,13,12,2,1,4,3,6,5,14,16,15,18,17,20,19)]#Re-arrange the data 
colnames(final700)<-c("HD1 no stim","HD1 stim","HD2 no stim","HD2 stim","HD3 no stim","HD4 no stim","HD4 stim","EC1 no stim","EC1 stim","EC2 no stim","EC2 stim","EC3 no stim","EC3 stim","PR1 no stim","PR2 no stim","PR2 stim","PR3 no stim","PR3 stim","PR4 no stim","PR4 stim")


#Make the heatmap for the processed data
library(pheatmap)
result <-pheatmap(final700, color = colorRampPalette (c("white","navy","firebrick3"))(100), cluster_rows = TRUE, cluster_cols = FALSE,clustering_distance_rows = "euclidean", clustering_method = "complete",gaps_col =c(2,4,5,7,9,11,13,14,16,18,20))




