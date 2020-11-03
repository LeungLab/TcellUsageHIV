#This code is written to find the CDR3b distribution as well as the length of CDR3b of the patient samples.
#We upload the excel data and then clean the data for the relevant information. Then we find the expanded clones (defined in terms of TCRusage and CDR3a and CDR3b sequences)
#For the expanded clones, we pool together the data for each different kind of patients
#For the expanded clones for the pooled patient data, we find the distribution of CDR3b sequences as well as the length of the CDR3b sequences

library(readxl)
my_data <- read_xlsx("HIV data_Pooled_Final.xlsx", col_types = "text", col_names = TRUE)#The data is read
new_data1<-my_data[,-c(1,3,4,6,8,11,12,16,17,33,34,35)]#The relevant columns are selected

library(dplyr)
new_data1<-new_data1%>%filter(TRBV != "NA", TRBJ != "NA", TRAV != "NA", TRAJ !="NA", CDR3a !="NA",  CDR3b !="NA", Stimulation != "BSV18")#Get rid of columns with NA values for TCRusage and BSV18 stimulation condition
new_data1<-new_data1%>%filter(TRAV !="TRAV1-1")#Get rid of data with TRAV1-1 because MAIT cells are defined as TRAV1-2
new_data1$TCRusage <-do.call(paste0, new_data1[c("TRAV","TRAJ","TRBV","TRBJ")])#Creating a new column which put together TRAV,TRAJ, TRBV, TRBJ into a new column called TCRusage


test<- count(new_data1, Subject, Stimulation,TCRusage,CDR3a,CDR3b)#Counting the number of cells with a given TCR usage and CDR3a and CDR3b sequences

#For expanded clones
library(reshape2)

final200<-test %>% filter(n != 1)# We are interested in the expanded clones so we get rid of cells that appear once with a given TCR usage and CDR3a and CDR3b sequences
a=grep("HIV", final200$Subject) #In the data set for expanded clones, any expression that has "HIV" in is chosen
final200$Subject[a]<-"HIV"#Any expression containing "HIV" is renamed as "HIV"

b=grep("EC", final200$Subject)#In the data set for expanded clones, any expression that has "EC" in is chosen
final200$Subject[b]<-"EC"#Any expression containing "EC" is renamed as "EC"

c=grep("HD", final200$Subject)#In the data set for expanded clones, any expression that has "HD" in is chosen
final200$Subject[c]<-"HD"#Any expression containing "HD" is renamed as "HD"

final300<-final200[,c(1,2,5,6)]#The relevant columns are selected
final600<-melt(final300)# This is done so that each row is a unique id-variable combination
final700<-final600[,c(1,2,3,5)]# The relevant columns are selected
final800<-final700[rep(row.names(final700), final700$value), 1:3]#Every CDR3b seqence that happens more then once gets a new row


write.csv(final800, "List of CDR3b_Pooled.csv")#Gives the required output in csv form
