#This code is written to find the CDR3a distribution as well as the length of CDR3a of the patient samples.
#We upload the excel data and then clean the data for the relevant information. Then we find the expanded clones (defined in terms of TCRusage and CDR3a and CDR3b sequences)
#For the expanded clones we find the distribution of CDR3a sequences as well as the length of the CDR3a sequences

library(readxl)
my_data <- read_xlsx("HIV data_Pooled_Final.xlsx", col_types = "text", col_names = TRUE)#The data is read
new_data1<-my_data[,-c(1,3,4,6,8,11,12,16,17,33,34,35)]#The relevant columns are selected

library(dplyr)
new_data1<-new_data1%>%filter(TRBV != "NA", TRBJ != "NA", TRAV != "NA", TRAJ !="NA", CDR3a !="NA",  CDR3b !="NA", Stimulation != "BSV18")#Get rid of columns with NA values for TCRusage and BSV18 stimulation condition
new_data1<-new_data1%>%filter(TRAV !="TRAV1-1") #Get rid of data with TRAV1-1 because MAIT cells are defined as TRAV1-2
new_data1$TCRusage <-do.call(paste0, new_data1[c("TRAV","TRAJ","TRBV","TRBJ")])#Create a new column called TCRusage which put together TRAV,TRAJ, TRBV, TRBJ

test<- count(new_data1, Subject, Stimulation,TCRusage,CDR3a,CDR3b) #Count the number of cells with a given TCR usage and CDR3a and CDR3b sequences

#For expanded clones
library(reshape2)

final200<-test %>% filter(n != 1)# We are interested in the expanded clones so we get rid of cells that appear once with a given TCR usage and CDR3a and CDR3b sequences


final300<-final200[,c(1,2,4,6)]#The relevant columns are selected

final600<-melt(final300)# This is done so that each row is a unique id-variable combination
final700<-final600[,c(1,2,3,5)]# The relevant columns are selected
final800<-final700[rep(row.names(final700), final700$value), 1:3]#Every CDR3a seqence that happens more then once gets a new row

final800$cdr3alength<-nchar(final800$CDR3a) #Add a new column which calculates the length of CDR3a

test100<- count(final800, Subject, Stimulation,cdr3alength) #This gives how many CDR3a are of a particular length for a given condition




write.csv(final800, "List of CDR3a_Individual donors.csv") #Gives the required output in csv form
write.csv(test100, "CDR3a length for individual donors.csv")#Gives the required output in csv form
