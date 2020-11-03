#This code gives the number of shared clones between the patient groups under stimulated conditions
#We upload the excel data and then clean the data for the relevant information. Then we find the expanded clones (defined in terms of TCRusage and CDR3a and CDR3b sequences)
#For the expanded clones, we pool together the data for each different kind of patients
#For the expanded clones of the pooled patient data, we find the number of shared clones between different patient group under stimulated conditions

library(readxl)
my_data <- read_xlsx("HIV data_Pooled_Final.xlsx", col_types = "text", col_names = TRUE)#The data is read
new_data1<-my_data[,-c(1,3,4,6,8,11,12,16,17,33,34,35)]#The relevant columns are selected

library(dplyr)
new_data1<-new_data1%>%filter(TRBV != "NA", TRBJ != "NA", TRAV != "NA", TRAJ !="NA", CDR3a !="NA",  CDR3b !="NA", Stimulation != "BSV18")#Get rid of columns with NA values for TCRusage and BSV18 stimulation condition
new_data1<-new_data1%>%filter(TRAV !="TRAV1-1")#Get rid of data with TRAV1-1 because MAIT cells are defined as TRAV1-2
new_data1$TCRusage <-do.call(paste0, new_data1[c("TRAV","TRAJ","TRBV","TRBJ")])#Create a new column called TCRusage which put together TRAV,TRAJ, TRBV, TRBJ

library(dplyr)
test<- count(new_data1, Subject, Stimulation,TCRusage,CDR3a,CDR3b)#Count the number of cells with a given TCR usage and CDR3a and CDR3b sequences

final3<-test %>% filter(n != 1)# We are interested in the expanded clones so we get rid of cells that appear once with a given TCR usage and CDR3a and CDR3b sequences

a=grep("HIV", final3$Subject)#In the data set for expanded clones, any expression that has "HIV" in it is chosen
final3$Subject[a]<-"HIV"#Any expression containing "HIV" is renamed as "HIV"

b=grep("EC", final3$Subject)#In the data set for expanded clones, any expression that has "EC" in it is chosen
final3$Subject[b]<-"EC"#Any expression containing "EC" is renamed as "EC"

c=grep("HD", final3$Subject)#In the data set for expanded clones, any expression that has "HD" in it is chosen
final3$Subject[c]<-"HD"#Any expression containing "HD" is renamed as "HD"



library(reshape2)

Elite_Healthy_stimulated<- final3 %>% filter( Subject !="HIV" & Stimulation =="stim")#Filter out the expanded clones for elite and healthy patients at stimulated condition
Elite_Healthy_stimulated1<-acast(Elite_Healthy_stimulated, TCRusage+CDR3a+CDR3b~Subject+Stimulation, value.var = "n", fun.aggregate = sum)#Calculate the number of clones with a ceratin TCR usage and CDR3a and CDR3b sequences
Elite_Healthy_stimulated2<-Elite_Healthy_stimulated1[!(apply(Elite_Healthy_stimulated1,1, function(y) any(y == 0))),]#Select only those clones that are shared between healthy and elite patients


Elite_HIV_stimulated <- final3 %>%filter(Subject !="HD" & Stimulation=="stim")#Filter out the expanded clones for elite and HIV patients at stimulated condition
Elite_HIV_stimulated1<-acast(Elite_HIV_stimulated, TCRusage+CDR3a+CDR3b~Subject+Stimulation, value.var = "n", fun.aggregate = sum)#Calculate the number of clones with a ceratin TCR usage and CDR3a and CDR3b sequences
Elite_HIV_stimulated2<-Elite_HIV_stimulated1[!(apply(Elite_HIV_stimulated1,1, function(y) any(y == 0))),]#Select only those clones that are shared between elite and HIV patients



HIV_Healthy_stimulated<-final3 %>% filter(Subject!="EC" & Stimulation=="stim")#Filter out the expanded clones for elite and healthy patients at un-stimulated condition
HIV_Healthy_stimulated1<-acast(HIV_Healthy_stimulated, TCRusage+CDR3a+CDR3b~Subject+Stimulation, value.var = "n", fun.aggregate = sum)#Calculate the number of clones with a ceratin TCR usage and CDR3a and CDR3b sequences
HIV_Healthy_stimulated2<-HIV_Healthy_stimulated1[!(apply(HIV_Healthy_stimulated1,1, function(y) any(y == 0))),]#Select only those clones that are shared between healthy and HIV patients


stimulated <- final3 %>% filter(Stimulation =="stim")#Filter out the expanded clones for stimulated consition for all the patients
stimulated1<-acast(stimulated, TCRusage+CDR3a+CDR3b~Subject+Stimulation, value.var = "n", fun.aggregate = sum)#Calculate the number of clones with a ceratin TCR usage and CDR3a and CDR3b sequences
stimulated2<-stimulated1[!(apply(stimulated1,1, function(y) any(y == 0))),]#Select the clones that are shared among the three patients


#Get the desired output in csv form
write.csv(Elite_Healthy_stimulated2, "Number_EC_HD_Shared clones_stimulated.csv")
write.csv(Elite_HIV_stimulated2, "Number_EC_HIV_shared clones_stimulated.csv")
write.csv(HIV_Healthy_stimulated2, "Number_HIV_HD_Shared clones_stimulated.csv")
write.csv(stimulated2, "Number_Shared clones among three donors_stimulated.csv")





