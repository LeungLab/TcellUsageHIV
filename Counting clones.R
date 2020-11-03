#This code gives the number of total clones, the number of unique clones, the number of clones that are present 2 times and the number of clones that are present more than 2 times
#We upload the excel data and then clean the data for the relevant information. Then we find the expanded clones (defined in terms of TCRusage and CDR3a and CDR3b sequences)
#For the expanded clones, we pool together the data for each different kind of patients
#Then we find the number of total clones, the number of unique clones, the number of clones that are present 2 times and the number of clones that are present more than 2 times 

library(readxl)
my_data <- read_xlsx("HIV data_Pooled_Final.xlsx", col_types = "text", col_names = TRUE)#The data is read
new_data1<-my_data[,-c(1,3,4,6,8,11,12,16,17,33,34,35)]#The relevant columns are selected

library(dplyr)
new_data1<-new_data1%>%filter(TRBV != "NA", TRBJ != "NA", TRAV != "NA", TRAJ !="NA", CDR3a !="NA",  CDR3b !="NA", Stimulation != "BSV18")#Get rid of columns with NA values for TCRusage and BSV18 stimulation condition
new_data1<-new_data1%>%filter(TRAV !="TRAV1-1")#Get rid of data with TRAV1-1 because MAIT cells are defined as TRAV1-2
new_data1$TCRusage <-do.call(paste0, new_data1[c("TRAV","TRAJ","TRBV","TRBJ")])#Create a new column called TCRusage which put together TRAV,TRAJ, TRBV, TRBJ

library(dplyr)
test<- count(new_data1, Subject, Stimulation,TCRusage,CDR3a,CDR3b)#Count the number of cells with a given TCR usage and CDR3a and CDR3b sequences

#For expanded clones

final200<-test %>% filter(n ==2)#Gives the clones that are present 2 times
final300<-test %>% filter(n==1) #Gives the unique clones
final400<-test%>%filter(n>2)# Gives the clones that are present more than 2 times
final500<-test%>%filter(n !=1) #Gives all the expanded clones

#Find any expression that contains "HIV" and replace them with "HIV"
a=grep("HIV", final200$Subject)
a1=grep("HIV", final300$Subject)
a2=grep("HIV", final400$Subject)
a3=grep("HIV", final500$Subject)

final200$Subject[a]<-"HIV"
final300$Subject[a1]<-"HIV"
final400$Subject[a2]<-"HIV"
final500$Subject[a3]<-"HIV"

#Find any expression that contains "EC" and replace them with "EC"

b=grep("EC", final200$Subject)
b1=grep("EC", final300$Subject)
b2=grep("EC", final400$Subject)
b3=grep("EC", final500$Subject)

final200$Subject[b]<-"EC"
final300$Subject[b1]<-"EC"
final400$Subject[b2]<-"EC"
final500$Subject[b3]<-"EC"

#Find any expression that contains "HD" and replace them with "HD"
c=grep("HD", final200$Subject)
c1=grep("HD", final300$Subject)
c2=grep("HD", final400$Subject)
c3=grep("HD", final500$Subject)

final200$Subject[c]<-"HD"
final300$Subject[c1]<-"HD"
final400$Subject[c2]<-"HD"
final500$Subject[c3]<-"HD"


final3000<-split(final200, list(final200$Subject, final200$Stimulation))#Splits the data set that give clones that occur 2 times in terms of subject and stimulation condition
Numerofclones_nequalto2<-sapply(final3000,function(x) sum(x$n))# This sums up the number of clones that occur 2 times in terms of subject and stimulation condition

final4000<-split(final300, list(final300$Subject, final300$Stimulation))#Splits the data set that give unique clones in terms of subject and stimulation condition
Number_uniqueclones<-sapply(final4000, function(x) sum(x$n))#This sums up the unique clones in terms of subject and stimulation condition

final5000<-split(final400,list(final400$Subject, final400$Stimulation))#Splits the data set that give clones that occur more than 2 times in terms of subject and stimulation condition
Numberofclones_ngreaterthan2<-sapply(final5000, function(x) sum(x$n))# This sums up the number of clones that occur more than 2 times in terms of subject and stimulation condition


result1000<-rbind(Numerofclones_nequalto2,Number_uniqueclones,Numberofclones_ngreaterthan2)#This put together the unique clones, the number of clones that occur 2 times and the number of clones that occur more than 2 times

final2018<-split(test,list(test$Subject, test$Stimulation))#Splits the data for expanded clones in terms of subject and stimulation condition
Totalclones<-sapply(final2018, function(x) sum(x$n))#This sums up the total number of expanded clones

#Get the desired outputs in csv form

write.csv(result1000, "Summary of unique and expanded clones.csv")
write.csv(Totalclones, "Total number of clones for all subjects.csv")
