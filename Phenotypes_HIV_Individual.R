#This code gives the phenotypic gene expressions for different patients under stimulated and un-stimulated conditions
#We upload the excel data and then clean the data for the relevant information. Then we find the expanded clones (defined in terms of TCRusage and CDR3a and CDR3b sequences)
#We denote gene expression as "1" and no expression as "0"
#The rows which have "NA" value for gene expressions are replaced with "0"
#The expression of each phenotypic gene is then found for the expanded clones

library(readxl)
my_data <- read_xlsx("HIV data_Pooled_Final.xlsx", col_types = "text", col_names = TRUE)#The data is read
new_data1<-my_data[,-c(1,3,4,6,8,11,12,16,17,33,34,35)]#The relevant columns are selected

library(dplyr)
new_data1<-new_data1%>%filter(TRBV != "NA", TRBJ != "NA", TRAV != "NA", TRAJ !="NA", CDR3a !="NA",  CDR3b !="NA", Stimulation != "BSV18")#Get rid of columns with NA values for TCRusage and BSV18 stimulation condition
new_data1<-new_data1%>%filter(TRAV !="TRAV1-1")#Get rid of data with TRAV1-1 because MAIT cells are defined as TRAV1-2
new_data1$TCRusage <-do.call(paste0, new_data1[c("TRAV","TRAJ","TRBV","TRBJ")])#Create a new column called TCRusage which put together TRAV,TRAJ, TRBV, TRBJ

#Denoting gene expression as "1" and no expression as "0"

a=grep("GATA3", new_data1$GATA3)
new_data1$GATA3[a]<-1
b=grep("FOXP3", new_data1$FOXP3)
new_data1$FOXP3[b]<-1
c=grep("GZMB", new_data1$GZMB)
new_data1$GZMB[c]<-1
d=grep("IFNG", new_data1$IFNG)
new_data1$IFNG[d]<-1
e=grep("IL10", new_data1$IL10)
new_data1$IL10[e]<-1
f=grep("IL12A", new_data1$IL12A)
new_data1$IL12A[f]<-1
g=grep("IL17A", new_data1$IL17A)
new_data1$IL17A[g]<-1
h=grep("PRF1", new_data1$PRF1)
new_data1$PRF1[h]<-1
i=grep("RORC", new_data1$RORC)
new_data1$RORC[i]<-1
j=grep("TBET", new_data1$TBET)
new_data1$TBET[j]<-1
k=grep("TGFB1 ", new_data1$TGFB1)
new_data1$TGFB1[k]<-1
l=grep("TNF", new_data1$TNF)
new_data1$TNF[l]<-1
m=grep("BCL6 ", new_data1$BCL6)
new_data1$BCL6 [m]<-1
n=grep("RUNX1", new_data1$RUNX1)
new_data1$RUNX1[n]<-1
o=grep("RUNX3", new_data1$RUNX3)
new_data1$RUNX3[o]<-1

#When the phenotypic gene expression is "NA", it is replaced with "0"
new_data1$GATA3[is.na(new_data1$GATA3)]<-0
new_data1$FOXP3[is.na(new_data1$FOXP3)]<-0
new_data1$GZMB[is.na(new_data1$GZMB)]<-0
new_data1$IFNG[is.na(new_data1$IFNG)]<-0
new_data1$IL10[is.na(new_data1$IL10)]<-0
new_data1$IL12A[is.na(new_data1$IL12A)]<-0
new_data1$IL17A[is.na(new_data1$IL17A)]<-0
new_data1$PRF1[is.na(new_data1$PRF1)]<-0
new_data1$RORC[is.na(new_data1$RORC)]<-0
new_data1$TBET[is.na(new_data1$TBET)]<-0
new_data1$TGFB1[is.na(new_data1$TGFB1)]<-0
new_data1$TNF[is.na(new_data1$TNF)]<-0
new_data1$BCL6[is.na(new_data1$BCL6)]<-0
new_data1$RUNX1[is.na(new_data1$RUNX1)]<-0
new_data1$RUNX3[is.na(new_data1$RUNX3)]<-0

new_data2<-new_data1[,-c(3,4,6,7)]


#For expanded clones

new_data100<-new_data2[,c(1,2,3,4,20)] #The relevant columns are selected


test100<- count(new_data100,Subject, Stimulation,TCRusage,CDR3a, CDR3b) #Count the number of cells with a given TCR usage and CDR3a and CDR3b sequences

final1000<-test100%>% filter(n != 1)# We are interested in the expanded clones, so we get rid of cells that appear once with a given TCR usage and CDR3a and CDR3b sequences
final2000<-final1000[,-6]
final2000<-final2000[,c(1,2,5,4,3)]


result1000<-merge(new_data2,final2000, by=c("Subject","Stimulation","CDR3a","CDR3b","TCRusage"))#The phenotypic gene expression for all the expanded clones are found

#For the expanded clones, the clones giving expression for a particular gene is found
new_data5<-result1000 %>% filter(GATA3==1)
new_data6<-result1000 %>% filter(FOXP3==1) 
new_data7<-result1000 %>% filter( GZMB==1)
new_data8<-result1000 %>% filter(IFNG==1)
new_data9<-result1000 %>% filter( IL10==1)
new_data10<-result1000 %>% filter(IL12A==1)
new_data11<-result1000 %>% filter(IL17A==1)
new_data12<-result1000 %>% filter( PRF1==1)
new_data13<-result1000%>% filter(RORC==1)
new_data14<-result1000 %>% filter(TBET==1)
new_data15<-result1000 %>% filter(TGFB1==1)
new_data16<-result1000 %>% filter(TNF==1)
new_data17<-result1000 %>% filter(BCL6==1)
new_data18<-result1000 %>% filter(RUNX1==1)
new_data19<-result1000 %>% filter(RUNX3==1)


#The desired outputs are written is csv files

write.csv(new_data5, " GATA3_Expanded cells_ID.csv")
write.csv(new_data6," FOXP3_Expanded cells_ID.csv")
write.csv(new_data7," GZMB_Expanded Cells_ID.csv")
write.csv(new_data8," IFNG_Expanded cells_ID.csv")
write.csv(new_data9, " IL10_Expanded cells_ID.csv")
write.csv(new_data10," IL12A_Expanded cells_ID.csv")
write.csv(new_data11, " IL17A_Expanded cells_ID.csv")
write.csv(new_data12, " PRF1_Expanded cells_ID.csv")
write.csv(new_data13, " RORC_Expanded cells_ID.csv")
write.csv(new_data14," TBET_Expanded cells_ID.csv")
write.csv(new_data15, "TGFB1_Expanded cells_ID.csv")
write.csv(new_data16, " TNF_Expanded cells_ID.csv")
write.csv(new_data17, " BCL6_Expanded cells_ID.csv")
write.csv(new_data18, " RUNX1_Expanded cells_ID.csv")
write.csv(new_data19, " RUNX3_Expanded cells_ID.csv")








