#This code is used to find out the phenotypic gene expressions for the top three dominant clones of the different patient types
#We upload the excel data and then clean the data for the relevant information. Then we find the expanded clones (defined in terms of TCRusage and CDR3a and CDR3b sequences)
#For the expanded clones, we pool together the data for each different kind of patients
#For the expanded clones of the pooled patient data, we find the top three dominant clones for each patient type
#We then find the phenotypic gene expressions for the dominant clones


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

#Any gene expression data with "NA" is converted to 0
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

new_data2<-new_data1[,-c(3,4,6,7)]#The relevant columns are selected

a=grep("HIV",new_data2$Subject)#In the data set, any expression that has "HIV" in it is chosen
new_data2$Subject[a]<-"HIV"#Any expression containing "HIV" is renamed as "HIV"

b=grep("EC", new_data2$Subject)#In the data, any expression that has "EC" in it is chosen
new_data2$Subject[b]<-"EC"#Any expression containing "EC" is renamed as "EC"


c=grep("HD", new_data2$Subject)#In the data set, any expression that has "HD" in it is chosen
new_data2$Subject[c]<-"HD"#Any expression containing "HD" is renamed as "HD"


#For expanded clones

new_data100<-new_data2[,c(1,2,3,4,20)]#The relevant columns are selected


test100<- count(new_data100,Subject, Stimulation,TCRusage,CDR3a, CDR3b) #Counting the number of cells with a given TCR usage and CDR3a and CDR3b sequences

final1000<-test100%>% filter(n != 1)# We are interested in the expanded clones so we get rid of cells that appear once with a given TCR usage and CDR3a and CDR3b sequences
final2000<-final1000[,-6]#The relevant columns are selected
final2000<-final2000[,c(1,2,5,4,3)]#The relevant columns are selected


a=grep("HIV",final2000$Subject)#In the data set for expanded clones, any expression that has "HIV" in it is chosen
final2000$Subject[a]<-"HIV"#Any expression containing "HIV" is renamed as "HIV"

b=grep("EC", final2000$Subject)#In the data set for expanded clones, any expression that has "EC" in it is chosen
final2000$Subject[b]<-"EC"#Any expression containing "EC" is renamed as "EC"

c=grep("HD", final2000$Subject)#In the data set for expanded clones, any expression that has "HD" in it is chosen
final2000$Subject[c]<-"HD"#Any expression containing "HD" is renamed as "HD"




result1000<-merge(new_data2,final2000, by=c("Subject","Stimulation","CDR3a","CDR3b","TCRusage"))#This gives the phenotypic gene expressions for the expanded clones


result3000<-final1000 %>%group_by(Subject,Stimulation) %>% top_n(3,n) #This gives the top 3 dominant clones out of all the expanded clones

result4000<-result3000[,c(1,2,4,5,3)]#The relevant columns are selected

new1<-merge(result1000,result4000, by=c("Subject","Stimulation","CDR3a","CDR3b","TCRusage"))#This gives the phenotypic gene expressions for the top 3 dominant clones


new1$Clonename <-do.call(paste0, new1[c("Subject","Stimulation","CDR3a","CDR3b","TCRusage")])#Create a new column called "Clonename" which put together subject. stimulation, CDR3a, CDR3b, TCRusage  
new1<-new1[,-c(1,2,3,4,5)]#The relevant columns are selected

#All the top dominant clones are named

a1=grep("ECno stimCAVRSSYNTDKLIFCASSQDTNTGELFFTRAV1-2TRAJ34TRBV7-2TRBJ2-2",new1$Clonename)
new1$Clonename[a1]<-"Clone 1"

a2=grep("ECno stimCAVTDSNYQLIWCSVEVGTAHSETQYFTRAV1-2TRAJ33TRBV29-1TRBJ2-5",new1$Clonename)
new1$Clonename[a2]<-"Clone 2"

a3=grep("ECno stimCAVSDSNYQLIWCSARSPGTHNEQFFTRAV1-2TRAJ33TRBV20-1TRBJ2-1",new1$Clonename)
new1$Clonename[a3]<-"Clone 3"

a4=grep("ECstimCAVRSSYNTDKLIFCASSQDTNTGELFFTRAV1-2TRAJ34TRBV7-2TRBJ2-2",new1$Clonename)
new1$Clonename[a4]<-"Clone 4"

a5=grep("ECstimCAVRDSNYQLIWCAWGQGGGHVGELFFTRAV1-2TRAJ33TRBV30TRBJ2-2",new1$Clonename)
new1$Clonename[a5]<-"Clone 5"

a6=grep("ECstimCAVRDSNYQLIWCASSQDTNTGELFFTRAV1-2TRAJ33TRBV7-2TRBJ2-2",new1$Clonename)
new1$Clonename[a6]<-"Clone 6"

a7=grep("HDno stimCAVLDSNYQLIWCSATRGPDFYEQYFTRAV1-2TRAJ33TRBV20-1TRBJ2-7",new1$Clonename)
new1$Clonename[a7]<-"Clone 7"

a8=grep("HDno stimCVPMDSNYQLIWCSARLGTPNQAGVQETQYFTRAV1-2TRAJ33TRBV20-1TRBJ2-5",new1$Clonename)
new1$Clonename[a8]<-"Clone 8"

a9=grep("HDno stimCAVLDSNYQLIWCSARDVAGDSYNEQFFTRAV1-2TRAJ33TRBV20-1TRBJ2-1",new1$Clonename)
new1$Clonename[a9]<-"Clone 9"

a10=grep("HDstimCVPMDSNYQLIWCSARLGTPNQAGVQETQYFTRAV1-2TRAJ33TRBV20-1TRBJ2-5",new1$Clonename)
new1$Clonename[a10]<-"Clone 10"

a11=grep("HDstimCAVLDSNYQLIWCSATRGPDFYEQYFTRAV1-2TRAJ33TRBV20-1TRBJ2-7",new1$Clonename)
new1$Clonename[a11]<-"Clone 11"

a12=grep("HDstimCAVRDSNYQLIWCSARLGTPNQAGVQETQYFTRAV1-2TRAJ33TRBV20-1TRBJ2-5",new1$Clonename)
new1$Clonename[a12]<-"Clone 12"

a13=grep("HDstimCAVRDSNYQLIWCASSETEDGANVLTFTRAV1-2TRAJ33TRBV10-2TRBJ2-6",new1$Clonename)
new1$Clonename[a13]<-"Clone 13"

a14=grep("HIVno stimCAVTNSNYQLIWCASSDGTGHTGELFFTRAV1-2TRAJ33TRBV6-1TRBJ2-2",new1$Clonename)
new1$Clonename[a14]<-"Clone 14"

a15=grep("HIVno stimCAGLDSNYQLIWCASSSNLGGGGDEQFFTRAV1-2TRAJ33TRBV6-2TRBJ2-1",new1$Clonename)
new1$Clonename[a15]<-"Clone 15"

a16=grep("HIVno stimCAVVLGDYKLSFCASSSTGEGNQPQHFTRAV1-2TRAJ20TRBV6-4TRBJ1-5",new1$Clonename)
new1$Clonename[a16]<-"Clone 16"

a17=grep("HIVstimCAGLDSNYQLIWCASSSNLGGGGDEQFFTRAV1-2TRAJ33TRBV6-2TRBJ2-1",new1$Clonename)
new1$Clonename[a17]<-"Clone 17"

a18=grep("HIVstimCAVTNSNYQLIWCASSDGTGHTGELFFTRAV1-2TRAJ33TRBV6-1TRBJ2-2",new1$Clonename)
new1$Clonename[a18]<-"Clone 18"

a19=grep("HIVstimCAVTDSNYQLIWCASSPLAGADTQYFTRAV1-2TRAJ33TRBV6-1TRBJ2-3",new1$Clonename)
new1$Clonename[a19]<-"Clone 19"

a20=grep("HIVstimCAVRDSNYQLIWCASSQEGASSYEQYFTRAV1-2TRAJ33TRBV4-2TRBJ2-7",new1$Clonename)
new1$Clonename[a20]<-"Clone 20"

a21=grep("HIVstimCAGLDSNYQLIWCASSIGLGSSYEQYVTRAV1-2TRAJ33TRBV6-1TRBJ2-7",new1$Clonename)
new1$Clonename[a21]<-"Clone 21"




write.csv(new1, "Phenotypic gene expressions for the dominant clones.csv") #Gives the required output in csv form







