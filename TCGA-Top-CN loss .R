install.packages("easypackages")
library(easypackages)    # Install -> install.packages("easypackages")
libraries("ggplot2", "tidyverse", "dplyr", "readxl")   # Upload packages
library(ggplot2)
library(tidyverse)
library(dplyr)
library(readxl)

path<-"/Users/parasts/Library/CloudStorage/GoogleDrive-p.shahrouzi@gmail.com/My Drive/Postdoc files/13q14.2 paper/Figure 1/TCGA_BRCA_CN_thresholded_by_gene.xlsx"  #1 copypaste my path

TCGA_data <- read_excel(path, sheet = "all_thresholded.by_genes")  #2 to read .xlsx documents and the sheet of intrest in excel

length(unique(TCGA_data$Cytoband))

TCGA_data_2 <- TCGA_data[,-c(1,2)]   #3 to remove the 2 firsts columns

TCGA_data_3 <- TCGA_data_2 %>%       #4 save the altered data.frame in TCGA_data_3
  pivot_longer(cols = !Cytoband,     # run pivot_longer function, conserve Cytoband column using !Cytoband (! means no)
               names_to = "Patient", # each column were 1 pacient and now each row will be 1 patient per cytoband
               values_to = "CNA") #  previous cell values are now the values of the column CNA

length(unique(TCGA_data_3$Patient))  #5 unique counts the number of patients, using length over unique I got the number
                                    
TCGA_no_duplicates<-TCGA_data_3  %>% distinct()    #6 or use the this function to easily remove duplicates (Paris discovered this haha!)


TCGA_no_duplicates <- unique( TCGA_data_3[ , 1:3 ] ) #7 unique counts unique patients of data.frame for the three columns

Plotting_data <- TCGA_no_duplicates %>%          #8 save the new data.frame in Plotting_data
    group_by(Cytoband, CNA) %>%       # you have to group the data by cytoband and CNA to count the number of patients
    summarise(Patient = n())                # counting the number of patinets



TCGA_Loss<- Plotting_data %>% 
  filter(CNA <= 0) #  < , <= , ==, >=, >


CNA_loss_only <- Plotting_data %>% 
  mutate_at(vars(CNA), factor) %>% #9 to select only -1 and -2
  filter(CNA %in% c(-1,-2))


ggplot(data = CNA_loss_only,          # data to plot
       mapping = aes(x = reorder(Cytoband, -Patient), y = Patient, fill = CNA)) +  # graph coordinates, x, y and group
  geom_bar(position="stack", stat = "identity") + 
  theme_bw() 

---------------

First_30 <- CNA_loss_only  %>% 
  group_by(CNA) %>% slice(which.max(Patient))   #this is to only find one max



CNA_loss_sorted <- CNA_loss_only[order(CNA_loss_only$Patient,decreasing=TRUE),] #this seems to give me the sorted ones but it is taking into consideration the CNA!!

filt50 <- CNA_loss_sorted[1:50,]

vector50 <- c(filt50$Cytoband)

prueba <- subset(CNA_loss_sorted, Cytoband %in% vector50)

new_data <- CNA_loss_sorted %>% filter(Cytoband == vector50)

library("plyr")
newdata <- arrange(CNA_loss_only,desc(Patient),CNA)

library("viridis")
ggplot(data = prueba,          # data to plot
       mapping = aes(x = reorder(Cytoband, -Patient), y = Patient, fill = CNA)) +  # graph coordinates, x, y and group
  geom_bar(position="stack", stat = "identity") + 
  theme_bw() +
  scale_fill_manual(values = c("green", "pink")) +
  # scale_fill_viridis_d() + 
  theme(axis.text.x = element_text(angle = 90, size = 12),
        legend.position = "bottom") +
  labs(x = "Cytoband", y ="Count of Pacients", fill = "Copy number alteration")

# finaldf <- CNA_loss_sorted %>%
#   group_by(Patient, Cytoband) %>%       # for each ID and Phase
#   top_n(100)

ggplot(data = CNA_loss_sorted,          # data to plot
       mapping = aes(x = reorder(Cytoband, -Patient), y = Patient, fill = CNA)) +  # graph coordinates, x, y and group
  geom_bar(position="stack", stat = "identity") + 
  theme_bw() 



first30 <- CNA_loss_only %>% 
  group_by(Cytoband) %>% top_n(n = 30, wt = Cytoband)

length(unique(first30$Cytoband))


ggplot(data = first30,          # data to plot
       mapping = aes(x = reorder(Cytoband, -Patient), y = Patient, fill = CNA)) +  # graph coordinates, x, y and group
  geom_bar(position="stack", stat = "identity") + 
  theme_bw() 




PatiensperCytoband <- CNA_loss_only%>% group_by(Cytoband) %>% 
  summarise(Patient = n())

df_2 <- merge(CNA_factor, PatiensperCytoband)

  ggplot(data = PatiensperCytoband, mapping = aes(x = Cytoband, y = Patient, fill = CNA)) +  # graph coordinates, x, y and group
  geom_bar(position="stack", stat = "identity") + 
  theme_bw()
