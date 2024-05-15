install.packages("easypackages") 
install.packages("tidyverse")
library(tidyverse)
library(readxl)

METABRIC_1 <- read_excel("/Users/parasts/Library/CloudStorage/GoogleDrive-p.shahrouzi@gmail.com/My Drive/Postdoc files/13q14.2 paper/untitled folder/Figure 2/METABRIC_data_cna_1.xlsx")
METABRIC_2 <- na.omit(METABRIC_1)  # To remove all the NA conditions

METABRIC_3 <- METABRIC_2[,-c(1,2)]

METABRIC_4 <- METABRIC_3 %>% distinct()



METABRIC_5 <- METABRIC_4 %>%    
  pivot_longer(
    cols = !Cytoband, 
    names_to = "Patient", 
    values_to = "CNA",
    values_transform = list(CNA = as.character))  # I had to add this code because I was getting this error:  Can't combine `MB-0000` <double> and `MB-4602` <character>.

METABRIC_6  <- filter(METABRIC_5 , CNA==c(-1,-2))

METABRIC_7 <- METABRIC_6  %>% distinct()

METABRIC_8 <- METABRIC_7  %>%    
  mutate_at(vars(CNA), factor) 

METABRIC_9 <- METABRIC_8 %>%   
  group_by(Cytoband, CNA) %>%       
  summarise(Patient = n()) 

METABRIC_10<- METABRIC_9[order(METABRIC_9$Patient,decreasing=TRUE),]

First_50 <- METABRIC_10 [1:50,]

vector50 <- c(First_50$Cytoband)

METABRIC_11<- subset(METABRIC_10, Cytoband %in% vector50)

METABRIC_12 <- METABRIC_11[-45,]

METABRIC_Final <- METABRIC_12[-c(1,50),]

METABRIC_Final_2 <- METABRIC_Final[-88,]

install.packages("viridis")  
library("viridis")

ggplot(data = METABRIC_Final_2 ,         
       mapping = aes(x = reorder(Cytoband, -Patient), y = Patient, fill = CNA)) +   # plot it awesomely ! :D
  geom_bar(position="stack", stat = "identity") + 
  theme_bw() +
  scale_fill_manual(values = c("pink", "blue")) +
  # scale_fill_viridis_d() + 
  theme(axis.text.x = element_text(angle = 90, size = 11),
        legend.position = "bottom") +
  labs(x = "Cytoband", y ="Patient number", fill = "Copy number alteration")

unique(METABRIC_Final_2 $CNA)


