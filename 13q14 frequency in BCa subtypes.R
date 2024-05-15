#########################################################
#########################################################
############13q14 CNA in breast cancer subtypes (Metabric)
#########################################################
########################################################
install.packages("tidyverse")
install.packages("easypackages")
library(tidyverse)
library(easypackages)
library(readxl)

path_1 <-"/Users/parasts/Library/CloudStorage/GoogleDrive-p.shahrouzi@gmail.com/My Drive/Postdoc files/13q14.2 paper/Figure 1/new analysis to replace pie charts/METABRIC_data_cna.xlsx"
path_2 <-"/Users/parasts/Library/CloudStorage/GoogleDrive-p.shahrouzi@gmail.com/My Drive/Postdoc files/13q14.2 paper/Figure 1/new analysis to replace pie charts/METABRIC_data_clinical.xlsx"

CNA <- read_excel(path_1)
clinical <- read_excel(path_2)


Metabric_CNA <- na.omit(CNA )
Metabric_CNA_2 <- Metabric_CNA[,-c(1,2)]

CNA13q14 <- Metabric_CNA_2$Cytoband == "13q14.2"
CNA13q14_df <- Metabric_CNA_2 %>% 
  filter(Cytoband == "13q14.2")

CNA13q14_df <- na.omit(CNA13q14_df )
CNA13q14  <- CNA13q14_df [,-1]


convert_to_numeric <- function(x) {
  as.numeric(as.character(x))
}
CNA13q14[] <- lapply(CNA13q14, convert_to_numeric) # to make the values numeric in each column


CNA13q_rounded <- data.frame (round (colMeans (CNA13q14))) # to round each culumn to one value
colnames (CNA13q_rounded) <- "av_13q"
CNA13q_rounded$Patient.ID <- gsub ("\\.", "-", rownames (CNA13q_rounded)) # to make a column. with patient ID


clinic_13q14 <- merge (clinical, CNA13q_rounded, by.x="Patient.ID" , by.y="Patient.ID", all.x=T)


install.packages("openxlsx")
library(openxlsx)

excel_file_path <- "/Users/parasts/Library/CloudStorage/GoogleDrive-p.shahrouzi@gmail.com/My Drive/Postdoc files/13q14.2 paper/Figure 1/new analysis to replace pie charts/Metabric_CNA-Clincal_merged.xlsx"
write.xlsx(clinic_13q14 , excel_file_path, rowNames = FALSE) #to save the merged CNA and clinical dataset in my folder



selected_columns <- clinic_13q14[, c("Pam50 + Claudin-low subtype", "av_13q")]
selected_columns <- na.omit(selected_columns)
colnames(selected_columns)[colnames(selected_columns) == "Pam50 + Claudin-low subtype"] <- "Pam50"
selected_columns$Pam50 <- factor(selected_columns$Pam50)
selected_columns <- na.omit(selected_columns)

excel_file_path <- "/Users/parasts/Library/CloudStorage/GoogleDrive-p.shahrouzi@gmail.com/My Drive/Postdoc files/13q14.2 paper/Figure 1/new analysis to replace pie charts/Metabric_CNA-Clincal_merged_2.xlsx"
write.xlsx(selected_columns , excel_file_path, rowNames = FALSE) #to save the merged CNA and clinical dataset in my folder

selected_columns$av_13q <- as.factor(selected_columns$av_13q) #covert the numeric CNAs to factor to be able to group them in ggplot

to_remove <- c("NA", "NC", "Normal")  # to remove those we dont want in the graph
filtered_data <- subset(selected_columns, !(Pam50 %in% to_remove))


ggplot(filtered_data, aes(x = Pam50, fill = av_13q)) +
  geom_bar() +
  labs( x = "PAM50 subtypes", y = "Patient count") +
  scale_fill_manual(values = c("blue", "pink", "#8da0cb", "#e78ac3", "#a6d854"), name = "13q14.2 CNA") +
  theme_minimal()


excel_file_path <- "/Users/parasts/Library/CloudStorage/GoogleDrive-p.shahrouzi@gmail.com/My Drive/Postdoc files/13q14.2 paper/Figure 1/new analysis to replace pie charts/Metabric_CNA-Clincal_merged_3.xlsx"
write.xlsx(selected_columns , excel_file_path, rowNames = FALSE) #to save the merged CNA and clinical dataset in my folder

ggplot(filtered_data, aes(x = Pam50, fill = av_13q)) +    # to have the percentage of subtypes
  geom_bar(position = "fill") +
  labs(x = "PAM50 subtypes", y = "Percentage of patients") +
  scale_fill_manual(values = c("blue", "pink", "#8da0cb", "#e78ac3", "#a6d854"), name = "13q14.2 CNA") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 12),  
        axis.title.x = element_text(size = 13),
        axis.text.y = element_text(size = 12),   
        axis.title.y = element_text(size = 13))  

summary_counts <- summary(filtered_data) # to count each variable
print(summary_counts)
contingency_table <- table(filtered_data$Pam50, filtered_data$av_13q)
print(contingency_table)



Counted_subtypes_CNA <- filtered_data %>%
  mutate(combined_variants = paste(Pam50, av_13q, sep = "_")) %>% # to chount CNAs in each subtypes
  count(combined_variants)





#########################################################
#########################################################
############13q14 CNA in breast cancer subtypes (TCGA)
#########################################################
########################################################


excel_file_path <- "/Users/parasts/Library/CloudStorage/GoogleDrive-p.shahrouzi@gmail.com/My Drive/Postdoc files/13q14.2 paper/Figure 1/new analysis to replace pie charts/TCGA 13q14.2 clinical and CNA data-merged.xlsx"
TCGA<- read_excel(excel_file_path)


TCGA_2 <- TCGA[, c("PAM50_mRNA_nature2012", "CNA 13q14.2")]
TCGA_2 <- na.omit(TCGA_2)
TCGA_2 $PAM50_mRNA_nature2012 <- factor(TCGA_2$PAM50_mRNA_nature2012)
colnames(TCGA_2)[colnames(TCGA_2) == "CNA 13q14.2"] <- "CNA_13q14.2"
TCGA_2$CNA_13q14.2 <- as.factor(TCGA_2$CNA_13q14.2)


ggplot(TCGA_2, aes(x = PAM50_mRNA_nature2012, fill = CNA_13q14.2)) +    # to have the percentage of subtypes
  geom_bar(position = "fill") +
  labs(x = "PAM50 subtypes", y = "Percentage of patients") +
  scale_fill_manual(values = c("blue", "pink", "#8da0cb", "#e78ac3", "#a6d854"), name = "13q14.2 CNA") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 12),  
        axis.title.x = element_text(size = 13),
        axis.text.y = element_text(size = 12),   
        axis.title.y = element_text(size = 13))  

summary_counts <- summary(TCGA_2) # to count each variable
print(summary_counts)
contingency_table <- table(TCGA_2$PAM50_mRNA_nature2012, TCGA_2$CNA_13q14.2)
print(contingency_table)



Counted_subtypes_CNA_TCGA <- TCGA_2 %>%
  mutate(combined_variants = paste(PAM50_mRNA_nature2012, CNA_13q14.2, sep = "_")) %>% # to chount CNAs in each subtypes
  count(combined_variants)




