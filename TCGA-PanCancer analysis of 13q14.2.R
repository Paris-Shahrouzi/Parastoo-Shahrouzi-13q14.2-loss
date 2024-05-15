#########################################################
#########################################################
############13q14 CNA in PAN CANCER (TCGA)
#########################################################
########################################################
install.packages("tidyverse")
install.packages("easypackages")
library(tidyverse)
library(easypackages)
library(readxl)

path_1 <-"/Users/parasts/Library/CloudStorage/GoogleDrive-p.shahrouzi@gmail.com/My Drive/Postdoc files/13q14.2 paper/Figure 1/PANCANCER analysis/TCGA-PAN CANCER-CN-13q14.2.xlsx"
path_2 <- "/Users/parasts/Library/CloudStorage/GoogleDrive-p.shahrouzi@gmail.com/My Drive/Postdoc files/13q14.2 paper/Figure 1/PANCANCER analysis/pancan_pcawg_2020_clinical_data.xlsx"


CNA_PanCancer_13q14.2 <- read_excel(path_1)
cancer_name <- read_excel(path_2, sheet="ID")
colnames(cancer_name)[which(names(cancer_name) == "Sample ID")] <- "Patient.ID"

CNA13q <- CNA_PanCancer_13q14.2[,-1]
CNA13q <- na.omit(CNA13q )
CNA13q <- apply(CNA13q, 2, as.numeric)
CNA13q <- data.frame (round (colMeans (CNA13q)))

CNA13q$Patient.ID <- gsub ("\\.", "-", rownames (CNA13q))

TCGA_PANCANCER_CNA13q <- merge(cancer_name, CNA13q, by = "Patient.ID", all = TRUE)

colnames(TCGA_PANCANCER_CNA13q)[which(names(TCGA_PANCANCER_CNA13q) == "round.colMeans.CNA13q..")] <- "13q14.2_CNA"
colnames(TCGA_PANCANCER_CNA13q)[which(names(TCGA_PANCANCER_CNA13q) == "Cancer Type")] <- "Cancer_type"

TCGA_PANCANCER_CNA13q$`13q14.2_CNA` <- as.factor(TCGA_PANCANCER_CNA13q$`13q14.2_CNA`)
TCGA_PANCANCER_CNA13q <- na.omit(TCGA_PANCANCER_CNA13q )

count_summary <- TCGA_PANCANCER_CNA13q %>%
  group_by(Cancer_type, `13q14.2_CNA`) %>%
  summarise(count = n()) %>%
  ungroup()

sorted_df <- TCGA_PANCANCER_CNA13q %>%
  arrange(desc(`13q14.2_CNA`))

patient_count <- TCGA_PANCANCER_CNA13q %>% #to sort from high to low
  count(Cancer_type) %>%
  arrange(desc(n))

TCGA_PANCANCER_CNA13q$Cancer_type <- factor(TCGA_PANCANCER_CNA13q$Cancer_type, levels = patient_count$Cancer_type)
TCGA_PANCANCER_CNA13q <- na.omit(TCGA_PANCANCER_CNA13q)


count_summary <- TCGA_PANCANCER_CNA13q %>%
  group_by(Cancer_type, `13q14.2_CNA`) %>%
  summarise(count = n()) %>%
  ungroup()

# Calculate the total count of 13q14.2_CNA for each cancer type
total_counts <- count_summary %>%
  group_by(Cancer_type) %>%
  summarise(total_count = sum(count))


# Merge the total counts with the count summary
count_summary <- merge(count_summary, total_counts, by = "Cancer_type")

# Calculate the percentage of each level of 13q14.2_CNA in each cancer type
count_summary <- count_summary %>%
  mutate(percentage = (count / total_count) * 100)


if (!requireNamespace("writexl", quietly = TRUE)) {
  install.packages("writexl")
}
library(writexl)

file_path <- "/Users/parasts/Library/CloudStorage/GoogleDrive-p.shahrouzi@gmail.com/My Drive/Postdoc files/13q14.2 paper/Figure 1/PANCANCER analysis/count_summary.xlsx"

# Write your dataset to an Excel file
write_xlsx(count_summary, file_path)



# Plot with reordered Cancer_type and percentage
ggplot(count_summary, aes(x = reorder(Cancer_type, -count), y = percentage, fill = factor(`13q14.2_CNA`))) +
  geom_bar(stat = "identity") +
  labs(x = "Cancer Type", y = "Percentage of 13q14.2_CNA", fill = "13q14.2 CNA") +
  scale_fill_manual(values = c("blue", "pink", "#8da0cb", "#e78ac3", "#a6d854"), name = "13q14.2 CNA") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))





ggplot(TCGA_PANCANCER_CNA13q, aes(x = reorder(Cancer_type, -`13q14.2_CNA`), fill = `13q14.2_CNA`)) +
  geom_bar() +
  labs(x = "Cancer Type", y = "Patient Number") +
  scale_fill_manual(values = c("blue", "pink", "#8da0cb", "#e78ac3", "#a6d854"), name = "13q14.2 CNA") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))


Counted_13q <- TCGA_PANCANCER_CNA13q %>%
  mutate(combined_variants = paste(Cancer_type, `13q14.2_CNA`, sep = "_")) %>% # to chount 13 q CNAs in each cancer type
  count(combined_variants)













