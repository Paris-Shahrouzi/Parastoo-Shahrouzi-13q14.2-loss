## To assess association between loss of 13q14.2 and
## cell fractions in TME calculated from spatial omics data
## in metabric dataset.  
## Youness Aizmzade, 23.05.2024
## younessazimzade@gmail.com

# Load required libraries
library(dplyr)
library(readr)
library(readxl)
library(stringi)
library(reshape2)
library(rstatix)
library(ggplot2)

# Set working directory (Modify this path for your local setup)
setwd("~/Data/MBRC")  

# Define constants and paths
path <- "MBRC"
readpath <- file.path("~/Data", path)

# Subgroups of interest
Subs <- c("all patients", "ER", "HER2", "TNBC")

# Read data
CellPop0 <- read.csv(file.path(readpath, "DataF.csv"))
raw_text <- readLines(file.path(readpath, "Frequencies.csv"), encoding = "UTF-8")
correct_text <- stri_encode(raw_text, from = "UTF-8", to = "UTF-8")

# Process and combine text data
combined_text <- paste(correct_text, collapse = "\n")
CellPop <- read_csv(combined_text)

# Data manipulation
CellPop <- CellPop %>%
  mutate(
    HER2 = CellPop0$HER2[match(metabric_id, CellPop0$Mixture)],
    PgR = CellPop0$PgR[match(metabric_id, CellPop0$Mixture)],
    TNBC = ifelse(ER == "Negative" & HER2 == "Negative" & PgR == "Negative", "TNBC", NA),
    `all patients` = "all patients"
  ) %>%
  rownames_to_column(var = "Mixture")

# Normalization
CellPop[2:33] <- sweep(CellPop[2:33], 1, rowSums(CellPop[2:33]), FUN = "/")

# Read and filter CNA data
CNA <- read_excel(file.path(readpath, "CNA.xlsx")) %>%
  filter(Cytoband == "13q14.2") %>%
  rename_with(~ gsub("-", ".", .))

# Analyze each marker
testAll <- data.frame()
for (targetMarker in Subs) {
  colmnames <- c(names(CellPop)[1:29], targetMarker)
  CellPop2 <- CellPop[colnames(CellPop) %in% colmnames]
  CNA2 <- CNA[colnames(CNA) %in% CellPop$Mixture]
  
  # Data reshaping and calculations
  CNA2 <- data.frame(t(CNA2)) %>%
    mutate_all(as.numeric) %>%
    mutate(Mean = rowMeans(.), 
           Mean = ifelse(Mean > 0, "WT", ifelse(Mean < 0.0, "Loss", "WT")))
  
  CNA2 <- cbind(rownames(CNA2), CNA2)
  CNA2$rownames.CNA2. <- gsub("-", ".", CNA2$rownames.CNA2.)
  
  CellPop2$LossStat <- CNA2$Mean[match(CellPop2$Mixture, CNA2$rownames.CNA2.)]
  
  # Statistical testing
  stat.test <- melt(CellPop2[,-1]) %>%
    na.omit() %>%
    group_by(CellType, targetMarker_stat) %>%
    t_test(Fraction ~ LossStat) %>%
    adjust_pvalue(method = "bonferroni") %>%
    add_significance("p.adj") %>%
    cbind(targetMarker, .)
  
  testAll <- rbind(testAll, stat.test)
}

# Output results
write.csv(testAll, file = file.path(readpath, "StatsSpTr.csv"), sep = "\t", row.names = FALSE, quote = FALSE)

# Additional processing for visualization (if needed)
# Visualizations and further analysis would follow here, similar to the ggplot2 usage and data transformation shown previously



 