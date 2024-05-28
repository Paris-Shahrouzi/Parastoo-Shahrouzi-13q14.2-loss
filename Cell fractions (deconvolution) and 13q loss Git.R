## To assess association between loss of 13q14.2 and
## cell fractions in TME calculated using deconvolution 
## in metabric dataset. Similar analysis was carried on TCGA
## Youness Aizmzade, 23.05.2024
## younessazimzade@gmail.com


# Load required libraries
library(reshape2)
library(dplyr)
library(readr)
library(readxl)
library(rstatix)

# Set relative paths and constants
base_path <- "~/Data/MBRC"
setwd(base_path)
path <- "MBRC"
readpath <- file.path(base_path, path)

# Define subgroups
Subs <- c("all patients", "ER", "HER2", "TNBC")

# Read data
CellPop <- read.csv(file.path(readpath, "DataF.csv"))
rownames(CellPop) <- CellPop$Mixture

# Manipulate data
CellPop <- CellPop %>%
  mutate(
    TNBC = ifelse(ER == "Negative" & HER2 == "Negative" & PgR == "Negative", "TNBC", NA),
    `all patients` = "all patients"
  )

# Read and filter CNA data
CNA <- read_excel(file.path(readpath, "CNA.xlsx")) %>%
  filter(Cytoband == "13q14.2") %>%
  rename_with(~gsub("-", ".", .))

# Initialize empty dataframe to accumulate results
testAll <- data.frame()

# Analyze each subgroup
for (targetMarker in Subs) {
  # Subset CellPop based on current targetMarker
  colmnames <- c(names(CellPop)[1:29], targetMarker)
  CellPop2 <- CellPop[, names(CellPop) %in% colmnames]
  
  # Subset and transform CNA data
  CNA2 <- CNA[, names(CNA) %in% CellPop$Mixture]
  CNA2 <- data.frame(t(CNA2)) %>%
    mutate_all(as.numeric) %>%
    mutate(Mean = rowMeans(.)) %>%
    mutate(Mean = ifelse(Mean > 0, "WT", ifelse(Mean < 0.0, "Loss", "WT")))
  
  # Combine data
  CNA2 <- cbind(rownames(CNA2), CNA2)
  rownames(CNA2) <- gsub("-", ".", rownames(CNA2))
  
  # Calculate Loss Status
  CellPop2$LossStat <- CNA2$Mean[match(CellPop2$Mixture, rownames(CNA2))]
  
  # Prepare for statistical test
  CellPop.m <- melt(CellPop2[,-1]) %>%
    rename_with(~c("targetMarker_stat", "LossStat", "CellType", "Fraction")) %>%
    na.omit()
  
  # Perform t-test
  stat.test <- CellPop.m %>%
    group_by(CellType, targetMarker_stat) %>%
    t_test(Fraction ~ LossStat) %>%
    adjust_pvalue(method = "bonferroni") %>%
    add_significance("p.adj") %>%
    cbind(targetMarker, .)
  
  # Append results
  testAll <- rbind(testAll, stat.test)
}

# Finalize data frame
testAll$Data <- "METABRIC"

# Write to CSV
write.csv(testAll, file = file.path(base_path, "Stats2.csv"), sep = "\t", row.names = FALSE, quote = FALSE)





 