## To assess association (co-occurance / exclussion) between loss of 13q14.2 and
## alterations other cytobands in metabric dataset. Similar analysis was carried on TCGA
## Youness Aizmzade, 23.05.2024
## younessazimzade@gmail.com




library(reshape2) 
library(reshape2)
library(dplyr) 
library(readr)
library(readxl)
topAltered <- 50 

CNA<-read_excel("~/MBRC/CNA.xlsx")
colnames(CNA) <- gsub("-", ".", colnames(CNA))
colnames(CNA)  <- sapply(colnames(CNA) , function(x) {
  n <- nchar(x)
  if (n > 14) {
    substr(x, 1, n - 16)
  } else {
    ""  # or return x if you want to keep the original string when it's shorter than 14 characters
  }
})


CNA.agg <- na.omit(CNA)

# Load the dplyr package
library(dplyr)

# Ensure all columns have valid, unique names
names(CNA) <- make.names(names(CNA), unique = TRUE)

# Now, try the aggregation again
CNA_aggregated <- CNA %>%
  group_by(V3 = CNA[,3]) %>%
  summarise(across(4:ncol(CNA), mean, na.rm = TRUE))

# Note: Replace `sum` with the aggregation function you need, and adjust na.rm as necessary

CNA_aggregated <- as.data.frame(t(CNA_aggregated))
colnames(CNA_aggregated) <- CNA_aggregated[1,]  
 CNA_aggregated  <- CNA_aggregated[-1,]  

CNA2 <- CNA_aggregated

CNA2<-as.data.frame(lapply(CNA2,as.numeric))
rownames(CNA2) <- rownames(CNA_aggregated)
colnames(CNA2) <- colnames(CNA_aggregated)

 
CNA3 <-  na.omit(CNA2 )
CNA3[CNA3>0] <- 0
 

 

CytobandAlter2 <-CNA3
 
 
CytobandAlter4 <- as.data.frame(t(rbind(colSums(CytobandAlter2),CytobandAlter2)))
CytobandAlter4 <- CytobandAlter4[order(-CytobandAlter4[,1]),]
CytobandAlter4 <- CytobandAlter4[1:topAltered,]
# Assuming CytobandAlter4 is your data frame
# Get the order of columns based on the first row values
 

CytobandAlter3 <- CytobandAlter2[, colnames(CytobandAlter2)%in% rownames(CytobandAlter4)]

 
library(dplyr)

# Initialize a data frame to store results
results <- CytobandAlter.frame(Gene = character(), P_Value = numeric(), Odds_Ratio = numeric(), stringsAsFactors = FALSE)

for(i in colnames(CytobandAlter3)) {  # Exclude the last column (13q_Loss) during iteration
  # Create contingency table for each gene vs. 13q_Loss
  contingency_table <- table(CytobandAlter3[[i]], CNA3Los$`13q14.2`)
  
  # Check if the contingency table has the required dimensions (2x2)
  if(all(dim(contingency_table) == c(2, 2))) {
    # Conduct Fisher's Exact Test
    test <- fisher.test(contingency_table, alternative = "two.sided")
    
    # Check if 'estimate' exists in the test result
    if("estimate" %in% names(test)) {
      odds_ratio = test$estimate["odds ratio"]
    } else {
      odds_ratio = NA  # Assign NA if 'estimate' does not exist
    }
    
    # Store results
    results <- rbind(results, data.frame(Gene = i, P_Value = test$p.value, Odds_Ratio = odds_ratio, stringsAsFactors = FALSE))
  } else {
    # Handle cases where the contingency table is not 2x2
    # For example, assign NA to p-value and odds ratio, or handle differently as required
    results <- rbind(results, data.frame(Gene = i, P_Value = NA, Odds_Ratio = NA, stringsAsFactors = FALSE))
  }
}

# View or return the results
results

results$Direction_Sign <- ifelse(results$Odds_Ratio > 1, 1, -1)

results$`-log10(p_value)*Direction` <- -log10(results$P_Value )*results$Direction_Sign

results <- results[results$P_Value <0.05,]
# View the updated results
results

results<-results[order(-results$`-log10(p_value)*Direction`),]

results$Gene<-as.character(results$Gene)
results$Gene<-factor(results$Gene,levels=results$Gene)

resultsLos <- results
 
#########  A similar analysis was carried out for gains and the results were put together.  
