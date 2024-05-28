## To assess association (co-occurance / exclussion) between loss of 13q14.2 and
## mutations in metabric dataset. Similar analysis was carried on TCGA
## Youness Aizmzade, 23.05.2024
## younessazimzade@gmail.com


library(reshape2)
library(dplyr) 
library(readr)
library(readxl)

Genes <- readxl::read_excel("~/Mutated genes in METABRIC.xlsx")
Genes$Freq <- gsub("%", "", Genes$Freq) # Remove all instances of '%'
Genes$Freq <- gsub("<", "", Genes$Freq) # Remove all instances of '<'
Genes$Freq <- as.numeric(Genes$Freq)

Genes2 <- Genes[Genes$Freq >2.0,]   ## Selecting genes with more than 2 % frequency

Data<-read.delim("~/metabric_mutation.txt")
colnames(Data) <- Data[1,]
Data <- Data[-1,]

Data2 <- dcast(Data, Hugo_Symbol ~ Tumor_Sample_Barcode )  
colnames(Data2) <- gsub("-", ".", colnames(Data2))

Data2 <- as.data.frame(t(Data2))

Data2 <-Data2[, Data2[1,]%in% Genes2$Gene ]

colnames(Data2) <- Data2[1,]; Data2 <- Data2[-1,]


CNA<-read_excel("~/MBRC/CNA.xlsx")
colnames(CNA) <- gsub("-", ".", colnames(CNA))


CNA2 <- as.data.frame(t(subset(CNA, CNA[,3]=="13q14.2")))

CNA2 <- CNA2[-1:-3,]
CNA3 <-CNA2
CNA2<-as.data.frame(lapply(CNA2,as.numeric))
CNA2$Mean<-rowMeans(CNA2)
rownames(CNA2) <- rownames(CNA3)

CNA2$Mean  <- ifelse(CNA2$Mean<0, 1,0)

Data2$`13q_Loss` <- CNA2$Mean[match( rownames(Data2),rownames(CNA2))]
Data3 <- Data2  


Data2<-as.data.frame(lapply(Data2,as.numeric))
Data2[Data2>1] <- 1 

rownames(Data2) <- rownames(Data3)
colnames(Data2) <- colnames(Data3)

Data2 <-na.omit(Data2)

# Assuming df is your dataframe and nameColumn is the column to check
Data2 <- Data2[!grepl("MTS", rownames(Data2), ignore.case = TRUE), ]


library(dplyr)

# Initialize a data frame to store results
results <- data.frame(Gene = character(), P_Value = numeric(), Odds_Ratio = numeric(), stringsAsFactors = FALSE)

for(i in colnames(Data2[, -ncol(Data2)])) {  # Exclude the last column (13q_Loss) during iteration
  # Create contingency table for each gene vs. 13q_Loss
  contingency_table <- table(Data2[[i]], Data2$`13q_Loss`)
  
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

 
results$Direction_Sign <- ifelse(results$Odds_Ratio > 1, 1, -1) 


results$`-log10(p_value)*Direction` <- -log10(results$P_Value )*results$Direction_Sign


results <- results[results$P_Value <0.05,]
# View the updated results
results

results<-results[order(-results$`-log10(p_value)*Direction`),]

results$Gene<-as.character(results$Gene)
results$Gene<-factor(results$Gene,levels=results$Gene)


library(ggplot2)
pdf("~/MBRC Mutation Co-occurance/Exclusion.pdf",width=14,height=7)
ggplot(results)+
  geom_bar(aes(x=`Gene`, y=`-log10(p_value)*Direction`, fill= as.factor(Direction_Sign)),stat="identity",width = 0.5)+ 
  theme_bw()+theme(text=element_text(size=30))+theme(axis.text.x=element_text(angle=90))+
  theme(legend.position="none")+ 
  scale_x_discrete(position="bottom") +
  ylim(-6,22)+geom_segment(x=0.55,y=20.5,xend=0.55,yend=-4.05,lineend="round", 
                           linejoin="round",size=1,arrow=arrow(length=unit(0.15,"inches"),ends="both"),
                           colour="black")+ annotate(geom="text", size =8,x=1.4,y=21.9,label="co-occurrence",color="black")+ 
  annotate(geom="text",size=8,x=1.6, y=-5.13,label="mutual exclusivity",color="black")+
  theme(axis.title.x=element_text(margin=margin(t=9)))
dev.off() 



