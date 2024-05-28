## To compare alteration in  13q14.2 and rest of cytobands in metabric dataset. Similar analysis was carried on TCGA
## Youness Aizmzade, 23.05.2024
## younessazimzade@gmail.com




library(reshape2)
library(dplyr) 
library(readr)
library(readxl)

path <- "MBRC"

 
readpath <-paste0("~/",path,"/")

 
CNA<-read_excel(paste0(readpath, "CNA.xlsx")) 
CNA2 <- subset(CNA, CNA$Cytoband=="13q14.2")
colnames(CNA2) <- gsub("-", ".", colnames(CNA2))


#CNA2 <- CNA 

CNA2<-data.frame(t(CNA2[,-1:-3]))
CNA3 <-CNA2

CNA2<-as.data.frame(lapply(CNA2,as.numeric))

CNA2$Mean<-rowMeans(CNA2)

hist(CNA2$Mean)

library(ggplot2)  

 



# Assuming your data is normally distributed and you want to overlay a normal distribution curve
CNAA<-read_excel(paste0(readpath, "CNA.xlsx"))
CNAA.m <- melt(lapply(CNAA[,4: ncol(CNAA)],as.numeric) )

CNAA.m <- na.omit(CNAA.m)
mean_val <- mean(CNAA.m$value)
sd_val <- sd(CNAA.m$value)

  

# Save the plot to a PDF
pdf("~/MBRC Selection.pdf", width = 7.5, height =5)
# Plotting
ggplot(CNA2, aes(x = Mean)) + 
  geom_histogram(aes(y = ..density.., fill = "Histogram"), binwidth = 1, color = "black", alpha = 0.5) +
  geom_line(aes(color = "Normal Distribution"), stat = "function", fun = dnorm, args = list(mean = mean_val, sd = sd_val), size = 1) +
  scale_color_manual(name = "", values = "red", labels = "Rest of the Genes") +
  scale_fill_manual(name = "", values = "blue", labels = "13q14.2 alterations") +
 #ggtitle("METABRIC") +
  xlab("CNA") +
  ylab("Frequency") +
  theme_bw()+xlim(-2.5,2.5) +  theme(text=element_text(size =18))+
  theme(legend.position = "right", legend.title = element_blank())
dev.off()


CNAA[,4: ncol(CNAA)] <- lapply(CNAA[,4: ncol(CNAA)],as.numeric)  

CNAA.m <- melt( CNAA[,3: ncol(CNAA)])

CNAA.m$Cytoband <- ifelse(CNAA.m$Cytoband=="13q14.2", "13q14.2", "Other")
colnames(CNAA.m)<- c("Cytoband", "Sample", "CNAs")

library(rstatix)

stat.test<-CNAA.m%>%
 # group_by(CellType,targetMarker_stat)%>%
  t_test(CNAs~Cytoband)%>%
  adjust_pvalue(method="bonferroni")%>%
  add_significance("p.adj")


