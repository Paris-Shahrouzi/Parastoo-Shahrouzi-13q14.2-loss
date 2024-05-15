############################
### METABRIC ###############
############################
setwd ("C:/Xavier/CNA_13q/files")

CNA <- data.frame(read.table("METABRIC_data_cna.txt", header=T, sep="\t", stringsAsFactors=F))
clinic <- data.frame(read.table("METABRIC_data_clinical.txt", header=T, sep="\t", stringsAsFactors=F))
CNA_13q14.2 <- CNA$Cytoband == "13q14.2"

CNA13q <- CNA[CNA_13q14.2, -c(1:3)]
CNA13q_av <- data.frame (round (colMeans (CNA13q)))
colnames (CNA13q_av) <- "av_13q"
CNA13q_av$Patient.ID <- gsub ("\\.", "-", rownames (CNA13q_av))

clinic_13q <- merge (clinic, CNA13q_av, by.x="Patient.ID" , by.y="Patient.ID", all.x=T)

clinic_13q$bin_13q <- gsub ("-2", "-1", clinic_13q$av_13q)
clinic_13q$bin_13q <- gsub ("2", "1", clinic_13q$bin_13q)
clinic_13q$TNBC <- ifelse (clinic_13q$HER2.Status=="Negative" & clinic_13q$ER.Status=="Negative" & clinic_13q$PR.Status=="Negative", 1, 0)

write.table (clinic_13q, "clinic_metabric_with_13q_xavier_09_11_22.txt",  sep = "\t ")


pdf ("Metabric_KM_13q.pdf")
library(survival)
#### ALL ###
time <- clinic_13q$Overall.Survival..Months
time <- time/12
timelim <- 15
time[time>timelim] <- timelim
status <- clinic_13q$OS== 1
status[time>=timelim] <- FALSE
surv <- Surv(time, status) 
fit <- survfit(surv ~  bin_13q, data=clinic_13q)
test <- survdiff(surv ~  bin_13q, data=clinic_13q)
p <- pchisq(test$chisq, length(test$n) - 1, lower.tail=FALSE)
col <- c("blue", "black", "red")
mylable <- substitute(italic(p) == MYVALUE, list(MYVALUE = format(p, digits= 3)))
main <- "KM 13q14.2 Metabric All samples" 
plot(fit, mark.time = T, xlab="Time (years)", ylab="OS survival", col=col, lwd=6, main=main, cex=1.5)
legend(x="bottomright",col=col,legend= c(paste("deletion (n=", test$n[1],  ")" , sep=""), paste("WT (n=", test$n[2], ")" , sep=""), paste("gain (n=", test$n[3],  ")" , sep="")),bty="n",lwd=2)
legend(x= "topright", legend= mylable,bty="n")

#### ALL  no gain###
clinic_13q_noamp <- clinic_13q[!clinic_13q$bin_13q==1,]
time <- clinic_13q_noamp$Overall.Survival..Months
time <- time/12
timelim <- 15
time[time>timelim] <- timelim
status <- clinic_13q_noamp$OS== 1
status[time>=timelim] <- FALSE
surv <- Surv(time, status) 
fit <- survfit(surv ~  bin_13q, data=clinic_13q_noamp)
test <- survdiff(surv ~  bin_13q, data=clinic_13q_noamp)
p <- pchisq(test$chisq, length(test$n) - 1, lower.tail=FALSE)
col <- c("blue", "black")
mylable <- substitute(italic(p) == MYVALUE, list(MYVALUE = format(p, digits= 3)))
main <- "KM 13q14.2 Metabric All samples  no gain" 
plot(fit, mark.time = T, xlab="Time (years)", ylab="OS survival", col=col, lwd=6, main=main, cex=1.5)
legend(x="bottomright",col=col,legend= c(paste("deletion (n=", test$n[1],  ")" , sep=""), paste("WT (n=", test$n[2], ")" , sep="")),bty="n",lwd=2)
legend(x= "topright", legend= mylable,bty="n")

#### ER+  no gain###
clinic_13q_noamp <- clinic_13q[clinic_13q$ER.Status=="Positive" & !clinic_13q$bin_13q==1,]
time <- clinic_13q_noamp$Overall.Survival..Months
time <- time/12
timelim <- 15
time[time>timelim] <- timelim
status <- clinic_13q_noamp$OS== 1
status[time>=timelim] <- FALSE
surv <- Surv(time, status) 
fit <- survfit(surv ~  bin_13q, data=clinic_13q_noamp)
test <- survdiff(surv ~  bin_13q, data=clinic_13q_noamp)
p <- pchisq(test$chisq, length(test$n) - 1, lower.tail=FALSE)
col <- c("blue", "black", "red")
mylable <- substitute(italic(p) == MYVALUE, list(MYVALUE = format(p, digits= 3)))
main <- "KM 13q14.2 Metabric ERpos samples  no gain" 
plot(fit, mark.time = T, xlab="Time (years)", ylab="OS survival", col=col, lwd=6, main=main, cex=1.5)
legend(x="bottomright",col=col,legend= c(paste("deletion (n=", test$n[1],  ")" , sep=""), paste("WT (n=", test$n[2], ")" , sep=""), paste("WT (n=", test$n[3], ")" , sep="")),bty="n",lwd=2)
legend(x= "topright", legend= mylable,bty="n")

#### ER-  no gain###
clinic_13q_noamp <- clinic_13q[clinic_13q$ER.Status=="Negative" & !clinic_13q$bin_13q==1,]
time <- clinic_13q_noamp$Overall.Survival..Months
time <- time/12
timelim <- 15
time[time>timelim] <- timelim
status <- clinic_13q_noamp$OS== 1
status[time>=timelim] <- FALSE
surv <- Surv(time, status) 
fit <- survfit(surv ~  bin_13q, data=clinic_13q_noamp)
test <- survdiff(surv ~  bin_13q, data=clinic_13q_noamp)
p <- pchisq(test$chisq, length(test$n) - 1, lower.tail=FALSE)
col <- c("blue", "black", "red")
mylable <- substitute(italic(p) == MYVALUE, list(MYVALUE = format(p, digits= 3)))
main <- "KM 13q14.2 Metabric ERneg samples  no gain" 
plot(fit, mark.time = T, xlab="Time (years)", ylab="OS survival", col=col, lwd=6, main=main, cex=1.5)
legend(x="bottomright",col=col,legend= c(paste("deletion (n=", test$n[1],  ")" , sep=""), paste("WT (n=", test$n[2], ")" , sep=""), paste("WT (n=", test$n[3], ")" , sep="")),bty="n",lwd=2)
legend(x= "topright", legend= mylable,bty="n")

#### HER2+  no gain###
clinic_13q_noamp <- clinic_13q[clinic_13q$HER2.Status=="Positive" & !clinic_13q$bin_13q==1,]
time <- clinic_13q_noamp$Overall.Survival..Months
time <- time/12
timelim <- 15
time[time>timelim] <- timelim
status <- clinic_13q_noamp$OS== 1
status[time>=timelim] <- FALSE
surv <- Surv(time, status) 
fit <- survfit(surv ~  bin_13q, data=clinic_13q_noamp)
test <- survdiff(surv ~  bin_13q, data=clinic_13q_noamp)
p <- pchisq(test$chisq, length(test$n) - 1, lower.tail=FALSE)
col <- c("blue", "black", "red")
mylable <- substitute(italic(p) == MYVALUE, list(MYVALUE = format(p, digits= 3)))
main <- "KM 13q14.2 Metabric HER2pos samples  no gain" 
plot(fit, mark.time = T, xlab="Time (years)", ylab="OS survival", col=col, lwd=6, main=main, cex=1.5)
legend(x="bottomright",col=col,legend= c(paste("deletion (n=", test$n[1],  ")" , sep=""), paste("WT (n=", test$n[2], ")" , sep=""), paste("WT (n=", test$n[3], ")" , sep="")),bty="n",lwd=2)
legend(x= "topright", legend= mylable,bty="n")

#### HER2-  no gain###
clinic_13q_noamp <- clinic_13q[clinic_13q$HER2.Status=="Negative" & !clinic_13q$bin_13q==1,]
time <- clinic_13q_noamp$Overall.Survival..Months
time <- time/12
timelim <- 15
time[time>timelim] <- timelim
status <- clinic_13q_noamp$OS== 1
status[time>=timelim] <- FALSE
surv <- Surv(time, status) 
fit <- survfit(surv ~  bin_13q, data=clinic_13q_noamp)
test <- survdiff(surv ~  bin_13q, data=clinic_13q_noamp)
p <- pchisq(test$chisq, length(test$n) - 1, lower.tail=FALSE)
col <- c("blue", "black", "red")
mylable <- substitute(italic(p) == MYVALUE, list(MYVALUE = format(p, digits= 3)))
main <- "KM 13q14.2 Metabric HER2neg samples  no gain" 
plot(fit, mark.time = T, xlab="Time (years)", ylab="OS survival", col=col, lwd=6, main=main, cex=1.5)
legend(x="bottomright",col=col,legend= c(paste("deletion (n=", test$n[1],  ")" , sep=""), paste("WT (n=", test$n[2], ")" , sep=""), paste("WT (n=", test$n[3], ")" , sep="")),bty="n",lwd=2)
legend(x= "topright", legend= mylable,bty="n")

#### Triple neg###
clinic_13q_noamp <- clinic_13q[clinic_13q$HER2.Status=="Negative" & clinic_13q$ER.Status=="Negative" & clinic_13q$PR.Status=="Negative" & !clinic_13q$bin_13q==1,]
time <- clinic_13q_noamp$Overall.Survival..Months
time <- time/12
timelim <- 15
time[time>timelim] <- timelim
status <- clinic_13q_noamp$OS== 1
status[time>=timelim] <- FALSE
surv <- Surv(time, status) 
fit <- survfit(surv ~  bin_13q, data=clinic_13q_noamp)
test <- survdiff(surv ~  bin_13q, data=clinic_13q_noamp)
p <- pchisq(test$chisq, length(test$n) - 1, lower.tail=FALSE)
col <- c("blue", "black", "red")
mylable <- substitute(italic(p) == MYVALUE, list(MYVALUE = format(p, digits= 3)))
main <- "KM 13q14.2 Metabric TNBC samples  no gain" 
plot(fit, mark.time = T, xlab="Time (years)", ylab="OS survival", col=col, lwd=6, main=main, cex=1.5)
legend(x="bottomright",col=col,legend= c(paste("deletion (n=", test$n[1],  ")" , sep=""), paste("WT (n=", test$n[2], ")" , sep=""), paste("WT (n=", test$n[3], ")" , sep="")),bty="n",lwd=2)
legend(x= "topright", legend= mylable,bty="n")


dev.off()

pdf("pie_charts_metabric.pdf")
pie((prop.table (table (clinic_13q$bin_13q))*100), labels=c("13q deletion", "WT", "13q amplification"), col=c("blue","black", "red"), main="ALL")
All <- prop.table (table (clinic_13q$bin_13q))*100
All_n <- table (clinic_13q$bin_13q)

pie((prop.table (table (clinic_13q$ER.Status, clinic_13q$bin_13q),1)*100)[1,], labels=c("13q deletion", "WT", "13q amplification"), col=c("blue","black", "red"), main="ERneg")
ERneg <- (prop.table (table (clinic_13q$ER.Status, clinic_13q$bin_13q),1)*100)[1,]
ERneg_n <- table (clinic_13q$ER.Status, clinic_13q$bin_13q)[1,]

pie((prop.table (table (clinic_13q$ER.Status, clinic_13q$bin_13q),1)*100)[2,], labels=c("13q deletion", "WT", "13q amplification"), col=c("blue","black", "red"), main="ERpos")
ERpos <- (prop.table (table (clinic_13q$ER.Status, clinic_13q$bin_13q),1)*100)[2,]
ERpos_n <- table (clinic_13q$ER.Status, clinic_13q$bin_13q)[2,]


pie((prop.table (table (clinic_13q$HER2.Status, clinic_13q$bin_13q),1)*100)[1,], labels=c("13q deletion", "WT", "13q amplification"), col=c("blue","black", "red"), main="HER2neg")
Her2neg <-  (prop.table (table (clinic_13q$HER2.Status, clinic_13q$bin_13q),1)*100)[1,]
Her2neg_n <- table (clinic_13q$HER2.Status, clinic_13q$bin_13q)[1,]


pie((prop.table (table (clinic_13q$HER2.Status, clinic_13q$bin_13q),1)*100)[2,], labels=c("13q deletion", "WT", "13q amplification"), col=c("blue","black", "red"), main="HER2neg")
Her2pos <-  (prop.table (table (clinic_13q$HER2.Status, clinic_13q$bin_13q),1)*100)[2,]
Her2pos_n <-  table (clinic_13q$HER2.Status, clinic_13q$bin_13q)[2,]


clinic_13q_noamp <- clinic_13q[clinic_13q$HER2.Status=="Negative" & clinic_13q$ER.Status=="Negative" & clinic_13q$PR.Status=="Negative",]
pie((prop.table (table (clinic_13q_noamp$bin_13q))*100), labels=c("13q deletion", "WT", "13q amplification"), col=c("blue","black", "red"), main="TNBC")
TNBC <- prop.table (table (clinic_13q_noamp$bin_13q))*100
TNBC_n <- table (clinic_13q_noamp$bin_13q)

mer <- rbind (All, All_n, ERpos, ERpos_n,  ERneg, ERneg_n, Her2neg, Her2neg_n, Her2pos, Her2pos_n, TNBC, TNBC_n)
dev.off()
write.table (mer, "pie_charts_underlying_data",  sep = "\t ")


#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
setwd ("C:/Xavier/CNA_13q/files")

CNA <- data.frame(read.table("TCGA_BRCA_CN_thresholded_by_gene.txt", header=T, sep="\t", stringsAsFactors=F))
CNA_13q14.2 <- CNA$Cytoband == "13q14.2"
CNA13q <- CNA[CNA_13q14.2, -c(1:3)]
CNA13q_av <- data.frame (round (colMeans (CNA13q)))
colnames (CNA13q_av) <- "av_13q"
CNA13q_av$Patient.ID <- gsub ("\\.", "-", rownames (CNA13q_av))
CNA13q_av$Patient.ID <- substr (CNA13q_av$Patient.ID, 1,12)

clinic1 <- data.frame(read.table("tcga 13q14.2 data-finalized.txt", header=T, sep="\t", stringsAsFactors=F))
clinic1$pat_id <- substr (clinic1$Patient.ID, 1,12)
clinic2 <- data.frame(read.table("TCGA_BRCA_clinical_data_tumor_only.txt", header=T, sep="\t", stringsAsFactors=F))
clinic <- merge (clinic1, clinic2, by.x='pat_id', by.y='X_PATIENT', all.x=T)

clinic_13q <- merge (clinic, CNA13q_av, by.x="pat_id" , by.y="Patient.ID", all.x=T)
clinic_13q$bin_13q <- gsub ("-2", "-1", clinic_13q$av_13q)
clinic_13q$bin_13q <- gsub ("2", "1", clinic_13q$bin_13q)



pdf ("TCGA_KM_13q.pdf")
library(survival)
#### ALL ###
time <- clinic_13q$OS.time
time <- time/365
timelim <- 15
time[time>timelim] <- timelim
status <- clinic_13q$OS== 1
status[time>=timelim] <- FALSE
surv <- Surv(time, status) 
fit <- survfit(surv ~  bin_13q, data=clinic_13q)
test <- survdiff(surv ~  bin_13q, data=clinic_13q)
p <- pchisq(test$chisq, length(test$n) - 1, lower.tail=FALSE)
col <- c("blue", "black", "red")
mylable <- substitute(italic(p) == MYVALUE, list(MYVALUE = format(p, digits= 3)))
main <- "KM 13q14.2 TCGA All samples" 
plot(fit, mark.time = T, xlab="Time (years)", ylab="OS survival", col=col, lwd=6, main=main, cex=1.5)
legend(x="bottomright",col=col,legend= c(paste("deletion (n=", test$n[1],  ")" , sep=""), paste("WT (n=", test$n[2], ")" , sep=""), paste("gain (n=", test$n[3],  ")" , sep="")),bty="n",lwd=2)
legend(x= "topright", legend= mylable,bty="n")

#### ALL  no gain###
clinic_13q_noamp <- clinic_13q[!clinic_13q$bin_13q==1,]
time <- clinic_13q_noamp$OS.time
time <- time/365
timelim <- 15
time[time>timelim] <- timelim
status <- clinic_13q_noamp$OS== 1
status[time>=timelim] <- FALSE
surv <- Surv(time, status) 
fit <- survfit(surv ~  bin_13q, data=clinic_13q_noamp)
test <- survdiff(surv ~  bin_13q, data=clinic_13q_noamp)
p <- pchisq(test$chisq, length(test$n) - 1, lower.tail=FALSE)
col <- c("blue", "black")
mylable <- substitute(italic(p) == MYVALUE, list(MYVALUE = format(p, digits= 3)))
main <- "KM 13q14.2 TCGA All samples  no gain" 
plot(fit, mark.time = T, xlab="Time (years)", ylab="OS survival", col=col, lwd=6, main=main, cex=1.5)
legend(x="bottomright",col=col,legend= c(paste("deletion (n=", test$n[1],  ")" , sep=""), paste("WT (n=", test$n[2], ")" , sep="")),bty="n",lwd=2)
legend(x= "topright", legend= mylable,bty="n")

#### ER+  no gain###
clinic_13q_noamp <- clinic_13q[clinic_13q$ER_up=="Positive" & !clinic_13q$bin_13q==1,]
time <- clinic_13q_noamp$OS.time
time <- time/365
timelim <- 15
time[time>timelim] <- timelim
status <- clinic_13q_noamp$OS== 1
status[time>=timelim] <- FALSE
surv <- Surv(time, status) 
fit <- survfit(surv ~  bin_13q, data=clinic_13q_noamp)
test <- survdiff(surv ~  bin_13q, data=clinic_13q_noamp)
p <- pchisq(test$chisq, length(test$n) - 1, lower.tail=FALSE)
col <- c("blue", "black", "red")
mylable <- substitute(italic(p) == MYVALUE, list(MYVALUE = format(p, digits= 3)))
main <- "KM 13q14.2 TCGA ERpos samples  no gain" 
plot(fit, mark.time = T, xlab="Time (years)", ylab="OS survival", col=col, lwd=6, main=main, cex=1.5)
legend(x="bottomright",col=col,legend= c(paste("deletion (n=", test$n[1],  ")" , sep=""), paste("WT (n=", test$n[2], ")" , sep="")),bty="n",lwd=2)
legend(x= "topright", legend= mylable,bty="n")

#### ER-  no gain###
clinic_13q_noamp <- clinic_13q[clinic_13q$ER_up=="Negative" & !clinic_13q$bin_13q==1,]
time <- clinic_13q_noamp$OS.time
time <- time/365
timelim <- 15
time[time>timelim] <- timelim
status <- clinic_13q_noamp$OS== 1
status[time>=timelim] <- FALSE
surv <- Surv(time, status) 
fit <- survfit(surv ~  bin_13q, data=clinic_13q_noamp)
test <- survdiff(surv ~  bin_13q, data=clinic_13q_noamp)
p <- pchisq(test$chisq, length(test$n) - 1, lower.tail=FALSE)
col <- c("blue", "black", "red")
mylable <- substitute(italic(p) == MYVALUE, list(MYVALUE = format(p, digits= 3)))
main <- "KM 13q14.2 TCGA ERneg samples  no gain" 
plot(fit, mark.time = T, xlab="Time (years)", ylab="OS survival", col=col, lwd=6, main=main, cex=1.5)
legend(x="bottomright",col=col,legend= c(paste("deletion (n=", test$n[1],  ")" , sep=""), paste("WT (n=", test$n[2], ")" , sep="")),bty="n",lwd=2)
legend(x= "topright", legend= mylable,bty="n")

#### HER2+  no gain###
clinic_13q_noamp <- clinic_13q[clinic_13q$HER2_up=="Positive" & !clinic_13q$bin_13q==1,]
time <- clinic_13q_noamp$OS.time
time <- time/365
timelim <- 15
time[time>timelim] <- timelim
status <- clinic_13q_noamp$OS== 1
status[time>=timelim] <- FALSE
surv <- Surv(time, status) 
fit <- survfit(surv ~  bin_13q, data=clinic_13q_noamp)
test <- survdiff(surv ~  bin_13q, data=clinic_13q_noamp)
p <- pchisq(test$chisq, length(test$n) - 1, lower.tail=FALSE)
col <- c("blue", "black", "red")
mylable <- substitute(italic(p) == MYVALUE, list(MYVALUE = format(p, digits= 3)))
main <- "KM 13q14.2 TCGA HER2pos samples  no gain" 
plot(fit, mark.time = T, xlab="Time (years)", ylab="OS survival", col=col, lwd=6, main=main, cex=1.5)
legend(x="bottomright",col=col,legend= c(paste("deletion (n=", test$n[1],  ")" , sep=""), paste("WT (n=", test$n[2], ")" , sep="")),bty="n",lwd=2)
legend(x= "topright", legend= mylable,bty="n")

#### HER2-  no gain###
clinic_13q_noamp <- clinic_13q[clinic_13q$HER2_up=="Negative" & !clinic_13q$bin_13q==1,]
time <- clinic_13q_noamp$OS.time
time <- time/365
timelim <- 15
time[time>timelim] <- timelim
status <- clinic_13q_noamp$OS== 1
status[time>=timelim] <- FALSE
surv <- Surv(time, status) 
fit <- survfit(surv ~  bin_13q, data=clinic_13q_noamp)
test <- survdiff(surv ~  bin_13q, data=clinic_13q_noamp)
p <- pchisq(test$chisq, length(test$n) - 1, lower.tail=FALSE)
col <- c("blue", "black", "red")
mylable <- substitute(italic(p) == MYVALUE, list(MYVALUE = format(p, digits= 3)))
main <- "KM 13q14.2 TCGA HER2neg samples  no gain" 
plot(fit, mark.time = T, xlab="Time (years)", ylab="OS survival", col=col, lwd=6, main=main, cex=1.5)
legend(x="bottomright",col=col,legend= c(paste("deletion (n=", test$n[1],  ")" , sep=""), paste("WT (n=", test$n[2], ")" , sep="")),bty="n",lwd=2)
legend(x= "topright", legend= mylable,bty="n")
dev.off()

pdf("pie_charts_TCGApdf")
pie((prop.table (table (clinic_13q$ER_up, clinic_13q$bin_13q),1)*100)[1,], labels=c("13q deletion", "WT", "13q amplification"), col=c("blue","black", "red"), main="ERneg")
pie((prop.table (table (clinic_13q$ER_up, clinic_13q$bin_13q),1)*100)[2,], labels=c("13q deletion", "WT", "13q amplification"), col=c("blue","black", "red"), main="ERpos")
pie((prop.table (table (clinic_13q$HER2_up, clinic_13q$bin_13q),1)*100)[1,], labels=c("13q deletion", "WT", "13q amplification"), col=c("blue","black", "red"), main="HER2neg")
pie((prop.table (table (clinic_13q$HER2_up, clinic_13q$bin_13q),1)*100)[2,], labels=c("13q deletion", "WT", "13q amplification"), col=c("blue","black", "red"), main="HER2pos")
pie((prop.table (table (clinic_13q$bin_13q))*100), labels=c("13q deletion", "WT", "13q amplification"), col=c("blue","black", "red"), main="all")
dev.off()


################################################
## Gene expression analysis - MEtabric #########
################################################
setwd ("C:/Xavier/CNA_13q/files")

CNA <- data.frame(read.table("METABRIC_data_cna.txt", header=T, sep="\t", stringsAsFactors=F))
clinic <- data.frame(read.table("METABRIC_data_clinical.txt", header=T, sep="\t", stringsAsFactors=F))
#CNA13q <- CNA [CNA$Cytoband == "13q14.2",]

CNA_13q14.2 <- CNA$Cytoband == "13q14.2"
CNA13q <- CNA[CNA_13q14.2, -c(1:3)]
CNA13q_av <- data.frame (round (colMeans (CNA13q)))

colname) <- "av_13q"


CNA13q_av$Patient.ID <- gsub ("\\.", "-", rownames (CNA13q_av))
clinic_13q <- merge (clinic, CNA13q_av, by.x="Patient.ID" , by.y=0, all.x=TRUE)

table (clinic$Patient.ID==rownames (clinic_13q))
clinic$av_13q <- clinic_13q$av_13q

clinic_13q$bin_13q <- gsub ("-2", "-1", clinic_13q$av_13q)
clinic_13q$bin_13q <- gsub ("2", "1", clinic_13q$bin_13q)
clinic$bin_13q <- clinic_13q$bin_13q

load("C:\\Xavier\\Nanostring_Immundata\\correl_maps\\heatmap_jan_2019\\expression.RData")

dim (Metabric)
clinic$Patient.ID <- gsub ( "-", ".", clinic$Patient.ID )

clinic_sub <- clinic[clinic$Patient.ID%in%colnames (Metabric),]
table (clinic_sub$Patient.ID==colnames (Metabric))
genes <- data.frame(read.table("genes.txt", header=F, sep="\t", stringsAsFactors=F))


wt <-  Metabric [rownames (Metabric) %in% genes$V1, clinic_sub$bin_13q==0]
del <- Metabric [rownames (Metabric) %in% genes$V1, clinic_sub$bin_13q==-1]
amp <- Metabric [rownames (Metabric) %in% genes$V1, clinic_sub$bin_13q==1]



genes_com <- genes$V1 [genes$V1%in%rownames (Metabric)]
res <- data.frame(type="_",WTvDEL=0, WTvAMP=0, DELvAMP=0,  Median_WT=0, Median_DEL=0, Median_AMP=0, stringsAsFactors=F)
k <- 0
for(probe in (genes_com)){
k=0

Median_WT <- median(wt[rownames(wt)%in%probe,])
Median_DEL <- median(del[rownames(del)%in%probe,])
Median_AMP <- median(amp[rownames(amp)%in%probe,])

WTvDEL <- wilcox.test(wt[rownames(wt)%in%probe,], del[rownames(del)%in%probe,], paired = F)[3]
WTvAMP <- wilcox.test(wt[rownames(wt)%in%probe,], amp[rownames(amp)%in%probe,], paired = F)[3]
DELvAMP <- wilcox.test(del[rownames(del)%in%probe,], amp[rownames(amp)%in%probe,], paired = F)[3]

res <- rbind(res,list(type=probe, WTvDEL=WTvDEL, WTvAMP=WTvAMP, DELvAMP=DELvAMP,  Median_WT=Median_WT, Median_DEL=Median_DEL, Median_AMP=Median_AMP))
k=k+1
}


res <- res[-1,]
res$fdr_WTvDEL <- p.adjust(res$WTvDEL, method = 'fdr', n = length(res$WTvDEL))
res$fdr_WTvAMP <- p.adjust(res$WTvAMP, method = 'fdr', n = length(res$WTvAMP))
res$fdr_DELvAMP <- p.adjust(res$DELvAMP, method = 'fdr', n = length(res$DELvAMP))

table (res$fdr_WTvDEL <0.05)
table (res$fdr_WTvAMP <0.05)
table (res$fdr_DELvAMP <0.05)



##########################
### Octobre 2022 #########
##########################

############################
### METABRIC ###############
############################
setwd ("C:/Xavier/CNA_13q/files")

CNA <- data.frame(read.table("METABRIC_data_cna.txt", header=T, sep="\t", stringsAsFactors=F))
clinic <- data.frame(read.table("METABRIC_data_clinical.txt", header=T, sep="\t", stringsAsFactors=F))
CNA_13q14.2 <- CNA$Cytoband == "13q14.2"

CNA13q <- CNA[CNA_13q14.2, -c(1:3)]
CNA13q_av <- data.frame (round (colMeans (CNA13q)))
colnames (CNA13q_av) <- "av_13q"
CNA13q_av$Patient.ID <- gsub ("\\.", "-", rownames (CNA13q_av))

clinic_13q <- merge (clinic, CNA13q_av, by.x="Patient.ID" , by.y="Patient.ID", all.x=T)

clinic_13q$bin_13q <- gsub ("-2", "-1", clinic_13q$av_13q)
clinic_13q$bin_13q <- gsub ("2", "1", clinic_13q$bin_13q)

library (survival)

#### ALL  no gain###
clinic_13q_noamp <- clinic_13q[!clinic_13q$bin_13q==1,]
time <- clinic_13q_noamp$Overall.Survival..Months
time <- time/12
timelim <- 15
time[time>timelim] <- timelim
status <- clinic_13q_noamp$OS== 1
status[time>=timelim] <- FALSE
surv <- Surv(time, status) 
fit <- survfit(surv ~  bin_13q, data=clinic_13q_noamp)
test <- survdiff(surv ~  bin_13q, data=clinic_13q_noamp)
p <- pchisq(test$chisq, length(test$n) - 1, lower.tail=FALSE)
col <- c("blue", "black")
mylable <- substitute(italic(p) == MYVALUE, list(MYVALUE = format(p, digits= 3)))
main <- "KM 13q14.2 Metabric All samples  no gain" 
plot(fit, mark.time = T, xlab="Time (years)", ylab="OS survival", col=col, lwd=6, main=main, cex=1.5)
legend(x="bottomright",col=col,legend= c(paste("deletion (n=", test$n[1],  ")" , sep=""), paste("WT (n=", test$n[2], ")" , sep="")),bty="n",lwd=2)
legend(x= "topright", legend= mylable,bty="n")

#### ER+  no gain###
clinic_13q_noamp <- clinic_13q[clinic_13q$ER.Status=="Positive" & !clinic_13q$bin_13q==1,]
time <- clinic_13q_noamp$Overall.Survival..Months
time <- time/12
timelim <- 15
time[time>timelim] <- timelim
status <- clinic_13q_noamp$OS== 1
status[time>=timelim] <- FALSE
surv <- Surv(time, status) 
fit <- survfit(surv ~  bin_13q, data=clinic_13q_noamp)
test <- survdiff(surv ~  bin_13q, data=clinic_13q_noamp)
p <- pchisq(test$chisq, length(test$n) - 1, lower.tail=FALSE)
col <- c("blue", "black")
mylable <- substitute(italic(p) == MYVALUE, list(MYVALUE = format(p, digits= 3)))
main <- "KM 13q14.2 Metabric ERpos samples  no gain" 
plot(fit, mark.time = T, xlab="Time (years)", ylab="OS survival", col=col, lwd=6, main=main, cex=1.5)
legend(x="bottomright",col=col,legend= c(paste("deletion (n=", test$n[1],  ")" , sep=""), paste("WT (n=", test$n[2], ")" , sep=""), paste("WT (n=", test$n[3], ")" , sep="")),bty="n",lwd=2)
legend(x= "topright", legend= mylable,bty="n")

#### ER-  no gain###
clinic_13q_noamp <- clinic_13q[clinic_13q$ER.Status=="Negative" & !clinic_13q$bin_13q==1,]
time <- clinic_13q_noamp$Overall.Survival..Months
time <- time/12
timelim <- 15
time[time>timelim] <- timelim
status <- clinic_13q_noamp$OS== 1
status[time>=timelim] <- FALSE
surv <- Surv(time, status) 
fit <- survfit(surv ~  bin_13q, data=clinic_13q_noamp)
test <- survdiff(surv ~  bin_13q, data=clinic_13q_noamp)
p <- pchisq(test$chisq, length(test$n) - 1, lower.tail=FALSE)
col <- c("blue", "black")
mylable <- substitute(italic(p) == MYVALUE, list(MYVALUE = format(p, digits= 3)))
main <- "KM 13q14.2 Metabric ERneg samples  no gain" 
plot(fit, mark.time = T, xlab="Time (years)", ylab="OS survival", col=col, lwd=6, main=main, cex=1.5)
legend(x="bottomright",col=col,legend= c(paste("deletion (n=", test$n[1],  ")" , sep=""), paste("WT (n=", test$n[2], ")" , sep=""), paste("WT (n=", test$n[3], ")" , sep="")),bty="n",lwd=2)
legend(x= "topright", legend= mylable,bty="n")

#### HER2+  no gain###
clinic_13q_noamp <- clinic_13q[clinic_13q$HER2.Status=="Positive" & !clinic_13q$bin_13q==1,]
time <- clinic_13q_noamp$Overall.Survival..Months
time <- time/12
timelim <- 15
time[time>timelim] <- timelim
status <- clinic_13q_noamp$OS== 1
status[time>=timelim] <- FALSE
surv <- Surv(time, status) 
fit <- survfit(surv ~  bin_13q, data=clinic_13q_noamp)
test <- survdiff(surv ~  bin_13q, data=clinic_13q_noamp)
p <- pchisq(test$chisq, length(test$n) - 1, lower.tail=FALSE)
col <- c("blue", "black")
mylable <- substitute(italic(p) == MYVALUE, list(MYVALUE = format(p, digits= 3)))
main <- "KM 13q14.2 Metabric HER2pos samples  no gain" 
plot(fit, mark.time = T, xlab="Time (years)", ylab="OS survival", col=col, lwd=6, main=main, cex=1.5)
legend(x="bottomright",col=col,legend= c(paste("deletion (n=", test$n[1],  ")" , sep=""), paste("WT (n=", test$n[2], ")" , sep=""), paste("WT (n=", test$n[3], ")" , sep="")),bty="n",lwd=2)
legend(x= "topright", legend= mylable,bty="n")

#### HER2-  no gain###
clinic_13q_noamp <- clinic_13q[clinic_13q$HER2.Status=="Negative" & !clinic_13q$bin_13q==1,]
time <- clinic_13q_noamp$Overall.Survival..Months
time <- time/12
timelim <- 15
time[time>timelim] <- timelim
status <- clinic_13q_noamp$OS== 1
status[time>=timelim] <- FALSE
surv <- Surv(time, status) 
fit <- survfit(surv ~  bin_13q, data=clinic_13q_noamp)
test <- survdiff(surv ~  bin_13q, data=clinic_13q_noamp)
p <- pchisq(test$chisq, length(test$n) - 1, lower.tail=FALSE)
col <- c("blue", "black")
mylable <- substitute(italic(p) == MYVALUE, list(MYVALUE = format(p, digits= 3)))
main <- "KM 13q14.2 Metabric HER2neg samples  no gain" 
plot(fit, mark.time = T, xlab="Time (years)", ylab="OS survival", col=col, lwd=6, main=main, cex=1.5)
legend(x="bottomright",col=col,legend= c(paste("deletion (n=", test$n[1],  ")" , sep=""), paste("WT (n=", test$n[2], ")" , sep=""), paste("WT (n=", test$n[3], ")" , sep="")),bty="n",lwd=2)
legend(x= "topright", legend= mylable,bty="n")

#### Triple neg###
clinic_13q_noamp <- clinic_13q[clinic_13q$HER2.Status=="Negative" & clinic_13q$ER.Status=="Negative" & clinic_13q$PR.Status=="Negative" & !clinic_13q$bin_13q==1,]
time <- clinic_13q_noamp$Overall.Survival..Months
time <- time/12
timelim <- 15
time[time>timelim] <- timelim
status <- clinic_13q_noamp$OS== 1
status[time>=timelim] <- FALSE
surv <- Surv(time, status) 
fit <- survfit(surv ~  bin_13q, data=clinic_13q_noamp)
test <- survdiff(surv ~  bin_13q, data=clinic_13q_noamp)
p <- pchisq(test$chisq, length(test$n) - 1, lower.tail=FALSE)
col <- c("blue", "black")
mylable <- substitute(italic(p) == MYVALUE, list(MYVALUE = format(p, digits= 3)))
main <- "KM 13q14.2 Metabric TNBC samples  no gain" 
plot(fit, mark.time = T, xlab="Time (years)", ylab="OS survival", col=col, lwd=6, main=main, cex=1.5)
legend(x="bottomright",col=col,legend= c(paste("deletion (n=", test$n[1],  ")" , sep=""), paste("WT (n=", test$n[2], ")" , sep=""), paste("WT (n=", test$n[3], ")" , sep="")),bty="n",lwd=2)
legend(x= "topright", legend= mylable,bty="n")


############################
### TCGA     ###############
############################

setwd ("C:/Xavier/CNA_13q/files")

CNA <- data.frame(read.table("TCGA_BRCA_CN_thresholded_by_gene.txt", header=T, sep="\t", stringsAsFactors=F))
CNA_13q14.2 <- CNA$Cytoband == "13q14.2"
CNA13q <- CNA[CNA_13q14.2, -c(1:3)]
CNA13q_av <- data.frame (round (colMeans (CNA13q)))
colnames (CNA13q_av) <- "av_13q"
CNA13q_av$Patient.ID <- gsub ("\\.", "-", rownames (CNA13q_av))
CNA13q_av$Patient.ID <- substr (CNA13q_av$Patient.ID, 1,12)

clinic1 <- data.frame(read.table("tcga 13q14.2 data-finalized.txt", header=T, sep="\t", stringsAsFactors=F))
clinic1$pat_id <- substr (clinic1$Patient.ID, 1,12)
clinic2 <- data.frame(read.table("TCGA_BRCA_clinical_data_tumor_only.txt", header=T, sep="\t", stringsAsFactors=F))
clinic <- merge (clinic1, clinic2, by.x='pat_id', by.y='X_PATIENT', all.x=T)

clinic_13q <- merge (clinic, CNA13q_av, by.x="pat_id" , by.y="Patient.ID", all.x=T)
clinic_13q$bin_13q <- gsub ("-2", "-1", clinic_13q$av_13q)
clinic_13q$bin_13q <- gsub ("2", "1", clinic_13q$bin_13q)

#### ALL  no gain###
clinic_13q_noamp <- clinic_13q[!clinic_13q$bin_13q==1,]
time <- clinic_13q_noamp$OS.time
time <- time/365
timelim <- 15
time[time>timelim] <- timelim
status <- clinic_13q_noamp$OS== 1
status[time>=timelim] <- FALSE
surv <- Surv(time, status) 
fit <- survfit(surv ~  bin_13q, data=clinic_13q_noamp)
test <- survdiff(surv ~  bin_13q, data=clinic_13q_noamp)
p <- pchisq(test$chisq, length(test$n) - 1, lower.tail=FALSE)
col <- c("blue", "black")
mylable <- substitute(italic(p) == MYVALUE, list(MYVALUE = format(p, digits= 3)))
main <- "KM 13q14.2 TCGA All samples  no gain" 
plot(fit, mark.time = T, xlab="Time (years)", ylab="OS survival", col=col, lwd=6, main=main, cex=1.5)
legend(x="bottomright",col=col,legend= c(paste("deletion (n=", test$n[1],  ")" , sep=""), paste("WT (n=", test$n[2], ")" , sep="")),bty="n",lwd=2)
legend(x= "topright", legend= mylable,bty="n")

#### ER+  no gain###
clinic_13q_noamp <- clinic_13q[clinic_13q$ER_up=="Positive" & !clinic_13q$bin_13q==1,]
time <- clinic_13q_noamp$OS.time
time <- time/365
timelim <- 15
time[time>timelim] <- timelim
status <- clinic_13q_noamp$OS== 1
status[time>=timelim] <- FALSE
surv <- Surv(time, status) 
fit <- survfit(surv ~  bin_13q, data=clinic_13q_noamp)
test <- survdiff(surv ~  bin_13q, data=clinic_13q_noamp)
p <- pchisq(test$chisq, length(test$n) - 1, lower.tail=FALSE)
col <- c("blue", "black", "red")
mylable <- substitute(italic(p) == MYVALUE, list(MYVALUE = format(p, digits= 3)))
main <- "KM 13q14.2 TCGA ERpos samples  no gain" 
plot(fit, mark.time = T, xlab="Time (years)", ylab="OS survival", col=col, lwd=6, main=main, cex=1.5)
legend(x="bottomright",col=col,legend= c(paste("deletion (n=", test$n[1],  ")" , sep=""), paste("WT (n=", test$n[2], ")" , sep="")),bty="n",lwd=2)
legend(x= "topright", legend= mylable,bty="n")

#### ER-  no gain###
clinic_13q_noamp <- clinic_13q[clinic_13q$ER_up=="Negative" & !clinic_13q$bin_13q==1,]
time <- clinic_13q_noamp$OS.time
time <- time/365
timelim <- 15
time[time>timelim] <- timelim
status <- clinic_13q_noamp$OS== 1
status[time>=timelim] <- FALSE
surv <- Surv(time, status) 
fit <- survfit(surv ~  bin_13q, data=clinic_13q_noamp)
test <- survdiff(surv ~  bin_13q, data=clinic_13q_noamp)
p <- pchisq(test$chisq, length(test$n) - 1, lower.tail=FALSE)
col <- c("blue", "black", "red")
mylable <- substitute(italic(p) == MYVALUE, list(MYVALUE = format(p, digits= 3)))
main <- "KM 13q14.2 TCGA ERneg samples  no gain" 
plot(fit, mark.time = T, xlab="Time (years)", ylab="OS survival", col=col, lwd=6, main=main, cex=1.5)
legend(x="bottomright",col=col,legend= c(paste("deletion (n=", test$n[1],  ")" , sep=""), paste("WT (n=", test$n[2], ")" , sep="")),bty="n",lwd=2)
legend(x= "topright", legend= mylable,bty="n")

#### HER2+  no gain###
clinic_13q_noamp <- clinic_13q[clinic_13q$HER2_up=="Positive" & !clinic_13q$bin_13q==1,]
time <- clinic_13q_noamp$OS.time
time <- time/365
timelim <- 15
time[time>timelim] <- timelim
status <- clinic_13q_noamp$OS== 1
status[time>=timelim] <- FALSE
surv <- Surv(time, status) 
fit <- survfit(surv ~  bin_13q, data=clinic_13q_noamp)
test <- survdiff(surv ~  bin_13q, data=clinic_13q_noamp)
p <- pchisq(test$chisq, length(test$n) - 1, lower.tail=FALSE)
col <- c("blue", "black", "red")
mylable <- substitute(italic(p) == MYVALUE, list(MYVALUE = format(p, digits= 3)))
main <- "KM 13q14.2 TCGA HER2pos samples  no gain" 
plot(fit, mark.time = T, xlab="Time (years)", ylab="OS survival", col=col, lwd=6, main=main, cex=1.5)
legend(x="bottomright",col=col,legend= c(paste("deletion (n=", test$n[1],  ")" , sep=""), paste("WT (n=", test$n[2], ")" , sep="")),bty="n",lwd=2)
legend(x= "topright", legend= mylable,bty="n")

#### HER2-  no gain###
clinic_13q_noamp <- clinic_13q[clinic_13q$HER2_up=="Negative" & !clinic_13q$bin_13q==1,]
time <- clinic_13q_noamp$OS.time
time <- time/365
timelim <- 15
time[time>timelim] <- timelim
status <- clinic_13q_noamp$OS== 1
status[time>=timelim] <- FALSE
surv <- Surv(time, status) 
fit <- survfit(surv ~  bin_13q, data=clinic_13q_noamp)
test <- survdiff(surv ~  bin_13q, data=clinic_13q_noamp)
p <- pchisq(test$chisq, length(test$n) - 1, lower.tail=FALSE)
col <- c("blue", "black", "red")
mylable <- substitute(italic(p) == MYVALUE, list(MYVALUE = format(p, digits= 3)))
main <- "KM 13q14.2 TCGA HER2neg samples  no gain" 
plot(fit, mark.time = T, xlab="Time (years)", ylab="OS survival", col=col, lwd=6, main=main, cex=1.5)
legend(x="bottomright",col=col,legend= c(paste("deletion (n=", test$n[1],  ")" , sep=""), paste("WT (n=", test$n[2], ")" , sep="")),bty="n",lwd=2)
legend(x= "topright", legend= mylable,bty="n")
dev.off()