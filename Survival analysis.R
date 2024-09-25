############################
### METABRIC ###############
############################
setwd ("/Users/parasts/Library/CloudStorage/GoogleDrive-p.shahrouzi@gmail.com/My Drive/Postdoc files/13q14.2 paper/Figure 2/Survival analysis")

CNA <- data.frame(read.table("METABRIC_cna.txt", header=T, sep="\t", stringsAsFactors=F))
clinic <- data.frame(read.table("METABRIC_clinical.txt", header=T, sep="\t", stringsAsFactors=F))
CNA_13q14.2 <- CNA$Cytoband == "13q14.2" #to pick 13q14.2 only from the dataset
clinic <- clinic[ , -ncol(clinic)]


CNA13q <- CNA[CNA_13q14.2, -c(1:3)]
CNA13q_av <- data.frame (round (colMeans (CNA13q))) #to get the average of 13q14.2 CNA for each patinet
colnames (CNA13q_av) <- "av_13q" # to name the column 
CNA13q_av$Patient.ID <- gsub ("\\.", "-", rownames (CNA13q_av)) # to replace periods (.) with dashes (-) in the rownames of your data frame CNA13q_av and store the result in a new column called Patient.ID.


clinic_13q <- merge (clinic, CNA13q_av, by.x="Patient.ID" , by.y="Patient.ID", all.x=T) #to merg CNA and clinic dataset


clinic_13q$bin_13q <- gsub ("-2", "-1", clinic_13q$av_13q) # we are replacing -2 with -1 to have all the losses in one group
clinic_13q$bin_13q <- gsub ("2", "1", clinic_13q$bin_13q) # we are replacing 2 with 1 to have all the gains in one group
clinic_13q$TNBC <- ifelse (clinic_13q$HER2.Status=="Negative" & clinic_13q$ER.Status=="Negative" & clinic_13q$PR.Status=="Negative", 1, 0) # to make a new subtype of Triple negative


library(survival)
#### ALL ###

# Rename the column
names(clinic_13q)[names(clinic_13q) == "Overall.Survival..Months."] <- "time"
clinic_13q$time <- gsub(",", ".", clinic_13q$time)


# Convert the 'time' column to numeric
clinic_13q$time <- as.numeric(as.character(clinic_13q$time))


# Remove rows with NA values in 'time'
clinic_13q <- clinic_13q[!is.na(clinic_13q$time), ]

# Convert time from months to years
clinic_13q$time <- clinic_13q$time / 12

# Cap the time variable at timelim
timelim <- 15
clinic_13q$time[clinic_13q$time > timelim] <- timelim


# Remove rows with NA values in 'bin_13q'
clinic_13q <- clinic_13q[!is.na(clinic_13q$bin_13q), ]

# Remove rows where bin_13q == "1"
clinic_13q <- clinic_13q[clinic_13q$bin_13q != "1", ]

# Ensure 'status' is computed after filtering 'clinic_13q'
clinic_13q <- clinic_13q[!is.na(clinic_13q$bin_13q), ]
clinic_13q <- clinic_13q[clinic_13q$bin_13q != "1", ]

# Recompute status
status <- clinic_13q$Overall.Survival.Status == "1:DECEASED"
status[clinic_13q$time >= timelim] <- FALSE
# Recreate the Surv object
surv <- Surv(clinic_13q$time, status)

# Perform survival analysis
fit <- survfit(surv ~ bin_13q, data = clinic_13q)
test <- survdiff(surv ~ bin_13q, data = clinic_13q)


p <- pchisq(test$chisq, length(test$n) - 1, lower.tail=FALSE) #calculate p value


install.packages("ggplot2")
install.packages("survminer")
library(ggplot2)
library(survminer)
install.packages("grid")

ggsurv <- ggsurvplot(
  fit, 
  data = clinic_13q, 
  risk.table = TRUE,       # Add a risk table below the plot
  pval = TRUE,             # Display the p-value of the log-rank test
  conf.int = TRUE,         # Add confidence intervals
  palette = c("blue", "black"), # Define colors for the groups
  xlab = "Time (years)",   # X-axis label
  ylab = "Overal Survival", # Y-axis label
  legend.title = "13q14.2", # Title for the legend
  legend.labs = c("Loss", "Wild Type"), # Labels for the legend
  ggtheme = theme_minimal() # Use a minimal theme for the plot
)

print(ggsurv$plot)

count_bin_13q_neg1 <- nrow(clinic_13q[clinic_13q$bin_13q == -1, ])
count_bin_13q_0 <- nrow(clinic_13q[clinic_13q$bin_13q == 0, ])
print(paste("Count of bin_13q == -1:", count_bin_13q_neg1))
print(paste("Count of bin_13q == 0:", count_bin_13q_0))

#### ER+  no gain###
clinic_13q_noamp <- clinic_13q[clinic_13q$ER.Status=="Positive" & !clinic_13q$bin_13q==1,]
timelim <- 15
clinic_13q_noamp$time[clinic_13q_noamp$time > timelim] <- timelim
status <- clinic_13q_noamp$Overall.Survival.Status == "1:DECEASED"
status[clinic_13q_noamp$time >= timelim] <- FALSE
surv <- Surv(clinic_13q_noamp$time, status)
fit <- survfit(surv ~ bin_13q, data = clinic_13q_noamp)
test <- survdiff(surv ~  bin_13q, data=clinic_13q_noamp)
p <- pchisq(test$chisq, length(test$n) - 1, lower.tail=FALSE)

ggsurv <- ggsurvplot(
  fit, 
  data = clinic_13q_noamp, 
  risk.table = TRUE,       # Add a risk table below the plot
  pval = TRUE,             # Display the p-value of the log-rank test
  conf.int = TRUE,         # Add confidence intervals
  palette = c("blue", "black"), # Define colors for the groups
  xlab = "Time (years)",   # X-axis label
  ylab = "Overal Survival", # Y-axis label
  legend.title = "13q14.2", # Title for the legend
  legend.labs = c("Loss", "Wild Type"), # Labels for the legend
  ggtheme = theme_minimal() # Use a minimal theme for the plot
)

print(ggsurv$plot)


count_bin_13q_neg1 <- nrow(clinic_13q_noamp[clinic_13q_noamp$bin_13q == -1, ])
count_bin_13q_0 <- nrow(clinic_13q_noamp[clinic_13q_noamp$bin_13q == 0, ])
print(paste("Count of bin_13q == -1:", count_bin_13q_neg1))
print(paste("Count of bin_13q == 0:", count_bin_13q_0))


#### ER-  no gain###
clinic_13q_ER_Negative <- clinic_13q[clinic_13q$ER.Status=="Negative" & !clinic_13q$bin_13q==1,]
timelim <- 15
clinic_13q_ER_Negative$time[clinic_13q_ER_Negative$time > timelim] <- timelim
status <- clinic_13q_ER_Negative$Overall.Survival.Status == "1:DECEASED"
status[clinic_13q_ER_Negative$time >= timelim] <- FALSE
surv <- Surv(clinic_13q_ER_Negative$time, status)
fit <- survfit(surv ~ bin_13q, data = clinic_13q_ER_Negative)
test <- survdiff(surv ~  bin_13q, data=clinic_13q_ER_Negative)
p <- pchisq(test$chisq, length(test$n) - 1, lower.tail=FALSE)

ggsurv <- ggsurvplot(
  fit, 
  data = clinic_13q_ER_Negative, 
  risk.table = TRUE,       # Add a risk table below the plot
  pval = TRUE,             # Display the p-value of the log-rank test
  conf.int = TRUE,         # Add confidence intervals
  palette = c("blue", "black"), # Define colors for the groups
  xlab = "Time (years)",   # X-axis label
  ylab = "Overal Survival", # Y-axis label
  legend.title = "13q14.2", # Title for the legend
  legend.labs = c("Loss", "Wild Type"), # Labels for the legend
  ggtheme = theme_minimal() # Use a minimal theme for the plot
)

print(ggsurv$plot)

#### HER2+  no gain###
clinic_13q_HER2_Positive <- clinic_13q[clinic_13q$HER2.Status=="Positive" & !clinic_13q$bin_13q==1,]
timelim <- 15
clinic_13q_HER2_Positive$time[clinic_13q_HER2_Positive$time > timelim] <- timelim
status <- clinic_13q_HER2_Positive$Overall.Survival.Status == "1:DECEASED"
status[clinic_13q_HER2_Positive$time >= timelim] <- FALSE
surv <- Surv(clinic_13q_HER2_Positive$time, status)
fit <- survfit(surv ~ bin_13q, data = clinic_13q_HER2_Positive)
test <- survdiff(surv ~  bin_13q, data=clinic_13q_HER2_Positive)
p <- pchisq(test$chisq, length(test$n) - 1, lower.tail=FALSE)

ggsurv <- ggsurvplot(
  fit, 
  data = clinic_13q_HER2_Positive, 
  risk.table = TRUE,       # Add a risk table below the plot
  pval = TRUE,             # Display the p-value of the log-rank test
  conf.int = TRUE,         # Add confidence intervals
  palette = c("blue", "black"), # Define colors for the groups
  xlab = "Time (years)",   # X-axis label
  ylab = "Overal Survival", # Y-axis label
  legend.title = "13q14.2", # Title for the legend
  legend.labs = c("Loss", "Wild Type"), # Labels for the legend
  ggtheme = theme_minimal() # Use a minimal theme for the plot
)

print(ggsurv$plot)

count_bin_13q_neg1 <- nrow(clinic_13q_HER2_Positive[clinic_13q_HER2_Positive$bin_13q == -1, ])
count_bin_13q_0 <- nrow(clinic_13q_HER2_Positive[clinic_13q_HER2_Positive$bin_13q == 0, ])
print(paste("Count of bin_13q == -1:", count_bin_13q_neg1))
print(paste("Count of bin_13q == 0:", count_bin_13q_0))

#### HER2-  no gain###
clinic_13q_HER2_Negative <- clinic_13q[clinic_13q$HER2.Status=="Negative" & !clinic_13q$bin_13q==1,]
timelim <- 15
clinic_13q_HER2_Negative$time[clinic_13q_HER2_Negative$time > timelim] <- timelim
status <- clinic_13q_HER2_Negative$Overall.Survival.Status == "1:DECEASED"
status[clinic_13q_HER2_Negative$time >= timelim] <- FALSE
surv <- Surv(clinic_13q_HER2_Negative$time, status)
fit <- survfit(surv ~ bin_13q, data = clinic_13q_HER2_Negative)
test <- survdiff(surv ~  bin_13q, data=clinic_13q_HER2_Negative)
p <- pchisq(test$chisq, length(test$n) - 1, lower.tail=FALSE)

ggsurv <- ggsurvplot(
  fit, 
  data = clinic_13q_HER2_Negative, 
  risk.table = TRUE,       # Add a risk table below the plot
  pval = TRUE,             # Display the p-value of the log-rank test
  conf.int = TRUE,         # Add confidence intervals
  palette = c("blue", "black"), # Define colors for the groups
  xlab = "Time (years)",   # X-axis label
  ylab = "Overal Survival", # Y-axis label
  legend.title = "13q14.2", # Title for the legend
  legend.labs = c("Loss", "Wild Type"), # Labels for the legend
  ggtheme = theme_minimal() # Use a minimal theme for the plot
)

print(ggsurv$plot)

#### Triple negative###
clinic_13q_Triple_Negative <- clinic_13q[clinic_13q$HER2.Status=="Negative" & clinic_13q$ER.Status=="Negative" & clinic_13q$PR.Status=="Negative" & !clinic_13q$bin_13q==1,]
timelim <- 15
clinic_13q_Triple_Negative$time[clinic_13q_Triple_Negative$time > timelim] <- timelim
status <- clinic_13q_Triple_Negative$Overall.Survival.Status == "1:DECEASED"
status[clinic_13q_Triple_Negative$time >= timelim] <- FALSE
surv <- Surv(clinic_13q_Triple_Negative$time, status)
fit <- survfit(surv ~ bin_13q, data = clinic_13q_Triple_Negative)
test <- survdiff(surv ~  bin_13q, data=clinic_13q_Triple_Negative)
p <- pchisq(test$chisq, length(test$n) - 1, lower.tail=FALSE)

ggsurv <- ggsurvplot(
  fit, 
  data = clinic_13q_Triple_Negative, 
  risk.table = TRUE,       # Add a risk table below the plot
  pval = TRUE,             # Display the p-value of the log-rank test
  conf.int = TRUE,         # Add confidence intervals
  palette = c("blue", "black"), # Define colors for the groups
  xlab = "Time (years)",   # X-axis label
  ylab = "Overal Survival", # Y-axis label
  legend.title = "13q14.2", # Title for the legend
  legend.labs = c("Loss", "Wild Type"), # Labels for the legend
  ggtheme = theme_minimal() # Use a minimal theme for the plot
)

print(ggsurv$plot)


#### PAM50_Basal###
clinic_13q_Basal <- clinic_13q[clinic_13q $Pam50...Claudin.low.subtype=="Basal"  & !clinic_13q$bin_13q==1,]
timelim <- 15
clinic_13q_Basal$time[clinic_13q_Basal$time > timelim] <- timelim
status <- clinic_13q_Basal$Overall.Survival.Status == "1:DECEASED"
status[clinic_13q_Basal$time >= timelim] <- FALSE
surv <- Surv(clinic_13q_Basal$time, status)
fit <- survfit(surv ~ bin_13q, data = clinic_13q_Basal)
test <- survdiff(surv ~  bin_13q, data=clinic_13q_Basal)
p <- pchisq(test$chisq, length(test$n) - 1, lower.tail=FALSE)

ggsurv <- ggsurvplot(
  fit, 
  data = clinic_13q_Basal, 
  risk.table = TRUE,       # Add a risk table below the plot
  pval = TRUE,             # Display the p-value of the log-rank test
  conf.int = TRUE,         # Add confidence intervals
  palette = c("blue", "black"), # Define colors for the groups
  xlab = "Time (years)",   # X-axis label
  ylab = "Overal Survival", # Y-axis label
  legend.title = "13q14.2", # Title for the legend
  legend.labs = c("Loss", "Wild Type"), # Labels for the legend
  ggtheme = theme_minimal() # Use a minimal theme for the plot
)

print(ggsurv$plot)

#### PAM50_Luminal AA###
clinic_13q_LuminalA <- clinic_13q[clinic_13q $Pam50...Claudin.low.subtype=="LumA"  & !clinic_13q$bin_13q==1,]
timelim <- 15
clinic_13q_LuminalA$time[clinic_13q_LuminalA$time > timelim] <- timelim
status <- clinic_13q_LuminalA$Overall.Survival.Status == "1:DECEASED"
status[clinic_13q_LuminalA$time >= timelim] <- FALSE
surv <- Surv(clinic_13q_LuminalA$time, status)
fit <- survfit(surv ~ bin_13q, data = clinic_13q_LuminalA)
test <- survdiff(surv ~  bin_13q, data=clinic_13q_LuminalA)
p <- pchisq(test$chisq, length(test$n) - 1, lower.tail=FALSE)

ggsurv <- ggsurvplot(
  fit, 
  data = clinic_13q_LuminalA, 
  risk.table = TRUE,       # Add a risk table below the plot
  pval = TRUE,             # Display the p-value of the log-rank test
  conf.int = TRUE,         # Add confidence intervals
  palette = c("blue", "black"), # Define colors for the groups
  xlab = "Time (years)",   # X-axis label
  ylab = "Overal Survival", # Y-axis label
  legend.title = "13q14.2", # Title for the legend
  legend.labs = c("Loss", "Wild Type"), # Labels for the legend
  ggtheme = theme_minimal() # Use a minimal theme for the plot
)

print(ggsurv$plot)

#### PAM50_Luminal B ###
clinic_13q_LuminalB <- clinic_13q[clinic_13q $Pam50...Claudin.low.subtype=="LumB"  & !clinic_13q$bin_13q==1,]
timelim <- 15
clinic_13q_LuminalB$time[clinic_13q_LuminalB$time > timelim] <- timelim
status <- clinic_13q_LuminalB$Overall.Survival.Status == "1:DECEASED"
status[clinic_13q_LuminalB$time >= timelim] <- FALSE
surv <- Surv(clinic_13q_LuminalB$time, status)
fit <- survfit(surv ~ bin_13q, data = clinic_13q_LuminalB)
test <- survdiff(surv ~  bin_13q, data=clinic_13q_LuminalB)
p <- pchisq(test$chisq, length(test$n) - 1, lower.tail=FALSE)

ggsurv <- ggsurvplot(
  fit, 
  data = clinic_13q_LuminalB, 
  risk.table = TRUE,       # Add a risk table below the plot
  pval = TRUE,             # Display the p-value of the log-rank test
  conf.int = TRUE,         # Add confidence intervals
  palette = c("blue", "black"), # Define colors for the groups
  xlab = "Time (years)",   # X-axis label
  ylab = "Overal Survival", # Y-axis label
  legend.title = "13q14.2", # Title for the legend
  legend.labs = c("Loss", "Wild Type"), # Labels for the legend
  ggtheme = theme_minimal() # Use a minimal theme for the plot
)

print(ggsurv$plot)

#### PAM50_Her2 enriched ###
clinic_13q_HER2_Enriched <- clinic_13q[clinic_13q $Pam50...Claudin.low.subtype=="Her2"  & !clinic_13q$bin_13q==1,]
timelim <- 15
clinic_13q_HER2_Enriched$time[clinic_13q_HER2_Enriched$time > timelim] <- timelim
status <- clinic_13q_HER2_Enriched$Overall.Survival.Status == "1:DECEASED"
status[clinic_13q_HER2_Enriched$time >= timelim] <- FALSE
surv <- Surv(clinic_13q_HER2_Enriched$time, status)
fit <- survfit(surv ~ bin_13q, data = clinic_13q_HER2_Enriched)
test <- survdiff(surv ~  bin_13q, data=clinic_13q_HER2_Enriched)
p <- pchisq(test$chisq, length(test$n) - 1, lower.tail=FALSE)

ggsurv <- ggsurvplot(
  fit, 
  data = clinic_13q_HER2_Enriched, 
  risk.table = TRUE,       # Add a risk table below the plot
  pval = TRUE,             # Display the p-value of the log-rank test
  conf.int = TRUE,         # Add confidence intervals
  palette = c("blue", "black"), # Define colors for the groups
  xlab = "Time (years)",   # X-axis label
  ylab = "Overal Survival", # Y-axis label
  legend.title = "13q14.2", # Title for the legend
  legend.labs = c("Loss", "Wild Type"), # Labels for the legend
  ggtheme = theme_minimal() # Use a minimal theme for the plot
)

print(ggsurv$plot)

#### ER POSITIVE/HER2 NEGATIVE ###
clinic_13q_ER_pos_Her2_neg <- clinic_13q[clinic_13q $ER.Status=="Positive"  & clinic_13q $HER2.Status=="Negative"& !clinic_13q$bin_13q==1,]
timelim <- 15
clinic_13q_ER_pos_Her2_neg$time[clinic_13q_ER_pos_Her2_neg$time > timelim] <- timelim
status <- clinic_13q_ER_pos_Her2_neg$Overall.Survival.Status == "1:DECEASED"
status[clinic_13q_ER_pos_Her2_neg$time >= timelim] <- FALSE
surv <- Surv(clinic_13q_ER_pos_Her2_neg$time, status)
fit <- survfit(surv ~ bin_13q, data = clinic_13q_ER_pos_Her2_neg)
test <- survdiff(surv ~  bin_13q, data=clinic_13q_ER_pos_Her2_neg)
p <- pchisq(test$chisq, length(test$n) - 1, lower.tail=FALSE)

ggsurv <- ggsurvplot(
  fit, 
  data = clinic_13q_ER_pos_Her2_neg, 
  risk.table = TRUE,       # Add a risk table below the plot
  pval = TRUE,             # Display the p-value of the log-rank test
  conf.int = TRUE,         # Add confidence intervals
  palette = c("blue", "black"), # Define colors for the groups
  xlab = "Time (years)",   # X-axis label
  ylab = "Overal Survival", # Y-axis label
  legend.title = "13q14.2", # Title for the legend
  legend.labs = c("Loss", "Wild Type"), # Labels for the legend
  ggtheme = theme_minimal() # Use a minimal theme for the plot
)

print(ggsurv$plot)

#To count 13q14.2 in this subtype

count_bin_13q_neg1 <- nrow(clinic_13q_ER_pos_Her2_neg[clinic_13q_ER_pos_Her2_neg$bin_13q == -1, ])
count_bin_13q_0 <- nrow(clinic_13q_ER_pos_Her2_neg[clinic_13q_ER_pos_Her2_neg$bin_13q == 0, ])
print(paste("Count of bin_13q == -1:", count_bin_13q_neg1))
print(paste("Count of bin_13q == 0:", count_bin_13q_0))

#### ER Negative/HER2 NEGATIVE ###
clinic_13q_ER_neg_Her2_neg <- clinic_13q[clinic_13q $ER.Status=="Negative"  & clinic_13q $HER2.Status=="Negative"& !clinic_13q$bin_13q==1,]
timelim <- 15
clinic_13q_ER_neg_Her2_neg$time[clinic_13q_ER_neg_Her2_neg$time > timelim] <- timelim
status <- clinic_13q_ER_neg_Her2_neg$Overall.Survival.Status == "1:DECEASED"
status[clinic_13q_ER_neg_Her2_neg$time >= timelim] <- FALSE
surv <- Surv(clinic_13q_ER_neg_Her2_neg$time, status)
fit <- survfit(surv ~ bin_13q, data = clinic_13q_ER_neg_Her2_neg)
test <- survdiff(surv ~  bin_13q, data=clinic_13q_ER_neg_Her2_neg)
p <- pchisq(test$chisq, length(test$n) - 1, lower.tail=FALSE)

ggsurv <- ggsurvplot(
  fit, 
  data = clinic_13q_ER_neg_Her2_neg, 
  risk.table = TRUE,       # Add a risk table below the plot
  pval = TRUE,             # Display the p-value of the log-rank test
  conf.int = TRUE,         # Add confidence intervals
  palette = c("blue", "black"), # Define colors for the groups
  xlab = "Time (years)",   # X-axis label
  ylab = "Overal Survival", # Y-axis label
  legend.title = "13q14.2", # Title for the legend
  legend.labs = c("Loss", "Wild Type"), # Labels for the legend
  ggtheme = theme_minimal() # Use a minimal theme for the plot
)

print(ggsurv$plot)

count_bin_13q_neg1 <- nrow(clinic_13q_ER_neg_Her2_neg[clinic_13q_ER_neg_Her2_neg$bin_13q == -1, ])
count_bin_13q_0 <- nrow(clinic_13q_ER_neg_Her2_neg[clinic_13q_ER_neg_Her2_neg$bin_13q == 0, ])
print(paste("Count of bin_13q == -1:", count_bin_13q_neg1))
print(paste("Count of bin_13q == 0:", count_bin_13q_0))

#############################################################################################################################
############################################################       TCGA      #######################################################
#############################################################################################################################

CNA_TCGA <- data.frame(read.table("TCGA_CNA.txt", header=T, sep="\t", stringsAsFactors=F))
CNA_13q14.2_TCGA <- CNA_TCGA$Cytoband == "13q14.2"
CNA13q <- CNA_TCGA[CNA_13q14.2_TCGA, -c(1:3)]
CNA13q_av_TCGA <- data.frame (round (colMeans (CNA13q)))
colnames (CNA13q_av_TCGA) <- "av_13q"
CNA13q_av_TCGA$Patient.ID <- gsub ("\\.", "-", rownames (CNA13q_av_TCGA))
CNA13q_av_TCGA$Patient.ID <- substr (CNA13q_av_TCGA$Patient.ID, 1,12)

TCGA_clinical <- data.frame(read.table("TCGA_CLINICAL.txt", header=T, sep="\t", stringsAsFactors=F))
TCGA_clinical$sampleID <- substr (TCGA_clinical$sampleID, 1,12)
names(TCGA_clinical)[names(TCGA_clinical) == "sampleID"] <- "Patient.ID"


TCGA_merged <- merge(CNA13q_av_TCGA, TCGA_clinical, by = "Patient.ID")

TCGA_merged$bin_13q <- gsub ("-2", "-1", TCGA_merged$av_13q)
TCGA_merged$bin_13q <- gsub ("2", "1", TCGA_merged$bin_13q)

#### all TCGA samples ####

clinic_13q_noamp <- TCGA_merged[!TCGA_merged$bin_13q==1,]
time <- clinic_13q_noamp $OS_Time_nature2012
time <- time/365
timelim <- 15
time[time>timelim] <- timelim
status <- clinic_13q_noamp $OS_event_nature2012== 1
status[time>=timelim] <- FALSE
surv <- Surv(time, status) 
fit <- survfit(surv ~  bin_13q, data=clinic_13q_noamp)
test <- survdiff(surv ~  bin_13q, data=clinic_13q_noamp)
p <- pchisq(test$chisq, length(test$n) - 1, lower.tail=FALSE)

ggsurv <- ggsurvplot(
  fit, 
  data = clinic_13q_noamp, 
  risk.table = TRUE,       
  pval = TRUE,             
  conf.int = TRUE,         
  palette = c("blue", "black"), 
  xlab = "Time (years)",   
  ylab = "Overal Survival", 
  legend.title = "13q14.2",
  legend.labs = c("Loss", "Wild Type"), 
  ggtheme = theme_minimal() 
)

print(ggsurv$plot)

#To count 13q14.2 in this subtype

count_bin_13q_neg1 <- nrow(clinic_13q_noamp[clinic_13q_noamp$bin_13q == -1, ])
count_bin_13q_0 <- nrow(clinic_13q_noamp[clinic_13q_noamp$bin_13q == 0, ])
print(paste("Count of bin_13q == -1:", count_bin_13q_neg1))
print(paste("Count of bin_13q == 0:", count_bin_13q_0))

#### TCGA-ER+ no gain ####

clinic_13q_ER_positive <- TCGA_merged[TCGA_merged $ER_Status_nature2012=="Positive" & !TCGA_merged$bin_13q==1,]
time <- clinic_13q_ER_positive $OS_Time_nature2012
time <- time/365
timelim <- 15
time[time>timelim] <- timelim
status <- clinic_13q_ER_positive $OS_event_nature2012== 1
status[time>=timelim] <- FALSE
surv <- Surv(time, status) 
fit <- survfit(surv ~  bin_13q, data=clinic_13q_ER_positive)
test <- survdiff(surv ~  bin_13q, data=clinic_13q_ER_positive)
p <- pchisq(test$chisq, length(test$n) - 1, lower.tail=FALSE)

ggsurv <- ggsurvplot(
  fit, 
  data = clinic_13q_ER_positive, 
  risk.table = TRUE,       
  pval = TRUE,             
  conf.int = TRUE,         
  palette = c("blue", "black"), 
  xlab = "Time (years)",   
  ylab = "Overal Survival", 
  legend.title = "13q14.2",
  legend.labs = c("Loss", "Wild Type"), 
  ggtheme = theme_minimal() 
)

print(ggsurv$plot)



#### TCGA-ER- no gain ####

clinic_13q_ER_negative <- TCGA_merged[TCGA_merged $ER_Status_nature2012=="Negative" & !TCGA_merged$bin_13q==1,]
time <- clinic_13q_ER_negative $OS_Time_nature2012
time <- time/365
timelim <- 15
time[time>timelim] <- timelim
status <- clinic_13q_ER_negative $OS_event_nature2012== 1
status[time>=timelim] <- FALSE
surv <- Surv(time, status) 
fit <- survfit(surv ~  bin_13q, data=clinic_13q_ER_negative)
test <- survdiff(surv ~  bin_13q, data=clinic_13q_ER_negative)
p <- pchisq(test$chisq, length(test$n) - 1, lower.tail=FALSE)

ggsurv <- ggsurvplot(
  fit, 
  data = clinic_13q_ER_negative, 
  risk.table = TRUE,       
  pval = TRUE,             
  conf.int = TRUE,         
  palette = c("blue", "black"), 
  xlab = "Time (years)",   
  ylab = "Overal Survival", 
  legend.title = "13q14.2",
  legend.labs = c("Loss", "Wild Type"), 
  ggtheme = theme_minimal() 
)

print(ggsurv$plot)

#### TCGA-HER2- no gain ####

clinic_13q_HER_negative <- TCGA_merged[TCGA_merged $HER2_Final_Status_nature2012=="Negative" & !TCGA_merged$bin_13q==1,]
time <- clinic_13q_HER_negative $OS_Time_nature2012
time <- time/365
timelim <- 15
time[time>timelim] <- timelim
status <- clinic_13q_HER_negative $OS_event_nature2012== 1
status[time>=timelim] <- FALSE
surv <- Surv(time, status) 
fit <- survfit(surv ~  bin_13q, data=clinic_13q_HER_negative)
test <- survdiff(surv ~  bin_13q, data=clinic_13q_HER_negative)
p <- pchisq(test$chisq, length(test$n) - 1, lower.tail=FALSE)

ggsurv <- ggsurvplot(
  fit, 
  data = clinic_13q_HER_negative, 
  risk.table = TRUE,       
  pval = TRUE,             
  conf.int = TRUE,         
  palette = c("blue", "black"), 
  xlab = "Time (years)",   
  ylab = "Overal Survival", 
  legend.title = "13q14.2",
  legend.labs = c("Loss", "Wild Type"), 
  ggtheme = theme_minimal() 
)

print(ggsurv$plot)

#### TCGA-HER2+ no gain ####

clinic_13q_HER_pos <- TCGA_merged[TCGA_merged $HER2_Final_Status_nature2012=="Positive" & !TCGA_merged$bin_13q==1,]
time <- clinic_13q_HER_pos $OS_Time_nature2012
time <- time/365
timelim <- 15
time[time>timelim] <- timelim
status <- clinic_13q_HER_pos $OS_event_nature2012== 1
status[time>=timelim] <- FALSE
surv <- Surv(time, status) 
fit <- survfit(surv ~  bin_13q, data=clinic_13q_HER_pos)
test <- survdiff(surv ~  bin_13q, data=clinic_13q_HER_pos)
p <- pchisq(test$chisq, length(test$n) - 1, lower.tail=FALSE)

ggsurv <- ggsurvplot(
  fit, 
  data = clinic_13q_HER_pos, 
  risk.table = TRUE,       
  pval = TRUE,             
  conf.int = TRUE,         
  palette = c("blue", "black"), 
  xlab = "Time (years)",   
  ylab = "Overal Survival", 
  legend.title = "13q14.2",
  legend.labs = c("Loss", "Wild Type"), 
  ggtheme = theme_minimal() 
)

print(ggsurv$plot)

#To count 13q14.2 in this subtype

count_bin_13q_neg1 <- nrow(clinic_13q_HER_pos[clinic_13q_HER_pos$bin_13q == -1, ])
count_bin_13q_0 <- nrow(clinic_13q_HER_pos[clinic_13q_HER_pos$bin_13q == 0, ])
print(paste("Count of bin_13q == -1:", count_bin_13q_neg1))
print(paste("Count of bin_13q == 0:", count_bin_13q_0))


#### TCGA-Triple negative###

TCGA_Triple_Negative <- TCGA_merged[TCGA_merged $HER2_Final_Status_nature2012=="Negative" & TCGA_merged $ER_Status_nature2012=="Negative" & TCGA_merged $PR_Status_nature2012=="Negative" & !TCGA_merged$bin_13q==1,]
time <- TCGA_Triple_Negative $OS_Time_nature2012
time <- time/365
timelim <- 15
time[time>timelim] <- timelim
status <- TCGA_Triple_Negative $OS_event_nature2012== 1
status[time>=timelim] <- FALSE
surv <- Surv(time, status) 
fit <- survfit(surv ~  bin_13q, data=TCGA_Triple_Negative)
test <- survdiff(surv ~  bin_13q, data=TCGA_Triple_Negative)
p <- pchisq(test$chisq, length(test$n) - 1, lower.tail=FALSE)


ggsurv <- ggsurvplot(
  fit, 
  data = TCGA_Triple_Negative, 
  risk.table = TRUE,       # Add a risk table below the plot
  pval = TRUE,             # Display the p-value of the log-rank test
  conf.int = TRUE,         # Add confidence intervals
  palette = c("blue", "black"), # Define colors for the groups
  xlab = "Time (years)",   # X-axis label
  ylab = "Overal Survival", # Y-axis label
  legend.title = "13q14.2", # Title for the legend
  legend.labs = c("Loss", "Wild Type"), # Labels for the legend
  ggtheme = theme_minimal() # Use a minimal theme for the plot
)

print(ggsurv$plot)

#To count 13q14.2 in this subtype

count_bin_13q_neg1 <- nrow(TCGA_Triple_Negative[TCGA_Triple_Negative$bin_13q == -1, ])
count_bin_13q_0 <- nrow(TCGA_Triple_Negative[TCGA_Triple_Negative$bin_13q == 0, ])
print(paste("Count of bin_13q == -1:", count_bin_13q_neg1))
print(paste("Count of bin_13q == 0:", count_bin_13q_0))


#### PAM50_Basal###
TCGA_AM50_Basal <- TCGA_merged[TCGA_merged $PAM50Call_RNAseq=="Basal" & !TCGA_merged$bin_13q==1,]
time <- TCGA_AM50_Basal $OS_Time_nature2012
time <- time/365
timelim <- 15
time[time>timelim] <- timelim
status <- TCGA_AM50_Basal $OS_event_nature2012== 1
status[time>=timelim] <- FALSE
surv <- Surv(time, status) 
fit <- survfit(surv ~  bin_13q, data=TCGA_AM50_Basal)
test <- survdiff(surv ~  bin_13q, data=TCGA_AM50_Basal)
p <- pchisq(test$chisq, length(test$n) - 1, lower.tail=FALSE)

ggsurv <- ggsurvplot(
  fit, 
  data = TCGA_AM50_Basal, 
  risk.table = TRUE,       
  pval = TRUE,             
  conf.int = TRUE,         
  palette = c("blue", "black"), 
  xlab = "Time (years)",   
  ylab = "Overal Survival", 
  legend.title = "13q14.2",
  legend.labs = c("Loss", "Wild Type"), 
  ggtheme = theme_minimal() 
)

print(ggsurv$plot)

#### PAM50_Luminal-A ###

TCGA_PAM50_LumA <- TCGA_merged[TCGA_merged $PAM50Call_RNAseq=="LumA" & !TCGA_merged$bin_13q==1,]
time <- TCGA_PAM50_LumA $OS_Time_nature2012
time <- time/365
timelim <- 15
time[time>timelim] <- timelim
status <- TCGA_PAM50_LumA $OS_event_nature2012== 1
status[time>=timelim] <- FALSE
surv <- Surv(time, status) 
fit <- survfit(surv ~  bin_13q, data=TCGA_PAM50_LumA)
test <- survdiff(surv ~  bin_13q, data=TCGA_PAM50_LumA)
p <- pchisq(test$chisq, length(test$n) - 1, lower.tail=FALSE)

ggsurv <- ggsurvplot(
  fit, 
  data = TCGA_PAM50_LumA, 
  risk.table = TRUE,       
  pval = TRUE,             
  conf.int = TRUE,         
  palette = c("blue", "black"), 
  xlab = "Time (years)",   
  ylab = "Overal Survival", 
  legend.title = "13q14.2",
  legend.labs = c("Loss", "Wild Type"), 
  ggtheme = theme_minimal() 
)

print(ggsurv$plot)

#### TCGA PAM50_Luminal B ###

TCGA_PAM50_LumB <- TCGA_merged[TCGA_merged $PAM50Call_RNAseq=="LumB" & !TCGA_merged$bin_13q==1,]
time <- TCGA_PAM50_LumB $OS_Time_nature2012
time <- time/365
timelim <- 15
time[time>timelim] <- timelim
status <- TCGA_PAM50_LumB $OS_event_nature2012== 1
status[time>=timelim] <- FALSE
surv <- Surv(time, status) 
fit <- survfit(surv ~  bin_13q, data=TCGA_PAM50_LumB)
test <- survdiff(surv ~  bin_13q, data=TCGA_PAM50_LumB)
p <- pchisq(test$chisq, length(test$n) - 1, lower.tail=FALSE)

ggsurv <- ggsurvplot(
  fit, 
  data = TCGA_PAM50_LumB, 
  risk.table = TRUE,       
  pval = TRUE,             
  conf.int = TRUE,         
  palette = c("blue", "black"), 
  xlab = "Time (years)",   
  ylab = "Overal Survival", 
  legend.title = "13q14.2",
  legend.labs = c("Loss", "Wild Type"), 
  ggtheme = theme_minimal() 
)

print(ggsurv$plot)

#### PAM50_Her2 enriched ###

TCGA_PAM50_HER2 <- TCGA_merged[TCGA_merged $PAM50Call_RNAseq=="Her2" & !TCGA_merged$bin_13q==1,]
time <- TCGA_PAM50_HER2 $OS_Time_nature2012
time <- time/365
timelim <- 15
time[time>timelim] <- timelim
status <- TCGA_PAM50_HER2 $OS_event_nature2012== 1
status[time>=timelim] <- FALSE
surv <- Surv(time, status) 
fit <- survfit(surv ~  bin_13q, data=TCGA_PAM50_HER2)
test <- survdiff(surv ~  bin_13q, data=TCGA_PAM50_HER2)
p <- pchisq(test$chisq, length(test$n) - 1, lower.tail=FALSE)

ggsurv <- ggsurvplot(
  fit, 
  data = TCGA_PAM50_HER2, 
  risk.table = TRUE,       
  pval = TRUE,             
  conf.int = TRUE,         
  palette = c("blue", "black"), 
  xlab = "Time (years)",   
  ylab = "Overal Survival", 
  legend.title = "13q14.2",
  legend.labs = c("Loss", "Wild Type"), 
  ggtheme = theme_minimal() 
)

print(ggsurv$plot)

#### ER POSITIVE/HER2 NEGATIVE ###

TCGA_ERpos_HER2neg <- TCGA_merged[TCGA_merged $ER_Status_nature2012=="Positive" & TCGA_merged $HER2_Final_Status_nature2012=="Negative" & !TCGA_merged$bin_13q==1,]
time <- TCGA_ERpos_HER2neg $OS_Time_nature2012
time <- time/365
timelim <- 15
time[time>timelim] <- timelim
status <- TCGA_ERpos_HER2neg $OS_event_nature2012== 1
status[time>=timelim] <- FALSE
surv <- Surv(time, status) 
fit <- survfit(surv ~  bin_13q, data=TCGA_ERpos_HER2neg)
test <- survdiff(surv ~  bin_13q, data=TCGA_ERpos_HER2neg)
p <- pchisq(test$chisq, length(test$n) - 1, lower.tail=FALSE)

ggsurv <- ggsurvplot(
  fit, 
  data = TCGA_ERpos_HER2neg, 
  risk.table = TRUE,       
  pval = TRUE,             
  conf.int = TRUE,         
  palette = c("blue", "black"), 
  xlab = "Time (years)",   
  ylab = "Overal Survival", 
  legend.title = "13q14.2",
  legend.labs = c("Loss", "Wild Type"), 
  ggtheme = theme_minimal() 
)

print(ggsurv$plot)

#To count 13q14.2 in this subtype

count_bin_13q_neg1 <- nrow(TCGA_ERpos_HER2neg[TCGA_ERpos_HER2neg$bin_13q == -1, ])
count_bin_13q_0 <- nrow(TCGA_ERpos_HER2neg[TCGA_ERpos_HER2neg$bin_13q == 0, ])
print(paste("Count of bin_13q == -1:", count_bin_13q_neg1))
print(paste("Count of bin_13q == 0:", count_bin_13q_0))



#### ER Negative/HER2 NEGATIVE ###

TCGA_ERneg_HER2neg <- TCGA_merged[TCGA_merged $ER_Status_nature2012=="Negative" & TCGA_merged $HER2_Final_Status_nature2012=="Negative" & !TCGA_merged$bin_13q==1,]
time <- TCGA_ERneg_HER2neg $OS_Time_nature2012
time <- time/365
timelim <- 15
time[time>timelim] <- timelim
status <- TCGA_ERneg_HER2neg $OS_event_nature2012== 1
status[time>=timelim] <- FALSE
surv <- Surv(time, status) 
fit <- survfit(surv ~  bin_13q, data=TCGA_ERneg_HER2neg)
test <- survdiff(surv ~  bin_13q, data=TCGA_ERneg_HER2neg)
p <- pchisq(test$chisq, length(test$n) - 1, lower.tail=FALSE)

ggsurv <- ggsurvplot(
  fit, 
  data = TCGA_ERneg_HER2neg, 
  risk.table = TRUE,       
  pval = TRUE,             
  conf.int = TRUE,         
  palette = c("blue", "black"), 
  xlab = "Time (years)",   
  ylab = "Overal Survival", 
  legend.title = "13q14.2",
  legend.labs = c("Loss", "Wild Type"), 
  ggtheme = theme_minimal() 
)

print(ggsurv$plot)

#To count 13q14.2 in this subtype

count_bin_13q_neg1 <- nrow(TCGA_ERneg_HER2neg[TCGA_ERneg_HER2neg$bin_13q == -1, ])
count_bin_13q_0 <- nrow(TCGA_ERneg_HER2neg[TCGA_ERneg_HER2neg$bin_13q == 0, ])
print(paste("Count of bin_13q == -1:", count_bin_13q_neg1))
print(paste("Count of bin_13q == 0:", count_bin_13q_0))
