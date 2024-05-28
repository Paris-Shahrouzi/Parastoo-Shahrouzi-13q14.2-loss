df1 <- read.csv("df1.csv", row.names=1, header=T, check.names=F)
df1 <- t(df1)
df1 <- as.data.frame(df1)
df1[, 2:ncol(df1)] <- lapply(2:ncol(df1), function(x) as.numeric(df1[[x]]))
df1_pval <- df1 %>% summarise(across(!group, ~wilcox.test(.x ~ group)$p.value), exact=NULL) %>%
  bind_rows(., p.adjust(., method = 'bonferroni')) %>%
  bind_rows(df1, .) %>%
  mutate(group=replace(group, is.na(group), c('p.values', 'adjusted_p.values')))
df1_scaled <- as.data.frame(scale(df1[,-1]))
df1_scaled <- cbind(df1$group, df1_scaled)
colnames(df1_scaled)[1] <- "group"
new <- sapply(df1_scaled[-1], \(x)diff(tapply(x, df1_scaled$group, median, na.rm=T)))
names(new) <- sub(".loss", "", names(new))
new <- cbind(data.frame(group = "median.diff"), t(new))
df1_scaled <- rbind(df1_scaled, new)
df1_final <- rbind(tail(df1_scaled, n=1), tail(df1_pval, n=2))
temp <- df1_final[,-1]
rownames(temp) <- df1_final[,1]
df1_final <- temp
df1_final <- as.data.frame(t(df1_final))
df1_final <- df1_final %>% arrange(adjusted_p.values, median.diff)
sessionInfo()
# R version 4.4.0 (2024-04-24)
# Platform: aarch64-apple-darwin20
# Running under: macOS Sonoma 14.4.1
#
# Matrix products: default
# BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
# LAPACK: /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.0
#
# locale:
# [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
#
# time zone: Australia/Adelaide
# tzcode source: internal
#
# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base     
#
# other attached packages:
# [1] dplyr_1.1.4
#
# loaded via a namespace (and not attached):
#  [1] utf8_1.2.4       R6_2.5.1         tidyselect_1.2.1 magrittr_2.0.3   glue_1.7.0      
#  [6] tibble_3.2.1     pkgconfig_2.0.3  generics_0.1.3   lifecycle_1.0.4  cli_3.6.2       
# [11] fansi_1.0.6      vctrs_0.6.5      withr_3.0.0      compiler_4.4.0   tools_4.4.0     
# [16] pillar_1.9.0     rlang_1.1.3     
