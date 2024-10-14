# R-script - Stacked bar graph 13q14.2 copy numberin TCGA PAN CANCER
library(ggplot2)
library(dplyr)
df <- read.csv("cn_13q14.2_TCGA_pancan.csv")
df_ordered <- df %>%
  group_by(type) %>%
  summarise(combined_fraction=sum(cn_13q14.2 %in% c("-2", "-1")) / n()) %>%
  arrange(desc(combined_fraction))
df <- df %>%
  mutate(type=factor(type, levels=df_ordered$type))
df_text <- df %>%
  group_by(type, cn_13q14.2) %>%
  summarise(count=n()) %>%
  mutate(percentage=count/sum(count)*100) %>%
  group_by(type) %>%
  mutate(cum_percentage=cumsum(percentage)) %>%
  left_join(df_ordered, by="type") %>%
  filter(cn_13q14.2=="-1")
df_text <- df_text %>%
  mutate(text_y_position=cum_percentage)
plot <- ggplot(df %>%
                 group_by(type, cn_13q14.2) %>%
                 summarise(count=n()) %>%
                 mutate(percentage=count/sum(count)*100), 
               aes(x=type, y=percentage, fill=factor(cn_13q14.2, levels=rev(c("-2", "-1", "0", "1", "2", "ambiguous"))))) + 
  geom_bar(stat="identity", position="stack") + 
  geom_text(
    data=df_text,
    aes(x=type, y=text_y_position, label=sprintf("%d%%", round(combined_fraction*100))),
    color="black",
    size=3,
    angle=90,
    hjust=0.5
  ) + 
  scale_y_continuous(
    breaks=c(0, 25, 50, 75, 100), 
    labels=c(0, 25, 50, 75, 100), 
    expand=c(0, 0)
  ) + 
  labs(x=NULL, y="Fraction of patients (%)", fill="13q14.2 cn") + 
  theme_minimal() +
  theme(
    axis.text.x=element_text(angle=90, hjust=1, color="black"), 
    axis.text.y=element_text(angle=90, hjust=0.5, color="black"), 
    axis.ticks.y=element_line(color="black", linewidth=0.25), 
    axis.ticks.x=element_line(color="black", linewidth=0.25), 
    axis.line.y=element_line(color="black", linewidth=0.25), 
    axis.line.x=element_line(color="black", linewidth=0.25), 
    panel.grid.major=element_blank(), 
    panel.grid.minor=element_blank(), 
    panel.border=element_rect(color="black", fill=NA, linewidth=0.5)
  ) + 
  scale_fill_manual(
    values=c(
      "-2"="#2a4d77", 
      "-1"="#4f8acb", 
      "0"="#d0d0d0", 
      "1"="#ff8070",  
      "2"="#d95252", 
      "ambiguous"="#fcfcfc"
    ), 
    labels=c(
      "ambiguous"="ambiguous", 
      "2"="amplification", 
      "1"="gain", 
      "0"="neutral", 
      "-1"="loss", 
      "-2"="deep loss"
    )
  ) + 
  guides(fill=guide_legend(reverse=T))
pdf("stacked_bar_plot.pdf", width=8, height=6)
plot
dev.off()

