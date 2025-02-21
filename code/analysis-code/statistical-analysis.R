###############################
# analysis script
#
#this script loads the processed, cleaned data, does a simple analysis
#and saves the results to the results folder

#load needed packages. make sure they are installed.
library(tidyverse) # for data analysis
library(gt) # for making tables
library(gtsummary) # for summary table
library(ggplot2) # for plotting
library(ggpubr) # for combining sub-figures
library(here) #for data loading/saving

#path to data
#note the use of the here() package and not absolute paths
data_location <- here::here("data","processed-data","processeddata.rds")

# Load data. 
data <- readRDS(data_location)
head(data)

# Define the order of sex
data$Sex <- factor(data$Sex, levels=c("M", "F"))

################################################################################
# Comparing body indices  between Male vs Female
################################################################################

# Make a summary table 
summary_tab <- data %>% select(-c(ID, CRP_mg_dl, INF_Gamma_pg_ml, TNF_Alpha_pg_ml, IL_10_pg_mL,
                                  CD4_plus, CD4_plus_pct, CD8_plus, CD8_plus_pct)) %>%
  tbl_summary(by=Sex, type=list(where(is.numeric) ~ "continuous"),
              statistic=list(all_continuous() ~ "{median} ({p25}, {p75})"),
              digits=all_continuous() ~ 2, missing="no",
              label=list(Age ~ "Age, years",
                         Height ~ "Height, cm",
                         Weight ~ "Weight, kg",
                         BMI ~ "Body mass index",
                         TBW_L ~ "Body water, L",
                         TBW_pct ~ "Percentage of body water, %",
                         Fat_kg ~ "Body fat, kg",
                         Fat_pct ~ "Percentage of body fat, %",
                         FMI ~ "Fat mass index",
                         LBM_kg ~ "Lean body mass, kg",
                         LBM_pct ~ "Percentage of lean body mass, %",
                         FFMI ~ "Fat-free mass index",
                         Leptin_ng_ml ~ "Letpin level, ng/ml")) %>%
  # add_p(test=all_continuous() ~ "wilcox.test",
  #       test.args=all_continuous() ~ list(exact=F),
  #       pvalue_fun=function(x) style_number(x, digits=3)) %>%
  modify_header(label="**Variable**, median(Q1, Q3)",
                stat_1="**Male**, N={n}",
                stat_2="**Female**, N={n}",
                p.value="*p*-value") %>%
  #modify_spanning_header(all_stat_cols() ~ "**Sex**") %>%
  modify_footnote(stat_1=NA, stat_2=NA) %>%
  as_gt() %>%
  #tab_header(title=md("**Table 1. Body Characteristics and Compositions by Sex**")) %>%
  tab_options(table_body.hlines.color="transparent", 
              table.border.top.color="transparent",
              table.border.bottom.color="black",
              table.border.bottom.width=3,
              table.width=pct(80)) %>%
  tab_style(style=cell_borders(sides="bottom", color="black", weight=px(2)), locations=cells_column_labels()) %>%
  tab_style(style=cell_borders(sides="top", color="black", weight=px(2)), locations=cells_column_labels()) %>%
  tab_style(style=cell_borders(sides="bottom", color="black", weight=px(3)), locations=cells_title()) %>%
  tab_style(style=cell_text(size=px(13)), locations=cells_body()) %>%
  tab_style(style=cell_text(size=px(13)), locations=cells_column_labels()) %>%
  tab_style(style=cell_text(size=px(13)), locations=cells_column_spanners()) %>%
  tab_style(style=cell_text(size=px(13)), locations=cells_title()) %>%
  tab_style(style=cell_text(size=px(14), weight="bold"), locations=cells_title(groups="title"))
summary_tab

# Save the table
summary_tab %>% gt::gtsave(filename="table1.png", path=here("results", "tables"))

################################################################################
# Comparing immune response between Male vs Female
################################################################################

# Define color scale and theme before making plots
colors_sex <- c("M"="steelblue2", "F"="orange2") # Colors for male and female
theme_MADA <- theme_bw()+
  theme(axis.text=element_text(size=10, color="black"),
        axis.title.x=element_text(size=10, face="bold"),
        axis.title.y=element_text(size=10, face="bold", margin=margin(r=0)),
        legend.position="none")

# Make a boxplot comparing CRP level by sex
boxplot_crp <- data %>% 
  ggplot(aes(x=Sex, y=CRP_mg_dl, fill=Sex, color=Sex))+
  geom_boxplot(linewidth=1, width=.4, alpha=.6)+
  scale_fill_manual(values=colors_sex)+
  scale_color_manual(values=colors_sex)+
  scale_x_discrete(breaks=c("M", "F"), labels=c("Male", "Female"))+
  scale_y_continuous(breaks=seq(0,16,2), limits=c(0,15))+
  labs(x="", y="CRP level, mg/dl")+
  theme_MADA
boxplot_crp

# Make a boxplot comparing INF-Gamma level by sex
boxplot_infgma <- data %>% 
  ggplot(aes(x=Sex, y=INF_Gamma_pg_ml, fill=Sex, color=Sex))+
  geom_boxplot(linewidth=1, width=.4, alpha=.6)+
  scale_fill_manual(values=colors_sex)+
  scale_color_manual(values=colors_sex)+
  scale_x_discrete(breaks=c("M", "F"), labels=c("Male", "Female"))+
  scale_y_continuous(breaks=seq(0,120,30), limits=c(0,120))+
  labs(x="", y="INF-gamma level, pg/ml")+
  theme_MADA
boxplot_infgma

# Make a boxplot comparing TNF-Alpha level by sex
boxplot_tnfaph <- data %>% 
  mutate(TNF_Alpha_pg_ml=TNF_Alpha_pg_ml/1000) %>% # Change the scale of TNF-alpha
  ggplot(aes(x=Sex, y=TNF_Alpha_pg_ml, fill=Sex, color=Sex))+
  geom_boxplot(linewidth=1, width=.4, alpha=.6)+
  scale_fill_manual(values=colors_sex)+
  scale_color_manual(values=colors_sex)+
  scale_x_discrete(breaks=c("M", "F"), labels=c("Male", "Female"))+
  scale_y_continuous(breaks=seq(1,3.5,.5), limits=c(1,3.2))+
  labs(x="", y="TNF-alpha level, 10^3 pg/ml")+
  theme_MADA
boxplot_tnfaph

# Make a boxplot comparing IL-10 level by sex
boxplot_il10 <- data %>% 
  ggplot(aes(x=Sex, y=IL_10_pg_mL, fill=Sex, color=Sex))+
  geom_boxplot(linewidth=1, width=.4, alpha=.6)+
  scale_fill_manual(values=colors_sex)+
  scale_color_manual(values=colors_sex)+
  scale_x_discrete(breaks=c("M", "F"), labels=c("Male", "Female"))+
  scale_y_continuous(breaks=seq(4,20,4), limits=c(4,20))+
  labs(x="", y="IL-10 level, pg/ml")+
  theme_MADA
boxplot_il10

# Make a boxplot comparing CD4+ count level by sex
boxplot_cd4c <- data %>% 
  ggplot(aes(x=Sex, y=CD4_plus, fill=Sex, color=Sex))+
  geom_boxplot(linewidth=1, width=.4, alpha=.6)+
  scale_fill_manual(values=colors_sex)+
  scale_color_manual(values=colors_sex)+
  scale_x_discrete(breaks=c("M", "F"), labels=c("Male", "Female"))+
  scale_y_continuous(breaks=seq(0,800,200), limits=c(0,800))+
  labs(x="", y="CD4+ count, /ul")+
  theme_MADA
boxplot_cd4c

# Make a boxplot comparing CD4+ percentage level by sex
boxplot_cd4p <- data %>% 
  ggplot(aes(x=Sex, y=CD4_plus_pct, fill=Sex, color=Sex))+
  geom_boxplot(linewidth=1, width=.4, alpha=.6)+
  scale_fill_manual(values=colors_sex)+
  scale_color_manual(values=colors_sex)+
  scale_x_discrete(breaks=c("M", "F"), labels=c("Male", "Female"))+
  scale_y_continuous(breaks=seq(0,60,15), limits=c(0,60))+
  labs(x="", y="CD4+ percentage, %")+
  theme_MADA
boxplot_cd4p

# Make a boxplot comparing CD8+ count level by sex
boxplot_cd8c <- data %>% 
  ggplot(aes(x=Sex, y=CD8_plus, fill=Sex, color=Sex))+
  geom_boxplot(linewidth=1, width=.4, alpha=.6)+
  scale_fill_manual(values=colors_sex)+
  scale_color_manual(values=colors_sex)+
  scale_x_discrete(breaks=c("M", "F"), labels=c("Male", "Female"))+
  scale_y_continuous(breaks=seq(0,800,200), limits=c(0,800))+
  labs(x="", y="CD8+ count, /ul")+
  theme_MADA
boxplot_cd8c

# Make a boxplot comparing CD8+ percentage level by sex
boxplot_cd8p <- data %>% 
  ggplot(aes(x=Sex, y=CD8_plus_pct, fill=Sex, color=Sex))+
  geom_boxplot(linewidth=1, width=.4, alpha=.6)+
  scale_fill_manual(values=colors_sex)+
  scale_color_manual(values=colors_sex)+
  scale_x_discrete(breaks=c("M", "F"), labels=c("Male", "Female"))+
  scale_y_continuous(breaks=seq(0,60,15), limits=c(0,60))+
  labs(x="", y="CD8+ percentage, %")+
  theme_MADA
boxplot_cd8p

# Combine the 8 boxplots for immune responses
boxplot_comb <- ggarrange(boxplot_crp, boxplot_infgma, boxplot_tnfaph, boxplot_il10,
                          boxplot_cd4c, boxplot_cd4p, boxplot_cd8c, boxplot_cd8p,
                          ncol=4, nrow=2, align="hv")
boxplot_comb

# Save the figure
figure1_file=here("results", "figures", "figure1.png")
ggsave(figure1_file, boxplot_comb, width=8, height=6, dpi=300)

################################################################################
# Correlation between leptin level and body composition
################################################################################

# Define a new color scale and theme
colors_sex <- c("M"="steelblue2", "F"="palevioletred2") # Colors for male and female
theme_MADA2 <- theme_bw()+
  theme(axis.text=element_text(size=12, color="black"),
        axis.title.x=element_text(size=12, face="bold"),
        axis.title.y=element_text(size=12, face="bold", margin=margin(r=5)),
        legend.position="top",
        legend.title=element_text(size=12, face="bold"),
        legend.text=element_text(size=12))

# Make a scatterplot for leptin ~ fat
cor_lpt_fat <- cor(data$Fat_kg, log10(data$Leptin_ng_ml), use="complete.obs", method="pearson")
sctplot_lpt_fat <- data %>% 
  ggplot(aes(x=Fat_kg, y=log10(Leptin_ng_ml), fill=Sex, size=Weight))+
  geom_point(alpha=.8, stroke=.5, color="black", shape=21)+
  scale_fill_manual(values=colors_sex, labels=c("M"="Male", "F"="Female"))+
  scale_size_continuous(range=c(2,8))+
  labs(x="Body fat, kg", y="Log leption level, ng/ml")+
  annotate("text", x=min(data$Fat_kg, na.rm=T), y=max(log10(data$Leptin_ng_ml), na.rm=T), 
           label=paste0("R = ", round(cor_lpt_fat, 2)), 
           hjust=0, vjust=1, fontface="bold", size=4)+
  theme_MADA2+
  guides(size="none", fill=guide_legend(override.aes=list(size=6)))
sctplot_lpt_fat

# Make a scatterplot for leptin ~ lean body mass
cor_lpt_lbm <- cor(data$LBM_kg, log10(data$Leptin_ng_ml), use="complete.obs", method="pearson")
sctplot_lpt_lbm <- data %>% 
  ggplot(aes(x=LBM_kg, y=log10(Leptin_ng_ml), fill=Sex, size=Weight))+
  geom_point(alpha=.8, stroke=.5, color="black", shape=21)+
  scale_fill_manual(values=colors_sex, labels=c("M"="Male", "F"="Female"))+
  scale_size_continuous(range=c(2,8))+
  labs(x="Lean body mass, kg", y="Log leption level, ng/ml")+
  annotate("text", x=min(data$LBM_kg, na.rm=T), y=max(log10(data$Leptin_ng_ml), na.rm=T), 
           label=paste0("R = ", round(cor_lpt_lbm, 2)), 
           hjust=0, vjust=1, fontface="bold", size=4)+
  theme_MADA2+
  guides(size="none", fill=guide_legend(override.aes=list(size=6)))
sctplot_lpt_lbm

# Make a scatterplot for CD4+ ~ leptin
cor_cd4_lpt <- cor(log10(data$Leptin_ng_ml), data$CD4_plus_pct, use="complete.obs", method="pearson")
sctplot_cd4_lpt <- data %>% 
  ggplot(aes(x=log10(Leptin_ng_ml), y=CD4_plus_pct, fill=Sex, size=Weight))+
  geom_point(alpha=.8, stroke=.5, color="black", shape=21)+
  scale_fill_manual(values=colors_sex, labels=c("M"="Male", "F"="Female"))+
  scale_size_continuous(range=c(2,8))+
  labs(x="Log leption level, ng/ml", y="Percentage of CD4+ cell, %")+
  annotate("text", x=min(log10(data$Leptin_ng_ml), na.rm=T), y=max(data$CD4_plus_pct, na.rm=T)-3, 
           label=paste0("R = ", round(cor_cd4_lpt, 2)), 
           hjust=0, vjust=1, fontface="bold", size=4)+
  theme_MADA2+
  guides(size="none", fill=guide_legend(override.aes=list(size=6)))
sctplot_cd4_lpt

# Make a scatterplot for CD8+ ~ leptin
cor_cd8_lpt <- cor(log10(data$Leptin_ng_ml), data$CD8_plus_pct, use="complete.obs", method="pearson")
sctplot_cd8_lpt <- data %>% 
  ggplot(aes(x=log10(Leptin_ng_ml), y=CD8_plus_pct, fill=Sex, size=Weight))+
  geom_point(alpha=.8, stroke=.5, color="black", shape=21)+
  scale_fill_manual(values=colors_sex, labels=c("M"="Male", "F"="Female"))+
  scale_size_continuous(range=c(2,8))+
  labs(x="Log leption level, ng/ml", y="Percentage of CD8+ cell, %")+
  annotate("text", x=min(log10(data$Leptin_ng_ml), na.rm=T), y=max(data$CD8_plus_pct, na.rm=T)-3, 
           label=paste0("R = ", round(cor_cd8_lpt, 2)), 
           hjust=0, vjust=1, fontface="bold", size=4)+
  theme_MADA2+
  guides(size="none", fill=guide_legend(override.aes=list(size=6)))
sctplot_cd8_lpt

# Combine the 4 scatterplots for correlation visualization
sctplot_comb <- ggarrange(sctplot_lpt_fat, sctplot_lpt_lbm, sctplot_cd4_lpt, sctplot_cd8_lpt,
                          ncol=2, nrow=2, align="hv", common.legend=T)
sctplot_comb

# Save the figure
figure2_file=here("results", "figures", "figure2.png")
ggsave(figure2_file, sctplot_comb, width=8, height=8, dpi=300)




















