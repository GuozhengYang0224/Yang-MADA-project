###############################
# analysis script
#
#this script loads the processed, cleaned data, does a simple analysis
#and saves the results to the results folder

rm(list = ls())

#load needed packages. make sure they are installed.
library(tidyverse) # for data analysis
library(gt) # for making tables
library(gtsummary) # for summary table
library(ggplot2) # for plotting
library(ggpubr) # for combining sub-figures
library(here) # for data loading/saving
library(boot) # for bootstrapping: estimate CI

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
summary_tab <- data %>% dplyr::select(-c(ID, CRP_mg_dl, INF_Gamma_pg_ml, TNF_Alpha_pg_ml, IL_10_pg_mL,
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
              table.width=pct(60)) %>%
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

################################################################################
# Mediation effect of Leptin using all observations
################################################################################

# Create variable: log-transformed leptin level
data$log_leptin <- log10(data$Leptin_ng_ml)

# To calculate standardized coefficient, need standard deviation of the variables.
sd_CD4 <- sd(log10(data$CD4_plus), na.rm=T)
sd_FMI <- sd(data$FMI, na.rm=T)
sd_leptin <- sd(data$log_leptin, na.rm=T)

# Model1: log_leptin ~ FMI (linear regression)
model1 <- lm(log_leptin ~ FMI, data=data)
summary(model1)

# Model2: CD4_plus ~ FMI + log_leptin (Poisson regression using a log link)
model2 <- glm(CD4_plus ~ FMI + log_leptin, family=poisson(link='log'), data=data)
summary(model2)

# Summary of estimated coefficients
b1 <- coef(model1)["FMI"]  # Effect of FMI on leptin
b2 <- coef(model2)["log_leptin"]  # Effect of leptin on CD4_plus (controlling for FMI)
b3 <- coef(model2)["FMI"]  # Direct effect of FMI on CD4_plus

# Standardize the coefficients
b1_std <- b1 * (sd_FMI / sd_leptin)
b2_std <- b2 * (sd_leptin / sd_CD4)
b3_std <- b3 * (sd_FMI / sd_CD4)

# Estimate direct effects
(eff_dir <- b3_std)

# Estimate indirect effects
(eff_ind <- b1_std * b2_std)

# Estimate total effects
(eff_tot <- eff_dir + eff_ind)

# Percentage of total effect mediated
(eff_pct <- 100 * eff_ind / eff_tot)

################################################################################
# Bootstrapping for 95% CI estimation
################################################################################

# A bootstrapping function for re-sampling and estimating mediation effects
bs_md <- function(data, indices) {
  d <- data[indices, ]
  
  sd_CD4 <- sd(log10(d$CD4_plus), na.rm=T)
  sd_FMI <- sd(d$FMI, na.rm=T)
  sd_leptin <- sd(d$log_leptin, na.rm=T)
  
  model1 <- lm(log_leptin ~ FMI, data=d)
  model2 <- glm(CD4_plus ~ FMI + log_leptin, family=poisson(link='log'), data=d)
  
  b1 <- coef(model1)["FMI"]
  b2 <- coef(model2)["log_leptin"]
  b3 <- coef(model2)["FMI"]
  
  b1_std <- b1 * (sd_FMI / sd_leptin)
  b2_std <- b2 * (sd_leptin / sd_CD4)
  b3_std <- b3 * (sd_FMI / sd_CD4)
  
  eff_dir <- b3_std
  eff_ind <- b1_std * b2_std
  eff_tot <- eff_dir + eff_ind
  
  return(c(Direct_Effect=eff_dir, 
           Indirect_Effect=eff_ind, 
           Total_Effect=eff_tot))
}

# Set seed and apply bootstrapping: generate 'rs' replicates
set.seed(123)
rs <- 500
bs_results <- boot(data, bs_md, R=rs)
(ci_dir <- boot.ci(bs_results, conf=.95, index=1, type="perc")$percent[4:5])
(ci_ind <- boot.ci(bs_results, conf=.95, index=2, type="perc")$percent[4:5])
(ci_tot <- boot.ci(bs_results, conf=.95, index=3, type="perc")$percent[4:5])

################################################################################
# Mediation effect of Leptin using males
################################################################################

# Subset the data set to males and calculate standard deviation of variables
data_m <- data[data$Sex=="M",]
sd_CD4_m <- sd(log10(data_m$CD4_plus), na.rm=T)
sd_FMI_m <- sd(data_m$FMI, na.rm=T)
sd_leptin_m <- sd(data_m$log_leptin, na.rm=T)

# Model1: log_leptin ~ FMI (linear regression)
model1_m <- lm(log_leptin ~ FMI, data=data_m)
summary(model1_m)

# Model2: CD4_plus ~ FMI + log_leptin (Poisson regression using a log link)
model2_m <- glm(CD4_plus ~ FMI + log_leptin, family=poisson(link='log'), data=data_m)
summary(model2_m)

# Summary of estimated coefficients
b1_m <- coef(model1_m)["FMI"]  # Effect of FMI on leptin
b2_m <- coef(model2_m)["log_leptin"]  # Effect of leptin on CD4_plus (controlling for FMI)
b3_m <- coef(model2_m)["FMI"]  # Direct effect of FMI on CD4_plus

# Standardize the coefficients
b1_std_m <- b1_m * (sd_FMI_m / sd_leptin_m)
b2_std_m <- b2_m * (sd_leptin_m / sd_CD4_m)
b3_std_m <- b3_m * (sd_FMI_m / sd_CD4_m)

# Estimate direct effects
(eff_dir_m <- b3_std_m)

# Estimate indirect effects
(eff_ind_m <- b1_std_m * b2_std_m)

# Estimate total effects
(eff_tot_m <- eff_dir_m + eff_ind_m)

# Percentage of total effect mediated
(eff_pct_m <- 100 * eff_ind_m / eff_tot_m)

# Estimate 95% CI using bootstrap
set.seed(123)
bs_results_m <- boot(data_m, bs_md, R=rs)
(ci_dir_m <- boot.ci(bs_results_m, conf=.95, index=1, type="perc")$percent[4:5])
(ci_ind_m <- boot.ci(bs_results_m, conf=.95, index=2, type="perc")$percent[4:5])
(ci_tot_m <- boot.ci(bs_results_m, conf=.95, index=3, type="perc")$percent[4:5])

################################################################################
# Mediation effect of Leptin using females
################################################################################

# Subset the data set to males and calculate standard deviation of variables
data_f <- data[data$Sex=="F",]
sd_CD4_f <- sd(log10(data_f$CD4_plus), na.rm=T)
sd_FMI_f <- sd(data_f$FMI, na.rm=T)
sd_leptin_f <- sd(data_f$log_leptin, na.rm=T)

# Model1: log_leptin ~ FMI (linear regression)
model1_f <- lm(log_leptin ~ FMI, data=data_f)
summary(model1_f)

# Model2: CD4_plus ~ FMI + log_leptin (Poisson regression using a log link)
model2_f <- glm(CD4_plus ~ FMI + log_leptin, family=poisson(link='log'), data=data_f)
summary(model2_f)

# Summary of estimated coefficients
b1_f <- coef(model1_f)["FMI"]  # Effect of FMI on leptin
b2_f <- coef(model2_f)["log_leptin"]  # Effect of leptin on CD4_plus (controlling for FMI)
b3_f <- coef(model2_f)["FMI"]  # Direct effect of FMI on CD4_plus

# Standardize the coefficients
b1_std_f <- b1_f * (sd_FMI_f / sd_leptin_f)
b2_std_f <- b2_f * (sd_leptin_f / sd_CD4_f)
b3_std_f <- b3_f * (sd_FMI_f / sd_CD4_f)

# Estimate direct effects
(eff_dir_f <- b3_std_f)

# Estimate indirect effects
(eff_ind_f <- b1_std_f * b2_std_f)

# Estimate total effects
(eff_tot_f <- eff_dir_f + eff_ind_f)

# Percentage of total effect mediated
(eff_pct_f <- 100 * eff_ind_f / eff_tot_f)

# Estimate 95% CI using bootstrap
set.seed(123)
bs_results_f <- boot(data_f, bs_md, R=rs)
(ci_dir_f <- boot.ci(bs_results_f, conf=.95, index=1, type="perc")$percent[4:5])
(ci_ind_f <- boot.ci(bs_results_f, conf=.95, index=2, type="perc")$percent[4:5])
(ci_tot_f <- boot.ci(bs_results_f, conf=.95, index=3, type="perc")$percent[4:5])

################################################################################
# Summarize coefficients from regression models into a table
################################################################################

# Confidence interval: 95%
z <- qnorm(1-0.05/2)

# Estimates and 95% CI from model1 - all observations:
est1 <- c(round(coef(model1)["FMI"],3), 
          round(coef(model1)["FMI"]-z*summary(model1)$coefficients["FMI",2],3),
          round(coef(model1)["FMI"]+z*summary(model1)$coefficients["FMI",2],3))
est1_str <- sprintf("%.3f (%.3f, %.3f)", est1[1], est1[2], est1[3])
# Estimates and 95% CI from model2 - all observations:
est2_1 <- c(round(coef(model2)["FMI"],3), 
            round(coef(model2)["FMI"]-z*summary(model2)$coefficients["FMI",2],3),
            round(coef(model2)["FMI"]+z*summary(model2)$coefficients["FMI",2],3))
est2_2 <- c(round(coef(model2)["log_leptin"],3), 
            round(coef(model2)["log_leptin"]-z*summary(model2)$coefficients["log_leptin",2],3),
            round(coef(model2)["log_leptin"]+z*summary(model2)$coefficients["log_leptin",2],3))
est2 <- matrix(c(est2_1, est2_2), nr=2, byrow=T)
est2_str <- c(sprintf("%.3f (%.3f, %.3f)", est2[1,1], est2[1,2], est2[1,3]),
              sprintf("%.3f (%.3f, %.3f)", est2[2,1], est2[2,2], est2[2,3]))

# Estimates and 95% CI from model1 - males:
est1_m <- c(round(coef(model1_m)["FMI"],3), 
            round(coef(model1_m)["FMI"]-z*summary(model1_m)$coefficients["FMI",2],3),
            round(coef(model1_m)["FMI"]+z*summary(model1_m)$coefficients["FMI",2],3))
est1_str_m <- sprintf("%.3f (%.3f, %.3f)", est1_m[1], est1_m[2], est1_m[3])
# Estimates and 95% CI from model2 - males:
est2_1_m <- c(round(coef(model2_m)["FMI"],3), 
              round(coef(model2_m)["FMI"]-z*summary(model2_m)$coefficients["FMI",2],3),
              round(coef(model2_m)["FMI"]+z*summary(model2_m)$coefficients["FMI",2],3))
est2_2_m <- c(round(coef(model2_m)["log_leptin"],3), 
              round(coef(model2_m)["log_leptin"]-z*summary(model2_m)$coefficients["log_leptin",2],3),
              round(coef(model2_m)["log_leptin"]+z*summary(model2_m)$coefficients["log_leptin",2],3))
est2_m <- matrix(c(est2_1_m, est2_2_m), nr=2, byrow=T)
est2_str_m <- c(sprintf("%.3f (%.3f, %.3f)", est2_m[1,1], est2_m[1,2], est2_m[1,3]),
                sprintf("%.3f (%.3f, %.3f)", est2_m[2,1], est2_m[2,2], est2_m[2,3]))

# Estimates and 95% CI from model1 - females:
est1_f <- c(round(coef(model1_f)["FMI"],3), 
            round(coef(model1_f)["FMI"]-z*summary(model1_f)$coefficients["FMI",2],3),
            round(coef(model1_f)["FMI"]+z*summary(model1_f)$coefficients["FMI",2],3))
est1_str_f <- sprintf("%.3f (%.3f, %.3f)", est1_f[1], est1_f[2], est1_f[3])
# Estimates and 95% CI from model2 - females:
est2_1_f <- c(round(coef(model2_f)["FMI"],3), 
              round(coef(model2_f)["FMI"]-z*summary(model2_f)$coefficients["FMI",2],3),
              round(coef(model2_f)["FMI"]+z*summary(model2_f)$coefficients["FMI",2],3))
est2_2_f <- c(round(coef(model2_f)["log_leptin"],3), 
              round(coef(model2_f)["log_leptin"]-z*summary(model2_f)$coefficients["log_leptin",2],3),
              round(coef(model2_f)["log_leptin"]+z*summary(model2_f)$coefficients["log_leptin",2],3))
est2_f <- matrix(c(est2_1_f, est2_2_f), nr=2, byrow=T)
est2_str_f <- c(sprintf("%.3f (%.3f, %.3f)", est2_f[1,1], est2_f[1,2], est2_f[1,3]),
                sprintf("%.3f (%.3f, %.3f)", est2_f[2,1], est2_f[2,2], est2_f[2,3]))

# Put estimated coefficients into a data frame
est_a <- data.frame(Group=c("Total", rep("",2), "Male", rep("",2), "Female", rep("",2)),
                    Variable=rep(c("", "Fat mass index, kg/m^2", "Log leption level, ng/ml"),3),
                    Model1=c("", est1_str, "-", "", est1_str_m, "-", "", est1_str_f, "-"),
                    Model2=c("", est2_str, "", est2_str_m, "", est2_str_f))

# Summarize coefficients into a table
est_a_tab <- est_a %>%
  gt() %>%
  #tab_header(title="Coefficient estimation (95% CI) from linear regression and Poisson regression") %>%
  cols_label(Group=md("**Sample**"),
             Variable=md("**Variable**"),
             Model1=md("**Linear regression**"),
             Model2=md("**Poisson regression**")) %>%
  tab_options(table_body.hlines.color="transparent", 
              table.border.top.color="transparent",
              table.border.bottom.color="black",
              table.border.bottom.width=3,
              table.width=pct(60)) %>%
  tab_style(style=cell_borders(sides="bottom", color="black", weight=px(2)), locations=cells_column_labels()) %>%
  tab_style(style=cell_borders(sides="top", color="black", weight=px(2)), locations=cells_column_labels()) %>%
  tab_style(style=cell_borders(sides="bottom", color="black", weight=px(3)), locations=cells_title()) %>%
  tab_style(style=cell_text(size=px(13)), locations=cells_body()) %>%
  tab_style(style=cell_text(size=px(13)), locations=cells_column_labels()) %>%
  tab_style(style=cell_text(size=px(13)), locations=cells_column_spanners()) %>%
  tab_style(style=cell_text(size=px(13)), locations=cells_title()) %>%
  tab_style(style=cell_text(size=px(14), weight="bold"), locations=cells_title(groups="title"))
est_a_tab

# Save the coefficient summary table
est_a_tab %>% gt::gtsave(filename="table2.png", path=here("results", "tables"))

################################################################################
# Summarize mediation effect estimation into a table
################################################################################

# Put estimated effects and CIs into a data frame:
eff_a <- data.frame(Group=c("Total", "Male", "Female"),
                    Total=c(sprintf("%.2f (%.2f, %.2f)", eff_tot, ci_tot[1], ci_tot[2]),
                            sprintf("%.2f (%.2f, %.2f)", eff_tot_m, ci_tot_m[1], ci_tot_m[2]),
                            sprintf("%.2f (%.2f, %.2f)", eff_tot_f, ci_tot_f[1], ci_tot_f[2])),
                    Direct=c(sprintf("%.2f (%.2f, %.2f)", eff_dir, ci_dir[1], ci_dir[2]),
                             sprintf("%.2f (%.2f, %.2f)", eff_dir_m, ci_dir_m[1], ci_dir_m[2]),
                             sprintf("%.2f (%.2f, %.2f)", eff_dir_f, ci_dir_f[1], ci_dir_f[2])),
                    Indirect=c(sprintf("%.2f (%.2f, %.2f)", eff_ind, ci_ind[1], ci_ind[2]),
                               sprintf("%.2f (%.2f, %.2f)", eff_ind_m, ci_ind_m[1], ci_ind_m[2]),
                               sprintf("%.2f (%.2f, %.2f)", eff_ind_f, ci_ind_f[1], ci_ind_f[2])),
                    Percent=as.integer(c(round(eff_pct,0), round(eff_pct_m,0), round(eff_pct_f,0))))

# Correct <0 or >1 percentages: 
eff_a$Percent <- ifelse((eff_a$Percent>0 & eff_a$Percent<100), eff_a$Percent, "-")

# Summarize all estimated effects into a table
eff_a_tab <- eff_a %>%
  gt() %>%
  #tab_header(title="Mediation analysis with estimated effects and 95% CI") %>%
  cols_label(Group=md("**Sample**"),
             Total=md("**Total effect**"),
             Direct=md("**Direct effect**"),
             Indirect=md("**Indirect effect**"),
             Percent=md("**Percent of effect mediated**"),) %>%
  tab_options(table_body.hlines.color="transparent", 
              table.border.top.color="transparent",
              table.border.bottom.color="black",
              table.border.bottom.width=3,
              table.width=pct(75)) %>%
  tab_style(style=cell_borders(sides="bottom", color="black", weight=px(2)), locations=cells_column_labels()) %>%
  tab_style(style=cell_borders(sides="top", color="black", weight=px(2)), locations=cells_column_labels()) %>%
  tab_style(style=cell_borders(sides="bottom", color="black", weight=px(3)), locations=cells_title()) %>%
  tab_style(style=cell_text(size=px(13)), locations=cells_body()) %>%
  tab_style(style=cell_text(size=px(13)), locations=cells_column_labels()) %>%
  tab_style(style=cell_text(size=px(13)), locations=cells_column_spanners()) %>%
  tab_style(style=cell_text(size=px(13)), locations=cells_title()) %>%
  tab_style(style=cell_text(size=px(14), weight="bold"), locations=cells_title(groups="title"))
eff_a_tab

# Save the mediation analysis table
eff_a_tab %>% gt::gtsave(filename="table3.png", path=here("results", "tables"))

