###############################
# Machine Learning Script
#
#this script applies ML methods to analyze the data set.
#Performances of regression tree and random forest are compared. 

rm(list = ls())

#load needed packages. make sure they are installed.
library(tidyverse) # for data analysis
library(gt) # for making tables
library(gtsummary) # for summary table
library(ggplot2) # for plotting
library(ggpubr) # for combining sub-figures
library(here) # for data loading/saving
library(tidymodels) # for ML models

#path to data
#note the use of the here() package and not absolute paths
data_location <- here::here("data","processed-data","processeddata.rds")

# Load data. 
data <- readRDS(data_location)
head(data)

# Define the sex variable: 1-Female
data$Female <- ifelse(data$Sex=="F", 1, 0)
# Convert scale of leptin
data$log_leptin <- log10(data$Leptin_ng_ml)

# Select predictors for ML models
data <- data %>% select(log_leptin, Female, Age, BMI, LBM_pct, Fat_pct)

# Set a seed
rngseed <- 1234

# 5-fold CV
data_split <- initial_split(data, prop=.8)
train_data <- training(data_split)
test_data  <- testing(data_split)
folds <- vfold_cv(train_data, v=5)

################################################################################
# Random Forest model
################################################################################

# Set seed
set.seed(rngseed)

# Random forest model
rcp <- recipe(log_leptin ~ Female + Age + BMI + LBM_pct + Fat_pct, data=train_data) %>%
  step_normalize(all_numeric_predictors())

# Add flexibility: tune the two parameters
rf_spec <- rand_forest(mtry=tune(), min_n=tune(), trees=10) %>%
  set_engine("ranger") %>%
  set_mode("regression")
rf_wf <- workflow() %>% add_model(rf_spec) %>% add_recipe(rcp)

# Tune the two parameters based on RMSE
rf_tune <- tune_grid(rf_wf, resamples=folds, grid=10, metrics=metric_set(rmse))

# Make a plot of tuning results
rf_tune_df <- as.data.frame(rf_tune$.metrics)
rf_tune_plot <- ggplot(rf_tune_df, aes(x=mtry, y=min_n, fill=.estimate))+
  geom_tile()+
  scale_fill_viridis_c(name="RMSE")+
  labs(x="Number of randomly sampled predictors", y="Minimum number of split points")+
  theme_bw()+
  theme(axis.text=element_text(size=12, color="black"),
        axis.title.x=element_text(size=12, face="bold"),
        axis.title.y=element_text(size=12, face="bold", margin=margin(r=5)),
        legend.position="top",
        legend.title=element_text(size=12, face="bold"),
        legend.text=element_text(size=12),
        legend.key.size=unit(1,"cm"))
figure_rf=here("results", "figures", "figure3.png")
ggsave(figure_rf, rf_tune_plot, width=8, height=6, dpi=300)

# Select the best fit model
best_rf <- select_best(rf_tune, metric="rmse")
final_rf <- finalize_workflow(rf_wf, best_rf)

# Try the model fit on the test set
rf_fit <- last_fit(final_rf, split=data_split)
rf_result <- collect_metrics(rf_fit)

################################################################################
# Regression tree
################################################################################

# Set seed
set.seed(rngseed)

# Regression tree model
tree_spec <- decision_tree(cost_complexity=tune(), tree_depth=tune()) %>%
  set_engine("rpart") %>%
  set_mode("regression")
tree_wf <- workflow() %>% add_model(tree_spec) %>% add_recipe(rcp)

# Tune the two parameters based on RMSE
tree_tune <- tune_grid(tree_wf, resamples=folds, grid=10, metrics=metric_set(rmse))

# Make a plot of tuning results
tree_tune_df <- as.data.frame(tree_tune$.metrics)
tree_tune_plot <- ggplot(tree_tune_df, aes(x=factor(tree_depth), y=factor(round(cost_complexity,5)), 
                                           fill=.estimate))+
  geom_tile()+
  scale_fill_viridis_c(name="RMSE")+
  labs(x="Tree depth", y="Cost complexity") +
  theme_bw()+
  theme(axis.text=element_text(size=12, color="black"),
        axis.title.x=element_text(size=12, face="bold"),
        axis.title.y=element_text(size=12, face="bold", margin=margin(r=5)),
        legend.position="top",
        legend.title=element_text(size=12, face="bold"),
        legend.text=element_text(size=12),
        legend.key.size=unit(1,"cm"))
figure_tree=here("results", "figures", "figure4.png")
ggsave(figure_tree, tree_tune_plot, width=8, height=6, dpi=300)

# Select the best fit model
best_tree <- select_best(tree_tune, metric="rmse")
final_tree <- finalize_workflow(tree_wf, best_tree)

# Try the model fit on the test set
tree_fit <- last_fit(final_tree, split=data_split)
tree_result <- collect_metrics(tree_fit)

################################################################################
# Compare: Random Forest vs Regression Tree
################################################################################

# RF: predicted value vs observed value
final_fit_rf <- fit(final_rf, data=data)
pred_rf <- predict(final_fit_rf, data) %>%
  bind_cols(data["log_leptin"])
colnames(pred_rf) <- c("pred", "Y")
rf_diagplot <- ggplot(pred_rf, aes(x=Y, y=pred))+
  geom_point(size=4, alpha=0.8, color="darkred")+
  geom_abline(intercept=0, slope=1, linetype="dashed", color="black", linewidth=1)+
  scale_x_continuous(limits=c(0, 1.8))+
  scale_y_continuous(limits=c(0, 1.8))+
  labs(x="Observed value", y="Predicted value", title="Random Forest")+
  theme_bw()+
  theme(axis.text=element_text(size=12, color="black"),
        axis.title.x=element_text(size=12, face="bold"),
        axis.title.y=element_text(size=12, face="bold", margin=margin(r=5)),
        title=element_text(size=14, color="black", face="bold"))

# RT: predicted value vs observed value
final_fit_tree <- fit(final_tree, data=data)
pred_tree <- predict(final_fit_tree, data) %>%
  bind_cols(data["log_leptin"])
colnames(pred_tree) <- c("pred", "Y")
tree_diagplot <- ggplot(pred_tree, aes(x=Y, y=pred))+
  geom_point(size=4, alpha=0.8, color="darkred")+
  geom_abline(intercept=0, slope=1, linetype="dashed", color="black", linewidth=1)+
  scale_x_continuous(limits=c(0, 1.8))+
  scale_y_continuous(limits=c(0, 1.8))+
  labs(x="Observed value", y="Predicted value", title="Regression Tree")+
  theme_bw()+
  theme(axis.text=element_text(size=12, color="black"),
        axis.title.x=element_text(size=12, face="bold"),
        axis.title.y=element_text(size=12, face="bold", margin=margin(r=5)),
        title=element_text(size=14, color="black", face="bold"))

# Combine the two plots and save
diag_comb <- ggarrange(rf_diagplot, tree_diagplot, ncol=2, nrow=1, align="h")
diag_comb
figure5_file=here("results", "figures", "figure5.png")
ggsave(figure5_file, diag_comb, width=12, height=6, dpi=300)





