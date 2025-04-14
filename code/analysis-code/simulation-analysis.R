###############################
# Simulation Script
#
#this script applies the RF model to simulated data set.

rm(list = ls())

#load needed packages. make sure they are installed.
library(tidyverse) # for data analysis
library(ggplot2) # for plotting
library(ggpubr) # for combining sub-figures
library(here) # for data loading/saving
library(tidymodels) # for ML models

# Simulate data set: N=1000
N <- 1000
# log_leptin ~ 8*Female + -0.1*Age^2 + 6*Age + exp(0.05*BMI) + error
# Female ~ Bernoulli(0.6)
# Age ~ N(35, 3)
# BMI ~ N(25, 1.5)
# error ~ N(0,1)

# Set a seed
rngseed <- 1234
set.seed(rngseed)

# Simulate a data set
data <- data.frame(Female=rbinom(N, size=1, prob=.6),
                   Age=rnorm(N, 35, 3),
                   BMI=rnorm(N, 25, 1.5),
                   error=rnorm(N, 0, 1))
data$log_leptin <- 8*data$Female + -0.1*data$Age^2 + 6*data$Age + exp(0.05*data$BMI) + data$error

# Select predictors for RF
data <- data %>% select(log_leptin, Female, Age, BMI)

# Set a seed
set.seed(rngseed)
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
rcp <- recipe(log_leptin ~ Female + Age + BMI, data=train_data) %>%
  step_normalize(all_numeric_predictors())

# Add flexibility: tune the two parameters
rf_spec <- rand_forest(mtry=tune(), min_n=tune(), trees=50) %>%
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
figure_rf2=here("results", "figures", "figure6.png")
ggsave(figure_rf2, rf_tune_plot, width=8, height=6, dpi=300)

# Select the best fit model
best_rf <- select_best(rf_tune, metric="rmse")
final_rf <- finalize_workflow(rf_wf, best_rf)

# Try the model fit on the test set
rf_fit <- last_fit(final_rf, split=data_split)
rf_result <- collect_metrics(rf_fit)

# RF: predicted value vs observed value
final_fit_rf <- fit(final_rf, data=data)
pred_rf <- predict(final_fit_rf, data) %>%
  bind_cols(data["log_leptin"])
colnames(pred_rf) <- c("pred", "Y")
rf_diagplot <- ggplot(pred_rf, aes(x=Y, y=pred))+
  geom_point(size=4, alpha=0.8, color="darkred")+
  geom_abline(intercept=0, slope=1, linetype="dashed", color="black", linewidth=1)+
  scale_x_continuous(limits=c(75, 105))+
  scale_y_continuous(limits=c(75, 105))+
  labs(x="Observed value", y="Predicted value", title="Random Forest")+
  theme_bw()+
  theme(axis.text=element_text(size=12, color="black"),
        axis.title.x=element_text(size=12, face="bold"),
        axis.title.y=element_text(size=12, face="bold", margin=margin(r=5)),
        title=element_text(size=14, color="black", face="bold"))

# Save the plot
figure7_file=here("results", "figures", "figure7.png")
ggsave(figure7_file, rf_diagplot, width=8, height=8, dpi=300)





