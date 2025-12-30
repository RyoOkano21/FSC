###-------------------------------------------------------------###
###             Analysis of service trade data                  ###
###-------------------------------------------------------------###

load("service.RData")
source("main_functions.R")
library(dplyr)
library(reshape2)
library(corrplot)
library(tidyverse)
library(Matrix)
library(expm)

## =======================
## Extract information of data
## =======================
func_vals_list <- covvec_list
N <- length(func_vals_list) # the number of units
T <- length(func_vals_list[[1]]) # the number of periods
T_0 <-  29
M <- length(func_vals_list[[1]][[1]]) # the number of grids
num_com <- 9 #the number of service categories
grids <- seq(0.01, 0.99, length=M) # the number of grids

## =======================
## Implement FSC
## =======================
# implement FSC method
weight_fscm <- FSCM(func_vals_list, T_0)
round(weight_fscm, 3)

# compute synthetic outcomes
scm_outcomes <- vector("list", T)
for(t in 1:T){
  control_matrix <- matrix(0, nrow = N-1, ncol=M)
  for(i in 2:N){
    control_matrix[i-1, ] <- func_vals_list[[i]][[t]]
  }
  scm_outcomes[[t]] <- t(weight_fscm) %*% control_matrix
}

# compute pre-treatment fit
fit_values_scm <- rep(0, T_0)
for(t in 1:T_0){
  diffs <- (func_vals_list[[1]][[t]] - scm_outcomes[[t]])^2
  fit_values_scm[t] <- sum(diffs[-1])
}
round(sqrt(sum(fit_values_scm)), 4)


## =======================
## Implement Augmented FSC
## =======================
# ## select hyperparameter by cross-validation
# obj_func <- function(lambda){
#   cross_val_covmat(lambda, func_vals_list, T_0, grids)
# }
#
# lambda_opt <- optimise(obj_func, interval = c(0, 1))[1] # 0.00186 is selected

# implement augmented FSC
weight_aug <- FSCM_aug_covmat(func_vals_list, T_0,
                              lambda = 0.00186, grids = grids)

# compute synthetic outcomes
ascm_outcomes <- vector("list", T)
for(t in 1:T){
  control_matrix <- matrix(0, nrow = N-1, ncol=M)
  for(i in 2:N){
    control_matrix[i-1, ] <- func_vals_list[[i]][[t]]
  }
  ascm_outcomes[[t]] <- t(weight_aug) %*% control_matrix
}

# compute pre-treatment fit
fit_values_aug <- rep(0, T_0)
for(t in 1:T_0){
  diffs <- (func_vals_list[[1]][[t]] - ascm_outcomes[[t]])^2
  fit_values_aug[t] <- sum(diffs)
}
sqrt(sum(fit_values_aug))

## =======================
## Convert synthetic vectors to matrices
## =======================
# FSC
scm_outcomes_mat <- vector("list", T)
for(t in 1:T){
  temp_mat <- matrix(0, nrow=num_com, ncol=num_com)
  temp_mat[upper.tri(temp_mat, diag = T)] <- scm_outcomes[[t]]
  temp_mat <- temp_mat + t(temp_mat)
  scm_outcomes_mat[[t]] <- temp_mat - diag(diag(temp_mat)/2)
}

# Augmented FSC
ascm_outcomes_mat <- vector("list", T)
for(t in 1:T){
  temp_mat <- matrix(0, nrow=num_com, ncol=num_com)
  temp_mat[upper.tri(temp_mat, diag = T)] <- ascm_outcomes[[t]]
  temp_mat <- temp_mat + t(temp_mat)
  temp_mat_sym <- temp_mat - diag(diag(temp_mat)/2)
  ascm_outcomes_mat[[t]] <- nearPD(temp_mat_sym, base.matrix=TRUE)$mat
}

## =======================
## Plot the differences between observed and FSC outcomes (2009Q1 to 2014Q4) 
## =======================
# create tibble for plot
service_label <- c("SC", "SD", "SE", "SG", "SH", "SI", "SJ", "SK", "SL")
list_cov <- vector("list", T_0-5)
for(s in 1:(T_0-5)){
  list_cov[[s]] <- covmat_list[[1]][[s]] -  scm_outcomes_mat[[s]]
}
df_cov <- map2_df(
  list_cov,
  period_names_post <- c("2009 Q1", "2009 Q2", "2009 Q3", "2009 Q4",
                         "2010 Q1", "2010 Q2", "2010 Q3", "2010 Q4",
                         "2011 Q1", "2011 Q2", "2011 Q3", "2011 Q4",
                         "2012 Q1", "2012 Q2", "2012 Q3", "2012 Q4",
                         "2013 Q1", "2013 Q2", "2013 Q3", "2013 Q4",
                         "2014 Q1", "2014 Q2", "2014 Q3", "2014 Q4"),
  ~{
    melt(.x) %>%
      mutate(matrix_id = .y)
  }
)

# plot
causal_est_plot <- ggplot(df_cov, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red", midpoint = 0, limits=c(-13, 13)
  ) +
  coord_fixed() +
  facet_wrap(~ matrix_id, ncol = 5) +
  labs(fill = "",
       x=NULL, y=NULL) +
  scale_x_continuous(breaks = 1:9, labels = service_label) +
  scale_y_continuous(breaks = 1:9, labels = service_label) +
  theme(text = element_text(size = 15),
        axis.text.x = element_text(
          angle = 90,
          hjust = 1,
          vjust = 0.5
        ))
causal_est_plot

## =======================
##  Plot the differences between observed and FSC outcomes (2015Q1 to 2017Q2) 
## =======================
# create tibble for plot
list_cov <- vector("list", 10)
for(s in 1:10){
  list_cov[[s]] <- covmat_list[[1]][[T_0+s-5]] -  scm_outcomes_mat[[T_0+s-5]]
}
df_cov <- map2_df(
  list_cov,
  period_names_post <- c("2015 Q1", "2015 Q2", "2015 Q3", "2015 Q4", "2016 Q1", "2016 Q2",
                         "2016 Q3", "2016 Q4", "2017 Q1", "2017 Q2"),
  ~{
    melt(.x) %>%
      mutate(matrix_id = .y)
  }
)

# lot
causal_est_plot <- ggplot(df_cov, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red", midpoint = 0, limits=c(-13, 13)
  ) +
  coord_fixed() +
  facet_wrap(~ matrix_id, ncol = 5) +
  labs(fill = "",
       x=NULL, y=NULL) +
  scale_x_continuous(breaks = 1:9, labels = service_label) +
  scale_y_continuous(breaks = 1:9, labels = service_label, ) +
  theme(text = element_text(size = 20),
        axis.text.x = element_text(
          angle = 90,
          hjust = 1,
          vjust = 0.5
        ))
causal_est_plot

## =======================
## plot the differences between observed and FSC outcomes (2017Q3 to 2018Q2)
## =======================
# create tibble for plot
list_cov <- vector("list", 4)
for(s in 1:4){
  list_cov[[s]] <- covmat_list[[1]][[T_0+5+s]] -  scm_outcomes_mat[[T_0+5+s]]
}
df_cov <- map2_df(
  list_cov,
  period_names_post <- c("2017 Q3", "2017 Q4", "2018 Q1", "2018 Q2"),
  ~{
    melt(.x) %>%
      mutate(matrix_id = .y)
  }
)

# plot
causal_est_plot <- ggplot(df_cov, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red", midpoint = 0, limits=c(-13, 13)
  ) +
  coord_fixed() +
  facet_wrap(~ matrix_id, nrow = 1) +
  labs(fill = "",
       x=NULL, y=NULL) +
  scale_x_continuous(breaks = 1:9, labels = service_label) +
  scale_y_continuous(breaks = 1:9, labels = service_label, ) +
  theme(text = element_text(size = 20),
        axis.text.x = element_text(
          angle = 90,
          hjust = 1,
          vjust = 0.5
        ))
causal_est_plot

## =======================
## plot the differences between observed and AFSC outcomes (2009Q1 to 2014Q4)
## =======================
# create tibble for plot
list_cov <- vector("list", T_0-5)
for(s in 1:(T_0-5)){
  list_cov[[s]] <- covmat_list[[1]][[s]] -  ascm_outcomes_mat[[s]]
}
df_cov <- map2_df(
  list_cov,
  period_names_post <- c("2009 Q1", "2009 Q2", "2009 Q3", "2009 Q4",
                         "2010 Q1", "2010 Q2", "2010 Q3", "2010 Q4",
                         "2011 Q1", "2011 Q2", "2011 Q3", "2011 Q4",
                         "2012 Q1", "2012 Q2", "2012 Q3", "2012 Q4",
                         "2013 Q1", "2013 Q2", "2013 Q3", "2013 Q4",
                         "2014 Q1", "2014 Q2", "2014 Q3", "2014 Q4"),
  ~{
    melt(.x) %>%
      mutate(matrix_id = .y)
  }
)

# plot
causal_est_plot <- ggplot(df_cov, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red", midpoint = 0, limits=c(-13, 13)
  ) +
  coord_fixed() +
  facet_wrap(~ matrix_id, ncol = 5) +
  labs(fill = "",
       x=NULL, y=NULL) +
  scale_x_continuous(breaks = 1:9, labels = service_label) +
  scale_y_continuous(breaks = 1:9, labels = service_label) +
  theme(text = element_text(size = 15),
        axis.text.x = element_text(
          angle = 90,
          hjust = 1,
          vjust = 0.5
        ))
causal_est_plot

## =======================
## plot the differences between observed and AFSC outcomes (2015Q1 to 2017Q2)
## =======================
# create tibble for plot
list_cov <- vector("list", 10)
for(s in 1:10){
  list_cov[[s]] <- covmat_list[[1]][[T_0+s-5]] -  ascm_outcomes_mat[[T_0+s-5]]
}
df_cov <- map2_df(
  list_cov,
  period_names_post <- c("2015 Q1", "2015 Q2", "2015 Q3", "2015 Q4", "2016 Q1", "2016 Q2",
                         "2016 Q3", "2016 Q4", "2017 Q1", "2017 Q2"),
  ~{
    melt(.x) %>%
      mutate(matrix_id = .y)
  }
)

# plot
causal_est_plot <- ggplot(df_cov, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red", midpoint = 0, limits=c(-13, 13)
  ) +
  coord_fixed() +
  facet_wrap(~ matrix_id, ncol = 5) +
  labs(fill = "",
       x=NULL, y=NULL) +
  scale_x_continuous(breaks = 1:9, labels = service_label) +
  scale_y_continuous(breaks = 1:9, labels = service_label) +
  theme(text = element_text(size = 20),
        axis.text.x = element_text(
          angle = 90,
          hjust = 1,
          vjust = 0.5
        ))
causal_est_plot

## =======================
## plot the differences between observed and AFSC outcomes (2017Q3 to 2018Q2)
## =======================
# create tibble for plot
list_cov <- vector("list", 4)
for(s in 1:4){
  list_cov[[s]] <- covmat_list[[1]][[T_0+5+s]] -  ascm_outcomes_mat[[T_0+5+s]]
}
df_cov <- map2_df(
  list_cov,
  period_names_post <- c("2017 Q3", "2017 Q4", "2018 Q1", "2018 Q2"),
  ~{
    melt(.x) %>%
      mutate(matrix_id = .y)
  }
)

# plot
causal_est_plot <- ggplot(df_cov, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red", midpoint = 0, limits=c(-13, 13)
  ) +
  coord_fixed() +
  facet_wrap(~ matrix_id, nrow = 1) +
  labs(fill = "",
       x=NULL, y=NULL) +
  scale_x_continuous(breaks = 1:9, labels = service_label) +
  scale_y_continuous(breaks = 1:9, labels = service_label) +
  theme(text = element_text(size = 20),
        axis.text.x = element_text(
          angle = 90,
          hjust = 1,
          vjust = 0.5
        ))
causal_est_plot

## =======================
## Plot prediction bands for causal effects
## =======================
# compute residuals and its quantiles
alpha <- 0.1
residual_mat <- matrix(0, T_0, M)
for(t in 1:T_0){
  residual_mat[t, ] <- abs(ascm_outcomes[[t]]  - func_vals_list[[1]][[t]])
}
quantile_vec <- rep(0, M)
for(x in 1:M){
  temp <- ((T_0+1)*alpha-1)/T_0
  resid_vec <- residual_mat[, x]
  quantile_vec[x] <- quantile(resid_vec, 1-temp)
}

# compute the estimates, lower bound and upper bound
causal_est_mat_aug <- matrix(0, nrow = T-T_0, M)
causal_lower_mat_aug <- matrix(0, nrow = T-T_0, M)
causal_upper_mat_aug <- matrix(0, nrow = T-T_0, M)
for(t in (T_0+1):T){
  causal_est_mat_aug[(t-T_0), ] <- func_vals_list[[1]][[t]] - ascm_outcomes[[t]]
  causal_lower_mat_aug[(t-T_0), ] <-  func_vals_list[[1]][[t]] - ascm_outcomes[[t]] - quantile_vec
  causal_upper_mat_aug[(t-T_0), ] <-  func_vals_list[[1]][[t]] - ascm_outcomes[[t]] + quantile_vec
}

# convert vectors to matrices
lower_covmat <- vector("list", T-T_0)
for(t in  (T_0+1):T){
  temp_mat <- matrix(0, nrow=num_com, ncol=num_com)
  temp_mat[upper.tri(temp_mat, diag = TRUE)] <- causal_lower_mat_aug[(t-T_0), ]
  temp_mat <- temp_mat + t(temp_mat)
  temp_mat_sym <- temp_mat - diag(diag(temp_mat)/2)
  lower_covmat[[t-T_0]] <- temp_mat_sym
}

upper_covmat <- vector("list", T-T_0)
for(t in  (T_0+1):T){
  temp_mat <- matrix(0, nrow=num_com, ncol=num_com)
  temp_mat[upper.tri(temp_mat, diag = TRUE)] <- causal_upper_mat_aug[(t-T_0), ]
  temp_mat <- temp_mat + t(temp_mat)
  temp_mat_sym <- temp_mat - diag(diag(temp_mat)/2)
  upper_covmat[[t-T_0]] <- temp_mat_sym
}

# create tibble for plot lower bounds
length_post <- T - T_0
list_cov <- lower_covmat
df_cov <- map2_df(
  list_cov,
  period_names_post <- c("2016 Q2", "2016 Q3", "2016 Q4", "2017 Q1", "2017 Q2", "2017 Q3",
                         "2017 Q4", "2018 Q1", "2018 Q2"),
  ~{
    melt(.x) %>%
      mutate(matrix_id = .y)
  }
)

# plot lower bounds
causal_lower_plot <- ggplot(df_cov, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red", midpoint = 0, limits=c(-13, 13)
  ) +
  coord_fixed() +
  facet_wrap(~ matrix_id, ncol = 5) +
  labs(fill = "",
       x=NULL, y=NULL) +
  scale_x_continuous(breaks = 1:9, labels = service_label) +
  scale_y_continuous(breaks = 1:9, labels = service_label) +
  theme(text = element_text(size = 20),
        axis.text.x = element_text(
          angle = 90,
          hjust = 1,
          vjust = 0.5
        ))
causal_lower_plot

# create tibble for upper bounds
length_post <- T - T_0
list_cov <- upper_covmat
df_cov <- map2_df(
  list_cov,
  period_names_post <- c("2016 Q2", "2016 Q3", "2016 Q4", "2017 Q1", "2017 Q2", "2017 Q3",
                         "2017 Q4", "2018 Q1", "2018 Q2"),
  ~{
    melt(.x) %>%
      mutate(matrix_id = .y)
  }
)

# plot upper bounds
causal_upper_plot <- ggplot(df_cov, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red", midpoint = 0, limits=c(-13, 13)
  ) +
  coord_fixed() +
  facet_wrap(~ matrix_id, ncol = 5) +
  labs(fill = "",
       x=NULL, y=NULL) +
  scale_x_continuous(breaks = 1:9, labels = service_label) +
  scale_y_continuous(breaks = 1:9, labels = service_label) +
  theme(text = element_text(size = 20),
        axis.text.x = element_text(
          angle = 90,
          hjust = 1,
          vjust = 0.5
        ))
causal_upper_plot

## =======================
## placebo permutation test (very slow)
## =======================
# implement the permutation test 
length_post <- T-T_0
mag_mat_service <- matrix(0, nrow = length_post, ncol=N)
for(t in 1:length_post){
  mag_mat_service[t, ] <- placebo_covmat(func_vals_list, T_0, grids,
                                     post_period = T_0+t, a=0, b=20)
}


# calclute p-values
p_vals <- rep(0, length_post)
for(t in 1:length_post){
  distances <- mag_mat_service[t, ]
  p_vals[t] <- (sum(distances>distances[1])+1)/N
}

# plot the results
post_year <- c("2016 Q2 (p-value = 0.130)", "2016 Q3 (p-value = 0.130)", "2016 Q4 (p-value = 0.130)", "2017 Q1 (p-value = 0.087)",
               "2017 Q2 (p-value = 0.130)", "2017 Q3 (p-value = 0.130)", "2017 Q4 (p-value = 0.130)", "2018 Q1 (p-value = 0.130)",
               "2018 Q2 (p-value = 0.130)")
placebo_df <- data.frame(id=1:N, as.vector(t(mag_mat_service)), treatment_indicator=c(1, rep(0, N-1)),
                         rep(post_year, each=N), rep(mag_mat_service[, 1], each=N))
placebo_tibble <- tibble(placebo_df)
names(placebo_tibble) <- c("id", "magnitude", "treatment_indicator", "year", "treat_magnitude")

placebo_plot <-
  ggplot(data = placebo_tibble, mapping = aes(x=magnitude)) +
  geom_histogram() +
  geom_vline(aes(xintercept = treat_magnitude), colour = "red", linetype=2, linewidth=1.5) +
  facet_wrap(~year) +
  theme(legend.position = "bottom", text = element_text(size = 20)) +
  labs(x="Magnitude", y="Count")
placebo_plot
