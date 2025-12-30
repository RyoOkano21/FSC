###-------------------------------------------------------------###
###             Simulation with autoregresssive model           ###
###-------------------------------------------------------------###

source("main_functions.R")
library(tidyverse)

## =======================
## Setting
## =======================
N <- 50
T_0 <- 9
T <- T_0+1
grids <- seq(0.01, 0.99, length=100)
M <- length(grids)

## =======================
## Generate pre-treatment outcomes
## =======================
set.seed(1)
func_vals_list <- vector("list", length = N)
for(i in 1:N){
  for(t in 1:T_0){
    func_vals_list[[i]][[t]] <- X_generate(grids, 10)
  }
}

## =======================
## Implement FSC
## =======================
weight_fscm <- FSCM(func_vals_list, T_0)

## =======================
## Implement augmented FSC
## =======================
K <- 20 # the number of basis functions

# select hyperparameter by cross-validation
obj_func <- function(lambda){
  cross_val(lambda, func_vals_list, T_0, K, grids)
}
lambda_opt <- as.numeric(optimise(obj_func, interval = c(0, 10000))[1])

# implement augmented FSC
weight_aug_1 <- FSCM_aug(func_vals_list, T_0, K = K, lambda = lambda_opt*100, grids = grids)
weight_aug_2 <- FSCM_aug(func_vals_list, T_0, K = K, lambda = lambda_opt, grids = grids)
weight_aug_3 <- FSCM_aug(func_vals_list, T_0, K = K, lambda = lambda_opt/100, grids = grids)

## =======================
## Monte Carlo simulation
## =======================
# low noise level
C <- 0.05
n_sim <- 500
error_box <- matrix(0, nrow = n_sim, ncol=5)
for(m in 1:n_sim){
  
  # generate post-treatment outcomes
  func_vals_list_post <- vector("list", N)
  set.seed(m)
  for(i in 1:N){
    func_vals_1 <- inner_vals(coef_func = coef_func_1, func_vals = func_vals_list[[i]][[T_0]], grids)
    func_vals_2 <- inner_vals(coef_func = coef_func_2, func_vals = func_vals_list[[i]][[T_0-1]], grids)
    func_vals_3 <- inner_vals(coef_func = coef_func_3, func_vals = func_vals_list[[i]][[T_0-2]], grids)
    error_vals <- error_generate_unif(grids, C)
    func_vals_list_post[[i]] <- func_vals_1 + func_vals_2 + func_vals_3 + error_vals
  }
  
  # compute synthetic outcomes
  control_matrix <- matrix(0, nrow = N-1, ncol=M)
  for(i in 2:N){
    control_matrix[i-1, ] <- func_vals_list_post[[i]]
  }
  predict_fscm <- t(weight_fscm) %*% control_matrix
  predict_ascm_1 <- t(weight_aug_1) %*% control_matrix
  predict_ascm_2 <- t(weight_aug_2) %*% control_matrix
  predict_ascm_3 <- t(weight_aug_3) %*% control_matrix
  predict_geodesic <- geodesic_OT(func_vals_list, func_vals_list_post, T_0 = T_0)
  
  # calculate errors
  error_box[m, 1]  <- sqrt(sum((predict_fscm - func_vals_list_post[[1]])^2))
  error_box[m, 2]  <- sqrt(sum((predict_ascm_1 - func_vals_list_post[[1]])^2))
  error_box[m, 3]  <- sqrt(sum((predict_ascm_2 - func_vals_list_post[[1]])^2))
  error_box[m, 4]  <- sqrt(sum((predict_ascm_3 - func_vals_list_post[[1]])^2))
  error_box[m, 5] <- sqrt(sum((predict_geodesic - func_vals_list_post[[1]])^2))
  print(m)
}
mat_0.05 <- error_box

# medium noise level
C <- 0.05
n_sim <- 500
error_box <- matrix(0, nrow = n_sim, ncol=5)
for(m in 1:n_sim){
  
  # generate post-treatment outcomes
  func_vals_list_post <- vector("list", N)
  set.seed(m)
  for(i in 1:N){
    func_vals_1 <- inner_vals(coef_func = coef_func_1, func_vals = func_vals_list[[i]][[T_0]], grids)
    func_vals_2 <- inner_vals(coef_func = coef_func_2, func_vals = func_vals_list[[i]][[T_0-1]], grids)
    func_vals_3 <- inner_vals(coef_func = coef_func_3, func_vals = func_vals_list[[i]][[T_0-2]], grids)
    error_vals <- error_generate_unif(grids, C)
    func_vals_list_post[[i]] <- func_vals_1 + func_vals_2 + func_vals_3 + error_vals
  }
  # compute synthetic outcomes
  control_matrix <- matrix(0, nrow = N-1, ncol=M)
  for(i in 2:N){
    control_matrix[i-1, ] <- func_vals_list_post[[i]]
  }
  predict_fscm <- t(weight_fscm) %*% control_matrix
  predict_ascm_1 <- t(weight_aug_1) %*% control_matrix
  predict_ascm_2 <- t(weight_aug_2) %*% control_matrix
  predict_ascm_3 <- t(weight_aug_3) %*% control_matrix
  predict_geodesic <- geodesic_OT(func_vals_list, func_vals_list_post, T_0 = T_0)
  
  # calculate errors
  error_box[m, 1]  <- sqrt(sum((predict_fscm - func_vals_list_post[[1]])^2))
  error_box[m, 2]  <- sqrt(sum((predict_ascm_1 - func_vals_list_post[[1]])^2))
  error_box[m, 3]  <- sqrt(sum((predict_ascm_2 - func_vals_list_post[[1]])^2))
  error_box[m, 4]  <- sqrt(sum((predict_ascm_3 - func_vals_list_post[[1]])^2))
  error_box[m, 5] <- sqrt(sum((predict_geodesic - func_vals_list_post[[1]])^2))
  print(m)
}
mat_0.2 <- error_box


# high noise level
C <- 1
n_sim <- 500
error_box <- matrix(0, nrow = n_sim, ncol=5)
for(m in 1:n_sim){
  
  # generate post-treatment outcomes
  func_vals_list_post <- vector("list", N)
  set.seed(m)
  for(i in 1:N){
    func_vals_1 <- inner_vals(coef_func = coef_func_1, func_vals = func_vals_list[[i]][[T_0]], grids)
    func_vals_2 <- inner_vals(coef_func = coef_func_2, func_vals = func_vals_list[[i]][[T_0-1]], grids)
    func_vals_3 <- inner_vals(coef_func = coef_func_3, func_vals = func_vals_list[[i]][[T_0-2]], grids)
    error_vals <- error_generate_unif(grids, C)
    func_vals_list_post[[i]] <- func_vals_1 + func_vals_2 + func_vals_3 + error_vals
  }
  # compute synthetic outcomes
  control_matrix <- matrix(0, nrow = N-1, ncol=M)
  for(i in 2:N){
    control_matrix[i-1, ] <- func_vals_list_post[[i]]
  }
  predict_fscm <- t(weight_fscm) %*% control_matrix
  predict_ascm_1 <- t(weight_aug_1) %*% control_matrix
  predict_ascm_2 <- t(weight_aug_2) %*% control_matrix
  predict_ascm_3 <- t(weight_aug_3) %*% control_matrix
  predict_geodesic <- geodesic_OT(func_vals_list, func_vals_list_post, T_0 = T_0)
  # calculate errors
  error_box[m, 1]  <- sqrt(sum((predict_fscm - func_vals_list_post[[1]])^2))
  error_box[m, 2]  <- sqrt(sum((predict_ascm_1 - func_vals_list_post[[1]])^2))
  error_box[m, 3]  <- sqrt(sum((predict_ascm_2 - func_vals_list_post[[1]])^2))
  error_box[m, 4]  <- sqrt(sum((predict_ascm_3 - func_vals_list_post[[1]])^2))
  error_box[m, 5] <- sqrt(sum((predict_geodesic - func_vals_list_post[[1]])^2))
  print(m)
}
mat_1.0 <- error_box

## =======================
##  Plot the results
## =======================
# create tibble for plot
mat_all <- rbind(mat_0.05, mat_0.2, mat_1.0)
tibble_all <- tibble(as.data.frame(mat_all))
noise_name_vec <- rep(c("Low noise level", "Medium noise level", "High noise level"), each=500)
tibble_all <- bind_cols(tibble_all, noise_name_vec)
names(tibble_all) <- c("FSC", "AFSC (i)", "AFSC (ii)",
                       "AFSC (iii)", "AGSC", "Noise")
tibble_long <- pivot_longer(tibble_all,
                            cols = c("FSC", "AFSC (i)", "AFSC (ii)",
                                     "AFSC (iii)", "AGSC"),
                            names_to = "method", values_to = "value")

# plot
tibble_long <- tibble_long %>%
  mutate(Noise = factor(Noise, levels = c("Low noise level", "Medium noise level", "High noise level")))
boxplot_auto <- ggplot(tibble_long, aes(x=factor(method), y=value)) +
  geom_boxplot() +
  facet_wrap(~ Noise) +
  scale_x_discrete(limits = rev.default(c("FSC", "AFSC (i)", "AFSC (ii)",
                                          "AFSC (iii)", "AGSC"))) +
  labs(x="Method", y="Estimation error") +
  theme(text = element_text(size = 20)) +
  scale_y_continuous(breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5), labels = c("0", "0.1", "0.2", "0.3", "0.4", "0.5")) +
  coord_flip()
boxplot_auto
