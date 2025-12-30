###-------------------------------------------------------------###
###                Simulation with factor model                 ###
###-------------------------------------------------------------###

source("main_functions.R")
library(tidyverse)

## =======================
## Settings
## =======================
N <- 50
T_0 <- 9
T <- T_0+1
grids <- seq(0.01, 0.99, length=100)
M <- length(grids)


## =======================
##  Generate factors and factor loadings
## =======================
J <- 5 # the number of factors

# generate factors
set.seed(2)
factor_vals_list <- vector("list", length = T)
for(t in 1:T){
  for(j in 1:J){
    factor_vals_list[[t]][[j]] <- mu_generate(grids, ind=t+j-1)
  }
}

# generate factor loadings
phi_mat <- matrix(rnorm(N*J), nrow= N, ncol=J)

## =======================
##  Monte Carlo simulations
## =======================
n_sim <- 500
error_box_list <- vector("list", 3)
for(c in c(1, 2, 3)){
  C <- c(0.02, 0.1, 0.5)[c] # noise level
  error_box <- matrix(0, nrow = n_sim, ncol=5)

  for(m in 1:n_sim){
    
    # generate outcomes
    func_vals_list <- vector("list", length = N)
    set.seed(m)
    for(i in 1:N){
      for(t in 1:T){
        factor_vals_mat <- matrix(0, nrow=J, ncol=M)
        for(j in 1:J){
          factor_vals_mat[j, ] <- factor_vals_list[[t]][[j]]
        }
        func_vals_list[[i]][[t]] <- (phi_mat[i, ] %*% factor_vals_mat)/100 + error_generate_unif(grids, C)
      }
    }

    # implement FSC
    weight_fscm <- FSCM(func_vals_list, T_0)

    # implement augmented FSC
    K <- 20
    obj_func <- function(lambda){
      cross_val(lambda, func_vals_list, T_0, K, grids)
    }
    lambda_opt <- as.numeric(optimise(obj_func, interval = c(0, 100000))[1])
    weight_aug_1 <- FSCM_aug(func_vals_list, T_0, K = K, lambda = lambda_opt*100, grids = grids) #cv optimalの100倍
    weight_aug_2 <- FSCM_aug(func_vals_list, T_0, K = K, lambda = lambda_opt, grids = grids) #cv optimal
    weight_aug_3 <- FSCM_aug(func_vals_list, T_0, K = K, lambda = lambda_opt/100, grids = grids) #cv optimalの0.01倍

    # divide outcomes into pre-treatment and post-treatment periods
    func_vals_list_pre <- vector("list", length = N)
    func_vals_list_post <- vector("list", length = N)
    for(i in 1:N){
      func_vals_list_pre[[i]] <- func_vals_list[[i]][1:T_0]
      func_vals_list_post[[i]] <- func_vals_list[[i]][[T]]
    }

    # compute synthetic outcomes
    control_matrix <- matrix(0, nrow = N-1, ncol=M)
    for(i in 2:N){
      control_matrix[i-1, ] <- func_vals_list[[i]][[T]]
    }
    predict_fscm <- t(weight_fscm) %*% control_matrix
    predict_ascm_1 <- t(weight_aug_1) %*% control_matrix
    predict_ascm_2 <- t(weight_aug_2) %*% control_matrix
    predict_ascm_3 <- t(weight_aug_3) %*% control_matrix
    predict_geodesic <- geodesic_OT(func_vals_list_pre, func_vals_list_post, T_0)

    # calculate prediction errors
    error_box[m, 1] <- sqrt(sum((predict_fscm -  func_vals_list[[1]][[T]])^2))
    error_box[m, 2] <- sqrt(sum((predict_ascm_1 - func_vals_list[[1]][[T]])^2))
    error_box[m, 3] <- sqrt(sum((predict_ascm_2 - func_vals_list[[1]][[T]])^2))
    error_box[m, 4] <- sqrt(sum((predict_ascm_3 - func_vals_list[[1]][[T]])^2))
    error_box[m, 5] <- sqrt(sum((predict_geodesic - func_vals_list[[1]][[T]])^2))
    print(m)
  }
  error_box_list[[c]] <- error_box
}

## =======================
## Plot the results
## =======================
# Create tibble
mat_1 <- error_box_list[[1]]
mat_2 <- error_box_list[[2]]
mat_3 <- error_box_list[[3]]
mat_all <- rbind(mat_1, mat_2, mat_3)
tibble_all <- tibble(as.data.frame(mat_all))
noise_name_vec <- rep(c("Low noise level", "Medium noise level", "High noise level"), each=n_sim)
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
boxplot_factor <- ggplot(tibble_long, aes(x=factor(method), y=value)) +
  geom_boxplot() +
  facet_wrap(~ Noise) +
  scale_x_discrete(limits = rev.default(c("FSC", "AFSC (i)", "AFSC (ii)",
                                          "AFSC (iii)", "AGSC"))) +
  labs(x="Method", y="Estimation error") +
  theme(text = element_text(size = 20)) +
  coord_flip()
boxplot_factor
