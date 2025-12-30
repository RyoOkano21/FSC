###-------------------------------------------------------------###
###          Main functions for the FSC methods                 ###
###-------------------------------------------------------------###

source("helper_functions.R")

## =======================
## Implement the FSC method
## =======================
FSCM <- function(func_vals_list, T_0){
  # input
  # func_vals_list: a list of outcomes of all units over all periods
  # T_0: the time period when the pre-treatment periods end

  # output
  # a vector of FSC weights

  # data info
  N <- length(func_vals_list) # the number of units
  N_0 <- N-1 # the number of control units
  M <- length(func_vals_list[[1]][[1]]) # the number of grids

  # list of outcomes in the pre-treatment periods
  pre_treat_vals_list <- vector("list", N)
  for(i in 1:N){
    pre_treat_vals_list[[i]] <- func_vals_list[[i]][1:T_0]
  }

  # matrices and vectors needed to conduct an optimization
  C_mat <- matrix(0, nrow = M*T_0, ncol = N_0)
  for(col in 1:N_0){
    i <- col+1
    unit_mat <- matrix(0, nrow=T_0, ncol=M)
    for(t in 1:T_0){
      unit_mat[t, ] <- pre_treat_vals_list[[i]][[t]]
    }
    C_mat[, col] <- as.vector(t(unit_mat))
  }

  unit_mat <- matrix(0, nrow=T_0, ncol=M)
  for(t in 1:T_0){
    unit_mat[t, ] <- pre_treat_vals_list[[1]][[t]]
  }
  d_vec <- as.vector(t(unit_mat))
  Aeq_mat <- matrix(1, nrow=1, ncol=N_0)
  beq_vec <- 1
  lb_vec <- rep(0, N_0)

  # conduct the optimization
  weight_scm <- my_lsqlincon(C = C_mat, d = d_vec, Aeq = Aeq_mat, beq = beq_vec, lb=lb_vec)
  return(round(weight_scm, 4))
}

## =======================
## Implement ridge augmented FSC with cubic B-splines
## =======================
FSCM_aug <- function(func_vals_list, T_0, K, lambda, grids){
  # input
  # func_vals_list: a list of outcomes of all units over all periods
  # T_0: the time period when the pre-treatment periods end
  # K: the number of orthogonal vectors
  # lambda: regularization parameter
  
  # output
  # a vector of augmented FSC weights
  
  # data info
  N <- length(func_vals_list)
  N_0 <- N-1
  T <- length(func_vals_list[[1]])
  M <- length(grids)

  # center the pre-treatment outcomes (demean across control units at each time)
  for(t in 1:T_0){
    time_func_mat <- matrix(0, nrow = N_0, ncol = M)
    for(i in 2:N){
      time_func_mat[i-1, ] <- func_vals_list[[i]][[t]]
    }
    mean_func_vec <- apply(X = time_func_mat, MARGIN = 2, FUN = mean)
    for(i in 1:N){
      func_vals_list[[i]][[t]] <- func_vals_list[[i]][[t]] - mean_func_vec
    }
  }

  # implement the FSC
  weight_fscm <- FSCM(func_vals_list, T_0)

  # compute inner products
  grids_width = 1/length(grids)
  Bsp <- Bsplines(x=grids,knots=seq(0,1,length=K-2))
  inner_prod_list <- vector("list", N)
  for(i in 1:N){
    for(t in 1:T){
      inner_prod <- (t(func_vals_list[[i]][[t]][-1]) %*% Bsp[-1, ])/grids_width
      inner_prod_list[[i]][[t]] <- inner_prod
    }
  }

  # compute the matrix r_{0\cdot} and vector r_1
  r_0dot <- matrix(0, nrow = N_0, ncol=K*T_0)
  for(i in 2:N){
    unit_inner_vec <- c()
    for(t in 1:T_0){
      unit_inner_vec <-  c(unit_inner_vec, inner_prod_list[[i]][[t]])
    }
    r_0dot[i-1, ] <- unit_inner_vec
  }
  r_1 <- numeric()
  for(t in 1:T_0){
    r_1 <-  c(r_1, inner_prod_list[[1]][[t]])
  }

  # compute the augmented weights
  weight_aug <- rep(0, N_0)
  vec_1 <- r_1 - t(r_0dot)%*%weight_fscm
  mat_1 <- solve(t(r_0dot)%*%r_0dot + lambda * diag(K*T_0))
  for(i in 2:N){
    vec_2 <- r_0dot[i-1, ]
    weight_aug[i-1] <- weight_fscm[i-1] + t(vec_1) %*%  mat_1  %*% vec_2
  }

  return(weight_aug)
}

## =======================
## Select the hyperparameter lambda by cross-validation (CV)
## =======================
cross_val <- function(lambda, func_vals_list, T_0, K, grids){
  # output
  # cross-validation error
  
  # data info
  N <- length(func_vals_list)

  # implement CV
  error_vals <- numeric()
  for(k in 1:T_0){

    # remove data of k-th time period
    list_leave <- vector("list", length = N)
    for(i in 1:N){
      list_leave[[i]] <- func_vals_list[[i]][-k]
    }

    # compute weights
    weight_aug_leave <- FSCM_aug(list_leave, (T_0-1), K,
                                 lambda = lambda, grids)

    # compute synthetic outcomes
    control_matrix <- matrix(0, nrow = N-1, ncol=length(grids))
    for(i in 2:N){
      control_matrix[i-1, ] <- func_vals_list[[i]][[k]]
    }
    ascm_outcome <- t(weight_aug_leave) %*% control_matrix

    # differences between the observations and synthetic outcomes
    error_vals[k] <- sum((func_vals_list[[1]][[k]] - ascm_outcome)^2)
  }

  # compute error
  sum(error_vals)
}

## =======================
## Implement the placebo permutation test using the augmented FSC
## =======================
placebo <- function(func_vals_list, T_0, K, grids, post_period, a, b){
  
  # data info
  N <- length(func_vals_list)
  
  d_vals <- rep(0, N)
  for(i in 1:N){
    # change the ordering of list
    unit_list <- func_vals_list[i]
    func_vals_list_new <- c(unit_list, func_vals_list)
    func_vals_list_new <- func_vals_list_new[-(i+1)]

    # select hyperparameter
    obj_func <- function(lambda){
      cross_val(lambda, func_vals_list_new, T_0, K, grids)
    }
    lambda_opt <- optimise(obj_func, interval = c(a, b))[1]

    # compute the augmented weights
    weight_aug <- FSCM_aug(func_vals_list_new, T_0, K,
                           lambda = as.numeric(lambda_opt), grids)

    # compute the synthetic outcomes
    control_matrix <- matrix(0, nrow = N-1, ncol=length(grids))
    for(j in 2:N){
      control_matrix[j-1, ] <- func_vals_list_new[[j]][[post_period]]
    }
    ascm_outcome <- t(weight_aug) %*% control_matrix

    # compute the magnitudes of causal effects
    d_vals[i] <-  sqrt(sum((func_vals_list_new[[1]][[post_period]] - ascm_outcome)^2))
  }

  d_vals
}

## =======================
## implement the augmented geodesic synthetic control
## =======================
geodesic_OT <- function(func_vals_list, func_vals_list_post, T_0){

  # matrix needed to compute the coefficients
  for(u in 1:length(grids)){
    mat_pre <- matrix(0, nrow = T_0, ncol = T_0)
    vec_pre <- rep(0, T_0)

    mat_pre_1 <- matrix(0, nrow = N-1, ncol = T_0)
    for(i in 2:N){
      for(j in 1:T_0){
        mat_pre_1[i-1, j] <- func_vals_list[[i]][[j]][u]
      }
    }
    mat_pre_2 <- t(mat_pre_1) %*% mat_pre_1
    mat_pre_3 <- rep(0, N-1)
    for(i in 2:N){
      mat_pre_3[i-1] <- func_vals_list_post[[i]][u]
    }
    mat_pre_4 <- t(mat_pre_1) %*% mat_pre_3

    mat_pre <- mat_pre + mat_pre_2
    vec_pre <- vec_pre + mat_pre_4
  }

  # compute coefficients
  coef_est_alpha <- solve(mat_pre, vec_pre)

  # compute model outcomes
  model_outcomes_list <- vector("list", length = N)
  for(i in 1:N){
    pre_vals_mat <- matrix(0, nrow = T_0, ncol = length(grids))
    for(t in 1:T_0){
      pre_vals_mat[t, ] <- func_vals_list[[i]][[t]]
    }
    model_outcomes_list[[i]] <- t(coef_est_alpha) %*% pre_vals_mat
  }

  # compute estimates of counterfactual outcomes
  control_matrix_model <- matrix(0, nrow = N-1, ncol=length(grids))
  for(i in 2:N){
    control_matrix_model[i-1, ] <- model_outcomes_list[[i]]
  }
  est <- t(weight_fscm) %*% control_matrix +
    model_outcomes_list[[1]] - t(weight_fscm) %*% control_matrix_model

  return(est)
}

## =======================
## ridge augmented FSC when the outcomes are covariance matrices
## =======================
FSCM_aug_covmat <- function(func_vals_list, T_0, lambda, grids){
  
  # data info
  N <- length(func_vals_list)
  N_0 <- N-1
  T <- length(func_vals_list[[1]])
  M <- length(grids)

  ## center the pre-treatment outcomes
  for(t in 1:T_0){
    time_func_mat <- matrix(0, nrow = N_0, ncol = M)
    for(i in 2:N){
      time_func_mat[i-1, ] <- func_vals_list[[i]][[t]]
    }
    mean_func_vec <- apply(X = time_func_mat, MARGIN = 2, FUN = mean)
    for(i in 1:N){
      func_vals_list[[i]][[t]] <- func_vals_list[[i]][[t]] - mean_func_vec
    }
  }

  # implement FSC
  weight_fscm <- FSCM(func_vals_list, T_0)

  # compute the matrix r_{0\cdot} and vector r_1
  r_0dot <- matrix(0, nrow = N_0, ncol=M*T_0)
  for(i in 2:N){
    unit_inner_vec <- c()
    for(t in 1:T_0){
      unit_inner_vec <-  c(unit_inner_vec, func_vals_list[[i]][[t]])
    }
    r_0dot[i-1, ] <- unit_inner_vec
  }
  r_1 <- numeric()
  for(t in 1:T_0){
    r_1 <-  c(r_1, func_vals_list[[1]][[t]])
  }

  # compute the augmented weights
  weight_aug <- rep(0, N_0)
  vec_1 <- r_1 - t(r_0dot)%*%weight_fscm
  mat_1 <- solve(t(r_0dot)%*%r_0dot + lambda * diag(M*T_0))
  for(i in 2:N){
    vec_2 <- r_0dot[i-1, ]
    weight_aug[i-1] <- weight_fscm[i-1] + t(vec_1) %*%  mat_1  %*% vec_2
  }
  return(weight_aug)
}

## =======================
## select hyperparameter in the ridge augmented FSC when outcomes are covariance matrices
## =======================
cross_val_covmat <- function(lambda, func_vals_list, T_0, grids){

  # data info
  N <- length(func_vals_list)

  # implement CV
  error_vals <- numeric()
  for(k in 1:T_0){
    # remove data of k-th time period
    list_leave <- vector("list", length = N)
    for(i in 1:N){
      list_leave[[i]] <- func_vals_list[[i]][-k]
    }

    # compute weights
    weight_aug_leave <- FSCM_aug_covmat(list_leave, (T_0-1),
                                        lambda = lambda, grids)

    # compute synthetic outcomes
    control_matrix <- matrix(0, nrow = N-1, ncol=length(grids))
    for(i in 2:N){
      control_matrix[i-1, ] <- func_vals_list[[i]][[k]]
    }
    ascm_outcome <- t(weight_aug_leave) %*% control_matrix

    # differences between the observations and synthetic outcomes
    error_vals[k] <- sum((func_vals_list[[1]][[k]] - ascm_outcome)^2)
  }

  # compute error
  sum(error_vals)
}

## =======================
## placebo permutation test using the augmented FSC when outcomes are covariance matrices
## =======================
placebo_covmat <- function(func_vals_list, T_0, grids, post_period, a, b){
  
  # data info
  N <- length(func_vals_list)
  
  d_vals <- rep(0, N)
  for(i in 1:N){
    # change the ordering of list
    unit_list <- func_vals_list[i]
    func_vals_list_new <- c(unit_list, func_vals_list)
    func_vals_list_new <- func_vals_list_new[-(i+1)]

    # select hyperparameter
    obj_func <- function(lambda){
      cross_val_covmat(lambda, func_vals_list_new, T_0, grids)
    }
    lambda_opt <- optimise(obj_func, interval = c(a, b))[1]

    # compute the augmented weights
    weight_aug <- FSCM_aug_covmat(func_vals_list_new, T_0,
                                  lambda = as.numeric(lambda_opt), grids)

    # compute the synthetic outcomes
    control_matrix <- matrix(0, nrow = N-1, ncol=length(grids))
    for(j in 2:N){
      control_matrix[j-1, ] <- func_vals_list_new[[j]][[post_period]]
    }
    ascm_outcome <- t(weight_aug) %*% control_matrix

    # Compute the magnitudes of causal effects
    d_vals[i] <-  sqrt(sum((func_vals_list_new[[1]][[post_period]] - ascm_outcome)^2))
    print(i)
  }
  
  d_vals
}


