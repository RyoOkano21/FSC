###-------------------------------------------------------------###
###                Analysis of fertility  data                  ###
###-------------------------------------------------------------###

load("asfr.RData")
source("main_functions.R")
library(tidyverse)


## =======================
## Extract information of data
## =======================
func_vals_list <- asfr_list
grids <- seq(0, 1, length=44)
N <- length(func_vals_list) # the number of units
T <- length(func_vals_list[[1]]) # the number of periods
M <- length(func_vals_list[[1]][[1]]) # the number of grids
T_0 <- 16 # the number of pre-treatment periods


## =======================
## Plot raw data 
## =======================
# create tibble for plot
func_vals_list <- asfr_list
N <- length(func_vals_list)
T <- length(func_vals_list[[1]])

visual_list <- asfr_list
visual_list[[N+1]] <- asfr_list[[1]] 
visual_list <-visual_list[2:(N+1)]
ids <- 1:N
times <- 1:T
age <- 12:55
overall_mat <- c()
for(i in 1:N){
  id <- i
  country_mat <- c()
  for(t in 1:T){
    year <- 1955+t
    asfr <- asfr_list[[i]][[t]]
    year_mat <- cbind(id, year, age, asfr)
    country_mat <- rbind(country_mat, year_mat)
  }
  overall_mat <- rbind(overall_mat, country_mat)
}
overall_mat <- cbind(overall_mat, overall_mat[, 1] != N) 
overall_df <- data.frame(overall_mat)
names(overall_df)[5] <- "treatment"
overall_tbl <- tibble(overall_df)

# Plot
my_plot <- ggplot(data =overall_tbl, mapping = aes(x = age, y=asfr, group = id)) +
  geom_line(mapping = aes(color=factor(treatment), size=factor(treatment))) +
  scale_size_manual(values = c(1, 0.5), labels=c("East Germany", "Control countries")) +
  scale_color_manual(values = c("#F8766D", "darkgray"), labels=c("East Germany", "Control countries")) +
  guides(size = guide_legend(title=NULL)) + #凡例タイトルを非表示
  guides(colour = guide_legend(title=NULL)) + #凡例タイトルを非表示
  facet_wrap(~year, nrow=5) +
  labs(x="Age", y="ASFR") +
  theme(legend.position = "top", text = element_text(size = 20))
my_plot


## =======================
## Implement FSC 
## =======================
# implement FSC method
weight_fscm <- FSCM(func_vals_list, T_0)

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
K <- 50 # the number of basis functions

# ## select hyperparameter by cross-validation
# obj_func <- function(lambda){
#   cross_val(lambda, func_vals_list, T_0, K, grids)
# }
# lambda_opt <- optimise(obj_func, interval = c(0, 10))[1] # 5.759146 is selected

# implement augmented FSC
weight_aug <- FSCM_aug(func_vals_list, T_0, K = K,
                       lambda = 5.759146, grids)

# compute synthetic outcomes
ascm_outcomes <- vector("list", T)
for(t in 1:T){
  control_matrix <- matrix(0, nrow = N-1, ncol=M)
  for(i in 2:N){
    control_matrix[i-1, ] <- func_vals_list[[i]][[t]]
  }
  ascm_outcomes[[t]] <- t(weight_aug) %*% control_matrix
}

# compute  pre-treatment fits
fit_values_aug <- rep(0, T_0)
for(t in 1:T_0){
  diffs <- (func_vals_list[[1]][[t]] - ascm_outcomes[[t]])^2
  fit_values_aug[t] <- sum(diffs)
}
sqrt(sum(fit_values_aug))


## =======================
## Plot the observed outcomes and synthetic outcomes (1956 to 1963)
## =======================
ascm_outcomes_2 <- ascm_outcomes
weight_aug_2 <- weight_aug

# create tibble for plot
overall_mat <- c()
age <- 12:55
for(t in 1:8){
  year <- 1955+t
  year_mat <- c()
  year_mat <- rbind(year_mat, cbind(year, age, as.vector(asfr_list[[1]][[t]]), 1))
  year_mat <- rbind(year_mat, cbind(year, age, as.vector(scm_outcomes[[t]]), 2))
  year_mat <- rbind(year_mat, cbind(year, age, as.vector(ascm_outcomes_2[[t]]), 3))
  overall_mat <- rbind(overall_mat, year_mat)
}
overall_df <- data.frame(overall_mat)
names(overall_df) <- c("year", "age", "asfr", "type")
overall_tib <- tibble(overall_df)

# plot
synthetic_plot <-
  ggplot(data = overall_tib, mapping = aes(x = age, y=asfr, group = type)) +
  geom_line(mapping = aes(color=factor(type), linetype = factor(type)), linewidth = 1.1) +
  scale_color_manual(values = c("#F8766D","#00BA38", "#619CFF"), labels=c("Observed", "FSC", "Augmented FSC")) +
  scale_linetype_manual(values = c("solid", "longdash", "dashed"), labels=c("Observed", "FSC", "Augmented FSC")) +
  guides(colour = guide_legend(title=NULL)) +
  guides(linetype = guide_legend(title=NULL)) +
  facet_wrap(~year, ncol=3) +
  labs(x="Age", y="ASFR") +
  theme(legend.position = "bottom", text = element_text(size = 20))
synthetic_plot


## =======================
## plot the observed outcomes and synthetic outcomes (1964 to 1975) 
## =======================
ascm_outcomes_2 <- ascm_outcomes
weight_aug_2 <- weight_aug

# create tibble for plot
overall_mat <- c()
age <- 12:55
for(t in c(T_0-7, T_0-6, T_0-5, T_0-4, T_0-3, T_0-2,
           T_0-1, T_0, T_0+1, T_0+2, T_0+3, T_0+4)){
  year <- 1955+t
  year_mat <- c()
  year_mat <- rbind(year_mat, cbind(year, age, as.vector(asfr_list[[1]][[t]]), 1))
  year_mat <- rbind(year_mat, cbind(year, age, as.vector(scm_outcomes[[t]]), 2))
  year_mat <- rbind(year_mat, cbind(year, age, as.vector(ascm_outcomes_2[[t]]), 3))
  overall_mat <- rbind(overall_mat, year_mat)
}
overall_df <- data.frame(overall_mat)
names(overall_df) <- c("year", "age", "asfr", "type")
overall_tib <- tibble(overall_df)

# plot
synthetic_plot <-
  ggplot(data = overall_tib, mapping = aes(x = age, y=asfr, group = type)) +
  geom_line(mapping = aes(color=factor(type), linetype = factor(type)), linewidth = 1.1) +
  scale_color_manual(values = c("#F8766D","#00BA38", "#619CFF"), labels=c("Observed", "FSC", "Augmented FSC")) +
  scale_linetype_manual(values = c("solid", "longdash", "dashed"), labels=c("Observed", "FSC", "Augmented FSC")) +
  guides(colour = guide_legend(title=NULL)) +
  guides(linetype = guide_legend(title=NULL)) +
  facet_wrap(~year, nrow=4) +
  labs(x="Age", y="ASFR") +
  theme(legend.position = "bottom", text = element_text(size = 20))
synthetic_plot


## =======================
## plot prediction bands for causal effects
## =======================
# compute residuals and its quantiles
alpha <- 0.1
residual_mat <- matrix(0, T_0, M)
for(t in 1:T_0){
  residual_mat[t, ] <- abs(ascm_outcomes_2[[t]] - func_vals_list[[1]][[t]])
}
quantile_vec <- rep(0, M)
for(x in 1:M){
  temp <- ((T_0+1)*alpha-1)/T_0
  resid_vec <- residual_mat[1:T_0, x]
  quantile_vec[x] <- quantile(resid_vec, 1-temp)
}

# compute the estimates, lower bound and upper bound
causal_est_mat_aug <- matrix(0, nrow = T-T_0, M)
causal_lower_mat_aug <- matrix(0, nrow = T-T_0, M)
causal_upper_mat_aug <- matrix(0, nrow = T-T_0, M)
for(t in (T_0+1):T){
  causal_est_mat_aug[(t-T_0), ] <- func_vals_list[[1]][[t]] - ascm_outcomes_2[[t]]
  causal_lower_mat_aug[(t-T_0), ] <-  func_vals_list[[1]][[t]] - ascm_outcomes_2[[t]] - quantile_vec
  causal_upper_mat_aug[(t-T_0), ] <-  func_vals_list[[1]][[t]] - ascm_outcomes_2[[t]] + quantile_vec
}

# creat tibble  for plot
prediction_df <- data.frame(rep(12:55, 4), as.vector(t(causal_est_mat_aug)),
                            as.vector(t(causal_lower_mat_aug)), as.vector(t(causal_upper_mat_aug)),rep(1972:1975, each=M))
names(prediction_df) <- c("x", "est", "lower", "upper", "year")
prediction_df <- rbind.data.frame(prediction_df)
prediction_tibble <- tibble(prediction_df)

# plot
prediction_plot <- ggplot(data=prediction_tibble, aes(x = x, y = est)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, color = NA, fill = "#619CFF") +
  geom_line(linewidth = 1.1, colour = "#619CFF") +
  facet_wrap(~year, nrow=2) +
  labs(
    x = "Age",
    y = "Causal Effect",
    color = "",
    fill = "",
    linetype=""
  )  +
  theme(legend.position = "none",
        text = element_text(size = 20)
  )
prediction_plot


## =======================
## placebo permutation test (very slow)
## =======================
# implement the permutation test
length_post <- T-T_0
mag_mat_asfr <- matrix(0, nrow = length_post, ncol=N)
for(t in 1:length_post){
  mag_mat_asfr[t, ] <- placebo(func_vals_list, T_0, K, grids, post_period = T_0+t, a=0, b=20)
}

# calclute p-values
p_vals <- rep(0, length_post)
for(t in 1:length_post){
  distances <- mag_mat_asfr[t, ]
  p_vals[t] <- (sum(distances>distances[1])+1)/N
}

# plot the results
my_vec <- mag_mat_asfr[, 1] + c(0.0005, 0, 0.0012, -0.0015)
post_year <- c("1972 (p-value = 0.095)", "1973 (p-value = 0.048)", "1974 (p-value = 0.048)", "1975 (p-value = 0.095)")
placebo_df <- data.frame(id=1:N, as.vector(t(mag_mat_asfr)), treatment_indicator=c(1, rep(0, N-1)),
                         rep(post_year, each=N), rep(my_vec, each=N))
placebo_tibble <- tibble(placebo_df)
names(placebo_tibble) <- c("id", "magnitude", "treatment_indicator", "year", "treat_magnitude")
country_names <- c("East Germany", control_country_names)

placebo_plot <-
  ggplot(data = placebo_tibble, mapping = aes(x=magnitude)) +
  geom_histogram() +
  geom_vline(aes(xintercept = treat_magnitude), colour = "red", linetype=2, linewidth=1.5) +
  facet_wrap(~year) +
  theme(legend.position = "bottom", text = element_text(size = 20)) +
  labs(x="Magnitude", y="Count")
placebo_plot
