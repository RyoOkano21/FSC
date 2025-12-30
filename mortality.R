###-------------------------------------------------------------###
###                Analysis of mortality  data                  ###
###-------------------------------------------------------------###


load("aad.RData")
source("main_functions.R")
library(tidyverse)
library(Rearrangement)
library(frechet)

Density_list <- AaD_Density_list
Quant_list <- AaD_Quant_list

## =======================
## Extract information of data
## =======================
func_vals_list <- Quant_list
grids <- seq(0.01, 0.99, length=100)
N <- length(func_vals_list) # the number of units
T <- length(func_vals_list[[1]]) # the number of periods
M <- length(func_vals_list[[1]][[1]]) # the number of grids
T_0 <- 21 # the number of pre-treatment periods

## =======================
## Plot raw data (quantile functions)
## =======================
# create tibble for plot
Rus_list <- Quant_list[[1]]
visual_list <- Quant_list
visual_list[[19]] <- Rus_list
visual_list <- visual_list[2:19]
overall_mat <- c()
grids <- seq(0.01, 0.99, length=100)
for(i in 1:18){
  id <- i
  country_mat <- c()
  for(t in 1:30){
    year <- 1969+t
    quants <- visual_list[[i]][[t]]
    year_mat <- cbind(id, year, grids, quants)
    country_mat <- rbind(country_mat, year_mat)
  }
  overall_mat <- rbind(overall_mat, country_mat)
}
overall_mat <- cbind(overall_mat, overall_mat[, 1] != 18)
overall_df <- data.frame(overall_mat)
names(overall_df)[5] <- "treatment"
overall_tbl <- tibble(overall_df)

# plot
my_plot <- ggplot(data = overall_tbl, mapping = aes(x = grids, y=quants, group = id)) +
  geom_line(mapping = aes(color=factor(treatment), size=factor(treatment))) +
  scale_size_manual(values = c(1, 0.5), labels=c("Russia", "Control countries")) +
  scale_color_manual(values = c("#F8766D", "darkgray"), labels=c("Russia", "Control countries")) +
  guides(size = guide_legend(title=NULL)) +
  guides(colour = guide_legend(title=NULL)) +
  facet_wrap(~year, ncol=5) +
  labs(x="Probability", y="Quantile") +
  theme(legend.position = "top", text = element_text(size = 20)) +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", "0.25", "0.5", "0.75", "1"))
my_plot

## =======================
## Plot raw data (density functions) 
## =======================
# create tibble for plot
Rus_list <- Density_list[[1]]
visual_list <- Density_list
visual_list[[19]] <- Rus_list
visual_list <- visual_list[2:19]
overall_mat <- c()
for(i in 1:18){
  id <- i
  country_mat <- c()
  for(t in 1:30){
    year <- 1969+t
    age <- visual_list[[i]][[t]]$x
    dens <- visual_list[[i]][[t]]$y
    year_mat <- cbind(id, year, age, dens)
    country_mat <- rbind(country_mat, year_mat)
  }
  overall_mat <- rbind(overall_mat, country_mat)
}
overall_mat <- cbind(overall_mat, overall_mat[, 1] != 18)
overall_df <- data.frame(overall_mat)
names(overall_df)[5] <- "treatment"
overall_tbl <- tibble(overall_df)

# plot
my_plot <- ggplot(data = overall_tbl, mapping = aes(x = age, y=dens, group = id)) +
  geom_line(mapping = aes(color=factor(treatment), size=factor(treatment))) +
  scale_size_manual(values = c(1, 0.5), labels=c("Russia", "Control countries")) +
  scale_color_manual(values = c("#F8766D", "darkgray"), labels=c("Russia", "Control countries")) +
  guides(size = guide_legend(title=NULL)) +
  guides(colour = guide_legend(title=NULL)) +
  facet_wrap(~year, ncol=5) +
  labs(x="Age (years)", y="Density") +
  theme(legend.position = "top", text = element_text(size = 20))
my_plot

## =======================
## Implement FSC
## =======================
# change scale of data to conduct stable optimizationz
for(i in 1:N){
  for(t in 1:T){
    func_vals_list[[i]][[t]] <- func_vals_list[[i]][[t]]/1000
  }
}

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
sqrt(sum(fit_values_scm))

## =======================
## Implement Augmented FSC
## =======================
K <- 50 # the number of basis functions

# ## select hyperparameter by cross-validation
# obj_func <- function(lambda){
#   cross_val(lambda, func_vals_list, T_0, K, grids)
# }
#
# lambda_opt <- optimise(obj_func, interval = c(0, 10))[1] # 5.889182 is selected

# implement augmented FSC
weight_aug <- FSCM_aug(func_vals_list, T_0, K = K,
                         lambda = 5.889182, grids = grids)

# compute synthetic outcomes
ascm_outcomes <- vector("list", T)
for(t in 1:T){
  control_matrix <- matrix(0, nrow = N-1, ncol=M)
  for(i in 2:N){
    control_matrix[i-1, ] <- func_vals_list[[i]][[t]]
  }
  ascm_outcomes[[t]] <- t(weight_aug) %*% control_matrix
  # modify the outcomes
  ascm_outcomes[[t]] <- modif(vals= ascm_outcomes[[t]], low=0, upp=110/1000, grids)
}

# compute pre-treatment fit
fit_values_aug <- rep(0, T_0)
for(t in 1:T_0){
  diffs <- (func_vals_list[[1]][[t]] - ascm_outcomes[[t]])^2
  fit_values_aug[t] <- sum(diffs)
}
sqrt(sum(fit_values_aug))

## =======================
## Plot the observed and synthetic quantile functions (1970 to 1984)
## =======================
# create tibble for plot
overall_mat <- c()
for(t in 1:15){
  year <- 1969+t
  year_mat <- c()
  year_mat <- rbind(year_mat, cbind(year, grids, as.vector(func_vals_list[[1]][[t]])*1000, 1))
  year_mat <- rbind(year_mat, cbind(year, grids, as.vector(scm_outcomes[[t]])*1000, 2))
  year_mat <- rbind(year_mat, cbind(year, grids, as.vector(ascm_outcomes[[t]])*1000, 3))
  overall_mat <- rbind(overall_mat, year_mat)
}
overall_df <- data.frame(overall_mat)
names(overall_df) <- c("year", "grids", "quantiles", "type")
overall_tib <- tibble(overall_df)

# plot
synthetic_plot <-
  ggplot(data = overall_tib, mapping = aes(x = grids, y=quantiles, group = type)) +
  geom_line(mapping = aes(color=factor(type), linetype = factor(type)), linewidth = 1.1) +
  scale_color_manual(values = c("#F8766D","#00BA38", "#619CFF"), labels=c("Observed", "FSC", "Augmented FSC")) +
  scale_linetype_manual(values = c("solid", "longdash", "dashed"), labels=c("Observed", "FSC", "Augmented FSC")) +
  guides(colour = guide_legend(title=NULL)) +
  guides(linetype = guide_legend(title=NULL)) +
  facet_wrap(~year, ncol=3) +
  labs(x="Probability", y="Quantile") +
  theme(legend.position = "bottom", text = element_text(size = 20)) +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", "0.25", "0.5", "0.75", "1"))
synthetic_plot

## =======================
## Plot the observed and synthetic quantile functions (1985 to 1999)
## =======================
# create tibble
overall_mat <- c()
for(t in 16:T){
  year <- 1969+t
  year_mat <- c()
  year_mat <- rbind(year_mat, cbind(year, grids, as.vector(func_vals_list[[1]][[t]])*1000, 1))
  year_mat <- rbind(year_mat, cbind(year, grids, as.vector(scm_outcomes[[t]])*1000, 2))
  year_mat <- rbind(year_mat, cbind(year, grids, as.vector(ascm_outcomes[[t]])*1000, 3))
  overall_mat <- rbind(overall_mat, year_mat)
}
overall_df <- data.frame(overall_mat)
names(overall_df) <- c("year", "grids", "quantiles", "type")
overall_tib <- tibble(overall_df)

# plot
synthetic_plot <-
  ggplot(data = overall_tib, mapping = aes(x = grids, y=quantiles, group = type)) +
  geom_line(mapping = aes(color=factor(type), linetype = factor(type)), linewidth = 1.1) +
  scale_color_manual(values = c("#F8766D","#00BA38", "#619CFF"), labels=c("Observed", "FSC", "Augmented FSC")) +
  scale_linetype_manual(values = c("solid", "longdash", "dashed"), labels=c("Observed", "FSC", "Augmented FSC")) +
  guides(colour = guide_legend(title=NULL)) +
  guides(linetype = guide_legend(title=NULL)) +
  facet_wrap(~year, ncol=3) +
  labs(x="Probability", y="Quantile") +
  theme(legend.position = "bottom", text = element_text(size = 20)) +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", "0.25", "0.5", "0.75", "1"))
synthetic_plot

## =======================
## Plot the observed and synthetic density functions
## =======================
# create tibble
grids <- seq(0, 1, length=100)
grid_density = seq(0, 110, length=300)
overall_mat <- c()
for(t in 1:T){
  year <- 1969+t
  year_mat <- c()
  dens_vals_obs <- quant_dens(func_vals_list[[1]][[t]]*1000, grids, 5000, grid_density)
  dens_vals_scm <- quant_dens(as.vector(scm_outcomes[[t]]*1000), grids, 5000, grid_density)
  dens_vals_ascm <- quant_dens(ascm_outcomes[[t]]*1000, grids, 5000, grid_density)
  year_mat <- rbind(year_mat, cbind(year, grid_density, dens_vals_obs, 1))
  year_mat <- rbind(year_mat, cbind(year, grid_density, dens_vals_scm, 2))
  year_mat <- rbind(year_mat, cbind(year, grid_density, dens_vals_ascm, 3))
  overall_mat <- rbind(overall_mat, year_mat)
}
overall_df <- data.frame(overall_mat)
names(overall_df) <- c("year", "grid_density", "Density", "type")
overall_tib <- tibble(overall_df)

# plot
synthetic_dens_plot <-
  ggplot(data = overall_tib, mapping = aes(x = grid_density, y=Density, group = type)) +
  geom_line(mapping = aes(color=factor(type), linetype = factor(type)), linewidth = 1.1) +
  scale_color_manual(values = c("#F8766D","#00BA38", "#619CFF"), labels=c("Observed", "FSC", "Augmented FSC")) +
  scale_linetype_manual(values = c("solid", "longdash", "dashed"), labels=c("Observed", "FSC", "Augmented FSC")) +
  guides(colour = guide_legend(title=NULL)) +
  guides(linetype = guide_legend(title=NULL)) +
  facet_wrap(~year, ncol=5) +
  labs(x="Age (years)", y="Density") +
  theme(legend.position = "bottom", text = element_text(size = 20))
synthetic_dens_plot

## =======================
## Plot prediction bands for causal effects
## =======================
# compute residuals and its quantiles
alpha <- 0.1
residual_mat <- matrix(0, T_0, M)
for(t in 1:T_0){
  residual_mat[t, ] <- abs(ascm_outcomes[[t]] - func_vals_list[[1]][[t]])
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

# creat tibble  for plot
prediction_df_aug <- data.frame(rep(grids, T-T_0), as.vector(t(causal_est_mat_aug*1000)),
                                as.vector(t(causal_lower_mat_aug*1000)), as.vector(t(causal_upper_mat_aug*1000)),rep(1991:1999, each=M), "Augmented FSC")
names(prediction_df_aug) <- c("x", "est", "lower", "upper", "year", "type")
prediction_tibble <- tibble(prediction_df_aug)

# plot
prediction_plot_aug <- ggplot(data=prediction_tibble, aes(x = x, y = est)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3,
              color = NA, fill = "#619CFF") +
  geom_line(linewidth = 1.1, colour = "#619CFF", ) +
  facet_wrap(~year, nrow=4) +
  labs(
    title = "",
    x = "Probability",
    y = "Causal Effect",
    color = "",
    fill = "",
    linetype=""
  ) +
  theme(legend.position = "none",
        text = element_text(size = 20)
  ) +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", "0.25", "0.5", "0.75", "1")) +
  scale_y_continuous(limits = c(-30, 10))
prediction_plot_aug

## =======================
## placebo permutation test (very slow)
## =======================
# implement the permutation test 
length_post <- T-T_0
mag_mat_aad <- matrix(0, nrow = length_post, ncol=N)
for(t in 1:length_post){
  mag_mat_aad[t, ] <- placebo(func_vals_list, T_0, K, grids, post_period = T_0+t, a=0, b=20)
}

# calclute p-values
p_vals <- rep(0, length_post)
for(t in 1:length_post){
  distances <- mag_mat_aad[t, ]
  p_vals[t] <- (sum(distances>distances[1])+1)/N
}

# plot the results
post_year <- c("1991 (p-value = 0.389)", "1992 (p-value = 0.056)", "1993 (p-value = 0.056)", "1994 (p-value = 0.056)",
               "1995 (p-value = 0.056)", "1996 (p-value = 0.056)", "1997 (p-value = 0.056)", "1998 (p-value = 0.056)",
               "1999 (p-value = 0.056)")
placebo_df <- data.frame(id=1:N, as.vector(t(mag_mat_aad)), treatment_indicator=c(1, rep(0, N-1)),
                         rep(post_year, each=N), rep(mag_mat_aad[, 1], each=N))
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
