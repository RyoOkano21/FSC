###-------------------------------------------------------------###
###                    Helper functions                         ###
###-------------------------------------------------------------###

library(pracma)
library(cubicBsplines)
library(Rearrangement)


## =======================
## Modified lsqlincon() to conduct stabilized optimizations
## =======================
my_lsqlincon <- function(C, d,                     # min ||C x - d||_2
                         A = NULL,   b = NULL,     # A x   <= b
                         Aeq = NULL, beq = NULL,   # Aeq x == beq
                         lb = NULL,  ub = NULL)    # lb <= x <= ub
{
  if (!requireNamespace("quadprog", quietly = TRUE)) {
    stop("quadprog needed for this function to work. Please install it.",
         call. = FALSE)
  }

  stopifnot(is.numeric(C), is.numeric(d))
  if (is.null(A) && !is.null(b) || !is.null(A) && is.null(b))
    stop("If any, both 'A' and 'b' must be NULL.")
  if (is.null(Aeq) && !is.null(beq) || !is.null(Aeq) && is.null(beq))
    stop("If any, both 'Aeq' and 'beq' must be NULL.")

  if (!is.matrix(C)) C <- matrix(C, 1)
  mc  <- nrow(C);   nc  <- ncol(C);  n <- nc
  if (length(d) != mc)
    stop("Dimensions of 'C' and 'd' do not fit.")
  if (is.null(A) && is.null(Aeq) && is.null(lb) && is.null(ub))
    return(qr.solve(C, d))

  if (!is.null(A)) {
    if (!is.matrix(A)) A <- matrix(A, 1)
    ma  <- nrow(A);   na  <- ncol(A)
    if (na != n)
      stop("Number of columns of 'A' does not fit with 'C'.")
    # ATTENTION: quadprog requires  A x >= b !
    A <- -A; b <- -b
  }
  if (is.null(Aeq)) {
    meq <- 0
  } else {
    if (!is.matrix(Aeq)) Aeq <- matrix(Aeq, 1)
    meq  <- nrow(Aeq);   neq  <- ncol(Aeq)
    if (neq != n)
      stop("Number of columns of 'Aeq' does not fit with 'C'.")
  }

  if (is.null(lb)) {
    diag_lb <- NULL
  } else {
    if (length(lb) == 1) {
      lb <- rep(lb, n)
    } else if (length(lb) != n) {
      stop("Length of 'lb' and dimensions of C do not fit.")
    }
    diag_lb <- diag(n)
  }
  if (is.null(ub)) {
    diag_ub <- NULL
  } else {
    if (length(ub) == 1) {
      ub <- rep(ub, n)
    } else if (length(ub) != n) {
      stop("Length of 'ub' and dimensions of C do not fit.")
    }
    # ATTENTION: quadprog requires -x >= -ub
    diag_ub <- -diag(n)
    ub <- -ub
  }

  Dmat <- t(C) %*% C
  dvec <- t(C) %*% d
  Dmat <- Dmat + + diag(1e-6, dim(Dmat)[1]) # ensure that Dmat is positive definite

  Amat <- rbind(Aeq, A, diag_lb, diag_ub)
  bvec <- c(beq, b, lb, ub)

  rslt <- quadprog::solve.QP(Dmat, dvec, t(Amat), bvec, meq=meq)
  rslt$solution
}


## =======================
## Modify a function so it becomes a valid (truncated + rearranged) quantile function
## =======================
modif <- function(vals, low, upp, grids){
  vals <- as.vector(vals)
  # truncation
  vals[which(vals < low)] <- low
  vals[which(vals > upp)] <- upp

  # rearrangement
  rearrangement(list(grids), vals)
}

## =======================
## Compute inner products used in simulations
## =======================
inner_vals <- function(coef_func, func_vals, grids){
  M <- length(grids)
  grids_width <- 1/M
  vals_vec_func <- rep(0, M)
  for(t in 1:M){
    x <- grids[t]
    vals_vec_point <- rep(0, M)
    for(s in 1:M){
      y <- grids[s]
      vals_vec_point[s] <- coef_func(x, y)*func_vals[s]*grids_width
    }
    vals_vec_func[t] <- sum(vals_vec_point)
  }
  return(vals_vec_func)
}


## =======================
## Generate error terms in simulations
## =======================
error_generate_unif <- function(grids, C){
  M <- length(grids)
  eps_vec <- runif(n=4, min = -C, max = C)
  error_vals <- rep(0, M)
  for(t in 1:M){
    x <- grids[t]
    error_vals[t] <- (eps_vec[1] + x^(1/2)*eps_vec[2]+ x^(1/3)*eps_vec[3] + x^(1/4)*eps_vec[4])/100
  }
  return(error_vals)
}


## =======================
## Transform quantile functions to density functions
## =======================
quant_dens <- function(quant_vals, grids, L, grid_density){
  unif_grids <- seq(0.001, 0.999, length=L)
  samples <- approx(x=grids, y=quant_vals, xout=unif_grids)
  CreateDensity(y = samples$y, optns = list(outputGrid=grid_density))$y
}


## =======================
## Basis functions used in simulations
## =======================
phi <- function(j, u){
  if(j == 1){
    1
  }
  else{
    sqrt(2)*cos((j-1)*pi*u)
  }
}


## =======================
## Generate values of a function for simulation
## =======================
X_generate <- function(grid, K){
  U <- runif(n = K, min = -3^(1/2), max = 3^(1/2))
  M <- length(grid)
  vals <- rep(0, M)
  for(m in 1:M){
    u <- grid[m]
    for(k in 1:K){
      vals[m] <- vals[m] + k^(-1.2) * U[k] * phi(k, u)
    }
  }
  return(vals/100)
}


## =======================
## Coefficient functions of autoregressive models used in simulations
## =======================
coef_func_1 <- function(x, y){
  dnorm(y, x, sd=0.1)*0.6
}
coef_func_2 <- function(x, y){
  dnorm(y, x, sd=0.1)*0.3
}
coef_func_3 <- function(x, y){
  dnorm(y, x, sd=0.1)*0.1
}


## =======================
## Generate values of a function for simulation
## =======================
mu_generate <- function(grids, ind){
  M <- length(grids)
  vals <- rep(0, M)
  for(m in 1:M){
    x <- grids[m]
    vals[m] <- phi(ind, x)
  }
  return(vals)
}
