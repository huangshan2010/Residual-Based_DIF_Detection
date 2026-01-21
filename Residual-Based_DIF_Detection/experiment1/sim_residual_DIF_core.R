############################################################
## Residual-Based DIF Simulation Core Script
## Benchmark design following Lim et al. (2024)
## Used for Simulation Study 1 
############################################################

############################################################
## 0) Setup
############################################################

## Required packages
req <- c("mirt","nnet","splines","dplyr","tidyr","purrr")

install_if_missing <- function(pkgs){
  to_install <- setdiff(pkgs, rownames(installed.packages()))
  if (length(to_install)) install.packages(to_install, dependencies = TRUE)
  invisible(lapply(pkgs, require, character.only = TRUE))
}
install_if_missing(req)

## Logistic link (2PL)
ilogit <- function(x) 1 / (1 + exp(-x))


############################################################
## 1) Reference Item Bank (Lim et al., 2024, Table A2)
############################################################

a_tableA2 <- c(
  1.67,2.18,1.95,1.50,2.77,1.09,1.57,1.27,1.29,2.09,
  1.27,1.60,1.47,1.72,1.35,3.21,1.57,1.66,1.33,1.13,
  1.74,1.52,1.30,1.79,2.18,1.97,2.70,1.84,1.21,1.81,
  1.06,1.21,2.15,1.28,2.02,1.42,2.27,1.36,1.86,1.58
)

b_tableA2 <- c(
  -1.10, .79,-1.33, .30,-.69, .28,-.11,1.32, .64,-1.05,
  .57,-.24,-.51, .77,-.08,-.32,-.18, .22,2.15,-1.92,
  .93,-.51, .50,-.05,-.42,1.55,-.52, .80,-1.62,-.22,
  .34,1.35,-.71,-.52, .04,-.20,-2.09,-.38,-.55, .55
)


############################################################
## 2) Sample-size schemes (Lim et al., Table 1)
############################################################

sample_scheme_12 <- tibble::tribble(
  ~name,         ~G, ~n_ref, ~n_focal,
  "UG3CN400",     3,   200,      100,
  "BG3CN450",     3,   150,      150,
  "UG5CN600",     5,   200,      100,
  "BG5CN750",     5,   150,      150,
  "UG8CN900",     8,   200,      100,
  "BG8CN1200",    8,   150,      150,
  "UG3CN800",     3,   400,      200,
  "BG3CN900",     3,   300,      300,
  "UG5CN1200",    5,   400,      200,
  "BG5CN1500",    5,   300,      300,
  "UG8CN1800",    8,   400,      200,
  "BG8CN2400",    8,   300,      300
)

decode_sizes <- function(G, n_ref, n_focal){
  c(n_ref, rep(n_focal, G - 1))
}


############################################################
## 3) DIF Parameter Generator (2PL)
############################################################

gen_2pl_par <- function(
  J = 40,
  DIF_type,
  DIF_pct,
  G,
  focal_groups,
  item_mode = c("tableA2","random"),
  delta_set = c(.3,.5,.7,.9)
){
  item_mode <- match.arg(item_mode)
  delta <- sample(delta_set, 1)

  ## Reference-group parameters
  if (item_mode == "tableA2") {
    if (J != 40) stop("Table A2 supports only J = 40.")
    a_ref <- a_tableA2
    b_ref <- b_tableA2
  } else {
    a_ref <- pmax(rnorm(J, 1, 0.2), .4)
    b_ref <- runif(J, -2, 2)
  }

  ## DIF item indices
  J_dif  <- round(J * DIF_pct / 100)
  dif_idx <- if (J_dif > 0) 1:J_dif else integer(0)

  ## Initialize group parameters
  a_list <- vector("list", G)
  b_list <- vector("list", G)
  for (g in 1:G){
    a_list[[g]] <- a_ref
    b_list[[g]] <- b_ref
  }

  ## Apply DIF to focal groups
  if (J_dif > 0){
    for (g in focal_groups){
      if (DIF_type %in% c("uniform","mixed"))
        b_list[[g]][dif_idx] <- b_ref[dif_idx] + delta
      if (DIF_type %in% c("nonuniform","mixed"))
        a_list[[g]][dif_idx] <- pmax(a_ref[dif_idx] - delta, .3)
    }
  }

  list(a_list = a_list,
       b_list = b_list,
       dif_idx = dif_idx)
}


############################################################
## 4) Ability distributions (impact manipulation)
############################################################

sim_abilities <- function(G, n_per_group, impact){
  if (impact == 0)
    mu <- rep(0, G)
  else
    mu <- c(0.5, rep(-0.5, G - 1))

  unlist(lapply(1:G, function(g)
    rnorm(n_per_group[g], mu[g], 1)))
}


############################################################
## 5) Multi-group response generation (2PL)
############################################################

sim_2pl_multi <- function(G, n_per_group, J, pars, impact){

  theta <- sim_abilities(G, n_per_group, impact)
  group <- rep(1:G, times = n_per_group)
  N <- sum(n_per_group)

  X <- matrix(0, N, J)
  for (j in 1:J){
    for (g in 1:G){
      idx <- which(group == g)
      p <- ilogit(pars$a_list[[g]][j] *
                   (theta[idx] - pars$b_list[[g]][j]))
      X[idx, j] <- rbinom(length(idx), 1, p)
    }
  }

  list(X = X, group = group, theta = theta)
}


############################################################
## 6) Pooled 2PL calibration and standardized residuals
############################################################

fit_2pl_mirt <- function(X){

  dat <- as.data.frame(X)
  colnames(dat) <- paste0("I", 1:ncol(X))

  fit <- mirt(dat, 1, itemtype = "2PL", verbose = FALSE)
  fs  <- fscores(fit, method = "EAP")

  pars <- coef(fit, IRTpars = TRUE, simplify = TRUE)$items
  a <- pars[, "a"]
  b <- pars[, "b"]
  theta <- as.numeric(fs)

  P_hat <- outer(theta, seq_along(a),
                 function(th, j) ilogit(a[j] * (th - b[j])))

  Z <- (X - P_hat) / sqrt(pmax(P_hat * (1 - P_hat), 1e-9))

  list(theta = theta,
       Z = Z)
}

############################################################
## End of core simulation script
############################################################
