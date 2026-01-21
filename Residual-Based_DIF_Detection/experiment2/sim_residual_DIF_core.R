############################################################
## Experiment 2 (Exp2): Simulation Core
##
## Targeted stress-test design following Section 4
## of the manuscript.
##
## This script defines data-generation utilities only.
## It does NOT execute DIF tests or produce results.
##
## To be sourced by:
##   experiment2/dif_tests_and_run.R
############################################################


############################################################
## 0) Package setup
############################################################
req <- c("mirt","nnet","splines","dplyr","tidyr","purrr","tibble")

install_if_missing <- function(pkgs){
  to_install <- setdiff(pkgs, rownames(installed.packages()))
  if (length(to_install))
    install.packages(to_install, dependencies = TRUE)
  invisible(lapply(pkgs, require, character.only = TRUE))
}
install_if_missing(req)


############################################################
## Logistic link (2PL)
############################################################
ilogit <- function(x) 1 / (1 + exp(-x))


############################################################
## Reference Item Bank (Lim et al., 2024, Table A2)
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
## Test forms / item subsets (Exp2)
##
## - full40  : all 40 items
## - target20: 20 items closest to b = 0 (by |b|)
## - hard20  : 20 hard items (excluding the most extreme 3)
############################################################
select_item_ids <- function(test_form){
  test_form <- match.arg(test_form,
                         c("full40","target20","hard20"))

  if (test_form == "full40"){
    return(1:40)
  }

  if (test_form == "target20"){
    return(order(abs(b_tableA2))[1:20])
  }

  if (test_form == "hard20"){
    ord <- order(b_tableA2, decreasing = TRUE)
    return(ord[4:23])  # exclude 3 most extreme items
  }
}

get_J_from_form <- function(test_form){
  if (test_form == "full40") 40 else 20
}


############################################################
## 1) DIF parameter generator (Exp2)
##
## - Item parameters are subsets of Table A2
## - DIF magnitude fixed at delta = 0.5
## - DIF proportion controlled externally (default: 10%)
## - DIF items are the first round(J * DIF_pct / 100)
############################################################
gen_2pl_par <- function(
  J,
  DIF_type,
  DIF_pct,
  G,
  focal_groups,
  item_ids,
  delta = 0.5
){

  if (length(item_ids) != J)
    stop("length(item_ids) must equal J.")

  ## Reference-group parameters
  a_ref <- a_tableA2[item_ids]
  b_ref <- b_tableA2[item_ids]

  ## DIF item indices
  J_dif  <- round(J * DIF_pct / 100)
  dif_idx <- if (J_dif > 0) 1:J_dif else integer(0)

  ## Initialize group-specific parameters
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

  list(
    a_list   = a_list,
    b_list   = b_list,
    dif_idx  = dif_idx,
    item_ids = item_ids
  )
}


############################################################
## 2) Ability distributions (Exp2: 4 impact levels)
##
## impact:
##   0 = no impact
##   1 = moderate (±0.5)
##   2 = severe   (±1.0)
##   3 = gradient / heterogeneous
############################################################
sim_abilities <- function(G, n_per_group, impact){

  if (impact == 0){
    mu <- rep(0, G)

  } else if (impact == 1){
    mu <- c(0.5, rep(-0.5, G - 1))

  } else if (impact == 2){
    mu <- c(1.0, rep(-1.0, G - 1))

  } else if (impact == 3){
    if (G == 8){
      mu <- c(1.0, 0.7, 0.4, 0.1,
              -0.1, -0.4, -0.7, -1.0)
    } else {
      mu <- seq(1, -1, length.out = G)
    }

  } else {
    stop("impact must be one of {0, 1, 2, 3}.")
  }

  unlist(lapply(
    1:G,
    function(g) rnorm(n_per_group[g], mu[g], 1)
  ))
}


############################################################
## 3) Multi-group response generation (2PL)
############################################################
sim_2pl_multi <- function(G, n_per_group, J, pars, impact){

  theta <- sim_abilities(G, n_per_group, impact)
  group <- rep(1:G, times = n_per_group)
  N <- sum(n_per_group)

  X <- matrix(0, N, J)
  for (j in 1:J){
    for (g in 1:G){
      idx <- which(group == g)
      p <- ilogit(
        pars$a_list[[g]][j] *
        (theta[idx] - pars$b_list[[g]][j])
      )
      X[idx, j] <- rbinom(length(idx), 1, p)
    }
  }

  list(X = X, group = group, theta = theta)
}


############################################################
## 4) Pooled 2PL calibration and standardized residuals
##
## - Uses mirt with a single pooled calibration
## - Robust to occasional convergence failures
############################################################
fit_2pl_mirt <- function(X){

  dat <- as.data.frame(X)
  colnames(dat) <- paste0("I", 1:ncol(X))

  fit <- tryCatch(
    mirt(dat, 1, itemtype = "2PL", verbose = FALSE),
    error = function(e) NULL
  )

  if (is.null(fit)){
    N <- nrow(X); J <- ncol(X)
    return(list(
      theta = rep(NA_real_, N),
      Z     = matrix(NA_real_, N, J)
    ))
  }

  fs <- tryCatch(
    fscores(fit, method = "EAP"),
    error = function(e) NULL
  )

  if (is.null(fs)){
    N <- nrow(theta); J <- ncol(X)
    return(list(
      theta = rep(NA_real_, N),
      Z     = matrix(NA_real_, N, J)
    ))
  }

  pars <- coef(fit, IRTpars = TRUE, simplify = TRUE)$items
  a <- pars[, "a"]
  b <- pars[, "b"]
  theta <- as.numeric(fs)

  P_hat <- outer(
    theta, seq_along(a),
    function(th, j) ilogit(a[j] * (th - b[j]))
  )

  Z <- (X - P_hat) /
       sqrt(pmax(P_hat * (1 - P_hat), 1e-9))

  list(theta = theta, Z = Z)
}


############################################################
## End of Experiment 2 simulation core
############################################################
