############################################################
## Experiment 2 (Exp2): DIF Tests and Batch Execution
##
## Targeted stress-test simulation corresponding to
## Section 4 of the manuscript.
##
## This script:
##  - applies DIF detection methods to standardized residuals
##  - executes the targeted Exp2 design grid
##  - returns item-level results (res2)
##
## IMPORTANT:
## This script depends on the simulation core functions.
## Please source the following file BEFORE running:
##
##   source("experiment2/sim_core.R")
##
## Required objects include:
##   - select_item_ids()
##   - get_J_from_form()
##   - gen_2pl_par()
##   - sim_2pl_multi()
##   - fit_2pl_mirt()
############################################################

## Safety check
if (!exists("gen_2pl_par") ||
    !exists("sim_2pl_multi") ||
    !exists("fit_2pl_mirt")) {
  stop(
    "Simulation core functions not found.\n",
    "Please run:\n",
    "  source('experiment2/sim_core.R')\n",
    "before executing this script."
  )
}

############################################################
## 5) DIF Tests (GRDIF, MI, Spline-LR)
##    with tryCatch guards for robustness
############################################################

############################################################
## 5.1 Parametric residual-moment DIF tests (GRDIF-type)
############################################################
grdif_tests <- function(zj, group, theta){

  tryCatch({
    df <- data.frame(
      z  = zj,
      g  = factor(group),
      th = scale(theta, scale = FALSE)
    )

    m_full <- lm(z ~ 0 + g + g:th, df)
    m_nR   <- lm(z ~ 0 + g:th, df)
    m_nS   <- lm(z ~ 0 + g,    df)
    m_null <- lm(z ~ 1,        df)

    c(
      p_GR_R  = anova(m_nR,   m_full)$`Pr(>F)`[2],
      p_GR_S  = anova(m_nS,   m_full)$`Pr(>F)`[2],
      p_GR_RS = anova(m_null, m_full)$`Pr(>F)`[2]
    )

  }, error = function(e){
    c(p_GR_R = NA_real_,
      p_GR_S = NA_real_,
      p_GR_RS = NA_real_)
  })
}


############################################################
## 5.2 Mutual-Information DIF (stratified G^2 aggregation)
############################################################
mi_dif_item <- function(
  xj, group, theta,
  K = 10,
  L = 4,
  min_per_bin = NULL,
  continuous_unique_ratio = 0.30,
  max_unique_discrete = 20
){

  tryCatch({

    stopifnot(length(xj) == length(group),
              length(theta) == length(group))

    ok <- complete.cases(xj, group, theta)
    xj <- xj[ok]; group <- group[ok]; theta <- theta[ok]

    group <- factor(group)
    if (nlevels(group) < 2)
      return(c(p_MI = NA_real_, G2 = NA_real_, df = NA_real_))

    ## Decide whether residuals are treated as continuous
    ux <- length(unique(xj))
    n  <- length(xj)
    unique_ratio <- ux / max(n, 1)

    is_numeric <- is.numeric(xj) || is.integer(xj)
    treat_as_continuous <- is_numeric &&
      (ux > max_unique_discrete ||
       unique_ratio > continuous_unique_ratio)

    ## Discretize residuals if needed
    if (treat_as_continuous){
      qx <- unique(as.numeric(
        stats::quantile(xj,
                         probs = seq(0, 1, length.out = L + 1),
                         na.rm = TRUE)
      ))
      if (length(qx) < 3)
        return(c(p_MI = NA_real_, G2 = NA_real_, df = NA_real_))

      xj_disc <- cut(xj, breaks = qx,
                     include.lowest = TRUE, labels = FALSE)
    } else {
      xj_disc <- xj
    }

    ## Ability stratification
    br <- unique(as.numeric(
      stats::quantile(theta,
                       probs = seq(0, 1, length.out = K + 1),
                       na.rm = TRUE)
    ))
    if (length(br) < 3){
      br <- c(min(theta) - 1,
              stats::median(theta),
              max(theta) + 1)
    }

    bin <- cut(theta, breaks = br,
               include.lowest = TRUE, labels = FALSE)

    ## Stratified G^2 aggregation
    G2_total <- 0
    df_total <- 0
    G <- nlevels(group)
    if (is.null(min_per_bin)) min_per_bin <- G + 2

    for (b in sort(unique(bin))){
      idx <- which(bin == b)
      if (length(idx) < min_per_bin) next

      tab <- table(group[idx], xj_disc[idx])
      if (nrow(tab) < 2 || ncol(tab) < 2) next

      rs <- rowSums(tab)
      cs <- colSums(tab)
      N  <- sum(tab)
      E  <- rs %*% t(cs) / N

      okcell <- (tab > 0) & (E > 0)
      G2_b <- 2 * sum(tab[okcell] *
                        log(tab[okcell] / E[okcell]))
      df_b <- (nrow(tab) - 1) * (ncol(tab) - 1)

      G2_total <- G2_total + G2_b
      df_total <- df_total + df_b
    }

    if (df_total <= 0)
      return(c(p_MI = NA_real_, G2 = G2_total, df = df_total))

    c(
      p_MI = stats::pchisq(G2_total, df_total, lower.tail = FALSE),
      G2   = G2_total,
      df   = df_total
    )

  }, error = function(e){
    c(p_MI = NA_real_, G2 = NA_real_, df = NA_real_)
  })
}


############################################################
## 5.3 Spline-Based Likelihood-Ratio DIF (Spline-LR)
############################################################
spline_lr_item <- function(zj, group, theta){

  tryCatch({

    df <- data.frame(
      g  = factor(group),
      z  = zj,
      th = theta
    )

    m0 <- nnet::multinom(
      g ~ splines::ns(th, 4),
      df, trace = FALSE
    )

    m1 <- nnet::multinom(
      g ~ splines::ns(th, 4) +
           splines::ns(z,  3) +
           splines::ns(th, 2):splines::ns(z, 3),
      df, trace = FALSE
    )

    LR   <- deviance(m0) - deviance(m1)
    dfLR <- max(length(coef(m1)) - length(coef(m0)), 1)

    c(p_SLR = pchisq(LR, dfLR, lower.tail = FALSE))

  }, error = function(e){
    c(p_SLR = NA_real_)
  })
}


############################################################
## 6) Single-condition runner (Exp2)
############################################################
run_one_condition_exp2 <- function(cond){

  G <- cond$G
  test_form <- cond$test_form
  item_ids  <- select_item_ids(test_form)
  J <- get_J_from_form(test_form)

  ## Balanced design: equal n per group (1:1)
  nvec <- rep(cond$n_per_group, G)

  focal_groups <- 2:G   # fixed: all focal groups

  pars <- gen_2pl_par(
    J = J,
    DIF_type = cond$DIF_type,
    DIF_pct  = cond$DIF_pct,
    G = G,
    focal_groups = focal_groups,
    item_ids = item_ids,
    delta = cond$delta
  )

  y_true <- as.integer(1:J %in% pars$dif_idx)

  sim <- sim_2pl_multi(G, nvec, J, pars, cond$impact)
  fit <- fit_2pl_mirt(sim$X)

  out <- lapply(1:J, function(j){

    zj <- fit$Z[, j]

    tibble::tibble(
      G = G,
      test_form  = test_form,
      J = J,
      n_per_group = cond$n_per_group,
      sample_ratio = "1:1",
      DIF_type = cond$DIF_type,
      DIF_pct  = cond$DIF_pct,
      delta    = cond$delta,
      impact   = cond$impact,
      rep      = cond$rep,
      item     = j,
      item_id  = pars$item_ids[j],
      y_true   = y_true[j],
      !!!as.list(c(
        grdif_tests(zj, sim$group, fit$theta),
        mi_dif_item(zj, sim$group, fit$theta),
        spline_lr_item(zj, sim$group, fit$theta)
      ))
    )
  })

  dplyr::bind_rows(out)
}


############################################################
## 7) Design grid for Experiment 2
############################################################
make_grid_exp2 <- function(reps_per = 1){

  tibble::tibble(
    G = 8,
    focal_mode = "all",
    DIF_type = "mixed",
    DIF_pct  = 10,
    delta    = 0.5
  ) %>%
    tidyr::expand_grid(
      impact      = c(0, 1, 2, 3),
      n_per_group = c(200, 600, 2000),
      test_form   = c("full40", "target20", "hard20"),
      rep         = 1:reps_per
    ) %>%
    dplyr::mutate(
      cell_id = dplyr::row_number(),
      seed    = cell_id + 200000
    ) %>%
    dplyr::group_by(
      G, focal_mode, DIF_type, DIF_pct, delta,
      impact, n_per_group, test_form
    ) %>%
    dplyr::mutate(
      cond_id = dplyr::cur_group_id()
    ) %>%
    dplyr::ungroup()
}


############################################################
## 8) Full batch execution (Exp2)
############################################################
run_full_grid_exp2 <- function(grid){

  total <- nrow(grid)

  purrr::map_dfr(seq_len(total), function(i){

    cond <- grid[i, ]

    cat(sprintf(
      "[%d / %d] cond_id=%d | form=%s | n=%d | impact=%d | rep=%d\n",
      i, total,
      cond$cond_id,
      cond$test_form,
      cond$n_per_group,
      cond$impact,
      cond$rep
    ))
    flush.console()

    set.seed(cond$seed)

    out <- run_one_condition_exp2(as.list(cond))

    dplyr::mutate(
      out,
      cell_id = cond$cell_id,
      cond_id = cond$cond_id,
      seed    = cond$seed
    )
  })
}


############################################################
## 9) Execute (Exp2)
##
## Example:
##   source("experiment2/sim_core.R")
##   source("experiment2/dif_tests_and_run.R")
##
##   grid2 <- make_grid_exp2(reps_per = 100)
##   res2  <- run_full_grid_exp2(grid2)
##
##   saveRDS(res2, file = "res2.rds")
############################################################
