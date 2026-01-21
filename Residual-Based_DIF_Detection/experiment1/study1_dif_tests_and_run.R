############################################################
## Study 1: DIF Testing and Batch Execution
##
## This script implements item-level DIF detection procedures
## and executes the full factorial design for Simulation Study 1.
##
## NOTE:
## This script depends on the simulation core functions.
## Please source the following file before running:
##
##   source("sim_residual_DIF_core.R")
############################################################




############################################################
## 6) DIF Tests
############################################################

############################################################
## 6.1 Parametric residual-moment DIF tests (GRDIF-type)
##
## These tests evaluate residual–group dependence using
## linear contrasts of residual moments, conditional on
## estimated ability.
############################################################
grdif_tests <- function(zj, group, theta){
  
  df <- data.frame(
    z  = zj,
    g  = factor(group),
    th = scale(theta, scale = FALSE)   # center ability
  )
  
  ## Full model: group main effects + group-by-ability interaction
  m_full <- lm(z ~ 0 + g + g:th, df)
  
  ## Reduced models
  m_nR   <- lm(z ~ 0 + g:th, df)   # remove group mean differences
  m_nS   <- lm(z ~ 0 + g,    df)   # remove interaction structure
  m_null <- lm(z ~ 1,        df)   # no group effects
  
  c(
    p_GR_R  = anova(m_nR,   m_full)$`Pr(>F)`[2],
    p_GR_S  = anova(m_nS,   m_full)$`Pr(>F)`[2],
    p_GR_RS = anova(m_null, m_full)$`Pr(>F)`[2]
  )
}


############################################################
## 6.2 Mutual-Information DIF (Stratified G^2 aggregation)
##
## Implements a conditional mutual-information test by:
## (a) discretizing ability into K strata,
## (b) discretizing residuals into L bins (if continuous),
## (c) aggregating likelihood-ratio statistics across strata.
##
## This provides a fully nonparametric assessment of
## residual–group dependence.
############################################################
mi_dif_item <- function(
    xj, group, theta,
    K = 10,                    # ability strata
    L = 4,                     # residual bins
    min_per_bin = NULL,
    continuous_unique_ratio = 0.30,
    max_unique_discrete = 20
){
  
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
    (ux > max_unique_discrete || unique_ratio > continuous_unique_ratio)
  
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
  if (length(br) < 3)
    br <- c(min(theta) - 1,
            stats::median(theta),
            max(theta) + 1)
  
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
}


############################################################
## 6.3 Spline-Based Likelihood-Ratio DIF (Spline-LR)
##
## Tests whether allowing residual-dependent smooth terms
## improves the prediction of group membership beyond ability.
############################################################
spline_lr_item <- function(zj, group, theta){
  
  df <- data.frame(
    g  = factor(group),
    z  = zj,
    th = theta
  )
  
  ## Baseline model: ability only
  m0 <- nnet::multinom(
    g ~ splines::ns(th, 4),
    df, trace = FALSE
  )
  
  ## Augmented model: residual smooth + interaction
  m1 <- nnet::multinom(
    g ~ splines::ns(th, 4) +
      splines::ns(z,  3) +
      splines::ns(th, 2):splines::ns(z, 3),
    df, trace = FALSE
  )
  
  LR   <- deviance(m0) - deviance(m1)
  dfLR <- max(length(coef(m1)) - length(coef(m0)), 1)
  
  c(p_SLR = pchisq(LR, dfLR, lower.tail = FALSE))
}


############################################################
## 7) Single-condition runner (Study 1)
##
## Executes one design condition × one replication:
## - simulate data (via sourced functions),
## - fit pooled 2PL,
## - compute residuals,
## - apply all DIF tests item-wise.
############################################################
run_one_condition <- function(cond){
  
  G  <- cond$G
  J  <- 40
  nvec <- decode_sizes(G, cond$n_ref, cond$n_focal)
  
  focal_groups <- if (cond$focal_mode == "one") 2 else 2:G
  
  pars <- gen_2pl_par(
    J = J,
    DIF_type = cond$DIF_type,
    DIF_pct  = cond$DIF_pct,
    G = G,
    focal_groups = focal_groups,
    item_mode = "tableA2"
  )
  
  y_true <- as.integer(1:J %in% pars$dif_idx)
  
  sim <- sim_2pl_multi(G, nvec, J, pars, cond$impact)
  fit <- fit_2pl_mirt(sim$X)
  
  out <- lapply(1:J, function(j){
    
    zj <- fit$Z[, j]
    xj <- sim$X[, j]
    
    tibble::tibble(
      name       = cond$name,
      G          = G,
      focal_mode = cond$focal_mode,
      n_ref      = cond$n_ref,
      n_focal    = cond$n_focal,
      DIF_type   = cond$DIF_type,
      DIF_pct    = cond$DIF_pct,
      impact     = cond$impact,
      rep        = cond$rep,
      item       = j,
      y_true     = y_true[j],
      !!!as.list(
        c(
          grdif_tests(zj, sim$group, fit$theta),
          mi_dif_item(zj, sim$group, fit$theta),
          spline_lr_item(zj, sim$group, fit$theta)
        )
      )
    )
  })
  
  dplyr::bind_rows(out)
}


############################################################
## 8) Design grid for Study 1
############################################################
make_grid_study1 <- function(reps_per = 1){
  
  grid <- sample_scheme_12 %>%
    tidyr::expand_grid(
      focal_mode = c("one", "all"),
      DIF_type   = c("uniform", "nonuniform", "mixed"),
      DIF_pct    = c(0, 10, 20),
      impact     = c(0, 1),
      rep        = 1:reps_per
    ) %>%
    dplyr::mutate(
      cell_id = dplyr::row_number(),
      seed    = cell_id + 100000
    ) %>%
    dplyr::group_by(
      name, G, n_ref, n_focal,
      focal_mode, DIF_type, DIF_pct, impact
    ) %>%
    dplyr::mutate(
      cond_id = dplyr::cur_group_id()
    ) %>%
    dplyr::ungroup()
  
  grid
}


############################################################
## 9) Full batch execution (Study 1)
############################################################
run_full_grid_study1 <- function(grid){
  
  total <- nrow(grid)
  
  purrr::map_dfr(seq_len(total), function(i){
    
    cond <- grid[i, ]
    
    cat(sprintf(
      "[%d / %d] cond_id=%d | scheme=%s | focal=%s | DIF=%s%d | impact=%d | rep=%d\n",
      i, total,
      cond$cond_id,
      cond$name,
      cond$focal_mode,
      cond$DIF_type,
      cond$DIF_pct,
      cond$impact,
      cond$rep
    ))
    flush.console()
    
    set.seed(cond$seed)
    
    out <- run_one_condition(as.list(cond))
    
    dplyr::mutate(
      out,
      cell_id = cond$cell_id,
      cond_id = cond$cond_id,
      seed    = cond$seed
    )
  })
}


############################################################
## 10) Execute (WARNING: long runtime)
############################################################
## Example:
##   source("sim_residual_DIF_core.R")
##   source("study1_dif_tests_and_run.R")
##
##   grid <- make_grid_study1(reps_per = 100)
##   res  <- run_full_grid_study1(grid)
##   save(res, file = "study1_results.RData")
############################################################
