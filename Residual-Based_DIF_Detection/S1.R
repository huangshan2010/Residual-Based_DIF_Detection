############################################################
## Residual-Based DIF Simulation Script
## Aligned with Lim et al. (2022)
##
## This script contains:
##   (1) data generation functions
##   (2) pooled IRT calibration and raw-residual extraction
##   (3) official RDIF analysis via irtQ::rdif()
##   (4) semi-parametric residual regression procedure
##   (5) MI-based residual contingency procedure
##   (6) purification wrappers
##   (7) single-condition and batch execution functions
##   (8) resume / checkpoint utilities
##
## Notes
## - The design is two-group only: reference vs focal.
## - RDIF uses the official irtQ implementation with ML scoring
##   and purification via purify.by = "rdifrs".
## - Semi and MI procedures are based on raw residuals:
##       r = x - P_hat
## - Person scoring is kept as ML.
## - ML finiteness is recorded rather than used as a stopping rule.
############################################################

############################################################
## 0) Setup
############################################################

req <- c("irtQ", "dplyr", "purrr", "tidyr", "splines", "tibble")

install_if_missing <- function(pkgs) {
  to_install <- setdiff(pkgs, rownames(installed.packages()))
  if (length(to_install)) install.packages(to_install, dependencies = TRUE)
  invisible(lapply(pkgs, require, character.only = TRUE))
}
install_if_missing(req)

`%||%` <- function(x, y) if (is.null(x)) y else x

############################################################
## 1) Small helpers
############################################################

ilogit <- function(x) 1 / (1 + exp(-x))

rtrunc_norm <- function(n, mean = 0, sd = 1, lower = -Inf, upper = Inf) {
  out <- numeric(0)
  while (length(out) < n) {
    cand <- rnorm(max(100, n - length(out)), mean = mean, sd = sd)
    cand <- cand[cand >= lower & cand <= upper]
    out <- c(out, cand)
  }
  out[seq_len(n)]
}

rtrunc_lnorm <- function(n, meanlog, sdlog, lower = 0, upper = Inf) {
  out <- numeric(0)
  while (length(out) < n) {
    cand <- rlnorm(max(100, n - length(out)), meanlog = meanlog, sdlog = sdlog)
    cand <- cand[cand >= lower & cand <= upper]
    out <- c(out, cand)
  }
  out[seq_len(n)]
}

rtrunc_beta <- function(n, shape1, shape2, lower = 0, upper = 1) {
  out <- numeric(0)
  while (length(out) < n) {
    cand <- rbeta(max(100, n - length(out)), shape1 = shape1, shape2 = shape2)
    cand <- cand[cand >= lower & cand <= upper]
    out <- c(out, cand)
  }
  out[seq_len(n)]
}

append_log <- function(msg, logfile = "simulation_log.txt") {
  cat(sprintf("[%s] %s\n", as.character(Sys.time()), msg),
      file = logfile, append = TRUE)
}

############################################################
## 2) Design factors
############################################################

sample_scheme <- tibble::tribble(
  ~scheme,       ~n_ref, ~n_focal,
  "250R_250F",     250,      250,
  "400R_100F",     400,      100,
  "1000R_1000F",  1000,     1000,
  "1600R_400F",   1600,      400
)

make_design <- function(reps_per = 100) {
  base <- tidyr::expand_grid(
    model      = c("2PLM", "3PLM"),
    scheme_tbl = seq_len(nrow(sample_scheme)),
    DIF_pct    = c(0, 10, 20, 40),
    DIF_type   = c("uniform", "nonuniform", "mixed"),
    impact     = c(0, 0.5),
    rep        = seq_len(reps_per)
  ) |>
    dplyr::mutate(
      scheme   = sample_scheme$scheme[scheme_tbl],
      n_ref    = sample_scheme$n_ref[scheme_tbl],
      n_focal  = sample_scheme$n_focal[scheme_tbl]
    ) |>
    dplyr::filter(!(DIF_pct == 0 & DIF_type != "uniform")) |>
    dplyr::select(-scheme_tbl)
  
  base |>
    dplyr::mutate(
      cond_id = dplyr::row_number(),
      seed    = 100000 + cond_id
    )
}

############################################################
## 3) Item generation
############################################################

## Lim et al. (2022):
## a ~ LN(.5, .3^2) truncated to [1.0, 3.5]
## b ~ N(0, 1) truncated to [-2.5, 2.5]
## g ~ Beta(5, 17) truncated to [0, .35] for 3PLM

gen_reference_item_bank <- function(J = 40, model = c("2PLM", "3PLM")) {
  model <- match.arg(model)
  
  a <- rtrunc_lnorm(J, meanlog = 0.5, sdlog = 0.3, lower = 1.0, upper = 3.5)
  b <- rtrunc_norm(J, mean = 0, sd = 1, lower = -2.5, upper = 2.5)
  g <- if (model == "3PLM") {
    rtrunc_beta(J, shape1 = 5, shape2 = 17, lower = 0, upper = 0.35)
  } else {
    rep(0, J)
  }
  
  list(a = a, b = b, g = g)
}

apply_dif <- function(ref_bank, DIF_type = c("uniform", "nonuniform", "mixed"),
                      DIF_pct = 0) {
  DIF_type <- match.arg(DIF_type)
  J <- length(ref_bank$a)
  j_dif <- round(J * DIF_pct / 100)
  dif_idx <- if (j_dif > 0) seq_len(j_dif) else integer(0)
  
  delta_cycle <- rep(c(0.3, 0.5, 0.7, 0.9), length.out = j_dif)
  
  focal_bank <- ref_bank
  
  if (j_dif > 0) {
    if (DIF_type %in% c("uniform", "mixed")) {
      focal_bank$b[dif_idx] <- ref_bank$b[dif_idx] + delta_cycle
    }
    if (DIF_type %in% c("nonuniform", "mixed")) {
      focal_bank$a[dif_idx] <- pmax(ref_bank$a[dif_idx] - delta_cycle, 0.05)
    }
  }
  
  list(ref = ref_bank, focal = focal_bank, dif_idx = dif_idx, delta = delta_cycle)
}

shape_item_df <- function(bank, model = c("2PLM", "3PLM")) {
  model <- match.arg(model)
  
  g_par <- if (!is.null(bank$g)) bank$g else rep(0, length(bank$a))
  
  irtQ::shape_df(
    par.drm = list(a = bank$a, b = bank$b, g = g_par),
    item.id = paste0("I", seq_along(bank$a)),
    cats = 2,
    model = model
  )
}

############################################################
## 4) Ability generation and response simulation
############################################################

sim_abilities <- function(n_ref, n_focal, impact = 0) {
  theta_ref   <- rnorm(n_ref, mean = 0, sd = 1)
  theta_focal <- rnorm(n_focal, mean = -impact, sd = 1)
  list(theta_ref = theta_ref, theta_focal = theta_focal)
}

simulate_two_group_data <- function(n_ref, n_focal, model = c("2PLM", "3PLM"),
                                    DIF_type = c("uniform", "nonuniform", "mixed"),
                                    DIF_pct = 0, impact = 0, D = 1) {
  model <- match.arg(model)
  DIF_type <- match.arg(DIF_type)
  
  ref_bank <- gen_reference_item_bank(J = 40, model = model)
  banks <- apply_dif(ref_bank, DIF_type = DIF_type, DIF_pct = DIF_pct)
  
  pars_ref   <- shape_item_df(banks$ref, model = model)
  pars_focal <- shape_item_df(banks$focal, model = model)
  
  th <- sim_abilities(n_ref = n_ref, n_focal = n_focal, impact = impact)
  resp_ref   <- irtQ::simdat(pars_ref, theta = th$theta_ref, D = D)
  resp_focal <- irtQ::simdat(pars_focal, theta = th$theta_focal, D = D)
  
  data  <- rbind(resp_ref, resp_focal)
  group <- c(rep(0, n_ref), rep(1, n_focal))
  
  list(
    data       = data,
    group      = group,
    dif_idx    = banks$dif_idx,
    ref_bank   = banks$ref,
    focal_bank = banks$focal,
    true_theta = c(th$theta_ref, th$theta_focal)
  )
}

############################################################
## 5) Pooled calibration + robust ML scoring via irtQ
############################################################

clean_theta_ml <- function(theta, data, range = c(-4, 4)) {
  raw_score <- rowSums(data, na.rm = TRUE)
  J <- ncol(data)
  
  theta_raw <- theta
  
  if (is.null(theta_raw) || length(theta_raw) != nrow(data)) {
    theta_raw <- rep(NA_real_, nrow(data))
  }
  
  theta_clean <- theta_raw
  
  ## 1. Infinite theta: cap to scoring range
  theta_clean[is.infinite(theta_clean) & theta_clean > 0] <- range[2]
  theta_clean[is.infinite(theta_clean) & theta_clean < 0] <- range[1]
  
  ## 2. NA theta for perfect / zero scores
  theta_clean[is.na(theta_clean) & raw_score == J] <- range[2]
  theta_clean[is.na(theta_clean) & raw_score == 0] <- range[1]
  
  ## 3. Remaining NA: replace by median of finite theta
  med_theta <- stats::median(theta_clean[is.finite(theta_clean)], na.rm = TRUE)
  if (!is.finite(med_theta)) med_theta <- 0
  
  theta_clean[is.na(theta_clean)] <- med_theta
  
  ## 4. Final cap
  theta_clean <- pmin(pmax(theta_clean, range[1]), range[2])
  
  list(
    theta_raw = theta_raw,
    theta_clean = theta_clean,
    
    n_na_theta_raw = sum(is.na(theta_raw)),
    n_infinite_theta_raw = sum(is.infinite(theta_raw), na.rm = TRUE),
    n_nonfinite_theta_raw = sum(!is.finite(theta_raw)),
    
    n_na_theta_clean = sum(is.na(theta_clean)),
    n_infinite_theta_clean = sum(is.infinite(theta_clean), na.rm = TRUE),
    n_nonfinite_theta_clean = sum(!is.finite(theta_clean)),
    
    all_theta_finite_raw = all(is.finite(theta_raw)),
    all_theta_finite_clean = all(is.finite(theta_clean))
  )
}


clean_p_value <- function(p) {
  p[is.na(p) | !is.finite(p)] <- 1
  p
}


clean_test_result <- function(x) {
  if (is.null(x) || length(x) == 0) {
    return(c(stat = NA_real_, df = NA_real_, p = 1))
  }
  
  if (is.null(x[["stat"]])) x[["stat"]] <- NA_real_
  if (is.null(x[["df"]]))   x[["df"]]   <- NA_real_
  if (is.null(x[["p"]]))    x[["p"]]    <- 1
  
  if (is.na(x[["p"]]) || !is.finite(x[["p"]])) {
    x[["p"]] <- 1
  }
  
  if (is.na(x[["stat"]]) || !is.finite(x[["stat"]])) {
    x[["stat"]] <- NA_real_
  }
  
  if (is.na(x[["df"]]) || !is.finite(x[["df"]])) {
    x[["df"]] <- NA_real_
  }
  
  c(
    stat = as.numeric(x[["stat"]]),
    df   = as.numeric(x[["df"]]),
    p    = as.numeric(x[["p"]])
  )
}


fit_pooled_irtq <- function(data, model = c("2PLM", "3PLM"), D = 1,
                            score_range = c(-4, 4),
                            item_skip = integer(0),
                            logfile = "simulation_log.txt") {
  model <- match.arg(model)
  
  J <- ncol(data)
  item_skip <- sort(unique(item_skip))
  item_skip <- item_skip[item_skip >= 1 & item_skip <= J]
  anchor_items <- setdiff(seq_len(J), item_skip)
  
  ## Safety fallback: if too few anchor items remain, use all items
  if (length(anchor_items) < 3) {
    append_log(
      sprintf(
        "Too few anchor items after purification: n_anchor=%s. Reverting to all items.",
        length(anchor_items)
      ),
      logfile = logfile
    )
    anchor_items <- seq_len(J)
    item_skip <- integer(0)
  }
  
  ############################################################
  ## 1. Estimate pooled item parameters using all items
  ##    These parameters are used to compute residuals for all items.
  ############################################################
  
  est_mod <- irtQ::est_irt(data = data, D = D, model = model)
  est_par <- est_mod$par.est
  
  ############################################################
  ## 2. Estimate ML theta using anchor items only
  ############################################################
  
  est_par_anchor <- est_par[anchor_items, , drop = FALSE]
  data_anchor <- data[, anchor_items, drop = FALSE]
  
  score_obj <- tryCatch(
    irtQ::est_score(
      x = est_par_anchor,
      data = data_anchor,
      method = "ML",
      range = score_range
    ),
    error = function(e) {
      append_log(sprintf("ML scoring error: %s", e$message), logfile = logfile)
      return(NULL)
    }
  )
  
  if (is.null(score_obj)) {
    score_raw <- rep(NA_real_, nrow(data))
    append_log(
      "ML scoring completely failed; raw theta set to NA before robust cleaning.",
      logfile = logfile
    )
  } else {
    score_raw <- score_obj$est.theta
    
    if (is.null(score_raw) || length(score_raw) != nrow(data)) {
      score_raw <- rep(NA_real_, nrow(data))
      append_log(
        "ML scoring returned invalid theta length; raw theta set to NA before robust cleaning.",
        logfile = logfile
      )
    }
  }
  
  theta_info <- clean_theta_ml(
    theta = score_raw,
    data = data_anchor,
    range = score_range
  )
  
  if (!theta_info$all_theta_finite_raw) {
    append_log(
      sprintf(
        paste0(
          "ML scoring produced non-finite raw theta values: ",
          "n_na_raw=%s, n_infinite_raw=%s, n_nonfinite_raw=%s, ",
          "n_na_clean=%s, n_infinite_clean=%s, n_nonfinite_clean=%s, N=%s"
        ),
        theta_info$n_na_theta_raw,
        theta_info$n_infinite_theta_raw,
        theta_info$n_nonfinite_theta_raw,
        theta_info$n_na_theta_clean,
        theta_info$n_infinite_theta_clean,
        theta_info$n_nonfinite_theta_clean,
        nrow(data)
      ),
      logfile = logfile
    )
  }
  
  list(
    est_mod = est_mod,
    est_par = est_par,
    
    score = theta_info$theta_clean,
    score_raw = theta_info$theta_raw,
    
    item_skip = item_skip,
    anchor_items = anchor_items,
    n_anchor_items = length(anchor_items),
    
    n_na_theta = theta_info$n_na_theta_raw,
    n_infinite_theta = theta_info$n_infinite_theta_raw,
    n_nonfinite_theta = theta_info$n_nonfinite_theta_raw,
    all_theta_finite = theta_info$all_theta_finite_raw,
    
    n_na_theta_clean = theta_info$n_na_theta_clean,
    n_infinite_theta_clean = theta_info$n_infinite_theta_clean,
    n_nonfinite_theta_clean = theta_info$n_nonfinite_theta_clean,
    all_theta_finite_clean = theta_info$all_theta_finite_clean
  )
}


extract_item_pars <- function(est_par) {
  out <- est_par[, c("id", "model", "par.1", "par.2", "par.3")]
  names(out) <- c("id", "model", "a", "b", "g")
  out$g[is.na(out$g)] <- 0
  out
}


compute_raw_residuals <- function(data, est_par, score, D = 1) {
  pars <- extract_item_pars(est_par)
  J <- nrow(pars)
  N <- nrow(data)
  
  P_hat <- matrix(NA_real_, nrow = N, ncol = J)
  
  for (j in seq_len(J)) {
    eta <- D * pars$a[j] * (score - pars$b[j])
    P_hat[, j] <- pars$g[j] + (1 - pars$g[j]) * ilogit(eta)
  }
  
  resid_raw <- data - P_hat
  colnames(resid_raw) <- pars$id
  
  list(P_hat = P_hat, resid_raw = resid_raw)
}


############################################################
## 6) Item-level DIF tests
############################################################

spline_lr_item <- function(rj, group, theta, focal_name = 1) {
  dat <- data.frame(
    g  = group,
    r  = rj,
    th = theta
  )
  
  dat <- dat[complete.cases(dat), ]
  
  if (nrow(dat) == 0 || length(unique(dat$g)) < 2) {
    return(c(stat = NA_real_, df = NA_real_, p = 1))
  }
  
  dat$g <- as.numeric(dat$g == focal_name)
  
  if (length(unique(dat$g)) < 2) {
    return(c(stat = NA_real_, df = NA_real_, p = 1))
  }
  
  m0 <- tryCatch(
    stats::glm(
      g ~ splines::ns(th, df = 4),
      data = dat,
      family = stats::binomial(link = "logit")
    ),
    error = function(e) NULL
  )
  
  m1 <- tryCatch(
    stats::glm(
      g ~ splines::ns(th, df = 4) +
        splines::ns(r, df = 3) +
        splines::ns(th, df = 2):splines::ns(r, df = 3),
      data = dat,
      family = stats::binomial(link = "logit")
    ),
    error = function(e) NULL
  )
  
  if (is.null(m0) || is.null(m1)) {
    return(c(stat = NA_real_, df = NA_real_, p = 1))
  }
  
  an <- tryCatch(
    stats::anova(m0, m1, test = "Chisq"),
    error = function(e) NULL
  )
  
  if (is.null(an)) {
    return(c(stat = NA_real_, df = NA_real_, p = 1))
  }
  
  out <- c(
    stat = unname(an$Deviance[2]),
    df   = unname(an$Df[2]),
    p    = unname(an$`Pr(>Chi)`[2])
  )
  
  clean_test_result(out)
}


mi_dif_item <- function(rj, group, theta,
                        K = 5,
                        L = 4,
                        min_per_bin = NULL) {
  
  ok <- complete.cases(rj, group, theta)
  rj <- rj[ok]
  group <- group[ok]
  theta <- theta[ok]
  
  group <- factor(group)
  
  if (length(rj) == 0 || nlevels(group) < 2) {
    return(c(stat = NA_real_, df = NA_real_, p = 1))
  }
  
  G <- nlevels(group)
  if (is.null(min_per_bin)) min_per_bin <- G + 2
  
  qx <- unique(as.numeric(
    stats::quantile(
      rj,
      probs = seq(0, 1, length.out = L + 1),
      na.rm = TRUE
    )
  ))
  
  if (length(qx) < 3) {
    return(c(stat = NA_real_, df = NA_real_, p = 1))
  }
  
  r_disc <- cut(
    rj,
    breaks = qx,
    include.lowest = TRUE,
    labels = FALSE
  )
  
  br <- unique(as.numeric(
    stats::quantile(
      theta,
      probs = seq(0, 1, length.out = K + 1),
      na.rm = TRUE
    )
  ))
  
  if (length(br) < 3) {
    br <- c(min(theta) - 1, stats::median(theta), max(theta) + 1)
  }
  
  bin <- cut(
    theta,
    breaks = br,
    include.lowest = TRUE,
    labels = FALSE
  )
  
  G2_total <- 0
  df_total <- 0
  
  for (b in sort(unique(bin))) {
    idx <- which(bin == b)
    
    if (length(idx) < min_per_bin) next
    
    tab <- table(droplevels(group[idx]), r_disc[idx])
    tab <- tab[rowSums(tab) > 0, colSums(tab) > 0, drop = FALSE]
    
    if (nrow(tab) < 2 || ncol(tab) < 2) next
    
    rs <- rowSums(tab)
    cs <- colSums(tab)
    N  <- sum(tab)
    
    E <- rs %*% t(cs) / N
    
    okcell <- (tab > 0) & (E > 0)
    
    G2_b <- 2 * sum(tab[okcell] * log(tab[okcell] / E[okcell]))
    df_b <- (nrow(tab) - 1) * (ncol(tab) - 1)
    
    G2_total <- G2_total + G2_b
    df_total <- df_total + df_b
  }
  
  if (df_total <= 0) {
    return(c(stat = G2_total, df = df_total, p = 1))
  }
  
  out <- c(
    stat = G2_total,
    df   = df_total,
    p    = stats::pchisq(G2_total, df_total, lower.tail = FALSE)
  )
  
  clean_test_result(out)
}


############################################################
## 7) Common purification wrapper
############################################################

flag_items_by_p <- function(pvals, alpha = 0.05) {
  pvals <- clean_p_value(pvals)
  which(pvals < alpha)
}


run_custom_purification <- function(data, group, model = c("2PLM", "3PLM"),
                                    test_fun,
                                    alpha = 0.05,
                                    purify = TRUE,
                                    max.iter = 20,
                                    D = 1,
                                    score_range = c(-4, 4),
                                    min.resp = NULL,
                                    logfile = "simulation_log.txt") {
  model <- match.arg(model)
  J <- ncol(data)
  
  item_skip <- integer(0)
  iter <- 1
  converged <- TRUE
  
  repeat {
    
    ############################################################
    ## 1. Fit pooled IRT model and estimate theta
    ##    If item_skip is not empty, theta is estimated using
    ##    anchor items only.
    ############################################################
    
    fit <- fit_pooled_irtq(
      data = data,
      model = model,
      D = D,
      score_range = score_range,
      item_skip = item_skip,
      logfile = logfile
    )
    
    score <- fit$score
    
    if (!is.null(min.resp)) {
      valid_n <- rowSums(!is.na(data[, fit$anchor_items, drop = FALSE]))
      score[valid_n < min.resp] <- NA_real_
    }
    
    if (any(!is.finite(score))) {
      score <- clean_theta_ml(
        theta = score,
        data = data[, fit$anchor_items, drop = FALSE],
        range = score_range
      )$theta_clean
    }
    
    ############################################################
    ## 2. Compute residuals for all items
    ############################################################
    
    raw <- compute_raw_residuals(
      data = data,
      est_par = fit$est_par,
      score = score,
      D = D
    )$resid_raw
    
    ############################################################
    ## 3. Test all items
    ##    Important: purification excludes suspected DIF items
    ##    from scoring, but all items are still tested.
    ############################################################
    
    pvals  <- rep(1, J)
    stats1 <- rep(NA_real_, J)
    dfs    <- rep(NA_real_, J)
    
    for (j in seq_len(J)) {
      tmp <- tryCatch(
        clean_test_result(test_fun(raw[, j], group, score)),
        error = function(e) {
          append_log(
            sprintf("Item-level test failed at item=%s: %s", j, e$message),
            logfile = logfile
          )
          c(stat = NA_real_, df = NA_real_, p = 1)
        }
      )
      
      stats1[j] <- tmp[["stat"]]
      dfs[j]    <- tmp[["df"]]
      pvals[j]  <- tmp[["p"]]
    }
    
    pvals <- clean_p_value(pvals)
    
    ############################################################
    ## 4. Identify suspected DIF items for the next iteration
    ##    For Semi: semi p-values are used.
    ##    For MI:   MI p-values are used.
    ############################################################
    
    flagged <- which(pvals < alpha)
    
    if (!purify) {
      break
    }
    
    new_skip <- sort(unique(flagged))
    
    ############################################################
    ## 5. Stop if the purification set is stable
    ############################################################
    
    if (identical(new_skip, sort(unique(item_skip)))) {
      break
    }
    
    item_skip <- new_skip
    iter <- iter + 1
    
    if (iter > max.iter) {
      converged <- FALSE
      break
    }
  }
  
  tibble::tibble(
    item = seq_len(J),
    stat = stats1,
    df   = dfs,
    p    = pvals,
    flagged = pvals < alpha,
    
    iter = iter,
    skipped_in_final_fit = seq_len(J) %in% item_skip,
    n_skipped_final = length(item_skip),
    n_anchor_final = J - length(item_skip),
    converged = converged,
    
    n_na_theta = fit$n_na_theta,
    n_infinite_theta = fit$n_infinite_theta,
    n_nonfinite_theta = fit$n_nonfinite_theta,
    all_theta_finite = fit$all_theta_finite,
    
    n_na_theta_clean = fit$n_na_theta_clean,
    n_infinite_theta_clean = fit$n_infinite_theta_clean,
    n_nonfinite_theta_clean = fit$n_nonfinite_theta_clean,
    all_theta_finite_clean = fit$all_theta_finite_clean
  )
}


############################################################
## 8) RDIF
############################################################

run_official_rdif <- function(data, group, model = c("2PLM", "3PLM"),
                              alpha = 0.05,
                              purify = TRUE,
                              purify.by = "rdifrs",
                              max.iter = 20,
                              D = 1,
                              score_range = c(-4, 4),
                              min.resp = NULL,
                              verbose = FALSE,
                              logfile = "simulation_log.txt") {
  model <- match.arg(model)
  
  fit <- fit_pooled_irtq(
    data = data,
    model = model,
    D = D,
    score_range = score_range,
    logfile = logfile
  )
  
  irtQ::rdif(
    x = fit$est_par,
    data = data,
    score = fit$score,
    group = group,
    focal.name = 1,
    D = D,
    alpha = alpha,
    purify = purify,
    purify.by = purify.by,
    max.iter = max.iter,
    min.resp = min.resp,
    method = "ML",
    range = score_range,
    verbose = verbose
  )
}


safe_col <- function(ds, pattern, default = NA_real_) {
  idx <- grep(pattern, names(ds))[1]
  if (is.na(idx)) return(rep(default, nrow(ds)))
  ds[[idx]]
}


rdif_to_tibble <- function(rdif_obj, purify = TRUE, alpha = 0.05) {
  block <- if (purify && !is.null(rdif_obj$with_purify)) {
    rdif_obj$with_purify
  } else {
    rdif_obj$no_purify
  }
  
  ds <- block$dif_stat
  names(ds) <- tolower(names(ds))
  
  item_col <- names(ds)[grepl("item", names(ds))][1]
  if (is.na(item_col)) item_col <- names(ds)[1]
  
  p_rdifr_raw  <- safe_col(ds, "p.*rdifr", NA_real_)
  p_rdifs_raw  <- safe_col(ds, "p.*rdifs", NA_real_)
  p_rdifrs_raw <- safe_col(ds, "p.*rdifrs", NA_real_)
  
  p_rdifr_clean  <- clean_p_value(p_rdifr_raw)
  p_rdifs_clean  <- clean_p_value(p_rdifs_raw)
  p_rdifrs_clean <- clean_p_value(p_rdifrs_raw)
  
  tibble::tibble(
    item           = seq_len(nrow(ds)),
    item_id        = ds[[item_col]],
    
    rdifr          = safe_col(ds, "^rdifr$", NA_real_),
    z_rdifr        = safe_col(ds, "std.*rdifr|z.*rdifr", NA_real_),
    rdifs          = safe_col(ds, "^rdifs$", NA_real_),
    z_rdifs        = safe_col(ds, "std.*rdifs|z.*rdifs", NA_real_),
    rdifrs         = safe_col(ds, "^rdifrs$", NA_real_),
    
    p_rdifr        = p_rdifr_clean,
    p_rdifs        = p_rdifs_clean,
    p_rdifrs       = p_rdifrs_clean,
    
    p_rdifr_raw    = p_rdifr_raw,
    p_rdifs_raw    = p_rdifs_raw,
    p_rdifrs_raw   = p_rdifrs_raw,
    
    n_ref          = safe_col(ds, "n.*ref|reference", NA_real_),
    n_focal        = safe_col(ds, "n.*foc|focal", NA_real_),
    n_total        = safe_col(ds, "n.*total|total", NA_real_),
    iter           = if ("iter" %in% names(ds)) ds[["iter"]] else 1L,
    
    flagged_rdifr  = p_rdifr_clean < alpha,
    flagged_rdifs  = p_rdifs_clean < alpha,
    flagged_rdifrs = p_rdifrs_clean < alpha,
    
    rdifr_p_was_na  = is.na(p_rdifr_raw)  | !is.finite(p_rdifr_raw),
    rdifs_p_was_na  = is.na(p_rdifs_raw)  | !is.finite(p_rdifs_raw),
    rdifrs_p_was_na = is.na(p_rdifrs_raw) | !is.finite(p_rdifrs_raw)
  )
}


############################################################
## 9) Semi procedure
############################################################

run_semi_dif <- function(data, group, model = c("2PLM", "3PLM"),
                         alpha = 0.05,
                         purify = TRUE,
                         max.iter = 10,
                         D = 1,
                         score_range = c(-4, 4),
                         min.resp = NULL,
                         logfile = "simulation_log.txt") {
  model <- match.arg(model)
  
  run_custom_purification(
    data = data,
    group = group,
    model = model,
    test_fun = spline_lr_item,
    alpha = alpha,
    purify = purify,
    max.iter = max.iter,
    D = D,
    score_range = score_range,
    min.resp = min.resp,
    logfile = logfile
  )
}




############################################################
## 10) MI procedure
############################################################

run_mi_dif <- function(data, group, model = c("2PLM", "3PLM"),
                       alpha = 0.05,
                       purify = TRUE,
                       max.iter = 10,
                       D = 1,
                       score_range = c(-4, 4),
                       min.resp = NULL,
                       logfile = "simulation_log.txt") {
  model <- match.arg(model)
  
  run_custom_purification(
    data = data,
    group = group,
    model = model,
    test_fun = mi_dif_item,
    alpha = alpha,
    purify = purify,
    max.iter = max.iter,
    D = D,
    score_range = score_range,
    min.resp = min.resp,
    logfile = logfile
  )
}


make_empty_custom_tab <- function(J, prefix) {
  out <- tibble::tibble(
    item = seq_len(J),
    stat = NA_real_,
    df = NA_real_,
    p = rep(1, J),
    flagged = FALSE,
    iter = NA_integer_,
    skipped_in_final_fit = NA,
    n_skipped_final = NA_integer_,
    n_anchor_final = NA_integer_,
    converged = FALSE,
    
    n_na_theta = NA_integer_,
    n_infinite_theta = NA_integer_,
    n_nonfinite_theta = NA_integer_,
    all_theta_finite = NA,
    
    n_na_theta_clean = NA_integer_,
    n_infinite_theta_clean = NA_integer_,
    n_nonfinite_theta_clean = NA_integer_,
    all_theta_finite_clean = NA
  )
  
  out |>
    dplyr::rename(
      "{prefix}_stat" := stat,
      "{prefix}_df" := df,
      "{prefix}_p" := p,
      "{prefix}_flagged" := flagged,
      "{prefix}_iter" := iter,
      "{prefix}_skipped_in_final_fit" := skipped_in_final_fit,
      "{prefix}_n_skipped_final" := n_skipped_final,
      "{prefix}_n_anchor_final" := n_anchor_final,
      "{prefix}_converged" := converged,
      
      "{prefix}_n_na_theta" := n_na_theta,
      "{prefix}_n_infinite_theta" := n_infinite_theta,
      "{prefix}_n_nonfinite_theta" := n_nonfinite_theta,
      "{prefix}_all_theta_finite" := all_theta_finite,
      
      "{prefix}_n_na_theta_clean" := n_na_theta_clean,
      "{prefix}_n_infinite_theta_clean" := n_infinite_theta_clean,
      "{prefix}_n_nonfinite_theta_clean" := n_nonfinite_theta_clean,
      "{prefix}_all_theta_finite_clean" := all_theta_finite_clean
    )
}


############################################################
## 12) Single-condition runner
############################################################

run_one_condition <- function(cond,
                              alpha = 0.05,
                              D = 1,
                              purify = TRUE,
                              purify.by = "rdifrs",
                              max.iter = 20,
                              score_range = c(-4, 4),
                              min.resp = NULL,
                              verbose = FALSE,
                              logfile = "simulation_log.txt") {
  
  sim <- simulate_two_group_data(
    n_ref    = cond$n_ref,
    n_focal  = cond$n_focal,
    model    = cond$model,
    DIF_type = cond$DIF_type,
    DIF_pct  = cond$DIF_pct,
    impact   = cond$impact,
    D        = D
  )
  
  J <- ncol(sim$data)
  y_true <- as.integer(seq_len(J) %in% sim$dif_idx)
  
  base_fit <- fit_pooled_irtq(
    data = sim$data,
    model = cond$model,
    D = D,
    score_range = score_range,
    logfile = logfile
  )
  
  rdif_error <- NA_character_
  semi_error <- NA_character_
  mi_error   <- NA_character_
  
  ############################################################
  ## RDIF
  ############################################################
  
  rdif_tab <- tryCatch(
    {
      rdif_obj <- run_official_rdif(
        data = sim$data,
        group = sim$group,
        model = cond$model,
        alpha = alpha,
        purify = purify,
        purify.by = purify.by,
        max.iter = max.iter,
        D = D,
        score_range = score_range,
        min.resp = min.resp,
        verbose = verbose,
        logfile = logfile
      )
      
      rdif_to_tibble(rdif_obj, purify = purify, alpha = alpha)
    },
    error = function(e) {
      rdif_error <<- e$message
      append_log(
        sprintf("RDIF failed at cond_id=%s: %s", cond$cond_id, e$message),
        logfile = logfile
      )
      make_empty_rdif_tab(J, alpha = alpha)
    }
  )
  
  ############################################################
  ## Semi
  ############################################################
  
  semi_res <- tryCatch(
    {
      run_semi_dif(
        data = sim$data,
        group = sim$group,
        model = cond$model,
        alpha = alpha,
        purify = purify,
        max.iter = max.iter,
        D = D,
        score_range = score_range,
        min.resp = min.resp,
        logfile = logfile
      ) |>
        dplyr::rename(
          semi_stat = stat,
          semi_df = df,
          semi_p = p,
          semi_flagged = flagged,
          semi_iter = iter,
          semi_skipped_in_final_fit = skipped_in_final_fit,
          semi_n_skipped_final = n_skipped_final,
          semi_n_anchor_final = n_anchor_final,
          semi_converged = converged,
          
          semi_n_na_theta = n_na_theta,
          semi_n_infinite_theta = n_infinite_theta,
          semi_n_nonfinite_theta = n_nonfinite_theta,
          semi_all_theta_finite = all_theta_finite,
          
          semi_n_na_theta_clean = n_na_theta_clean,
          semi_n_infinite_theta_clean = n_infinite_theta_clean,
          semi_n_nonfinite_theta_clean = n_nonfinite_theta_clean,
          semi_all_theta_finite_clean = all_theta_finite_clean
        )
    },
    error = function(e) {
      semi_error <<- e$message
      append_log(
        sprintf("Semi failed at cond_id=%s: %s", cond$cond_id, e$message),
        logfile = logfile
      )
      make_empty_custom_tab(J, "semi")
    }
  )
  
  ############################################################
  ## MI
  ############################################################
  
  mi_res <- tryCatch(
    {
      run_mi_dif(
        data = sim$data,
        group = sim$group,
        model = cond$model,
        alpha = alpha,
        purify = purify,
        max.iter = max.iter,
        D = D,
        score_range = score_range,
        min.resp = min.resp,
        logfile = logfile
      ) |>
        dplyr::rename(
          mi_stat = stat,
          mi_df = df,
          mi_p = p,
          mi_flagged = flagged,
          mi_iter = iter,
          mi_skipped_in_final_fit = skipped_in_final_fit,
          mi_n_skipped_final = n_skipped_final,
          mi_n_anchor_final = n_anchor_final,
          mi_converged = converged,
          
          mi_n_na_theta = n_na_theta,
          mi_n_infinite_theta = n_infinite_theta,
          mi_n_nonfinite_theta = n_nonfinite_theta,
          mi_all_theta_finite = all_theta_finite,
          
          mi_n_na_theta_clean = n_na_theta_clean,
          mi_n_infinite_theta_clean = n_infinite_theta_clean,
          mi_n_nonfinite_theta_clean = n_nonfinite_theta_clean,
          mi_all_theta_finite_clean = all_theta_finite_clean
        )
    },
    error = function(e) {
      mi_error <<- e$message
      append_log(
        sprintf("MI failed at cond_id=%s: %s", cond$cond_id, e$message),
        logfile = logfile
      )
      make_empty_custom_tab(J, "mi")
    }
  )
  
  method_status <- dplyr::case_when(
    is.na(rdif_error) & is.na(semi_error) & is.na(mi_error) ~ "ok",
    TRUE ~ "partial"
  )
  
  method_error_message <- paste(
    na.omit(c(
      if (!is.na(rdif_error)) paste0("RDIF: ", rdif_error),
      if (!is.na(semi_error)) paste0("Semi: ", semi_error),
      if (!is.na(mi_error)) paste0("MI: ", mi_error)
    )),
    collapse = " | "
  )
  
  if (identical(method_error_message, "")) {
    method_error_message <- NA_character_
  }
  
  out <- rdif_tab |>
    dplyr::left_join(semi_res, by = "item") |>
    dplyr::left_join(mi_res, by = "item") |>
    dplyr::mutate(
      y_true = y_true,
      
      model = cond$model,
      scheme = cond$scheme,
      n_ref = cond$n_ref,
      n_focal = cond$n_focal,
      DIF_pct = cond$DIF_pct,
      DIF_type = cond$DIF_type,
      impact = cond$impact,
      rep = cond$rep,
      cond_id = cond$cond_id,
      seed = cond$seed,
      
      status = method_status,
      error_message = method_error_message,
      rdif_error = rdif_error,
      semi_error = semi_error,
      mi_error = mi_error,
      
      ml_all_theta_finite = base_fit$all_theta_finite,
      ml_n_na_theta = base_fit$n_na_theta,
      ml_n_infinite_theta = base_fit$n_infinite_theta,
      ml_n_nonfinite_theta = base_fit$n_nonfinite_theta,
      
      ml_all_theta_finite_clean = base_fit$all_theta_finite_clean,
      ml_n_na_theta_clean = base_fit$n_na_theta_clean,
      ml_n_infinite_theta_clean = base_fit$n_infinite_theta_clean,
      ml_n_nonfinite_theta_clean = base_fit$n_nonfinite_theta_clean
    )
  
  out
}


############################################################
## 13) Error object for completely failed conditions
############################################################

make_error_result <- function(cond, error_message) {
  tibble::tibble(
    item = NA_integer_,
    item_id = NA_character_,
    
    rdifr = NA_real_,
    z_rdifr = NA_real_,
    rdifs = NA_real_,
    z_rdifs = NA_real_,
    rdifrs = NA_real_,
    
    p_rdifr = NA_real_,
    p_rdifs = NA_real_,
    p_rdifrs = NA_real_,
    
    p_rdifr_raw = NA_real_,
    p_rdifs_raw = NA_real_,
    p_rdifrs_raw = NA_real_,
    
    n_ref = cond$n_ref,
    n_focal = cond$n_focal,
    n_total = cond$n_ref + cond$n_focal,
    iter = NA_integer_,
    
    flagged_rdifr = NA,
    flagged_rdifs = NA,
    flagged_rdifrs = NA,
    
    rdifr_p_was_na = NA,
    rdifs_p_was_na = NA,
    rdifrs_p_was_na = NA,
    
    semi_stat = NA_real_,
    semi_df = NA_real_,
    semi_p = NA_real_,
    semi_flagged = NA,
    semi_iter = NA_integer_,
    semi_skipped_in_final_fit = NA,
    semi_converged = NA,
    semi_n_na_theta = NA_integer_,
    semi_n_infinite_theta = NA_integer_,
    semi_n_nonfinite_theta = NA_integer_,
    semi_all_theta_finite = NA,
    semi_n_na_theta_clean = NA_integer_,
    semi_n_infinite_theta_clean = NA_integer_,
    semi_n_nonfinite_theta_clean = NA_integer_,
    semi_all_theta_finite_clean = NA,
    
    mi_stat = NA_real_,
    mi_df = NA_real_,
    mi_p = NA_real_,
    mi_flagged = NA,
    mi_iter = NA_integer_,
    mi_skipped_in_final_fit = NA,
    mi_converged = NA,
    mi_n_na_theta = NA_integer_,
    mi_n_infinite_theta = NA_integer_,
    mi_n_nonfinite_theta = NA_integer_,
    mi_all_theta_finite = NA,
    mi_n_na_theta_clean = NA_integer_,
    mi_n_infinite_theta_clean = NA_integer_,
    mi_n_nonfinite_theta_clean = NA_integer_,
    mi_all_theta_finite_clean = NA,
    
    y_true = NA_integer_,
    
    model = cond$model,
    scheme = cond$scheme,
    DIF_pct = cond$DIF_pct,
    DIF_type = cond$DIF_type,
    impact = cond$impact,
    rep = cond$rep,
    cond_id = cond$cond_id,
    seed = cond$seed,
    
    status = "error",
    error_message = as.character(error_message),
    rdif_error = NA_character_,
    semi_error = NA_character_,
    mi_error = NA_character_,
    
    ml_all_theta_finite = NA,
    ml_n_na_theta = NA_integer_,
    ml_n_infinite_theta = NA_integer_,
    ml_n_nonfinite_theta = NA_integer_,
    ml_all_theta_finite_clean = NA,
    ml_n_na_theta_clean = NA_integer_,
    ml_n_infinite_theta_clean = NA_integer_,
    ml_n_nonfinite_theta_clean = NA_integer_
  )
}


############################################################
## 14) Batch execution
############################################################

run_full_grid <- function(grid,
                          alpha = 0.05,
                          D = 1,
                          purify = TRUE,
                          purify.by = "rdifrs",
                          max.iter = 20,
                          score_range = c(-4, 4),
                          min.resp = NULL,
                          verbose = FALSE,
                          logfile = "simulation_log.txt") {
  
  purrr::map_dfr(seq_len(nrow(grid)), function(i) {
    cond <- as.list(grid[i, ])
    
    cat(sprintf(
      "[%d/%d] cond_id=%s | model=%s | scheme=%s | DIF=%s %s%% | impact=%.1f | rep=%s\n",
      i, nrow(grid), cond$cond_id, cond$model, cond$scheme,
      cond$DIF_type, cond$DIF_pct, cond$impact, cond$rep
    ))
    flush.console()
    
    set.seed(cond$seed)
    
    tryCatch(
      {
        run_one_condition(
          cond = cond,
          alpha = alpha,
          D = D,
          purify = purify,
          purify.by = purify.by,
          max.iter = max.iter,
          score_range = score_range,
          min.resp = min.resp,
          verbose = verbose,
          logfile = logfile
        )
      },
      error = function(e) {
        append_log(
          sprintf("run_full_grid failed at cond_id=%s: %s", cond$cond_id, e$message),
          logfile = logfile
        )
        make_error_result(cond, e$message)
      }
    )
  })
}


############################################################
## 15) Resume / checkpoint execution
############################################################

run_full_grid_resume <- function(grid,
                                 out_dir = "sim_chunks_robust",
                                 alpha = 0.05,
                                 D = 1,
                                 purify = TRUE,
                                 purify.by = "rdifrs",
                                 max.iter = 20,
                                 score_range = c(-4, 4),
                                 min.resp = NULL,
                                 verbose = FALSE,
                                 logfile = file.path(out_dir, "simulation_log.txt")) {
  
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  get_outfile <- function(cond_id) {
    file.path(out_dir, sprintf("cond_%05d.rds", cond_id))
  }
  
  finished <- vapply(grid$cond_id, function(id) file.exists(get_outfile(id)), logical(1))
  grid_todo <- grid[!finished, , drop = FALSE]
  
  cat(sprintf("Total conditions: %d\n", nrow(grid)))
  cat(sprintf("Already finished: %d\n", sum(finished)))
  cat(sprintf("Remaining: %d\n", nrow(grid_todo)))
  
  append_log(
    sprintf(
      "Resume run started. Total=%d, finished=%d, remaining=%d",
      nrow(grid), sum(finished), nrow(grid_todo)
    ),
    logfile = logfile
  )
  
  if (nrow(grid_todo) == 0) {
    message("All conditions are already completed.")
    return(invisible(NULL))
  }
  
  for (i in seq_len(nrow(grid_todo))) {
    cond <- as.list(grid_todo[i, ])
    outfile <- get_outfile(cond$cond_id)
    
    cat(sprintf(
      "[%d/%d remaining] cond_id=%s | model=%s | scheme=%s | DIF=%s %s%% | impact=%.1f | rep=%s\n",
      i, nrow(grid_todo), cond$cond_id, cond$model, cond$scheme,
      cond$DIF_type, cond$DIF_pct, cond$impact, cond$rep
    ))
    flush.console()
    
    append_log(
      sprintf(
        "Starting cond_id=%s | model=%s | scheme=%s | DIF=%s %s%% | impact=%.1f | rep=%s",
        cond$cond_id, cond$model, cond$scheme,
        cond$DIF_type, cond$DIF_pct, cond$impact, cond$rep
      ),
      logfile = logfile
    )
    
    set.seed(cond$seed)
    
    res_i <- tryCatch(
      {
        run_one_condition(
          cond = cond,
          alpha = alpha,
          D = D,
          purify = purify,
          purify.by = purify.by,
          max.iter = max.iter,
          score_range = score_range,
          min.resp = min.resp,
          verbose = verbose,
          logfile = logfile
        )
      },
      error = function(e) {
        append_log(
          sprintf("Condition completely failed cond_id=%s: %s", cond$cond_id, e$message),
          logfile = logfile
        )
        make_error_result(cond, e$message)
      }
    )
    
    saveRDS(res_i, outfile)
    append_log(sprintf("Saved %s", outfile), logfile = logfile)
    
    rm(res_i)
    invisible(gc())
  }
  
  invisible(NULL)
}


############################################################
## 16) Collect / progress utilities
############################################################

collect_results <- function(out_dir = "sim_chunks_robust") {
  files <- list.files(out_dir, pattern = "^cond_\\d+\\.rds$", full.names = TRUE)
  if (length(files) == 0) stop("No result files found in out_dir.")
  purrr::map_dfr(files, readRDS)
}


check_progress <- function(grid, out_dir = "sim_chunks_robust") {
  get_outfile <- function(cond_id) file.path(out_dir, sprintf("cond_%05d.rds", cond_id))
  finished <- vapply(grid$cond_id, function(id) file.exists(get_outfile(id)), logical(1))
  
  data.frame(
    total = nrow(grid),
    finished = sum(finished),
    remaining = sum(!finished),
    pct_finished = round(mean(finished) * 100, 2)
  )
}


collect_errors <- function(res) {
  res |>
    dplyr::filter(status == "error") |>
    dplyr::distinct(
      cond_id, model, scheme, DIF_pct, DIF_type, impact, rep,
      error_message
    )
}


collect_partial <- function(res) {
  res |>
    dplyr::filter(status == "partial") |>
    dplyr::distinct(
      cond_id, model, scheme, DIF_pct, DIF_type, impact, rep,
      rdif_error, semi_error, mi_error, error_message
    )
}


############################################################
## 17) ML finiteness summary
############################################################

summarise_ml_finiteness <- function(res) {
  res |>
    dplyr::distinct(
      cond_id, model, scheme, DIF_pct, DIF_type, impact, rep,
      status,
      ml_all_theta_finite,
      ml_n_na_theta,
      ml_n_infinite_theta,
      ml_n_nonfinite_theta,
      ml_all_theta_finite_clean,
      ml_n_na_theta_clean,
      ml_n_infinite_theta_clean,
      ml_n_nonfinite_theta_clean
    ) |>
    dplyr::summarise(
      n_conditions_total = dplyr::n(),
      n_conditions_ok = sum(status == "ok", na.rm = TRUE),
      n_conditions_partial = sum(status == "partial", na.rm = TRUE),
      n_conditions_error = sum(status == "error", na.rm = TRUE),
      
      n_conditions_raw_all_theta_finite =
        sum(ml_all_theta_finite %in% TRUE, na.rm = TRUE),
      n_conditions_raw_not_all_theta_finite =
        sum(ml_all_theta_finite %in% FALSE, na.rm = TRUE),
      
      n_conditions_clean_all_theta_finite =
        sum(ml_all_theta_finite_clean %in% TRUE, na.rm = TRUE),
      n_conditions_clean_not_all_theta_finite =
        sum(ml_all_theta_finite_clean %in% FALSE, na.rm = TRUE),
      
      total_raw_na_theta = sum(ml_n_na_theta, na.rm = TRUE),
      total_raw_infinite_theta = sum(ml_n_infinite_theta, na.rm = TRUE),
      total_raw_nonfinite_theta = sum(ml_n_nonfinite_theta, na.rm = TRUE),
      
      total_clean_na_theta = sum(ml_n_na_theta_clean, na.rm = TRUE),
      total_clean_infinite_theta = sum(ml_n_infinite_theta_clean, na.rm = TRUE),
      total_clean_nonfinite_theta = sum(ml_n_nonfinite_theta_clean, na.rm = TRUE)
    )
}


extract_nonfinite_ml_conditions <- function(res) {
  res |>
    dplyr::distinct(
      cond_id, model, scheme, DIF_pct, DIF_type, impact, rep,
      status,
      ml_all_theta_finite,
      ml_n_na_theta,
      ml_n_infinite_theta,
      ml_n_nonfinite_theta,
      ml_all_theta_finite_clean,
      ml_n_na_theta_clean,
      ml_n_infinite_theta_clean,
      ml_n_nonfinite_theta_clean
    ) |>
    dplyr::filter(
      ml_all_theta_finite %in% FALSE |
        ml_all_theta_finite_clean %in% FALSE
    )
}


############################################################
## 18) Evaluation utilities
############################################################

compute_metrics <- function(flagged, truth) {
  valid <- !is.na(flagged) & !is.na(truth)
  flagged <- flagged[valid]
  truth <- truth[valid]
  
  if (length(flagged) == 0 || length(truth) == 0) {
    return(
      tibble::tibble(
        TP = NA_integer_,
        FP = NA_integer_,
        TN = NA_integer_,
        FN = NA_integer_,
        sensitivity = NA_real_,
        specificity = NA_real_,
        precision   = NA_real_,
        npv         = NA_real_,
        type1       = NA_real_
      )
    )
  }
  
  tp <- sum(flagged == 1 & truth == 1)
  fp <- sum(flagged == 1 & truth == 0)
  tn <- sum(flagged == 0 & truth == 0)
  fn <- sum(flagged == 0 & truth == 1)
  
  tibble::tibble(
    TP = tp,
    FP = fp,
    TN = tn,
    FN = fn,
    sensitivity = if ((tp + fn) == 0) NA_real_ else tp / (tp + fn),
    specificity = if ((tn + fp) == 0) NA_real_ else tn / (tn + fp),
    precision   = if ((tp + fp) == 0) NA_real_ else tp / (tp + fp),
    npv         = if ((tn + fn) == 0) NA_real_ else tn / (tn + fn),
    type1       = if ((tn + fp) == 0) NA_real_ else fp / (tn + fp)
  )
}


summarise_by_replication <- function(res, grid = NULL) {
  
  res_cond <- res |>
    dplyr::distinct(
      cond_id, model, scheme, n_ref, n_focal,
      DIF_pct, DIF_type, impact, rep,
      status, error_message,
      rdif_error, semi_error, mi_error,
      ml_all_theta_finite, ml_n_na_theta, ml_n_infinite_theta, ml_n_nonfinite_theta,
      ml_all_theta_finite_clean, ml_n_na_theta_clean,
      ml_n_infinite_theta_clean, ml_n_nonfinite_theta_clean
    )
  
  res_usable <- res |>
    dplyr::filter(status %in% c("ok", "partial"))
  
  perf_usable <- res_usable |>
    dplyr::group_by(
      cond_id, model, scheme, n_ref, n_focal,
      DIF_pct, DIF_type, impact, rep
    ) |>
    dplyr::group_modify(~ {
      dplyr::bind_cols(
        compute_metrics(as.integer(.x$flagged_rdifr), .x$y_true) |>
          dplyr::rename_with(~ paste0("rdifr_", .x), everything()),
        
        compute_metrics(as.integer(.x$flagged_rdifs), .x$y_true) |>
          dplyr::rename_with(~ paste0("rdifs_", .x), everything()),
        
        compute_metrics(as.integer(.x$flagged_rdifrs), .x$y_true) |>
          dplyr::rename_with(~ paste0("rdifrs_", .x), everything()),
        
        compute_metrics(as.integer(.x$semi_flagged), .x$y_true) |>
          dplyr::rename_with(~ paste0("semi_", .x), everything()),
        
        compute_metrics(as.integer(.x$mi_flagged), .x$y_true) |>
          dplyr::rename_with(~ paste0("mi_", .x), everything())
      )
    }) |>
    dplyr::ungroup()
  
  if (is.null(grid)) {
    out <- res_cond |>
      dplyr::left_join(
        perf_usable,
        by = c(
          "cond_id", "model", "scheme", "n_ref", "n_focal",
          "DIF_pct", "DIF_type", "impact", "rep"
        )
      )
  } else {
    out <- grid |>
      dplyr::select(
        cond_id, model, scheme, n_ref, n_focal,
        DIF_pct, DIF_type, impact, rep
      ) |>
      dplyr::left_join(
        res_cond |>
          dplyr::select(
            cond_id, status, error_message,
            rdif_error, semi_error, mi_error,
            ml_all_theta_finite, ml_n_na_theta,
            ml_n_infinite_theta, ml_n_nonfinite_theta,
            ml_all_theta_finite_clean, ml_n_na_theta_clean,
            ml_n_infinite_theta_clean, ml_n_nonfinite_theta_clean
          ),
        by = "cond_id"
      ) |>
      dplyr::left_join(
        perf_usable,
        by = c(
          "cond_id", "model", "scheme", "n_ref", "n_focal",
          "DIF_pct", "DIF_type", "impact", "rep"
        )
      )
  }
  
  out |>
    dplyr::arrange(cond_id)
}


############################################################
## 19) Optional final summary across replications
############################################################

summarise_perf_mean <- function(perf_all) {
  perf_all |>
    dplyr::group_by(model, scheme, n_ref, n_focal, DIF_pct, DIF_type, impact) |>
    dplyr::summarise(
      n_conditions = dplyr::n(),
      n_ok = sum(status == "ok", na.rm = TRUE),
      n_partial = sum(status == "partial", na.rm = TRUE),
      n_error = sum(status == "error", na.rm = TRUE),
      
      dplyr::across(
        where(is.numeric) & 
          !dplyr::all_of(c(
            "cond_id", "n_ref", "n_focal", "DIF_pct", "impact", "rep",
            "ml_n_na_theta", "ml_n_infinite_theta", "ml_n_nonfinite_theta",
            "ml_n_na_theta_clean", "ml_n_infinite_theta_clean", "ml_n_nonfinite_theta_clean"
          )),
        ~ mean(.x, na.rm = TRUE),
        .names = "mean_{.col}"
      ),
      .groups = "drop"
    )
}


############################################################
## 20) Example execution
############################################################


grid_full <- make_design(reps_per = 100)

## Use a new folder to avoid mixing old failed chunks with robust chunks
out_dir <- "sim_chunks_purified_v5"

run_full_grid_resume(
  grid = grid_full,
  out_dir = out_dir,
  alpha = 0.05,
  D = 1,
  purify = TRUE,
  purify.by = "rdifrs",
  max.iter = 20,
  score_range = c(-4, 4),
  min.resp = NULL,
  verbose = FALSE
)

check_progress(grid_full, out_dir)

res_all <- collect_results(out_dir)

error_conditions <- collect_errors(res_all)
partial_conditions <- collect_partial(res_all)

ml_summary <- summarise_ml_finiteness(res_all)
ml_problem_conditions <- extract_nonfinite_ml_conditions(res_all)

perf_all <- summarise_by_replication(res_all, grid = grid_full)

## Check whether all 160 conditions are retained
nrow(perf_all)
table(perf_all$status, useNA = "ifany")
 
