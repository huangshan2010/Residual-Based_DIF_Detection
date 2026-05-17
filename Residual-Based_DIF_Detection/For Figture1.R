############################################################
## Final conceptual figure for residual-based DIF procedures
## A. Raw residual pattern
## B. RDIF
## C. Spline-LR
## D. MI
############################################################

library(dplyr)
library(ggplot2)
library(patchwork)
library(mgcv)

make_residual_method_figure_final <- function(
    out_prefix = "Figure_residual_methods_illustration_final",
    theta_bins_a = 12,
    theta_bins_mi = 5,
    residual_bins_mi = 4,
    y_lim_resid = c(-0.12, 0.12),
    interval_digits = 2
) {
  
  ############################################################
  ## 0. Helper: pretty interval labels
  ############################################################
  make_interval_labels <- function(breaks, digits = 2) {
    fmt <- function(x) format(round(x, digits), nsmall = digits, trim = TRUE)
    labs <- character(length(breaks) - 1)
    for (i in seq_len(length(breaks) - 1)) {
      left_bracket  <- if (i == 1) "[" else "("
      right_bracket <- "]"
      labs[i] <- paste0(left_bracket, fmt(breaks[i]), ", ", fmt(breaks[i + 1]), right_bracket)
    }
    labs
  }
  
  ############################################################
  ## 1. Check required functions from the Study 2 script
  ############################################################
  needed_objs <- c(
    "make_design",
    "simulate_two_group_data",
    "fit_pooled_irtq",
    "compute_raw_residuals"
  )
  
  missing_objs <- needed_objs[!vapply(needed_objs, exists, logical(1))]
  if (length(missing_objs) > 0) {
    stop(
      "Please source the Study 2 simulation script first. Missing: ",
      paste(missing_objs, collapse = ", ")
    )
  }
  
  ############################################################
  ## 2. Pick one representative Study 2 condition
  ############################################################
  grid_tmp <- if (exists("grid_full")) {
    grid_full
  } else {
    make_design(reps_per = 100)
  }
  
  cond <- grid_tmp %>%
    filter(
      scheme == "1000R_1000F",
      item_profile == "centered",
      impact == "no_impact",
      DIF_magnitude == "large",
      rep == 1
    ) %>%
    slice(1)
  
  if (nrow(cond) == 0) {
    message("Target condition not found. Using first available condition.")
    cond <- make_design(reps_per = 1) %>% slice(1)
  }
  
  cond <- as.list(cond)
  
  ############################################################
  ## 3. Regenerate the representative dataset
  ############################################################
  set.seed(cond$seed)
  
  sim <- simulate_two_group_data(
    n_ref      = cond$n_ref,
    n_focal    = cond$n_focal,
    model      = cond$model,
    n_items    = cond$n_items,
    b_mean     = cond$b_mean,
    DIF_pct    = cond$DIF_pct,
    DIF_type   = cond$DIF_type,
    impact     = cond$impact,
    DIF_gamma  = cond$DIF_gamma,
    D          = 1,
    item_seed  = cond$item_seed,
    sim_seed   = cond$sim_seed
  )
  
  fit <- fit_pooled_irtq(
    data = sim$data,
    model = cond$model,
    D = 1,
    score_range = c(-4, 4)
  )
  
  raw_out <- compute_raw_residuals(
    data = sim$data,
    est_par = fit$est_par,
    score = fit$score,
    D = 1
  )
  
  resid_raw <- raw_out$resid_raw
  
  ## Use the first true DIF item for illustration
  item_j <- sim$dif_idx[1]
  
  dat <- data.frame(
    theta = fit$score,
    residual = resid_raw[, item_j],
    group01 = sim$group,
    group = factor(
      sim$group,
      levels = c(0, 1),
      labels = c("Reference", "Focal")
    )
  ) %>%
    filter(is.finite(theta), is.finite(residual)) %>%
    mutate(
      sq_residual = residual^2
    )
  
  theta_min <- min(dat$theta, na.rm = TRUE)
  theta_max <- max(dat$theta, na.rm = TRUE)
  theta_grid <- seq(theta_min, theta_max, length.out = 200)
  
  ############################################################
  ## A. Raw residual pattern
  ############################################################
  theta_breaks_a <- quantile(
    dat$theta,
    probs = seq(0, 1, length.out = theta_bins_a + 1),
    na.rm = TRUE
  )
  theta_breaks_a <- unique(as.numeric(theta_breaks_a))
  
  if (length(theta_breaks_a) < 4) {
    theta_breaks_a <- seq(theta_min, theta_max, length.out = theta_bins_a + 1)
  }
  
  theta_centers_a <- (head(theta_breaks_a, -1) + tail(theta_breaks_a, -1)) / 2
  
  theta_center_df_a <- data.frame(
    theta_bin = seq_along(theta_centers_a),
    theta_center = theta_centers_a
  )
  
  raw_pat <- dat %>%
    mutate(
      theta_bin = cut(
        theta,
        breaks = theta_breaks_a,
        include.lowest = TRUE,
        labels = FALSE
      )
    ) %>%
    filter(!is.na(theta_bin)) %>%
    group_by(group, theta_bin) %>%
    summarise(
      mean_residual = mean(residual),
      .groups = "drop"
    ) %>%
    left_join(theta_center_df_a, by = "theta_bin")
  
  p_a <- ggplot(
    raw_pat,
    aes(x = theta_center, y = mean_residual, linetype = group)
  ) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2) +
    scale_linetype_manual(values = c("solid", "dashed")) +
    coord_cartesian(ylim = y_lim_resid) +
    labs(
      title = "A. Raw residual pattern",
      x = expression(hat(theta)),
      y = "Residual",
      linetype = "Group"
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold")
    )
  
  ############################################################
  ## B. RDIF
  ## Schematic near-linear residual-moment comparison
  ############################################################
  rdif_tab <- dat %>%
    group_by(group) %>%
    summarise(
      mean_residual = mean(residual),
      mean_sq_residual = mean(sq_residual),
      .groups = "drop"
    )
  
  mean_r_ref <- rdif_tab$mean_residual[rdif_tab$group == "Reference"]
  mean_r_foc <- rdif_tab$mean_residual[rdif_tab$group == "Focal"]
  
  rdif_r <- mean_r_foc - mean_r_ref
  
  mean_s_ref <- rdif_tab$mean_sq_residual[rdif_tab$group == "Reference"]
  mean_s_foc <- rdif_tab$mean_sq_residual[rdif_tab$group == "Focal"]
  
  rdif_s <- mean_s_foc - mean_s_ref
  
  theta_mid <- mean(c(theta_min, theta_max))
  
  ## Schematic lines only
  small_slope_ref <- 0.006
  small_slope_foc <- 0.004
  
  rdif_line_dat <- bind_rows(
    data.frame(
      theta = theta_grid,
      fitted = mean_r_ref + small_slope_ref * (theta_grid - theta_mid),
      group = "Reference"
    ),
    data.frame(
      theta = theta_grid,
      fitted = mean_r_foc + small_slope_foc * (theta_grid - theta_mid),
      group = "Focal"
    )
  ) %>%
    mutate(group = factor(group, levels = c("Reference", "Focal")))
  
  p_b <- ggplot(
    rdif_line_dat,
    aes(x = theta, y = fitted, linetype = group)
  ) +
    geom_line(linewidth = 1) +
    scale_linetype_manual(values = c("solid", "dashed")) +
    coord_cartesian(ylim = y_lim_resid) +
    labs(
      title = "B. RDIF",
      x = expression(hat(theta)),
      y = "Residual",
      linetype = "Group"
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold")
    )
  
  ############################################################
  ## C. Spline-LR
  ## Group-specific nonlinear fitted residual curves
  ############################################################
  m_spline_vis <- mgcv::gam(
    residual ~ group + s(theta, by = group, k = 5),
    data = dat,
    method = "REML"
  )
  
  ## Dense grid for fitted lines
  theta_dense <- seq(theta_min, theta_max, length.out = 200)
  
  pred_c_line <- expand.grid(
    theta = theta_dense,
    group = levels(dat$group)
  )
  pred_c_line$group <- factor(pred_c_line$group, levels = levels(dat$group))
  pred_c_line$fitted_residual <- predict(m_spline_vis, newdata = pred_c_line)
  
  ## Sparse key points only
  theta_key <- quantile(
    dat$theta,
    probs = c(0.10, 0.35, 0.65, 0.90),
    na.rm = TRUE
  )
  theta_key <- unique(as.numeric(theta_key))
  
  pred_c_pts <- expand.grid(
    theta = theta_key,
    group = levels(dat$group)
  )
  pred_c_pts$group <- factor(pred_c_pts$group, levels = levels(dat$group))
  pred_c_pts$fitted_residual <- predict(m_spline_vis, newdata = pred_c_pts)
  
  p_c <- ggplot(
    pred_c_line,
    aes(x = theta, y = fitted_residual, linetype = group)
  ) +
    geom_line(linewidth = 1) +
    geom_point(
      data = pred_c_pts,
      size = 2
    ) +
    scale_linetype_manual(values = c("solid", "dashed")) +
    coord_cartesian(ylim = y_lim_resid) +
    labs(
      title = "C. Spline-LR",
      x = expression(hat(theta)),
      y = "Residual",
      linetype = "Group"
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold")
    )
  
  ############################################################
  ## D. MI
  ## Ability strata by residual category
  ## Fill = deviation from 0.50
  ############################################################
  theta_breaks_mi <- quantile(
    dat$theta,
    probs = seq(0, 1, length.out = theta_bins_mi + 1),
    na.rm = TRUE
  )
  theta_breaks_mi <- unique(as.numeric(theta_breaks_mi))
  
  if (length(theta_breaks_mi) < 4) {
    theta_breaks_mi <- seq(theta_min, theta_max, length.out = theta_bins_mi + 1)
  }
  
  resid_breaks_mi <- quantile(
    dat$residual,
    probs = seq(0, 1, length.out = residual_bins_mi + 1),
    na.rm = TRUE
  )
  resid_breaks_mi <- unique(as.numeric(resid_breaks_mi))
  
  if (length(resid_breaks_mi) < 4) {
    resid_breaks_mi <- seq(
      min(dat$residual, na.rm = TRUE),
      max(dat$residual, na.rm = TRUE),
      length.out = residual_bins_mi + 1
    )
  }
  
  theta_labels_mi <- make_interval_labels(theta_breaks_mi, digits = interval_digits)
  resid_labels_mi <- make_interval_labels(resid_breaks_mi, digits = interval_digits)
  
  dat_mi <- dat %>%
    mutate(
      theta_stratum = cut(
        theta,
        breaks = theta_breaks_mi,
        include.lowest = TRUE,
        labels = theta_labels_mi
      ),
      residual_cat = cut(
        residual,
        breaks = resid_breaks_mi,
        include.lowest = TRUE,
        labels = resid_labels_mi
      )
    ) %>%
    filter(!is.na(theta_stratum), !is.na(residual_cat))
  
  mi_tab <- dat_mi %>%
    group_by(theta_stratum, residual_cat) %>%
    summarise(
      focal_prop = mean(group01),
      n = n(),
      .groups = "drop"
    ) %>%
    mutate(
      dev50 = abs(focal_prop - 0.50),
      text_col = ifelse(dev50 >= 0.12, "white", "black")
    )
  
  ## Keep interval order
  mi_tab$theta_stratum <- factor(mi_tab$theta_stratum, levels = theta_labels_mi)
  mi_tab$residual_cat  <- factor(mi_tab$residual_cat, levels = resid_labels_mi)
  
  p_d <- ggplot(
    mi_tab,
    aes(x = theta_stratum, y = residual_cat, fill = dev50)
  ) +
    geom_tile(color = "black", linewidth = 0.25) +
    geom_text(
      aes(
        label = paste0("p=", round(focal_prop, 2), "\n", "n=", n),
        color = text_col
      ),
      size = 2.8
    ) +
    scale_fill_gradient(
      low = "white",
      high = "grey30",
      name = "|p - 0.50|"
    ) +
    scale_color_identity() +
    labs(
      title = "D. MI",
      x = "Ability stratum",
      y = "Residual category"
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold"),
      axis.text.x = element_text(angle = 25, hjust = 1),
      axis.text.y = element_text(size = 9)
    )
  
  ############################################################
  ## Combine panels
  ############################################################
  fig <- (p_a | p_b) / (p_c | p_d) +
    plot_annotation(
      title = paste0(
        "Illustration of Residual-Based DIF Procedures"
      )
    )
  
  ############################################################
  ## Save output
  ############################################################
  ggsave(
    filename = paste0(out_prefix, ".png"),
    plot = fig,
    width = 13.5,
    height = 8.7,
    dpi = 300
  )
  
  ggsave(
    filename = paste0(out_prefix, ".pdf"),
    plot = fig,
    width = 13.5,
    height = 8.7
  )
  
  return(list(
    figure = fig,
    condition = cond,
    illustrated_item = item_j,
    rdif_r = rdif_r,
    rdif_s = rdif_s,
    theta_breaks_mi = theta_breaks_mi,
    resid_breaks_mi = resid_breaks_mi
  ))
}

## Run
fig_out_final <- make_residual_method_figure_final()
fig_out_final$figure