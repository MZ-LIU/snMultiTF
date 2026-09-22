#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# Helper script for TF driver identification.
# Uses mgcv::gam for ActDyn and ProgramCouple.

suppressPackageStartupMessages({
  library(mgcv)
  library(dplyr)
  library(readr)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 6) {
  stop(
    "Usage: Rscript run_mgcv_tf_driver.R <mode> <input_csv> <output_csv> ",
    "<k_value> <epsilon> <grid_size>\n",
    "mode must be 'actdyn' or 'program'"
  )
}

mode <- args[1]
input_csv <- args[2]
output_csv <- args[3]
k_value <- as.integer(args[4])
epsilon <- as.numeric(args[5])
grid_size <- as.integer(args[6])

rank_scaled_high_better <- function(x) {
  out <- rep(0, length(x))
  ok <- is.finite(x)
  n <- sum(ok)
  if (n == 0) {
    return(out)
  }
  if (n == 1) {
    out[ok] <- 1
    return(out)
  }
  r <- rank(-x[ok], ties.method = "average")
  out[ok] <- 1 - ((r - 1) / (n - 1))
  out
}

fit_one_gam <- function(df, y_col, k_value, grid_size = 100) {
  df <- df %>%
    filter(is.finite(.data[[y_col]]), is.finite(pseudotime))

  n <- nrow(df)

  if (n < max(20, k_value + 5)) {
    return(list(
      pvalue = 1,
      fitted = rep(NA_real_, n),
      resid = rep(NA_real_, n),
      pred_grid = rep(NA_real_, grid_size),
      amp = NA_real_,
      n_used = n,
      edf = NA_real_,
      deviance_explained = NA_real_
    ))
  }

  if (sd(df[[y_col]], na.rm = TRUE) == 0) {
    return(list(
      pvalue = 1,
      fitted = rep(NA_real_, n),
      resid = rep(NA_real_, n),
      pred_grid = rep(NA_real_, grid_size),
      amp = NA_real_,
      n_used = n,
      edf = NA_real_,
      deviance_explained = NA_real_
    ))
  }

  form <- as.formula(paste0(y_col, " ~ s(pseudotime, k = ", k_value, ")"))

  fit <- tryCatch(
    gam(form, data = df, family = gaussian(), method = "REML"),
    error = function(e) NULL
  )

  if (is.null(fit)) {
    return(list(
      pvalue = 1,
      fitted = rep(NA_real_, n),
      resid = rep(NA_real_, n),
      pred_grid = rep(NA_real_, grid_size),
      amp = NA_real_,
      n_used = n,
      edf = NA_real_,
      deviance_explained = NA_real_
    ))
  }

  sm <- summary(fit)

  pval <- tryCatch(
    as.numeric(sm$s.table[1, "p-value"]),
    error = function(e) 1
  )

  edf <- tryCatch(
    as.numeric(sm$s.table[1, "edf"]),
    error = function(e) NA_real_
  )

  deviance_explained <- tryCatch(
    as.numeric(sm$dev.expl),
    error = function(e) NA_real_
  )
  if (is.finite(deviance_explained)) {
    deviance_explained <- max(0, min(1, deviance_explained))
  }

  fitted_vals <- as.numeric(fitted(fit))
  resid_vals <- as.numeric(residuals(fit, type = "response"))

  grid <- data.frame(
    pseudotime = seq(
      min(df$pseudotime, na.rm = TRUE),
      max(df$pseudotime, na.rm = TRUE),
      length.out = grid_size
    )
  )

  pred_grid <- tryCatch(
    as.numeric(predict(fit, newdata = grid, type = "response")),
    error = function(e) rep(NA_real_, grid_size)
  )

  amp <- sd(pred_grid, na.rm = TRUE)

  list(
    pvalue = pval,
    fitted = fitted_vals,
    resid = resid_vals,
    pred_grid = pred_grid,
    amp = amp,
    n_used = n,
    edf = edf,
    deviance_explained = deviance_explained
  )
}

if (mode == "actdyn") {
  dat <- read_csv(input_csv, show_col_types = FALSE)

  res <- dat %>%
    group_by(TF) %>%
    group_modify(function(.x, .y) {
      fit <- fit_one_gam(.x, "activity", k_value, grid_size)

      tibble(
        gam_pvalue = fit$pvalue,
        Amp = fit$amp,
        activity_std = sd(.x$activity, na.rm = TRUE),
        gam_n_used = fit$n_used,
        gam_edf = fit$edf,
        deviance_explained = fit$deviance_explained
      )
    }) %>%
    ungroup()

  res <- res %>%
    mutate(
      FDR_GAM = p.adjust(gam_pvalue, method = "BH"),
      Amp_rank_scaled = rank_scaled_high_better(Amp),
      DevExpl_rank_scaled = rank_scaled_high_better(deviance_explained),
      ActDyn = Amp_rank_scaled * DevExpl_rank_scaled
    ) %>%
    arrange(desc(ActDyn), desc(deviance_explained), desc(Amp)) %>%
    mutate(ActDyn_rank = row_number()) %>%
    select(ActDyn_rank, everything())

  write_csv(res, output_csv)
} else if (mode == "program") {
  dat <- read_csv(input_csv, show_col_types = FALSE)

  res <- dat %>%
    group_by(TF) %>%
    group_modify(function(.x, .y) {
      n_targets <- suppressWarnings(as.integer(unique(.x$n_targets)[1]))

      if (!is.finite(n_targets) || n_targets < 10) {
        return(tibble(
          ProgramCouple = NA_real_,
          program_couple_pvalue = NA_real_,
          target_program_n_targets = n_targets,
          program_couple_spline_k = k_value,
          program_couple_grid_size = grid_size,
          program_couple_activity_edf = NA_real_,
          program_couple_program_edf = NA_real_,
          status = "too_few_targets"
        ))
      }

      .x <- .x %>%
        filter(is.finite(activity), is.finite(program_score), is.finite(pseudotime))

      if (nrow(.x) < max(20, k_value + 5)) {
        return(tibble(
          ProgramCouple = NA_real_,
          program_couple_pvalue = NA_real_,
          target_program_n_targets = n_targets,
          program_couple_spline_k = k_value,
          program_couple_grid_size = grid_size,
          program_couple_activity_edf = NA_real_,
          program_couple_program_edf = NA_real_,
          status = "too_few_cells"
        ))
      }

      fit_a <- fit_one_gam(.x, "activity", k_value, grid_size)
      fit_s <- fit_one_gam(.x, "program_score", k_value, grid_size)

      ok <- is.finite(fit_a$pred_grid) & is.finite(fit_s$pred_grid)

      if (sum(ok) < 10 ||
          sd(fit_a$pred_grid[ok], na.rm = TRUE) == 0 ||
          sd(fit_s$pred_grid[ok], na.rm = TRUE) == 0) {
        corr <- NA_real_
        pval <- NA_real_
        status <- "failed"
      } else {
        ct <- cor.test(fit_a$pred_grid[ok], fit_s$pred_grid[ok], method = "pearson")
        corr <- abs(as.numeric(ct$estimate))
        pval <- as.numeric(ct$p.value)
        status <- "ok"
      }

      tibble(
        ProgramCouple = corr,
        program_couple_pvalue = pval,
        target_program_n_targets = n_targets,
        program_couple_spline_k = k_value,
        program_couple_grid_size = grid_size,
        program_couple_activity_edf = fit_a$edf,
        program_couple_program_edf = fit_s$edf,
        status = status
      )
    }) %>%
    ungroup() %>%
    arrange(desc(ProgramCouple), desc(target_program_n_targets)) %>%
    mutate(ProgramCouple_rank = row_number()) %>%
    select(ProgramCouple_rank, everything())

  write_csv(res, output_csv)
} else {
  stop("Unknown mode: ", mode, ". Expected 'actdyn' or 'program'.")
}
