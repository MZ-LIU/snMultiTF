#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# Helper script for TF driver identification.
# Uses RobustRankAggreg::aggregateRanks for rank integration.

suppressPackageStartupMessages({
  library(RobustRankAggreg)
  library(readr)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 5) {
  stop(
    "Usage: Rscript run_rra_tf_driver.R <network_csv> <actdyn_csv> ",
    "<program_csv> <universe_csv> <output_csv>"
  )
}

network_csv <- args[1]
actdyn_csv <- args[2]
program_csv <- args[3]
universe_csv <- args[4]
output_csv <- args[5]

network <- read_csv(network_csv, show_col_types = FALSE)
actdyn <- read_csv(actdyn_csv, show_col_types = FALSE)
program <- read_csv(program_csv, show_col_types = FALSE)
universe <- read_csv(universe_csv, show_col_types = FALSE)

universe_tfs <- unique(as.character(universe$TF))
N <- length(universe_tfs)

if (N == 0) {
  stop("Universe is empty.")
}

network_list <- network %>%
  filter(TF %in% universe_tfs) %>%
  arrange(rank) %>%
  pull(TF) %>%
  as.character() %>%
  unique()

actdyn_list <- actdyn %>%
  filter(TF %in% universe_tfs) %>%
  arrange(rank) %>%
  pull(TF) %>%
  as.character() %>%
  unique()

program_list <- program %>%
  filter(TF %in% universe_tfs, is.finite(ProgramCouple)) %>%
  arrange(rank) %>%
  pull(TF) %>%
  as.character() %>%
  unique()

glist <- list(
  NetworkSupport = network_list,
  ActDyn = actdyn_list,
  ProgramCouple = program_list
)

rra <- aggregateRanks(
  glist = glist,
  N = N,
  method = "RRA",
  full = FALSE
)

rra <- as.data.frame(rra)

# RobustRankAggreg usually returns columns Name and Score.
if ("Name" %in% colnames(rra)) {
  rra <- rra %>% rename(TF = Name)
} else {
  colnames(rra)[1] <- "TF"
}

if ("Score" %in% colnames(rra)) {
  rra <- rra %>% rename(RRA_pvalue = Score)
} else {
  colnames(rra)[2] <- "RRA_pvalue"
}

rra <- rra %>%
  mutate(
    TF = as.character(TF),
    RRA_pvalue = as.numeric(RRA_pvalue),
    RRA_FDR = p.adjust(RRA_pvalue, method = "BH")
  ) %>%
  arrange(RRA_pvalue) %>%
  mutate(rank = row_number()) %>%
  select(rank, TF, RRA_pvalue, RRA_FDR)

write_csv(rra, output_csv)
