args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Usage: Rscript run_archr_trajectory.R <config.json>")
}

if (!requireNamespace("jsonlite", quietly = TRUE)) {
  stop("The R package 'jsonlite' is required.")
}
if (!requireNamespace("nabor", quietly = TRUE)) {
  stop("The R package 'nabor' is required. Install with: install.packages('nabor')")
}

library(magrittr)
library(nabor)

config <- jsonlite::fromJSON(args[[1]], simplifyVector = TRUE)
set.seed(as.integer(config$seed))

# Helper function (same as scMEGA)
getQuantiles <- function(x) {
  return(round(rank(x) / length(x), digits = 5))
}

# Read input CSVs
metadata <- read.csv(config$metadata_csv, stringsAsFactors = FALSE, check.names = FALSE)
embedding <- read.csv(config$embedding_csv, stringsAsFactors = FALSE, check.names = FALSE)

if (!"mdata_cell_id" %in% colnames(metadata)) {
  stop("metadata_csv must contain 'mdata_cell_id'.")
}
if (!config$group_by %in% colnames(metadata)) {
  stop(sprintf("metadata_csv must contain group column '%s'.", config$group_by))
}
if (!"mdata_cell_id" %in% colnames(embedding)) {
  stop("embedding_csv must contain 'mdata_cell_id'.")
}

# Merge metadata and embedding
input_df <- merge(metadata, embedding, by = "mdata_cell_id", all = FALSE, sort = FALSE)
if (nrow(input_df) == 0) {
  stop("No overlapping cells between metadata_csv and embedding_csv.")
}

# Get trajectory groups
trajectory <- as.character(config$trajectory)
group_by <- config$group_by

# Filter to trajectory cells only (same as scMEGA)
df.group <- data.frame(group = input_df[[group_by]], row.names = input_df$mdata_cell_id)
df.group <- df.group[df.group$group %in% trajectory, , drop = FALSE]

if (nrow(df.group) == 0) {
  stop("No cells found in the specified trajectory groups.")
}

# Get embedding dimensions
dim_cols <- setdiff(colnames(embedding), "mdata_cell_id")
if (length(dim_cols) < 2) {
  stop("embedding_csv must contain at least 2 dimensions.")
}

# Extract coordinates for all cells
data.use.all <- as.matrix(input_df[, dim_cols, drop = FALSE])
rownames(data.use.all) <- input_df$mdata_cell_id

# Filter to trajectory cells
data.use <- data.use.all[rownames(df.group), , drop = FALSE]

# Pre-filter outliers (same as scMEGA)
filterObj <- lapply(seq_along(trajectory), function(x) {
  groupsx <- rownames(df.group)[df.group$group == trajectory[x]]
  matx <- data.use[groupsx, , drop = FALSE]
  matMeanx <- colMeans(matx)
  diffx <- sqrt(colSums((t(matx) - matMeanx)^2))
  idxKeep <- which(diffx <= quantile(diffx, config$pre_filter_quantile))
  list(mat = matx[idxKeep, , drop = FALSE], groups = groupsx[idxKeep])
})

matList <- lapply(seq_along(filterObj), function(x) filterObj[[x]]$mat)
matFilter <- Reduce("rbind", matList)

groupList <- lapply(seq_along(filterObj), function(x) filterObj[[x]]$groups)
allGroups <- Reduce("c", groupList)
groupsFilter <- df.group[allGroups, , drop = FALSE]

# Compute initial time based on centroid distances (same as scMEGA)
initialTime <- lapply(seq_along(trajectory), function(x) {
  groupsx <- rownames(groupsFilter)[groupsFilter$group == trajectory[x]]
  matx <- matFilter[groupsx, , drop = FALSE]
  if (x != length(trajectory)) {
    groupsxp1 <- rownames(groupsFilter)[groupsFilter$group == trajectory[x + 1]]
    meanx <- colMeans(matFilter[groupsxp1, , drop = FALSE])
    diffx <- sqrt(colSums((t(matx) - meanx)^2))
    timex <- (1 - getQuantiles(diffx)) + x
  } else {
    groupsxm1 <- rownames(groupsFilter)[groupsFilter$group == trajectory[x - 1]]
    meanx <- colMeans(matFilter[groupsxm1, , drop = FALSE])
    diffx <- sqrt(colSums((t(matx) - meanx)^2))
    timex <- getQuantiles(diffx) + x
  }
  timex
}) %>% unlist

# Fit smooth splines (same as scMEGA)
dof <- as.integer(config$dof)
spar <- as.numeric(config$spar)

matSpline <- lapply(seq_len(ncol(matFilter)), function(x) {
  stats::smooth.spline(x = initialTime, y = matFilter[names(initialTime), x], df = dof, spar = spar)[[2]]
}) %>% Reduce("cbind", .) %>% data.frame()

# KNN projection (same as scMEGA)
knnObj <- nabor::knn(data = matSpline, query = data.use, k = 3)
knnIdx <- knnObj[[1]]
knnDist <- knnObj[[2]]
knnDiff <- ifelse(knnIdx[, 2] > knnIdx[, 3], 1, -1)
knnDistQ <- getQuantiles(knnDist[, 1])

# Post-filter (same as scMEGA)
idxKeep <- which(knnDist[, 1] <= quantile(knnDist[, 1], config$post_filter_quantile))
dfTrajectory <- data.frame(
  row.names = rownames(data.use),
  Distance = knnDist[, 1],
  Trajectory = knnIdx[, 1] + knnDiff * knnDistQ
)[idxKeep, , drop = FALSE]

# Scale pseudotime to [0, 100] (same as scMEGA)
pseudotime_raw <- dfTrajectory$Trajectory
names(pseudotime_raw) <- rownames(dfTrajectory)
pseudotime <- 100 * getQuantiles(pseudotime_raw)

# Create output dataframe with all cells
output_df <- data.frame(
  mdata_cell_id = input_df$mdata_cell_id,
  stringsAsFactors = FALSE
)
output_df[[config$trajectory_name]] <- NA_real_

# Map pseudotime back to cells
matched_cells <- intersect(names(pseudotime), input_df$mdata_cell_id)
output_df[[config$trajectory_name]][match(matched_cells, input_df$mdata_cell_id)] <- pseudotime[matched_cells]

# Write output
write.csv(output_df, file = config$output_csv, row.names = FALSE, quote = TRUE)
cat(sprintf("Trajectory saved to: %s\n", config$output_csv))
cat(sprintf("Number of trajectory cells: %d\n", length(matched_cells)))
cat(sprintf("Pseudotime range: [%.3f, %.3f]\n", min(pseudotime, na.rm = TRUE), max(pseudotime, na.rm = TRUE)))
cat("Done!\n")
