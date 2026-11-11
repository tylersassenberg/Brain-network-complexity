#data_sorting.R
#compiles per-subject results from rqa_mrn_metrics.R into sample-level datasets
#and computes summary statistics for MVMD mode central frequencies

#user settings
output_dir <- "output"
casnet_dir <- file.path(output_dir, "casnet_results")
modes      <- 2:5

net_names <- c("VisCent", "VisPeri", "SomMotA", "SomMotB", "DorsAttnA", "DorsAttnB",
               "SalVentAttnA", "SalVentAttnB", "LimbicB", "LimbicA", "ContA", "ContB",
               "ContC", "DefaultA", "DefaultB", "DefaultC", "TempPar")

#rqa and recurrence network measures kept from the layer results
rqa_vars <- c("DET", "MEAN_dl", "ENT_dl", "LAM_vl", "TT_vl", "ENT_vl",
              "GlobalClustering", "NetworkTransitivity", "GlobalEfficiency")

#list of subject ids
subject_ids <- read.csv(file.path(output_dir, "subject_ids.csv"), colClasses = "character")$subject_id

#upper triangle of an MRN interlayer MI matrix, read row-wise
read_upper_tri <- function(file_name) {
  mrn_data <- as.matrix(read.csv(file_name, row.names = 1))
  t(mrn_data)[lower.tri(mrn_data)]
}

#join and save rqa and inter-network mrn results for each mode
net_pairings <- apply(combn(net_names, 2), 2, paste, collapse = "_")

for (mode in modes) {
  mode_dir <- file.path(casnet_dir, paste0("mode", mode))

  #one row per subject; columns = all networks for metric 1, then metric 2, ...
  results <- t(sapply(subject_ids, function(id) {
    layer_data <- read.csv(file.path(mode_dir, sprintf("%s_layer_results_mode%d.csv", id, mode)))
    unlist(layer_data[, rqa_vars])
  }))
  combos <- expand.grid(net_names, rqa_vars)
  colnames(results) <- paste0(combos$Var1, "_", combos$Var2, "_mode", mode)

  mrn_results <- t(sapply(subject_ids, function(id) {
    read_upper_tri(file.path(mode_dir, sprintf("%s_MRN_results_mode%d.csv", id, mode)))
  }))
  colnames(mrn_results) <- paste0(net_pairings, "_MI_mode", mode)

  joined <- data.frame(subject_id = subject_ids, results, mrn_results, check.names = FALSE)
  write.csv(joined, file.path(output_dir, sprintf("MVMD_casnet_data_mode%d.csv", mode)), row.names = FALSE)
}

#sort intermode MRN data and catch missing subjects
mode_names <- paste0("mode", modes)
mode_pairings <- apply(combn(mode_names, 2), 2, paste, collapse = "_")
n_pairs <- length(mode_pairings)

#rows = subjects, columns = (mode pairings * 17 networks)
intermode_results <- matrix(NA, length(subject_ids), n_pairs * length(net_names))
missing_subjects <- c() #list to store IDs of participants with missing data

for (i in seq_along(subject_ids)) {
  id <- subject_ids[i]

  for (j in seq_along(net_names)) {
    net <- net_names[j]
    file_name <- file.path(casnet_dir, "intermode_MRN", sprintf("%s_MRN_%s_results.csv", id, net))

    #use tryCatch to handle missing files
    tryCatch({
      cols <- ((j - 1) * n_pairs + 1):(j * n_pairs)
      intermode_results[i, cols] <- read_upper_tri(file_name)
    }, error = function(e) {
      message(paste("Error: Could not read file for ID:", id, "Network:", net))
      if (!(id %in% missing_subjects)) {
        missing_subjects <<- c(missing_subjects, id)
      }
    })
  }
}

#full column names (mode pairing + network name)
all_colnames <- expand.grid(mode_pairings, net_names)
colnames(intermode_results) <- paste0(all_colnames$Var1, "_MI_", all_colnames$Var2)

intermode_joined <- data.frame(subject_id = subject_ids, intermode_results, check.names = FALSE)
write.csv(intermode_joined, file.path(output_dir, "MVMD_intermode_MRN.csv"), row.names = FALSE)

#view any missing subjects
print(missing_subjects)

#MVMD central frequency summary statistics
freq <- read.csv(file.path(output_dir, "MVMD_central_frequencies.csv"), header = TRUE)
freq_modes <- freq[, grep("^Mode\\d+$", colnames(freq))] * 1000

freq_summary <- data.frame(
  Mode     = colnames(freq_modes),
  Mean     = sapply(freq_modes, mean),
  SD       = sapply(freq_modes, sd),
  Skew     = sapply(freq_modes, psych::skew),
  Kurtosis = sapply(freq_modes, e1071::kurtosis),
  row.names = NULL
)
print(freq_summary)
write.csv(freq_summary, file.path(output_dir, "MVMD_central_frequency_summary.csv"), row.names = FALSE)
