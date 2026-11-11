#rqa_mrn_metrics.R
#recurrence quantification analysis and multiplex recurrence network
#analysis of the network-averaged MVMD mode timeseries created by run_MVMD.m
#  part 1 - for each mode: RQA and recurrence network measures for each of the
#           17 networks, plus an inter-network MRN (layers = networks)
#  part 2 - for each network: an inter-mode MRN (layers = modes)

#load packages
library(casnet)
library(igraph)

#user settings
output_dir <- "output"                                   #folder used by run_MVMD.m
net_ts_dir <- file.path(output_dir, "network_timeseries")
casnet_dir <- file.path(output_dir, "casnet_results")
modes      <- 2:5                                        #modes to analyze
target_rr  <- .05                                        #radius based on 5% recurrence

net_names <- c("VisCent", "VisPeri", "SomMotA", "SomMotB", "DorsAttnA", "DorsAttnB",
               "SalVentAttnA", "SalVentAttnB", "LimbicB", "LimbicA", "ContA", "ContB",
               "ContC", "DefaultA", "DefaultB", "DefaultC", "TempPar")

#list of subject ids
subject_ids <- read.csv(file.path(output_dir, "subject_ids.csv"), colClasses = "character")$subject_id

#helper functions
read_mode_ts <- function(id, mode) {
  mode_data <- read.csv(file.path(net_ts_dir, sprintf("%s_mode%d.csv", id, mode)), header = FALSE)
  colnames(mode_data) <- net_names
  mode_data
}

fit_rn <- function(y, emDim, emLag) {
  casnet::rn(y1 = y, emDim = emDim, emLag = emLag, emRad = NA, targetValue = target_rr)
}

rn_graph <- function(RN) {
  igraph::graph_from_adjacency_matrix(RN, mode = "undirected", diag = FALSE)
}

#part 1: recurrence quantification analysis and inter-network MRN
for (mode in modes) {
  mode_dir <- file.path(casnet_dir, paste0("mode", mode))
  dir.create(mode_dir, recursive = TRUE, showWarnings = FALSE)

  for (id in subject_ids) {
    mode_data <- read_mode_ts(id, mode)

    #matrix to save results
    network_results <- matrix(NA, length(net_names), 48, dimnames = list(net_names, NULL))
    optimdims <- numeric(length(net_names))
    optimlags <- numeric(length(net_names))

    for (net in seq_along(net_names)) {
      #estimate parameters
      network <- est_parameters(y = mode_data[, net], silent = TRUE)
      optimdims[net] <- network$optimDim
      optimlags[net] <- network$optimLag

      #fit recurrence matrix
      RN_network <- fit_rn(mode_data[, net], network$optimDim, network$optimLag)
      RQA_network <- casnet::rp_measures(RN_network, silent = TRUE)
      network_results[net, 1:42] <- as.matrix(RQA_network)

      #recurrence network graph
      g_network <- rn_graph(RN_network)
      igraph::V(g_network)$size <- igraph::degree(g_network)

      #recurrence network measures
      network_measures <- casnet::rn_measures(g_network, silent = TRUE)
      network_results[net, 43:48] <- as.matrix(network_measures$graph_measures)
    }
    colnames(network_results) <- c(colnames(RQA_network), colnames(network_measures$graph_measures))
    write.csv(network_results, file.path(mode_dir, sprintf("%s_layer_results_mode%d.csv", id, mode)))

    #save embedding parameters (reused for the inter-mode MRN in part 2)
    write.csv(data.frame(network = net_names, emDim = optimdims, emLag = optimlags),
              file.path(mode_dir, sprintf("%s_embedding_mode%d.csv", id, mode)), row.names = FALSE)

    #multiplex recurrence network (across networks)
    emDim <- floor(median(optimdims))
    emLag <- floor(median(optimlags))

    layers <- lapply(net_names, function(net) rn_graph(fit_rn(mode_data[[net]], emDim, emLag)))
    names(layers) <- paste0("g", seq_along(layers))
    MRN <- mrn(layers = layers)
    MRN_results <- as.matrix(MRN$interlayerMI[[1]])
    write.csv(MRN_results, file.path(mode_dir, sprintf("%s_MRN_results_mode%d.csv", id, mode)))
  }
}

#part 2: inter-mode MRN (across modes, within each network)
intermode_dir <- file.path(casnet_dir, "intermode_MRN")
dir.create(intermode_dir, recursive = TRUE, showWarnings = FALSE)

for (id in subject_ids) {
  mode_ts <- lapply(modes, function(mode) read_mode_ts(id, mode))
  embedding <- lapply(modes, function(mode) {
    read.csv(file.path(casnet_dir, paste0("mode", mode), sprintf("%s_embedding_mode%d.csv", id, mode)))
  })

  for (j in seq_along(net_names)) {
    net <- net_names[j]

    #embedding parameters: median of this network's optimal parameters across modes
    emDim <- floor(median(sapply(embedding, function(e) e$emDim[j])))
    emLag <- floor(median(sapply(embedding, function(e) e$emLag[j])))

    layers <- lapply(mode_ts, function(ts) rn_graph(fit_rn(ts[[net]], emDim, emLag)))
    names(layers) <- paste0("mode", modes)
    MRN <- mrn(layers = layers)
    MRN_results <- as.matrix(MRN$interlayerMI[[1]])
    write.csv(MRN_results, file.path(intermode_dir, sprintf("%s_MRN_%s_results.csv", id, net)))
  }
}
