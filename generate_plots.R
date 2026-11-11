#generate_plots.R
#creates boxplots and ggseg renders of mean differences across recurrence quantification analysis,
#multiplex recurrence network, and relative energy metrics

#load packages
library(tidyverse)
library(patchwork)
library(ggseg)
library(ggsegYeo2011)

#folders for input data (outputs of run_MVMD.m and data_sorting.R) and figures
data_dir <- "output"
fig_dir  <- "figures"
dir.create(fig_dir, showWarnings = FALSE)

#publication theme
theme_pub <- function(base_size = 12) {
  theme_classic(base_size = base_size) +
    theme(
      legend.position = "none",
      axis.title = element_text(face = "bold"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold")
    )
}

#colors for 17-network schaefer atlas
schaefer_colors <- c(
  "#E69423","#87324B","#778CB1","#FFFF01","#CD3E50","#000083",
  "#4A9B3D","#00760C","#DCF8A2","#7A8731","#C43AFB","#FF99D6",
  "#4681B4","#2BCCA0","#0929FA","#7C1383","#FA0000"
)

#load data across modes
modes <- 2:5
all_data <- map_dfr(modes, function(m) {
  file_path <- file.path(data_dir, sprintf("MVMD_casnet_data_mode%d.csv", m))
  read.csv(file_path) %>%
    select(2:290) %>%
    mutate(ID = row_number(), Mode = sprintf("%02d", m))
})

#load intermode IMI
intermode <- read.csv(file.path(data_dir, "MVMD_intermode_MRN.csv"))

net_names <- c("VisCent", "VisPeri", "SomMotA", "SomMotB", "DorsAttnA", "DorsAttnB", 
               "SalVentAttnA", "SalVentAttnB", "LimbicB", "LimbicA", "ContA", 
               "ContB", "ContC", "DefaultA", "DefaultB", "DefaultC", "TempPar")

#helper to pivot and compute means across network groups
inter_long <- intermode %>%
  mutate(ID = row_number()) %>%
  pivot_longer(cols = 2:103, names_to = "Variable", values_to = "Value") %>%
  mutate(Network = factor(rep(net_names, each = 6, length.out = n()))) %>%
  group_by(ID, Network) %>%
  summarise(Average_IMI = mean(Value, na.rm = TRUE), .groups = "drop")

#plot average IMI across networks
p_network_imi <- ggplot(inter_long, aes(x = reorder(Network, Average_IMI), y = Average_IMI, fill = Network)) +
  geom_boxplot(outlier.size = 0.5, alpha = 0.8) +
  scale_fill_manual(values = schaefer_colors) +
  labs(x = "Network", y = "Average IMI") +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#average IMI across modes
all_nets <- map_dfr(modes, function(m) {
  data_m <- read.csv(file.path(data_dir, sprintf("MVMD_casnet_data_mode%d.csv", m)), header = T)[, 2:290]
  data.frame(Mode = sprintf("%02d", m), value = rowMeans(data_m[, 154:289], na.rm = TRUE))
}) %>%
  mutate(Mode = factor(Mode))

p_mode_imi <- ggplot(all_nets, aes(x = Mode, y = value, fill = Mode)) +
  geom_boxplot(width = 0.5, outlier.size = .5, alpha = 0.85) +
  scale_fill_brewer(palette = "Set2") +
  labs(x = "Mode", y = "Average IMI") +
  coord_cartesian(ylim = range(all_nets$value, na.rm = TRUE)) + 
  theme_pub()

#combine plot panels
IMI_plot <- p_network_imi + p_mode_imi + plot_layout(
    widths = c(2.5, 1), 
    guides = "keep"
  ) +
  plot_annotation(
    tag_levels = 'A' 
  ) & 
  theme(
    plot.tag = element_text(size = 14, face = "bold")
  )
IMI_plot
ggsave(file.path(fig_dir, "fig7_IMI_plot.png"), IMI_plot, width = 8, height = 6, units = "in", dpi = 300)


#define rqa metrics titles
rqa_vars <- c("DET", "MEAN_dl", "LAM_vl", "TT_vl", "ENT_dl", "ENT_vl")

rqa_labels <- c(
  "DET"     = "Determinism",
  "MEAN_dl" = "Mean Diagonal Line Length",
  "LAM_vl"  = "Laminarity",
  "TT_vl"   = "Trapping Time",
  "ENT_dl"  = "Diagonal Line Entropy",
  "ENT_vl"  = "Vertical Line Entropy"
)

rqa_long <- all_data %>%
  pivot_longer(
    cols = matches("_(DET|MEAN_dl|LAM_vl|TT_vl|ENT_dl|ENT_vl)_mode"),
    names_to = c("Network", "Metric", "Mode_Num"),
    names_pattern = "^(.*)_(DET|MEAN_dl|LAM_vl|TT_vl|ENT_dl|ENT_vl)_mode(\\d+)$",
    values_to = "Value"
  ) %>%
  group_by(ID, Mode, Metric) %>%
  summarise(Mean_Value = mean(Value, na.rm = TRUE), .groups = "drop")

#facetted 2x3 rqa grid plot
p_rqa_grid <- ggplot(rqa_long, aes(x = Mode, y = Mean_Value, fill = Mode)) +
  geom_boxplot(outlier.size = 0.5, alpha = 0.85) +
  facet_wrap(
    ~ Metric, 
    ncol = 2, 
    scales = "free_y", 
    labeller = labeller(Metric = rqa_labels)
  ) +
  scale_fill_brewer(palette = "Set2") +
  labs(x = "Mode", y = "Value") +
  theme_pub() +
  theme(strip.text = element_text(size = 11, face = "bold"))

p_rqa_grid
ggsave(file.path(fig_dir, "fig5_rqa_grid.png"), p_rqa_grid, width = 8, height = 6, units = "in", dpi = 300)


#define recurrence network titles
net_vars <- c("GlobalClustering", "GlobalEfficiency", "NetworkTransitivity")

net_labels <- c(
  "GlobalClustering"   = "Global Clustering",
  "GlobalEfficiency"   = "Global Efficiency",
  "NetworkTransitivity" = "Network Transitivity"
)

net_long <- all_data %>%
  pivot_longer(
    cols = matches("_(GlobalClustering|GlobalEfficiency|NetworkTransitivity)_mode"),
    names_to = c("Network", "Metric", "Mode_Num"),
    names_pattern = "^(.*)_(GlobalClustering|GlobalEfficiency|NetworkTransitivity)_mode(\\d+)$",
    values_to = "Value"
  ) %>%
  group_by(ID, Mode, Metric) %>%
  summarise(Mean_Value = mean(Value, na.rm = TRUE), .groups = "drop")

#facetted recurrence network plots
p_net_grid <- ggplot(net_long, aes(x = Mode, y = Mean_Value, fill = Mode)) +
  geom_boxplot(outlier.size = 0.5, alpha = 0.85, width = 0.5) +
  facet_wrap(
    ~ Metric, 
    ncol = 1, 
    scales = "free_y", 
    labeller = labeller(Metric = net_labels)
  ) +
  scale_x_discrete(expand = expansion(mult = c(0.3, 0.3))) +
  scale_fill_brewer(palette = "Set2") +
  labs(x = "Mode", y = "Value") +
  theme_pub() +
  theme(
    strip.text = element_text(size = 11, face = "bold"),
    panel.spacing = unit(1, "lines") 
  ) + 
  theme(plot.margin = margin(5,75,5,75))
p_net_grid
ggsave(file.path(fig_dir, "fig6_net_grid.png"), p_net_grid, width = 8, height = 8, units = "in", dpi = 300)

#ggseg relative energy plots
#load relative energy dataset
energy <- read.csv(file.path(data_dir, "MVMD_network_relative_energy.csv"), header = TRUE)

#yeo17 atlas region name for each network
#NOTE: check these against unique(as.data.frame(yeo17())$region) for your ggsegYeo2011 version
yeo17_regions <- c(
  VisCent = "Visual A", VisPeri = "Visual B",
  SomMotA = "Somatomotor A", SomMotB = "Somatomotor B",
  DorsAttnA = "Dorsal Attention A", DorsAttnB = "Dorsal Attention B",
  SalVentAttnA = "Salience/VentAttn A", SalVentAttnB = "Salience/VentAttn B",
  LimbicA = "Limibic A", LimbicB = "Limibic B",
  ContA = "Control A", ContB = "Control B", ContC = "Control C",
  TempPar = "Temporal Parietal",
  DefaultA = "Default A", DefaultB = "Default B", DefaultC = "Default C"
)

#mean relative energy of each network for each analysis mode (columns mode2, mode3, ...)
df_brain <- energy %>%
  pivot_longer(
    cols = matches("_RE_mode\\d+$"),
    names_to = c("net_name", "Mode_Num"),
    names_pattern = "^(.*)_RE_mode(\\d+)$",
    values_to = "Value"
  ) %>%
  filter(as.numeric(Mode_Num) %in% modes) %>%
  group_by(net_name, Mode_Num) %>%
  summarise(Value = mean(Value, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Mode_Num, values_from = Value, names_prefix = "mode") %>%
  mutate(region = yeo17_regions[net_name], net_name = factor(net_name))

#atlas reference plot (top panel)
p_atlas <- ggplot(df_brain) +
  geom_brain(
    atlas = yeo17(), 
    hemi = "left",
    view = c("medial", "lateral"),
    aes(fill = as.factor(net_name)), 
    color = "black", 
    size = .75
  ) + 
  scale_fill_manual(values = schaefer_colors, na.value = "darkgray") +
  labs(title = "Yeo 17-Network Atlas") +
  theme_void() +
  theme(
    plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
    legend.position = "none"
  )

#helper to generate plots of relative energy
plot_brain_mode <- function(data, fill_var, title) {
  ggplot(data) +
    geom_brain(
      atlas = yeo17(),
      hemi = "left",
      view = c("lateral", "medial"),
      aes(fill = .data[[fill_var]]),
      color = "black", 
      size = .75
    ) +
    scale_fill_viridis_c(option = "magma", name = "Relative\nEnergy") +
    labs(title = title) +
    theme_void() +
    theme(
      plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
      legend.key.height = unit(0.3, "cm"),
      legend.key.width = unit(0.25, "cm"),
      legend.title = element_text(size = 8)
    )
}

#create list of energy plots and combine
mode_plots <- map(2:5, ~plot_brain_mode(df_brain, paste0("mode", .x), paste("Mode", .x)))
combined_brain_plots <- wrap_plots(c(list(p_atlas), mode_plots), ncol = 1)
combined_brain_plots
ggsave(file.path(fig_dir, "fig4_combined_brain_plots.png"), combined_brain_plots, width = 8, height = 8, units = "in", dpi = 300)


#calculate cortical means across modes 1 to 10
energy_long <- energy %>%
  mutate(ID = row_number()) %>%
  pivot_longer(
    cols = matches("_RE_mode\\d+$"),
    names_to = c("Network", "Mode_Num"),
    names_pattern = "^(.*)_RE_mode(\\d+)$",
    values_to = "Value"
  ) %>%
  mutate(
    Mode = sprintf("%02d", as.numeric(Mode_Num)) # Format as 01, 02, ..., 10
  ) %>%
  group_by(ID, Mode) %>%
  summarise(Mean_Energy = mean(Value, na.rm = TRUE), .groups = "drop")

#create boxplot
p_energy_all <- ggplot(energy_long, aes(x = Mode, y = Mean_Energy, fill = Mode)) +
  geom_boxplot(width = 0.5, outlier.size = 0.5, alpha = 0.85) +
  scale_fill_viridis_d(option = "turbo", direction = -1) +
  scale_x_discrete(expand = expansion(mult = c(0.05, 0.05))) +
  labs(
    x = "Mode",
    y = "Relative Energy"
  ) +
  theme_pub()

p_energy_all
ggsave(file.path(fig_dir, "fig3_energy_plots.png"), p_energy_all, width = 8, height = 6, units = "in", dpi = 300)

