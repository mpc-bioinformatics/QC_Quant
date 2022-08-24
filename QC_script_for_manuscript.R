
library(openxlsx)
library(ggpubr)     ### for arranging plots
library(cowplot)    ### for aligning plots
library(tidyverse)
library(ggplot2)
library(ComplexHeatmap)

path <- ""  ## add path to data set here
setwd(path)

data_path <- "proteinGroups.txt"
output_path <- "graphics/"
RScript_path <-"/QC_scripts/"


D <- read.table(data_path, header = TRUE, sep = "\t")
### remove decoys and if necessary contaminants and "only identified by site":
D <- D[D$Reverse == "", ]
#D <- D[D$Potential.contaminant == "", ]
D <- D[D$Only.identified.by.site == "", ]

id <- D$Protein.IDs
D_LFQ <- D[grepl("LFQ.", colnames(D))]
#D_LFQ <- D_LFQ[, c(1:4, 20, 5:19)]
colnames(D_LFQ) <- c(paste0("Std_in_sol_", 1:5), paste0("Rapid_in_sol_", 1:5), paste0("Std_FASP_", 1:5), paste0("Rapid_FASP_", 1:5))

### raw intensities
D_RI <- D[, c(122:141)]
#D_RI <- D_RI[, c(1:4, 20, 5:19)]
colnames(D_RI) <- c(paste0("Std_in_sol_", 1:5), paste0("Rapid_in_sol_", 1:5), paste0("Std_FASP_", 1:5), paste0("Rapid_FASP_", 1:5))


### check if output path exists
if (!dir.exists(paste0(path, output_path))){
  dir.create(paste0(path, output_path))
}
###
log_data = TRUE

### choose normalization: "nonorm" = without normalization, "median" = Median normalization,
###                       "loess" = LOESS normalization, "quantile" = Quantile normalization
normalization = "nonorm" #"nonorm
suffix <- normalization

group <- c(rep("Std_in_sol", 5), rep("Rapid_in_sol", 5), rep("Std_FASP", 5), rep("Rapid_FASP", 5))
group <- factor(group, levels = c("Std_in_sol", "Rapid_in_sol", "Std_FASP", "Rapid_FASP"))
group_colours <- c("#00BFC4", "deeppink", "darkgoldenrod1", "#7CAE00")


################################################################################

groupvar_name = "Group"

plot_device = "tiff"

plot_dpi = 600


plot_height_validvalueplot = 10
plot_width_validvalueplot  = 15
plot_height_boxplots = 10
plot_width_boxplots  = 15
plot_height_PCA = 15
plot_width_PCA  = 25
plot_height_MAPlot = 15
plot_width_MAPlot  = 20
plot_height_hist = 35
plot_width_hist  = 38


plot_xlim_PCA <- NULL
plot_ylim_PCA <- NULL
plot_xlim_hist <- NULL
plot_ylim_hist <- NULL

sample_filter <- NULL

zero_to_NA = TRUE

log_base = 2

################################################################################

library(openxlsx)   ### for reading and writing xlsx files
library(limma)      ### e.g. for normalization
library(tidyverse)  ### tidyverse functionalities + ggplot2
library(scales)     ### for colours
library(matrixStats)### e.g. for rowVars used in PCA plot function
library(cowplot)
library(affy)       ### for MA-Plots
library(beepr)      ### for sound after completion of MA-Plots

setwd(path)

source(paste0(RScript_path, "PCA_plot_v1_2.R"))
source(paste0(RScript_path, "MA_Plots_v1_2.R"))
source(paste0(RScript_path, "ValidValue_Plot_v1_2.R"))
source(paste0(RScript_path, "Boxplot_v1_2.R"))
source(paste0(RScript_path, "automatedNormalization_v1_2.R"))

################################################################################

output_path <- paste0(path, output_path)


if(zero_to_NA) {
  D_LFQ[D_LFQ == 0] <- NA
  D_RI[D_RI == 0] <- NA
}

if(log_data) {
  D_LFQ <- log(D_LFQ, base = log_base)
  D_RI <- log(D_RI, base = log_base)
}

nr_groups <- length(levels(group))

### set group colours if not previously specified
if (is.null(group_colours) & nr_groups >= 1) group_colours <- scales::hue_pal()(nr_groups)


### normalize data:
D_loess_gw <- automatedNormalization(DATA = D_RI, method = "loess", log = FALSE,
                                  id_columns = id, output_path = output_path, suffix = "loess_proteins_gw",
                                  groupwise = TRUE, group = group)
D_loess_gw <- data.frame(D_loess_gw)
D_median_gw <- automatedNormalization(DATA = D_RI, method = "median", log = FALSE,
                                     id_columns = id, output_path = output_path, suffix = "median_proteins_gw",
                                     groupwise = TRUE, group = group)
D_median_gw <- data.frame(D_median_gw)
D_quantile_gw <- automatedNormalization(DATA = D_RI, method = "quantile", log = FALSE,
                                      id_columns = id, output_path = output_path, suffix = "quantile_proteins_gw",
                                      groupwise = TRUE, group = group)
D_quantile_gw <- data.frame(D_quantile_gw)


### convert data to long format
to_long_format <- function(D) {
  D_long <- pivot_longer(data = D, cols = 1:ncol(D))
  D_long$group <- factor(substr(D_long$name, 1, nchar(D_long$name)-2),
                         levels = c("Std_in_sol","Rapid_in_sol", "Std_FASP", "Rapid_FASP"))
  D_long$name <- factor(D_long$name,
                        levels = c(paste0("Std_in_sol_", 1:5), paste0("Rapid_in_sol_", 1:5), paste0("Std_FASP_", 1:5), paste0("Rapid_FASP_", 1:5)))
  return(D_long)
}

D_LFQ_long <- to_long_format(D_LFQ)
D_RI_long <- to_long_format(D_RI)
D_loess_gw_long <- to_long_format(D_loess_gw)
D_median_gw_long <- to_long_format(D_median_gw)
D_quantile_gw_long <- to_long_format(D_quantile_gw)

################################################################################
################################################################################
##### PCAs

PCA_RI_nonorm <- PCA_Plot(X = D_RI, groupvar1 = group, groupvar2 = NULL,
         groupvar1_name = groupvar_name, groupvar2_name = NULL,
         plot_device = plot_device, group_colours = group_colours,
         plot_height = plot_height_PCA, plot_width = plot_width_PCA,
         log_data = FALSE, log_base = log_base,  scale. = TRUE,
         impute = FALSE, impute_method = "mean", propNA = 0,
         point.size = 4, base_size = 20,
         returnPCA = FALSE, title = NULL,
         output_path = output_path, suffix = "nonorm_proteins",
         ylim = plot_ylim_PCA, xlim = plot_xlim_PCA)

PCA_LFQ <- PCA_Plot(X = D_LFQ, groupvar1 = group, groupvar2 = NULL,
                          groupvar1_name = groupvar_name, groupvar2_name = NULL,
                          plot_device = plot_device, group_colours = group_colours,
                          plot_height = plot_height_PCA, plot_width = plot_width_PCA,
                          log_data = FALSE, log_base = log_base,  scale. = TRUE,
                          impute = FALSE, impute_method = "mean", propNA = 0,
                          point.size = 4, base_size = 20,
                          returnPCA = FALSE, title = NULL,
                          output_path = output_path, suffix = "LFQ_proteins",
                          ylim = plot_ylim_PCA, xlim = plot_xlim_PCA)


PCA_loess_gw <- PCA_Plot(X = D_loess_gw, groupvar1 = group, groupvar2 = NULL,
                      groupvar1_name = groupvar_name, groupvar2_name = NULL,
                      plot_device = plot_device, group_colours = group_colours,
                      plot_height = plot_height_PCA, plot_width = plot_width_PCA,
                      log_data = FALSE, log_base = log_base,  scale. = TRUE,
                      impute = FALSE, impute_method = "mean", propNA = 0,
                      point.size = 4, base_size = 20,
                      returnPCA = FALSE, title = NULL,
                      output_path = output_path, suffix = "loess_proteins_gw",
                      ylim = plot_ylim_PCA, xlim = plot_xlim_PCA)

PCA_median_gw <- PCA_Plot(X = D_median_gw, groupvar1 = group, groupvar2 = NULL,
                         groupvar1_name = groupvar_name, groupvar2_name = NULL,
                         plot_device = plot_device, group_colours = group_colours,
                         plot_height = plot_height_PCA, plot_width = plot_width_PCA,
                         log_data = FALSE, log_base = log_base,  scale. = TRUE,
                         impute = FALSE, impute_method = "mean", propNA = 0,
                         point.size = 4, base_size = 20,
                         returnPCA = FALSE, title = NULL,
                         output_path = output_path, suffix = "median_proteins_gw",
                         ylim = plot_ylim_PCA, xlim = plot_xlim_PCA)

PCA_quantile_gw <- PCA_Plot(X = D_quantile_gw, groupvar1 = group, groupvar2 = NULL,
                          groupvar1_name = groupvar_name, groupvar2_name = NULL,
                          plot_device = plot_device, group_colours = group_colours,
                          plot_height = plot_height_PCA, plot_width = plot_width_PCA,
                          log_data = FALSE, log_base = log_base,  scale. = TRUE,
                          impute = FALSE, impute_method = "mean", propNA = 0,
                          point.size = 4, base_size = 20,
                          returnPCA = FALSE, title = NULL,
                          output_path = output_path, suffix = "quantile_proteins_gw",
                          ylim = plot_ylim_PCA, xlim = plot_xlim_PCA)


################################################################################
## valid value plots

vvpl_RI <- ValidValuePlot(X = D_RI_long, groupvar_name = groupvar_name, sample_filter = sample_filter,
               plot_device = plot_device, group_colours = group_colours,
               plot_height = plot_height_validvalueplot, plot_width = plot_width_validvalueplot,
               plot_dpi = plot_dpi, suffix = "nonorm_proteins", output_path = output_path,
               ylim = c(0, 600), title = "Raw Intensities")

vvpl_LFQ <- ValidValuePlot(X = D_LFQ_long, groupvar_name = groupvar_name, sample_filter = sample_filter,
                           plot_device = plot_device, group_colours = group_colours,
                           plot_height = plot_height_validvalueplot, plot_width = plot_width_validvalueplot,
                           plot_dpi = plot_dpi, suffix = "LFQ_proteins", output_path = output_path,
                           ylim = c(0, 600), title = "LFQ Values")


vvpl_panel <- ggpubr::ggarrange(vvpl_RI$pl_valid_values, vvpl_LFQ$pl_valid_values,
                                common.legend = FALSE, legend = "none",
                              labels = "AUTO",
                              font.label = list(size = 20, color = "black", face = "bold", family = NULL))
print(vvpl_panel)

ggsave(paste0(output_path, "validvalues_proteins_panel.tif"), plot = vvpl_panel, width = 30, height = 13,
              device = "tiff", units = "cm", dpi = 600)


################################################################################
### boxplots

bp_RI <- Boxplots(X = D_RI_long, groupvar_name = groupvar_name, sample_filter = sample_filter,
         plot_device = plot_device, group_colours = group_colours,
         plot_height = plot_height_boxplots, plot_width = plot_width_boxplots,
         plot_dpi = plot_dpi, log_data = FALSE, log_base = log_base,
         suffix = "nonorm_proteins", method = "boxplot", base_size = 25, legend.key.size = unit(1.5, 'lines'))

bp_LFQ <- Boxplots(X = D_LFQ_long, groupvar_name = groupvar_name, sample_filter = sample_filter,
                   plot_device = plot_device, group_colours = group_colours,
                   plot_height = plot_height_boxplots, plot_width = plot_width_boxplots,
                   plot_dpi = plot_dpi, log_data = FALSE, log_base = log_base,
                   suffix = "LFQ_proteins", method = "boxplot", base_size = 25, legend.key.size = unit(1.5, 'lines'))

bp_loess_gw <- Boxplots(X = D_loess_gw_long, groupvar_name = groupvar_name, sample_filter = sample_filter,
                  plot_device = plot_device, group_colours = group_colours,
                  plot_height = plot_height_boxplots, plot_width = plot_width_boxplots,
                  plot_dpi = plot_dpi, log_data = FALSE, log_base = log_base,
                  suffix = "loess_proteins_gw", method = "boxplot", base_size = 25, legend.key.size = unit(1.5, 'lines'))

bp_median_gw <- Boxplots(X = D_median_gw_long, groupvar_name = groupvar_name, sample_filter = sample_filter,
                        plot_device = plot_device, group_colours = group_colours,
                        plot_height = plot_height_boxplots, plot_width = plot_width_boxplots,
                        plot_dpi = plot_dpi, log_data = FALSE, log_base = log_base,
                        suffix = "median_proteins_gw", method = "boxplot", base_size = 25, legend.key.size = unit(1.5, 'lines'))

bp_quantile_gw <- Boxplots(X = D_quantile_gw_long, groupvar_name = groupvar_name, sample_filter = sample_filter,
                        plot_device = plot_device, group_colours = group_colours,
                        plot_height = plot_height_boxplots, plot_width = plot_width_boxplots,
                        plot_dpi = plot_dpi, log_data = FALSE, log_base = log_base,
                        suffix = "quantile_proteins_gw", method = "boxplot", base_size = 25, legend.key.size = unit(1.5, 'lines'))

################################################################################
### MA-Plots

MAPlots(X = D_RI, log = TRUE, alpha = FALSE, output_path = output_path,
        suffix = "RawIntensities_proteins")

MAPlots(X = D_LFQ, log = TRUE, alpha = FALSE, output_path = output_path,
        suffix = "LFQ_proteins")

MAPlots(X = D_loess_gw, log = TRUE, alpha = FALSE, output_path = output_path,
        suffix = "loess_proteins_gw")

MAPlots(X = D_median_gw, log = TRUE, alpha = FALSE, output_path = output_path,
        suffix = "median_proteins_gw")

MAPlots(X = D_quantile_gw, log = TRUE, alpha = FALSE, output_path = output_path,
        suffix = "quantile_proteins_gw")


################################################################################
################################################################################
################################################################################
################################################################################
### Input of peptide data

D_pep <- read.table("peptides.txt", header = TRUE, sep = "\t")
id_pep <- D_pep$Sequence
D_pep_LFQ <- D_pep[,118:137]  # LFQ normalized
D_pep_LFQ[D_pep_LFQ == 0] <- NA
#D_pep_LFQ <- D_pep_LFQ[, c(1:4, 20, 5:19)]
colnames(D_pep_LFQ) <- c(paste0("Std_in_sol_", 1:5), paste0("Rapid_in_sol_", 1:5), paste0("Std_FASP_", 1:5), paste0("Rapid_FASP_", 1:5))

D_pep_RI <- D_pep[,87:106]   # Raw Intensities
D_pep_RI[D_pep_RI == 0] <- NA
#D_pep_RI <- D_pep_RI[, c(1:4, 20, 5:19)]
colnames(D_pep_RI) <- c(paste0("Std_in_sol_", 1:5), paste0("Rapid_in_sol_", 1:5), paste0("Std_FASP_", 1:5), paste0("Rapid_FASP_", 1:5))

if(log_data) {
  D_pep_LFQ <- log(D_pep_LFQ, base = log_base)
  D_pep_RI <- log(D_pep_RI, base = log_base)
}


D_pep_RI_long <- to_long_format(D_pep_RI)


################################################################################
#### valid value plot peptides


vvpl_RI <- ValidValuePlot(X = D_pep_RI_long, groupvar_name = groupvar_name, sample_filter = sample_filter,
                          plot_device = plot_device, group_colours = group_colours,
                          plot_height = plot_height_validvalueplot, plot_width = plot_width_validvalueplot,
                          plot_dpi = plot_dpi, suffix = "nonorm_peptides", output_path = output_path,
                          title = "Raw Intensities")

################################################################################
### correlation heatmap

library(circlize)

K_pep <- cor(D_pep_RI, use = "pairwise.complete.obs", method = "pearson")
K_pep_mat <- as.matrix(K_pep)

K_prot <- cor(D_RI, use = "pairwise.complete.obs", method = "pearson")
K_prot_mat <- as.matrix(K_prot)

col_fun <- colorRamp2(c(0.7, 0.85, 1), c("blue", "white", "red"))
lgd = Legend(col_fun = col_fun, title = "Pearson corr.", title_position = "leftcenter-rot")

ht <<- Heatmap(K_pep_mat,
               name= "Pearson corr.",
               cluster_rows = TRUE,
               clustering_distance_rows = "pearson", clustering_method_rows = "complete",
               cluster_columns = TRUE, clustering_distance_columns = "pearson",
               clustering_method_columns = "complete",
               # bottom_annotation = HeatmapAnnotation(Group = as.matrix(group), col = bottom_groupcolours,
               #                                       annotation_legend_param = list(nrow = 1)),
              # row_labels = row_labels,
               col = c("blue", "white", "red"),
              heatmap_legend_param = list(direction = "vertical", title_position = "leftcenter-rot"),
              show_heatmap_legend = FALSE)

tiff(paste0(output_path, "Corr_Heatmap_peptides.tif"), width = 15, height = 15, units = "cm", res = 600)
draw(ht)
draw(lgd, x = unit(0.9, "npc"), y = unit(0.1, "npc"))
dev.off()


ht2 <<- Heatmap(K_prot_mat,
               name= "Pearson corr.",
               cluster_rows = TRUE,
               clustering_distance_rows = "pearson", clustering_method_rows = "complete",
               cluster_columns = TRUE, clustering_distance_columns = "pearson",
               clustering_method_columns = "complete",
               # bottom_annotation = HeatmapAnnotation(Group = as.matrix(group), col = bottom_groupcolours,
               #                                       annotation_legend_param = list(nrow = 1)),
               # row_labels = row_labels,
               col = c("blue", "white", "red"),
               heatmap_legend_param = list(direction = "vertical", title_position = "leftcenter-rot"),
               show_heatmap_legend = FALSE)
ht2
draw(lgd, x = unit(0.9, "npc"), y = unit(0.15, "npc"))

tiff(paste0(output_path, "Corr_Heatmap_proteins.tif"), width = 15, height = 15, units = "cm", res = 600)
draw(ht2)
draw(lgd, x = unit(0.9, "npc"), y = unit(0.1, "npc"))
dev.off()


write.xlsx(as.data.frame(K_pep), "graphics/Correlations_peptides_new.xlsx", rowNames = TRUE)
write.xlsx(as.data.frame(K_prot), "graphics/Correlations_proteins_new.xlsx", rowNames = TRUE)
