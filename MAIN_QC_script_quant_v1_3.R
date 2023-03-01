################################################################################
#### R script for QC of quantitative proteomics data version 1.3

## !IMPORTANT!:
## Before using this script, make sure that your R and RStudio are updated to
## the newest version. We tested this script with R version 4.2.1 and RStudio
## version 2022.12.0 and cannot guarantee that this script will work with older
## versions.
## Important: R itself will not update if you use the "Check for updates"
## functionality of RStudio (this will only update RStudio).
## Please update R by downloading and installing the newest version from
## https://cran.r-project.org/

## Also, additional R packages need to be installed in their newest version.
## To install and/or update these packages, please run the following script
## (adjust the path to your computer, the script can be found in the R_source_code folder).
## If you are asked to restart R between,
## please answer "yes". It may then be necessary to run the package installation
## script again, until all packages are installed/updated.

source("package_installation_v1_3.R")

#### Instructions:
## Data in form of an xlsx file can be imported (rows = proteins/peptides, columns = samples).
## This can be the raw intensities or already normalized data.

## Requirements on the data file:
## The data file may contain columns that are not intensity columns (e.g. protein accession, gene name etc.),
## these can be set as "id_columns" later and will be skipped for the processing and plots.
## All other columns must belong to samples and contain intensities. All other columns
## that may be present in the original file must be either deleted or specified as "id_columns".

## Peptides or proteins that you do not want to be used for the diagnostic graphics
## (e.g. contaminants), have to be removed in Excel before loading the data into R.

## The columns describing the samples need to have column names of the following form:
## groupname_samplenumber
## e.g.: control_1, control_2, control_3, ...., patient_1, patient_2, patient_3, ...
## (!IMPORTANT!: The column names must not contain any further underscore (_)
## nor blanks or other special characters except dots).
## If you have more than 10 or 100 samples in at least one group, you need to
## add leading zeros to the samplenumbers, otherwise R will not be able to
## correctly sort the samples for the graphics.

## e.g. control_01, control_02, control_03, .... control10, control11 for >= 10 samples
## control_001, control_002, ..., control_010, control_011, ..., control_100, control_101 for >= 100 samples

## The analysis of data with more than two groups is possible without any problems.

## The number of MA-Plots that are generated can be very large, which leads to
## a long run time and a large output file. If more than 1000 plots will be
## generated, you are asked if you want to continue or skipt the MA-Plots.

## The following user setting can be changed and adjusted to the data.
## The basic settings have to be adjusted in almost every case.
## For the advanced settings, the default values will work in most cases,
## but can still be adjusted if you like, e.g. the colours used for the plots.


################################################################################
### Basic settings:


### set path to working directory (all other paths can be defined relatively to this)
path <- ""
setwd(path)
### set path to data file (relativ to the working directory specified above)
data_path <- "testdata.xlsx"
### set output_path were all created graphics are saved (relativ to the working directory specified above)
output_path <- "QC_results/"
### set RScripts path were the necessary RScripts are located (e.g. for PCA and MA plots)
### (relativ to the working directory specified above)
RScript_path <- "QC_scripts/"

### check if output path exists
if (!dir.exists(paste0(path, output_path))){
  dir.create(paste0(path, output_path))
}
###


### Specify which columns contain the peptide or protein intensities.
### All other columns will be deleted before generating the plots.
### e.g. 3:29 indicates all columns from 3 to 29. Columns 1 and 2 will be removed.
#intensity_columns = 3:29
intensity_columns = 5:54

### log_data: TRUE, if data should be log-transformed. FALSE if not (e.g. if data is already on log-level)
log_data = TRUE

### choose normalization: "nonorm" = without normalization, "median" = Median normalization,
###                       "loess" = LOESS normalization, "quantile" = Quantile normalization
normalization = "nonorm"

### set use_groups to FALSE, if you have no groups in the data
use_groups = TRUE

################################################################################
### Advanced settings:

### Here you can set colours for the different experimental groups, which are used
### for the graphics. This has to be in form of a vector with colour names
### or hexadecimal codes, one colour for each group (groups will be sorted alphabetically)
### Example:
### group_colours <- c("red", "green", "blue")
### named colours in R can be found here: http://www.stat.columbia.edu/~tzheng/files/Rcolor.pdf
### If you leave this as NULL, the default ggplot2 colours will be used.
group_colours <- NULL

### Name of the group variable that will show up as the legend title in the Plots
### (e.g. "Treatment", "Group", "Condition", ...)
groupvar_name = "Group"

### Here you can define the plot device that will be used to save the graphics
### possible values are "pdf" (default), "png", "svg", "tiff"
### Please note that the histograms and MA-Plots will always be saved as
### pdf because multiple pages are needed.
plot_device = "pdf"

### Resolution (in dots per inch) for the plots. This only does have an effect for
### png or tiff.
plot_dpi = 300

### Here you can adjust the heights and widths for the different plot types in cm,
### e.g. if you have many samples and the plots are too small.
plot_height_validvalueplot = 10
plot_width_validvalueplot  = 15
plot_height_boxplots = 10
plot_width_boxplots  = 15
plot_height_PCA = 15
plot_width_PCA  = 20
plot_height_MAPlot = 15
plot_width_MAPlot  = 15
plot_height_hist = 35
plot_width_hist  = 38

### Adjust x and y axis limits for the PCA plots and Histograms, e.g. plot_xlim_PCA = c(-5, 5).
### The default value NULL will calculate appropriate limits from the data.
plot_xlim_PCA <- NULL
plot_ylim_PCA <- NULL
plot_xlim_hist <- NULL
plot_ylim_hist <- NULL

### Change number of bins in histogram (if NULL, Sturges rule is applied, i.e. ceiling(1+log2(length(n))))
hist_numbins <- NULL

### Here you can specify samples that should be used for plotting by giving a vector
### of column names, e.g.
### sample_filter <- c("control_1", "control_3", "patient_2")
### If you leave this as NULL, all samples will be plotted.
### This will currently be used for the boxplots and Valid values plots ONLY!
sample_filter <- NULL

### Set symbols, that should be recognized as a missing value (do not include "0", this will be treated by zero_to_NA)
na_strings = c("NA", "NaN", "Filtered","#NV")
### TRUE, if 0 should be treated as a missing value. FALSE, if 0 should be left as a 0 (does not work in combination with log_data = TRUE).
zero_to_NA = TRUE

### Base of the logarithm used
log_base = 2

### Suffix to add to the file names (usually the normalization type, but can be changed, e.g. to "LFQ")
suffix <- normalization



################################################################################
################################################################################
################################################################################
#### start of QC Script!
#### only change if you know what you are doing! ;-)

################################################################################
### loading required R packages
library(openxlsx)   ### for reading and writing xlsx files
library(limma)      ### e.g. for normalization
library(tidyverse)  ### tidyverse functionalities + ggplot2
library(scales)     ### for colours
library(matrixStats)### e.g. for rowVars used in PCA plot function
library(cowplot)
library(affy)       ### for MA-Plots
#library(ggplus)     ### used in histogram function
library(beepr)      ### for sound after completion of MA-Plots
library(ggrepel)

setwd(path)

source(paste0(RScript_path, "PCA_plot_v1_3.R"))
source(paste0(RScript_path, "MA_Plots_v1_3.R"))
source(paste0(RScript_path, "ValidValue_Plot_v1_3.R"))
source(paste0(RScript_path, "automatedNormalization_v1_3.R"))


################################################################################
#### read in data file

output_path <- paste0(path, output_path)

D <- read.xlsx(data_path, na.strings = na_strings)

id <- D[, -intensity_columns]

D <- D[, intensity_columns]

if (use_groups) {
  group <- factor(limma::strsplit2(colnames(D), "_")[,1])
} else {
  group <- NULL
}

if(zero_to_NA) {
  D[D == 0] <- NA
}

if(log_data) {
  D <- log(D, base = log_base)

}

nr_groups <- length(levels(group))

### set group colours if not previously specified
if (is.null(group_colours) & nr_groups >= 1) group_colours <- scales::hue_pal()(nr_groups)


### normalize data:
D <- automatedNormalization(DATA = D, method = normalization, log = FALSE,
                            id_columns = id, output_path = output_path)

### convert data to long format
D_long <- pivot_longer(data = D, cols = 1:ncol(D))
if (use_groups) {
  D_long$group <- factor(strsplit2(D_long$name, "_")[,1])
} else {
  D_long$group <- NA
}



################################################################################
#### Valid values plots

ValidValuePlot(X = D_long, groupvar_name = groupvar_name, sample_filter = sample_filter,
               plot_device = plot_device, group_colours = group_colours,
               plot_height = plot_height_validvalueplot, plot_width = plot_width_validvalueplot,
               plot_dpi = plot_dpi, suffix = suffix, output_path = output_path)


################################################################################
#### Boxplots and Violin Plots

Boxplots(X = D_long, groupvar_name = groupvar_name, sample_filter = sample_filter,
         plot_device = plot_device, group_colours = group_colours,
         plot_height = plot_height_boxplots, plot_width = plot_width_boxplots,
         plot_dpi = plot_dpi, log_data = FALSE, log_base = log_base,
         suffix = suffix, method = "boxplot")

Boxplots(X = D_long, groupvar_name = groupvar_name, sample_filter = sample_filter,
         plot_device = plot_device, group_colours = group_colours,
         plot_height = plot_height_boxplots, plot_width = plot_width_boxplots,
         plot_dpi = plot_dpi, log_data = FALSE, log_base = log_base,
         suffix = suffix, method = "violinplot")


################################################################################
#### PCA plots

PCA_Plot(X = D, id = id, groupvar1 = group, groupvar2 = NULL, groupvar1_name = groupvar_name, groupvar2_name = NULL,
         plot_device = plot_device, group_colours = group_colours,
         plot_height = plot_height_PCA, plot_width = plot_width_PCA,
         log_data = FALSE, log_base = log_base,  scale. = TRUE,
         impute = FALSE, impute_method = "mean", propNA = 0,
         point.size = 4, base_size = 20,
         returnPCA = FALSE, title = NULL,
         output_path = output_path, suffix = suffix,
         ylim = plot_ylim_PCA, xlim = plot_xlim_PCA, PCx = 1, PCy = 2)

PCA_Plot(X = D, id = id, groupvar1 = group, groupvar2 = NULL, groupvar1_name = groupvar_name, groupvar2_name = NULL,
         plot_device = plot_device, group_colours = group_colours,
         plot_height = plot_height_PCA, plot_width = plot_width_PCA,
         log_data = FALSE, log_base = log_base,  scale. = TRUE,
         impute = FALSE, impute_method = "mean", propNA = 0,
         point.size = 4, base_size = 20,
         returnPCA = FALSE, title = NULL,
         output_path = output_path, suffix = paste0(suffix, "_labelled"),
         ylim = plot_ylim_PCA, xlim = plot_xlim_PCA, label = TRUE, PCx = 1, PCy = 2,
         label_size = 4)


PCA_Plot(X = D, id = id, groupvar1 = group, groupvar2 = NULL,
         groupvar1_name = groupvar_name, groupvar2_name = NULL,
         plot_device = plot_device, group_colours = group_colours,
         plot_height = plot_height_PCA, plot_width = plot_width_PCA,
         log_data = FALSE, log_base = log_base,  scale. = TRUE,
         impute = TRUE, impute_method = "mean", propNA = 0.5,
         point.size = 4, base_size = 20,
         returnPCA = FALSE, title = NULL,
         output_path = output_path, suffix = paste0(suffix, "_imputed_labelled"),
         ylim = plot_ylim_PCA, xlim = plot_xlim_PCA, label = TRUE, PCx = 1, PCy = 2,
         label_size = 4)

################################################################################
### MA-Plots

MAPlots(X = D, log = FALSE, alpha = FALSE,
        plot_height = plot_height_MAPlot, plot_width = plot_width_MAPlot,
        output_path = output_path, suffix = suffix)

################################################################################
################################################################################
################################################################################
#### End of QC Script ##########################################################
