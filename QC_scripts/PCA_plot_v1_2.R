
# X: data set (may contain missing values) in wide format
# log_data: logical(1) should data be log-transformed?
# impute: logical(1) should missing values be imputed?
# impute_method: imputation method ("mean" or "median")
# propNA: proportion of missing values that are allowed for a protein
#       - proteins with more NAs will be discarded
#       - for all others, missing values will be imputed by mean or median of the respective protein
# scale.: should data be scaled before computing PCA? (argument of prcomp())
# groupvar1: variable used for colouring
# groupvar2: variable used for shape
# groupvar1_name, groupvar2_name: Titles of legends for colour and shape
# point.size size of points in the plot
# base_size: base size
# group_colours: colours for the groups
# return_PCA: if TRUE, function returns PCA values (coordinates)
# title: title for the plot
# plot_device: default is pdf
# plot_height: in cm
# plot_width: in cm
# plot_dpi
# output_path
# normalization: determines suffix for file name
# xlim, ylim: optional x and y axis limits


PCA_Plot <- function(X, id, log_data = TRUE, log_base = 2,
                     impute, impute_method = "mean", propNA = 0,
                     scale. = TRUE,
                     groupvar1, groupvar2 = NULL, groupvar1_name, groupvar2_name = NULL,
                     point.size = 4, base_size = 11,
                     group_colours = NULL, returnPCA = FALSE, title = NULL,
                     plot_device = "pdf", plot_height = 10, plot_width = 10, plot_dpi = 300,
                     output_path = "", suffix = "nonorm", xlim = NULL, ylim = NULL,
                     label = FALSE, PCx = 1, PCy = 2) {

  require(ggplot2)
  require(matrixStats)
  require(cowplot)
  require(ggrepel)

  use_groups <- !is.null(groupvar1)

  if (log_data) {
    X <- log(X, base = log_base)
  }


  mean_NA <- apply(X, 1, function(x) mean(is.na(x)))

  ### remove rows with too many missing values
  X_2 <- X[mean_NA <= propNA, ]

  ## perform imputation
  if (impute) {
    X_3 <- as.data.frame(t(apply(X_2, 1, function(x) {
      if(anyNA(x)) {
        x[is.na(x)] <- switch(impute_method, mean = mean(x, na.rm = TRUE),
                              median = median(x, na.rm = TRUE))
      }
      return(x)
    })))
  } else {
    X_3 <- na.omit(X_2)
  }

  X_3 <<- X_3

  ### remove proteins/peptides with a constant value (variance near zero)
  v <- matrixStats::rowVars(as.matrix(X_3))
  ind_zeroVar <- which(v < 1e-25)
  if (length(ind_zeroVar) > 0) X_3 <- X_3[-ind_zeroVar,]

  print(paste0(nrow(X_3), " of ", nrow(X), " rows are used for PCA."))

  ### calculate PCA
  pca <- prcomp(t(X_3), scale. = scale.)
  pred <- predict(pca, t(X_3))
  summ <- summary(pca)

  var50 <- which(summ$importance[3,] >= 0.5)[1]
  print(paste0("50% explained variance is reached with ", var50, " principle components."))

  ### version with colour and shape
  if (!is.null(groupvar1) & !is.null(groupvar2)) {
    D_PCA <- data.frame(pred[,c(PCx,PCy)], groupvar1 = groupvar1, groupvar2 = groupvar2)
    colnames(D_PCA)[1:2] <- c("PCx", "PCy")
    pl <- ggplot(data = D_PCA, aes(x=PCx, y=PCy)) +
      geom_point(aes(colour = groupvar1, shape = groupvar2), size = point.size)
    pl <- pl + labs(colour = groupvar1_name, shape = groupvar2_name)
    if (!is.null(group_colours)) pl <- pl + scale_colour_manual(values = group_colours)
    ### more than 6 different shapes will otherwise give an error message:
    if (nlevels(D_PCA$groupvar2) > 6) pl <- pl + scale_shape_manual(values = 1:nlevels(D$groupvar2))
  }

  ### version with only colour
  if (!is.null(groupvar1) & is.null(groupvar2)) {
    D_PCA <- data.frame(pred[,c(PCx,PCy)], groupvar1 = groupvar1, label = colnames(X_3))
    colnames(D_PCA)[1:2] <- c("PCx", "PCy")
    pl <- ggplot(data = D_PCA, aes(x=PCx, y=PCy)) +
      geom_point(aes(colour = groupvar1), size = point.size)
    pl <- pl + labs(colour = groupvar1_name)
    if (!is.null(group_colours)) pl <- pl + scale_colour_manual(values = group_colours)
    if(label) pl <- pl + geom_text_repel(aes(x=PCx, y=PCy, label = label, colour = groupvar1))
  }


  ### version without colour or shape
  if (is.null(groupvar1) & is.null(groupvar2)) {
    D_PCA <- data.frame(pred[,c(PCx,PCy)], label = colnames(X_3))
    colnames(D_PCA)[1:2] <- c("PCx", "PCy")
    pl <- ggplot(data = D_PCA, aes(x=PCx, y=PCy)) +
      geom_point(size = point.size)
    if(label) pl <- pl + geom_text_repel(aes(x=PCx, y=PCy, label = label))
  }

  pl <- pl + theme_bw(base_size = base_size) +
    xlab(paste0("PC", PCx, " (", round(100*summ$importance[2,PCx], 1), "%)")) +
    ylab(paste0("PC", PCy, " (", round(100*summ$importance[2,PCy], 1), "%)"))

  ### add plot title
  if (!is.null(title)) {
    pl <- pl + ggtitle(title)
  }

  ### add xlim and ylim if specified
  if (!is.null(xlim)) pl <- pl + xlim(xlim)
  if (!is.null(ylim)) pl <- pl + xlim(ylim)

  #print(pl)

  ggsave(paste0(output_path, "PCA_plot_", suffix, ".", plot_device),
         device = plot_device, plot = pl, dpi = plot_dpi, height = plot_height,
         width = plot_width, units = "cm")

  if (returnPCA) {
    return(list(pl = pl, D_PCA_plot = cbind(D_PCA, Sample = colnames(X)), pca = pca))#,
    # D_PCA = cbind(id3, X_3)))
  } else {
    return(pl)
  }

}


################################################################################
### Changelog:

### Version 1.0 (2021-07-06):
# minor changes, e.g. saving of the PCA plot

### Version 1.1 (2021-10-18):
# add option to set x and y axis limits

### Version 1.2
# add label option
