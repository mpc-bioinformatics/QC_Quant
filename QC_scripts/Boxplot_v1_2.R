
## X: Data in long format

Boxplots <- function(X, groupvar_name = "Group", sample_filter = NULL,
                     plot_device = "pdf", suffix = "nonorm",
                     plot_height = 10, plot_width = 15, plot_dpi = 300,
                     log_data = FALSE, log_base = 2, group_colours = NULL,
                     method = "boxplot", base_size = 20, ...){

  require(ggplot2)
  require(tidyverse)

  use_groups <- !all(is.na(X$group))

  if(log_data) {
    X$value <- log(X$value, base = log_base)
  }

  ## filtering of samples
  if (!is.null(sample_filter)) {
  X <- X %>% filter(name %in% sample_filter)
  }

  if (use_groups) {
    pl_boxplot <- ggplot(X, aes(x = name, y = value, fill = group)) +
      labs(fill = groupvar_name)
    if (!is.null(group_colours)) pl_boxplot <- pl_boxplot + scale_fill_manual(values = group_colours)
  } else {
    pl_boxplot <- ggplot(X, aes(x = name, y = value))
  }


  pl_boxplot <- pl_boxplot +
    theme_bw(base_size = base_size) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), ...) +
    ylab("Log intensity") + xlab("Sample")

  if (method != "boxplot"){
    pl_boxplot <- pl_boxplot + geom_violin()
  }else{
    pl_boxplot <- pl_boxplot + geom_boxplot()
  }


  ggsave(paste0(output_path, method, "_", suffix, ".", plot_device), device = plot_device, height = plot_height,
         width = plot_width, plot = pl_boxplot, dpi = plot_dpi)

  return(pl_boxplot)
}


################################################################################
### Changelog:

### Version 1.0 (2021-07-06):
# moved code for boxplot into separate function

### Version 1.1 (2022-03-15):
# add option for a violin plot
# make script work if no groups are defined

### Version 1.2:
# add option to change base size

