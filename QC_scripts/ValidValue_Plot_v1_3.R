
## X: Data in long format

ValidValuePlot <- function(X, groupvar_name = "Group", sample_filter=NULL,
                           suffix = "nonorm", plot_device = "pdf",
                           group_colours = NULL, plot_height = 10,
                           plot_width = 15, plot_dpi = 300, ylim = NULL, title = NULL,
                           output_path = "",  ...){

  require(ggplot2)
  require(tidyverse)

  use_groups <- !all(is.na(X$group))


  # filter data set for samples defined in sample_filter
  if (!is.null(sample_filter)) {
    X <- X %>% filter(name %in% sample_filter)
  }

  X <- X %>% group_by(name, group) %>% summarize(nrvalid = sum(!is.na(value)), meanvalid = mean(!is.na(value)), .groups = 'drop')

  write.xlsx(x = X, file = paste0(output_path, "validvalues_", suffix, ".xlsx"), overwrite = TRUE, keepNA = TRUE)


  pl_valid_values <- ggplot(X) +
    theme_bw(base_size = 15) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
          legend.position = "none", plot.title = element_text(hjust = 0.5)) +
    ylab("Number of valid values") + xlab("Sample")


  if (use_groups) {
    pl_valid_values <- pl_valid_values +
      geom_bar(stat = "identity",aes(x=name,y= nrvalid, fill=group)) +
      labs(fill = groupvar_name)
    if (!is.null(group_colours)) pl_valid_values <- pl_valid_values + scale_fill_manual(values = group_colours)
  } else {
    pl_valid_values <- pl_valid_values +
      geom_bar(stat = "identity",aes(x=name,y= nrvalid))
  }

  if (!is.null(ylim)) pl_valid_values <- pl_valid_values + ylim(ylim)
  if (!is.null(title)) pl_valid_values <- pl_valid_values + ggtitle(title)

  ggsave(paste0(output_path,"valid_value_plot_", suffix,".",plot_device),
         plot = pl_valid_values, device = plot_device,
         height = plot_height, width = plot_width, dpi = plot_dpi)


  if (use_groups) {
    pl_valid_values_rel <-  ggplot(X, aes(x = name, y = meanvalid*100, fill = group)) +
      labs(fill = groupvar_name)
    if (!is.null(group_colours)) pl_valid_values_rel <- pl_valid_values_rel + scale_fill_manual(values = group_colours)
  } else {
    pl_valid_values_rel <-  ggplot(X, aes(x = name, y = meanvalid*100))
  }
  pl_valid_values_rel <- pl_valid_values_rel +
    geom_bar(stat = "identity") +
    theme_bw(base_size = 15) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          legend.position = "none", plot.title = element_text(hjust = 0.5)) +
    ylab("Percentage of valid values") + xlab("Sample")

  if (!is.null(title)) pl_valid_values <- pl_valid_values + ggtitle(title)

  ggsave(paste0(output_path, "valid_values_percentage_plot_",suffix,".",plot_device),
         device = plot_device,height = plot_height,width = plot_width,
         plot = pl_valid_values_rel, dpi = plot_dpi)

  return(list(pl_valid_values = pl_valid_values, pl_valid_values_rel = pl_valid_values_rel))
}


################################################################################
### Changelog:

### Version 1.0 (2021-07-06):
# moved code for valid value plot to separate function

### Version 1.1 (2022-03-15):
# add white background (theme_bw())
# make script work if no groups are defined

