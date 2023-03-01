

#### Function for a single MA-Plot:
# x1: Sample 1
# x2: Sample 2
# log: Should data be log-transformed?
#      TRUE, if not already log-transformed, FALSE, if already log-transformed
# alpha: Should points be transparent?
# col: colours of the data points
# ...: further arguments for ma.plot
MAPlot_single <- function(x1, x2, log = TRUE, alpha = FALSE, col = "black", ...) {

  if(log) {
    x1 <- log2(x1)
    x2 <- log2(x2)
  }
  if(alpha) {
    col = alpha(col, 0.5)
  }

  M <- na.omit(x1 - x2)
  A <- na.omit((x1 + x2)/2)

  if (length(col) > 1) {
    na.ind <- attr(M, "na.action")
    col <- col[-na.ind]
  }


  ma.plot(A = A, M = M, pch = 16, cex = 0.7, col = col, show.statistics = FALSE, ...)
}



# function to check if user is sure to plot more than 1000 plots --> if yes it return 1 and the plots will be created
MAPlots_check <- function(X, maxPlots, ...){
  number_states <- max(as.integer(as.factor(colnames(X))))
  number_plots <- choose(number_states,2)
  return_value <- 2

  if(number_plots >= maxPlots){
    beepr::beep(sound = 10)
    user_input <- readline(prompt = paste("Are you sure to plot",number_plots,"MA_plots? [yes/no]"))
    if (user_input == "yes"){
      return_value <- 1
    }else return_value <- 0
  }else return_value <- 1

  return(return_value)
}



### main function for MA-Plots
# X: Data in wide format
# labels: labels of the samples for the title of the MA-Plot
# labels2: second line in title, e.g. group membership
MAPlots <- function(X, log = TRUE, alpha = FALSE, suffix="nonorm",
                    labels = 1:ncol(X), labels2 = colnames(X), maxPlots = 5000,
                    plot_height=15, plot_width=15, output_path = "", ...) {

  require(limma)
  require(affy)
  require(scales)
  require(beepr)

  number_states <- max(as.integer(as.factor(colnames(X))))
  number_plots <- choose(number_states,2)

  if(MAPlots_check(X, maxPlots) == 1){

  num <- 0

  print("Generating MA-Plots ...")

  pb <- txtProgressBar(min = 0,max = number_plots,char = "#",style = 3)

  pdf(paste0(output_path, "MA_Plots_", suffix, ".pdf"), height = plot_height/2.54, width = plot_width/2.54)

  for(i in 1:(ncol(X)-1)) {
    for (j in (i + 1):ncol(X)) {

      if (is.null(labels2)) {
        main = paste(labels[i], labels[j])
      } else  {
        main = paste(labels[i], labels[j], "\n", labels2[i], labels2[j])
      }

      num <- num + 1
      setTxtProgressBar(pb, num)

      MAPlot_single(X[,i], X[, j], log = log, main = main, ...)
    }
  }
  # sound chosen, "treasure", "facebook" also cool :)
  beepr::beep("coin")
  close(pb)
  print("MA-Plots finished!")

  dev.off()

  }

}

################################################################################
### Changelog:

### Version 1.0 (2021-07-06):
# added progress bar
# added handling of situations where more than 1000 plots will be generated

### Version 1.1 (2021-10-18)
# no changes
## 2022-01-27
# print number_plots "are u sure to plot number_plots?"

