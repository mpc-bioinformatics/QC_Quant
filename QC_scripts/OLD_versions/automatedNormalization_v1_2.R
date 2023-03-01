
automatedNormalization <- function(DATA, DATA.name = deparse(substitute(DATA)),
                                   method = "loess", suffix = method, log = TRUE, id_columns = NULL,
                                   output_path = "", groupwise = FALSE, group = NULL) {

  require(limma)
  require(beepr)
  require(openxlsx)

  if(method == "loess" | method == "quantile" | method == "median"){

    if(log) {
      DATA <- log2(DATA)
    }

    #### choose normalization function
    fun <- switch(method,
                  "loess" = limma::normalizeBetweenArrays,
                  "quantile" = limma::normalizeBetweenArrays,
                  "median" = limma::normalizeBetweenArrays)

    ### choose arguments for normalization function
    args <- switch(method,
                   "loess" = list(object = DATA, method = "cyclicloess"),
                   "quantile" = list(object = DATA, method = "quantile"),
                   "median" = list(object = DATA, method = "scale"))

    if (!groupwise) {
      DATA_norm <- do.call(fun, args)
      DATA_norm <- as.data.frame(DATA_norm)
    } else {
      DATA_split <- split.default(DATA, group)
      DATA_split_norm <- lapply(DATA_split, limma::normalizeBetweenArrays, method = args$method)
      DATA_norm <- do.call(cbind, DATA_split_norm)
    }



    if (length(id_columns) >= 1 ) {
      DATA_norm_2 <- cbind(id_columns, DATA_norm)
    }

    tryCatch(expr = {
      DATA_norm_2 <- as.data.frame(DATA_norm_2)
      write.xlsx(x = DATA_norm_2, file = paste0(output_path, DATA.name, "_", suffix, ".xlsx"), keepNA = TRUE)
      message("Normalized data successfully saved!")},
      error = function(err) {
        # error handler picks up where error was generated
        print(paste("MY_ERROR:  ",err))
        beepr::beep(sound = 10)
        user_input <- readline(prompt = paste0("+++ Do you want to overwrite ", paste0(DATA.name,"_",method,".xlsx"), "? +++ [yes/no] "))
        if(user_input == "yes"){
          DATA_norm_2 <- as.data.frame(DATA_norm_2)
          write.xlsx(x = DATA_norm_2, file = paste0(output_path, DATA.name,"_",suffix,".xlsx"), overwrite = TRUE, keepNA = TRUE)
          message("Normalized data successfully saved!")
        } else {
          message("Overwriting of normalized data failed. Please allow overwriting, remove the data file or choose different normalization method!")
        }
      })
  }else{ # if method == "nonorm"
    DATA_norm <- DATA
    cat("No normalization applied.")
    if (length(id_columns) >= 1 ) {
      DATA_norm_2 <- cbind(id_columns, DATA_norm)
    }
    DATA_norm_2 <- as.data.frame(DATA_norm_2)

    write.xlsx(x = DATA_norm_2, file = paste0(output_path, DATA.name, "_", suffix, ".xlsx"), keepNA = TRUE)
  }

  return(DATA_norm)
}


################################################################################
### Changelog:

### Version 1.0 (2021-07-06):
# added function for normalization of data
# added handling for overwriting existing normalized data

### Version 1.1 (2021-10-18):
# No changes.

### Version 1.2:
# allow group-wise normalization (in general not recommended, only use when you know what you are doing!)

