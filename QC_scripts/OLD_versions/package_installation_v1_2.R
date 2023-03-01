#checking for necessary packages, installing missing ones and updating old ones


update_available <- old.packages()[,1]


if("affy" %in% update_available | !require(limma)){
  if(!requireNamespace("BiocManager", quietly = TRUE)){
    install.packages("BiocManager")
  }
  BiocManager::install("affy", update = TRUE, ask = FALSE, force = TRUE)
}

if("beepr" %in% update_available | !require(beepr)) {install.packages("beepr", dependencies = TRUE, quiet = TRUE)}

if("cowplot" %in% update_available | !require(cowplot)) {install.packages("cowplot", dependencies = TRUE, quiet = TRUE)}

if("devtools" %in% update_available | !require(devtools)) {install.packages("devtools", dependencies = TRUE, quiet = TRUE)}

if("ggplot2" %in% update_available | !require(ggplot2)) {install.packages("ggplot2", dependencies = TRUE, quiet = TRUE)}

if("ggplus" %in% update_available | !require(ggplus)) {devtools::install_github("guiastrennec/ggplus", dependencies = TRUE, quiet = TRUE)}

if("ggpubr" %in% update_available | !require(ggpubr)) {install.packages("ggpubr", dependencies = TRUE, quiet = TRUE)}

if("limma" %in% update_available | !require(limma)){
  if(!requireNamespace("BiocManager", quietly = TRUE)){
    install.packages("BiocManager")
  }
  BiocManager::install("limma", update = TRUE, ask = FALSE, force = TRUE)
}

if("matrixStats" %in% update_available | !require(matrixStats)) {install.packages("matrixStats", dependencies = TRUE, quiet = TRUE)}

if("openxlsx" %in% update_available | !require(openxlsx)) {install.packages("openxlsx", dependencies = TRUE, quiet = TRUE)}

if("scales" %in% update_available | !require(scales)) {install.packages("scales", dependencies = TRUE, quiet = TRUE)}

if("tidyverse" %in% update_available | !require(tidyverse)) {install.packages("tidyverse", dependencies = TRUE, quiet = TRUE)}

if("ComplexHeatmap" %in% update_available | !require(limma)){
  if(!requireNamespace("BiocManager", quietly = TRUE)){
    install.packages("BiocManager")
  }
  BiocManager::install("ComplexHeatmap", update = TRUE, ask = FALSE, force = TRUE)
}


################################################################################
### Changelog:

### Version 1.0 (2021-07-06):
# added possibility to install/update all necessary R packages

### Version 1.0 (2021-10-18):
# changed for limma and affy package the option to ask = FALSE
#  (User will not be asked on the console or just once)

### Version 1.2
# add ComplexHeatmap package

