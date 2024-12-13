### ROC analysis script ####
### Example data to run the code can be found in example_files. Full analysis can be run on Dataset S1 (to be published).
### ROC analysis on following metrics with Category as the factor variable:
# `average AFM pDockQ`
# `average AFM model confidence`
# `average AFM ipTM`
# `average AFM pDockq2`
# `average AFM LIS`
# `average AF3 LIS`
# `average AF3 ipTM`
# `average AF3 model confidence`
# `best AFM pDockQ`
# `best AFM model confidence`
# `best AFM ipTM`
# `best AFM pDockq2`
# `best AFM LIS`
# `best AF3 LIS`
# `best AF3 ipTM`
# `best AF3 model confidence`
# `Mean AFM score`
# `Combined AFM & AF3 score`
### therefore pick columns 6:22

library(data.table)
library(pROC)


dt <- fread("path/to/example_data_3.csv")

## Pick relevant columns
ROC_analysis_data <- dt[, c(7:11,13:17,19:21,23:28)]
names(ROC_analysis_data)



## Set the levels of the factor variable explicity.
## "AtPRS-v1" = interactor
## "AtRRS-v1" = non-interactor
category_levels <- c("AtRRS-v1", "AtPRS-v1")
ROC_analysis_data$Category <- factor(ROC_analysis_data$Category, levels = category_levels)

## ROC analysis of different metrics
ROCs <- list()
for(icol in 1:18){
  ROCs[[icol]] <- roc(ROC_analysis_data$Category, ROC_analysis_data[[icol]], quiet = TRUE)
}
names(ROCs) <- colnames(ROC_analysis_data)[1:18]



nROCs <- length(ROCs)

first <- TRUE

for(icol in 1:nROCs){
  plot(ROCs[[icol]], add = !first,
       col = icol, lty = icol, 
       xaxs = 'i', yaxs = 'i') 
  first <- FALSE
}

ROC_labels <- names(ROCs)

## get the AUC values
AUC_values <- data.frame(method = ROC_labels,
                         AUC = sapply(ROCs, function(x) as.numeric(auc(x))))


## sort by the value of the AUC
AUC_values <- sort_by(AUC_values, ~AUC, decreasing = TRUE)
AUC_values$method <- paste0(AUC_values$method," (AUC = ",round(AUC_values$AUC,4),")")

## Legend 
legend('bottomright',
       AUC_values$method,
       col = 1:nROCs, lty = 1:nROCs, lwd = 2, ncol = 1, bg = 'white',cex=0.8)

## Print the values
for(i in 1:nrow(AUC_values)){
  cat(AUC_values$method[i],"\n")
