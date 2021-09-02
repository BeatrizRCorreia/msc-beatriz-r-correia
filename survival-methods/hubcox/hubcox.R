# INPUT PARAMETERS ---------------------------------------------------
# Set seed to "1997" for reproducible results
params <- list(seed = 1997)
set.seed(params$seed)

# Number of successful iterations that we want
nr.of.iterations.to.perform <- 25

# Loss measure used for cross-validation
measures.list <- c("deviance")
# measures.list <- c("C")

# Alpha for the elastic net mixing parameter (alpha=0 (ridge), alpha=1 (lasso))
alpha.values.list <- c(0.1)

# Min.degree
min.degree.v <- 0.6

# Pdf for the plots
pdf('~/Documents/msc-beatriz-r-correia/survival-methods/hubcox/results/alpha0.1_test.pdf')
# -------------------------------------------------------------------

library(dplyr)
library(ggplot2)
library(survival)
library(loose.rock)
library(futile.logger)
library(curatedTCGAData)
library(TCGAutils)
library(R.utils)
library(glmSparseNet)
library(caret)

# Import functions from another file
setwd("~/Documents/msc-beatriz-r-correia/other-functions/")
source("auxiliary_functions.R")

# Download a dataset from the package curatedTCGAData by the string code passed as argument
# source("extract-dataset.R")
# extract_data("BRCA")
# extract_primary_solid_tumor_data("BRCA")

# To run on my pc
setwd("~/Documents/msc-beatriz-r-correia/data/BRCA_primary_solid_tumor/")

# To generate the transposed file available in the data folder 
# data <- read.csv("data2_01_BRCA_RNASeq2GeneNorm-20160128.csv", check.names=FALSE)
# transposeddata = setNames(data.frame(t(data[,-1])), data[,1])
# write.csv(transposeddata, file="data2_01_BRCA_RNASeq2GeneNorm-20160128_transposed.csv")

transposeddata <- read.csv("data2_01_BRCA_RNASeq2GeneNorm-20160128_transposed.csv", check.names=FALSE)
coldata <- read.csv("data2_colData.csv", check.names=FALSE)

coldata <- coldata %>% as.data.frame %>%
  select(patientID, vital_status, Days.to.date.of.Death, Days.to.Date.of.Last.Contact, days_to_death, days_to_last_followup) %>%
  # Convert days to integer
  mutate(Days.to.date.of.Death = as.integer(Days.to.date.of.Death)) %>%
  mutate(Days.to.Date.of.Last.Contact  = as.integer(Days.to.Date.of.Last.Contact)) %>%
  mutate(days_to_death  = as.integer(days_to_death)) %>%
  mutate(days_to_last_followup  = as.integer(days_to_last_followup)) %>%
  # Find max time between all days (ignoring missings)
  rowwise %>%
  mutate(time = max(days_to_last_followup, Days.to.date.of.Death, Days.to.Date.of.Last.Contact, days_to_death, na.rm = TRUE)) %>%
  # Keep only survival variables and codes
  select(patientID, status = vital_status, time) %>% 
  # Discard individuals with survival time less or equal to 0
  filter(!is.na(time) & time > 0) %>% as.data.frame
print("---------------------------------------------------------------------------------------")

rownames(transposeddata) <- transposeddata[,1]
colnames(transposeddata)[1]<-"patient_colname"
transposeddata$patient_colname <- NULL

# Set index as the patientID
rownames(coldata) <- coldata$patientID

# Get matches between survival and assay data
transposeddata <- transposeddata[TCGAbarcode(rownames(transposeddata)) %in% rownames(coldata),]

censored = subset(coldata, coldata$status == 0)
not_censored = subset(coldata, coldata$status == 1)
# Splits the censored data randomly into 80% and 20%
censored80 = sort(sample(nrow(censored), nrow(censored)*.8))
ytrain_censored <- censored[censored80,]
ytest_censored <- censored[-censored80,]
# Splits the not_censored data randomly into 80% and 20%
not_censored80 = sort(sample(nrow(not_censored), nrow(not_censored)*.2))
ytest_notcensored <- not_censored[not_censored80,]
ytrain_notcensored <- not_censored[-not_censored80,]
# Merge the censored and not censored parts
ytrain = rbind(ytrain_censored, ytrain_notcensored)
ytest = rbind(ytest_censored, ytest_notcensored)

transposeddata000 <- transposeddata
new_column <- TCGAbarcode(rownames(transposeddata000))
rownames(transposeddata000) <- new_column
xtrain <- transposeddata000[rownames(ytrain),]

# Scaling xtrain
preprocess_train <- preProcess(xtrain, method = c("center", "scale"))
xtrain <- predict(preprocess_train, xtrain)
variables_with_variance_zero_or_close <- nearZeroVar(xtrain)
xtrain <- xtrain[,-variables_with_variance_zero_or_close]
xtrain <- as.matrix(xtrain)

xtest <- transposeddata000[rownames(ytest),]

# Scaling xtest
xtest <- predict(preprocess_train, xtest)
xtest <- xtest[,-variables_with_variance_zero_or_close]
xtest <- as.matrix(xtest)

# Using only a subset of genes previously selected to keep this short example.
# set.seed(params$seed)
small.subset <- c('CD5', 'CSF2RB', 'IRGC', 'NEUROG2', 'NLRC4', 'PDE11A',
                  'PTEN', 'TP53', 'BRAF',
                  'PIK3CB', 'QARS', 'RFC3', 'RPGRIP1L', 'SDC1', 'TMEM31',
                  'YME1L1', 'ZBTB11', sample(colnames(transposeddata), 100)) %>%
  unique

# Use only a few selected genes from the dataset
# xdata <- transposeddata[, small.subset[small.subset %in% colnames(transposeddata)]]

# Use all genes in the dataset divided in train and test sets
# xtrain
# xtest
ytrain <- ytrain %>% select(time, status)
ytest <- ytest %>% select(time, status)

run.counter <- 1
all.outputs.info <- list()
start.counter <- 1
nr.of.null.runs <- 0
x <- 1
start_time <- Sys.time()

while (x <= nr.of.iterations.to.perform) {
  for (k in 1:length(measures.list)) {
    list.of.fits <- list()
    for (i in 1:length(alpha.values.list)) {
      tryCatch({
        print(run.counter)
        start_time2 <- Sys.time()
        output.info <- vector("list", length = 9)
        names(output.info) <- c("min.degree", "alpha", "nr.selected.genes", "list.of.genes", "pvalue.train", "pvalue.test", "iteration.nr", "cindex.train", "cindex.test")

        fit.name <- paste0("BRCA_", measures.list[k])
        fit.name <- paste0(fit.name, "_alpha")
        fit.name <- paste0(fit.name, alpha.values.list[i])
        pdf.name <- paste(fit.name, "pdf", sep=".")
        txt.name <- paste(fit.name, "txt", sep=".")
        
        output.info[1] <- min.degree.v
        output.info[2] <- alpha.values.list[i]
        list.of.fits[[fit.name]] <- cv.glmHub(xtrain, Surv(ytrain$time, ytrain$status),
                                              family  = 'cox',
                                              nlambda = 100,
                                              nfolds = 10,
                                              alpha = alpha.values.list[i],
                                              network = "correlation",
                                              network.options = networkOptions(cutoff = 0.005, min.degree = min.degree.v))
        coefs.v <- coef(list.of.fits[[fit.name]], s = 'lambda.min')[,1] %>% { .[. != 0]}
        output.info[3] <- length(names(coefs.v))
        if (output.info$nr.selected.genes == 0) {
          x <- x
          nr.of.null.runs <- nr.of.null.runs + 1
        }
        else {
          x <- x + 1
          print("mais uma run successful :)")
        }
        output.info[4] <- list(names(coefs.v))
        if (output.info$nr.selected.genes != 0) {
          output_of_separated1 <- separate2GroupsCox(as.vector(coefs.v),
                                                     xtrain[, names(coefs.v)],
                                                     ytrain,
                                                     xlab = "Time in days",
                                                     plot.title = 'Train set', legend.outside = FALSE,
                                                     risk.table = TRUE,
                                                     tables.height = 0.2,
                                                     surv.median.line = "hv", # adds the median survival pointer
                                                     cumcensor = TRUE)
          output.info[5] <- output_of_separated1$pvalue
          print(output_of_separated1$plot)
        }
        else {
          output.info[5] <- NA
        }
        if (output.info$nr.selected.genes != 0) {
          output_of_separated2 <- separate2GroupsCox(as.vector(coefs.v),
                                                     xtest[, names(coefs.v)],
                                                     ytest,
                                                     xlab = "Time in days",
                                                     plot.title = 'Test set', legend.outside = FALSE,
                                                     risk.table = TRUE,
                                                     tables.height = 0.2,
                                                     surv.median.line = "hv", # adds the median survival pointer
                                                     cumcensor = TRUE)
          output.info[6] <- output_of_separated2$pvalue
          print(output_of_separated2$plot)
        }
        else {
          output.info[6] <- NA
        }
        output.info[7] <- run.counter
        
        if (output.info$nr.selected.genes > 0) {
          c.index.train <- c.index.calculator(coefs.v, xtrain, ytrain)
          output.info[8] <- c.index.train
          c.index.test <- c.index.calculator(coefs.v, xtest, ytest)
          output.info[9] <- c.index.test
        }
        else {
          output.info[8] <- NA
          output.info[9] <- NA
        }
        all.outputs.info[start.counter] <- list(output.info)
        start.counter <- start.counter + 1
        run.counter <- run.counter + 1
        end_time2 <- Sys.time()
        total_time2 <- end_time2 - start_time2
        print(total_time2)
      }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    }
  }
}

print("nr of iterations to perform:")
print(nr.of.iterations.to.perform)

print("nr of runs executed:")
print(length(all.outputs.info))

print("nr of unsuccessful runs:")
nr.runs.executed <- length(all.outputs.info) - nr.of.iterations.to.perform
print(nr.runs.executed)

genes.selected <- get_genes_selected(nr.of.iterations.to.perform, all.outputs.info)

print("genes selected with the nr of iterations they were selected:")
print(genes.selected)

print("---------------------------------------------------------")

print("average nr of genes selected in successful runs:")
average.genes.selected <- average_number_of_genes_selected(nr.of.iterations.to.perform, all.outputs.info)
print(average.genes.selected)

print("standard deviation of the average nr of genes selected in successful runs:")
standard.deviation.average.genes <- calculate_standard_deviation("genes.selected", average.genes.selected, nr.of.iterations.to.perform, all.outputs.info)
print(standard.deviation.average.genes)

print("----- PVALUE TRAIN --------------------------------------")

print("average pvalue train in successful runs:")
average.pvalue.train <- average_pvalue_train(nr.of.iterations.to.perform, all.outputs.info)
print(average.pvalue.train)

print("standard deviation of the average pvalue train:")
standard.deviation.average.pvalue.train <- calculate_standard_deviation("pvalue.train", average.pvalue.train, nr.of.iterations.to.perform, all.outputs.info)
print(standard.deviation.average.pvalue.train)

print("---------------------------------------------------------")

print("average pvalue train in statistically significant successful runs:")
average.pvalue.train.stsg.runs <- average_pvalue_train_significant(all.outputs.info, "average")
print(average.pvalue.train.stsg.runs)

print("nr of statistically significant pvalue train successful runs:")
nr.of.stsg.runs.train <- average_pvalue_train_significant(all.outputs.info, "counter")
print(nr.of.stsg.runs.train)

print("---------------------------------------------------------")

print("average pvalue train in NON-statistically significant successful runs:")
average.pvalue.train.non.stsg.runs <- average_pvalue_train_non_significant(all.outputs.info, "average")
print(average.pvalue.train.non.stsg.runs)

print("nr of NON-statistically significant pvalue train successful runs:")
nr.of.non.stsg.runs.train <- average_pvalue_train_non_significant(all.outputs.info, "counter")
print(nr.of.non.stsg.runs.train)

print("----- PVALUE TEST --------------------------------------")

print("average pvalue test in successful runs:")
average.pvalue.test <- average_pvalue_test(nr.of.iterations.to.perform, all.outputs.info)
print(average.pvalue.test)

print("standard deviation of the average pvalue test:")
standard.deviation.average.pvalue.test <- calculate_standard_deviation("pvalue.test", average.pvalue.test, nr.of.iterations.to.perform, all.outputs.info)
print(standard.deviation.average.pvalue.test)

print("---------------------------------------------------------")

print("average pvalue test in statistically significant successful runs:")
average.pvalue.test.stsg.runs <- average_pvalue_test_significant(all.outputs.info, "average")
print(average.pvalue.test.stsg.runs)

print("nr of statistically significant pvalue test successful runs:")
nr.of.stsg.runs.test <- average_pvalue_test_significant(all.outputs.info, "counter")
print(nr.of.stsg.runs.test)

print("---------------------------------------------------------")

print("average pvalue test in NON-statistically significant successful runs:")
average.pvalue.test.non.stsg.runs <- average_pvalue_test_non_significant(all.outputs.info, "average")
print(average.pvalue.test.non.stsg.runs)

print("nr of NON-statistically significant pvalue test successful runs:")
nr.of.non.stsg.runs.test <- average_pvalue_test_non_significant(all.outputs.info, "counter")
print(nr.of.non.stsg.runs.test)

print("---------------------------------------------------------")

print("iteration with the lowest pvalue train:")
iteration.with.lowest.pvalue.train <- run_with_lowest_pvalue_train(all.outputs.info)
print(iteration.with.lowest.pvalue.train)

print("---------------------------------------------------------")

print("iteration with the lowest pvalue test:")
iteration.with.lowest.pvalue.test <- run_with_lowest_pvalue_test(all.outputs.info)
print(iteration.with.lowest.pvalue.test)

print("----- C-INDEX TRAIN --------------------------------------")

print("average c-index train in successful runs:")
average.c.index.train <- c.index.average.train(nr.of.iterations.to.perform, all.outputs.info, "average")
print(average.c.index.train)

print("iterations accounted for the average calculation (not NA):")
not.NA.iterations.c.index.train <- c.index.average.train(nr.of.iterations.to.perform, all.outputs.info, "iterations")
print(not.NA.iterations.c.index.train)

print("standard deviation of the average c-index train:")
standard.deviation.average.c.index.train <- calculate_standard_deviation("cindex.train", average.c.index.train, not.NA.iterations.c.index.train, all.outputs.info)
print(standard.deviation.average.c.index.train)

print("----- C-INDEX TEST --------------------------------------")

print("average c-index test in successful runs:")
average.c.index.test <- c.index.average.test(nr.of.iterations.to.perform, all.outputs.info, "average")
print(average.c.index.test)

print("iterations accounted for the average calculation (not NA):")
not.NA.iterations.c.index.test <- c.index.average.test(nr.of.iterations.to.perform, all.outputs.info, "iterations")
print(not.NA.iterations.c.index.test)

print("standard deviation of the average c-index test:")
standard.deviation.average.c.index.test <- calculate_standard_deviation("cindex.test", average.c.index.test, not.NA.iterations.c.index.test, all.outputs.info)
print(standard.deviation.average.c.index.test)

print("---------------------------------------------------------")

end_time <- Sys.time()

print(end_time - start_time)

dev.off()
