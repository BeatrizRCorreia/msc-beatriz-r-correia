# INPUT PARAMETERS ---------------------------------------------------
# Set seed to "1997" for reproducible results
params <- list(seed = 1997)
set.seed(params$seed)

# Number of successful iterations that we want
nr.of.iterations.to.perform <- 25

# Tumoral weights to use
# direct weights (w)
weights.used <- c("w")
# inverse of the weights (1/w)
# weights.used <- c("1/w")

# Alpha for the elastic net mixing parameter (alpha=0 (ridge), alpha=1 (lasso))
alpha.values.list <- c(0.1)

# Fit to which data
# BRCA (breast)
fit.to.dataset <- "BRCA"
# PRAD (prostate)
# fit.to.dataset <- "PRAD"

# Pdf for the plots
# Server
pdf('~/msc-beatriz-correia/results/BRCA_alpha0.1_w_test.pdf')
# My pc
# pdf('~/Documents/msc-beatriz-r-correia/survival-methods/tcox/results/BRCA_alpha0.1_w_test.pdf')
# -------------------------------------------------------------------

suppressMessages(library(dplyr))
suppressMessages(library(DT))
suppressMessages(library(tidyverse))
suppressMessages(library(futile.logger))
suppressMessages(library(loose.rock))
suppressMessages(library(survival))
suppressMessages(library(GGally))
suppressMessages(library(ggplot2))
suppressMessages(library(glmnet))
suppressMessages(library(glmSparseNet))
suppressMessages(library(gplots))
suppressMessages(library(readr))
suppressMessages(library(maftools))
suppressMessages(library(glmSparseNet))
suppressMessages(library(TCGAutils))
suppressMessages(library(caret))

# Import functions from another file
# Server
source("auxiliary_functions.R")
# My pc
# setwd("~/Documents/msc-beatriz-r-correia/other-functions/")
# source("auxiliary_functions.R")

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# FUNCTIONS
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

twiner_glmnet <- function (x, y, alpha, w) {
  cv.glmnet(x, y, family = 'cox', penalty.factor = w, alpha = alpha, nlambda = 100, nfolds= 10)
}

twiner_coefs <- function (cv.fit){
  coef(cv.fit, s = 'lambda.min')[,1] %>% { .[. != 0]}
}

twiner_genes <- function (cv.fit) {
  b <- which(cv.fit$glmnet.fit$beta[,which(cv.fit$cvm == min(cv.fit$cvm))] != 0)
  return (names(b))
}

start_time111 <- Sys.time()

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# DATASET: BRCA (breast)
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# To run on server
setwd("~/msc-beatriz-correia/data/BRCA_primary_solid_tumor/")
# To run on my pc
# setwd("~/Documents/msc-beatriz-r-correia/data/BRCA_primary_solid_tumor/")

# To generate the transposed file available in the data folder 
# data <- read.csv("data2_01_BRCA_RNASeq2GeneNorm-20160128.csv", check.names=FALSE)
# transposeddata = setNames(data.frame(t(data[,-1])), data[,1])
# write.csv(transposeddata, file="data2_01_BRCA_RNASeq2GeneNorm-20160128_transposed.csv")

transposeddata_BRCA <- read.csv("data2_01_BRCA_RNASeq2GeneNorm-20160128_transposed.csv", check.names=FALSE)
coldata_BRCA <- read.csv("data2_colData.csv", check.names=FALSE)

coldata_BRCA <- coldata_BRCA %>% as.data.frame %>%
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

rownames(transposeddata_BRCA) <- transposeddata_BRCA[,1]
colnames(transposeddata_BRCA)[1] <- "patient_colname"
transposeddata_BRCA$patient_colname <- NULL

# Set index as the patientID
rownames(coldata_BRCA) <- coldata_BRCA$patientID

# Get matches between survival and assay data
transposeddata_BRCA <- transposeddata_BRCA[TCGAbarcode(rownames(transposeddata_BRCA)) %in% rownames(coldata_BRCA),]

censored_BRCA = subset(coldata_BRCA, coldata_BRCA$status == 0)
not_censored_BRCA = subset(coldata_BRCA, coldata_BRCA$status == 1)
# Splits the censored data randomly into 80% and 20%
censored80_BRCA = sort(sample(nrow(censored_BRCA), nrow(censored_BRCA)*.8))
ytrain_censored_BRCA <- censored_BRCA[censored80_BRCA,]
ytest_censored_BRCA <- censored_BRCA[-censored80_BRCA,]
# Splits the not_censored data randomly into 80% and 20%
not_censored80_BRCA = sort(sample(nrow(not_censored_BRCA), nrow(not_censored_BRCA)*.2))
ytest_notcensored_BRCA <- not_censored_BRCA[not_censored80_BRCA,]
ytrain_notcensored_BRCA <- not_censored_BRCA[-not_censored80_BRCA,]
# Merge the censored and not censored parts
ytrain_BRCA = rbind(ytrain_censored_BRCA, ytrain_notcensored_BRCA)
ytest_BRCA = rbind(ytest_censored_BRCA, ytest_notcensored_BRCA)

transposeddata000_BRCA <- transposeddata_BRCA
new_column_BRCA <- TCGAbarcode(rownames(transposeddata000_BRCA))
rownames(transposeddata000_BRCA) <- new_column_BRCA
xtrain_BRCA <- transposeddata000_BRCA[rownames(ytrain_BRCA),]

# Scaling xtrain
preprocess_train_BRCA <- preProcess(xtrain_BRCA, method = c("center", "scale"))
xtrain_BRCA <- predict(preprocess_train_BRCA, xtrain_BRCA)
variables_with_variance_zero_or_close_BRCA <- nearZeroVar(xtrain_BRCA)
xtrain_BRCA <- xtrain_BRCA[,-variables_with_variance_zero_or_close_BRCA]
xtrain_BRCA <- as.matrix(xtrain_BRCA)

xtest_BRCA <- transposeddata000_BRCA[rownames(ytest_BRCA),]

# Scaling xtest
xtest_BRCA <- predict(preprocess_train_BRCA, xtest_BRCA)
xtest_BRCA <- xtest_BRCA[,-variables_with_variance_zero_or_close_BRCA]
xtest_BRCA <- as.matrix(xtest_BRCA)

# Use all genes in the dataset divided in train and test sets
# xtrain
# xtest
ytrain_BRCA <- ytrain_BRCA %>% select(time, status)
ytest_BRCA <- ytest_BRCA %>% select(time, status)

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# DATASET: PRAD (prostate)
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# To run on server
setwd("~/msc-beatriz-correia/data/PRAD_primary_solid_tumor/")
# To run on my pc
# setwd("~/Documents/msc-beatriz-correia/data/PRAD_primary_solid_tumor/")

# To generate the transposed file available in the data folder 
# data <- read.csv("data2_01_PRAD_RNASeq2GeneNorm-20160128.csv", check.names=FALSE)
# transposeddata = setNames(data.frame(t(data[,-1])), data[,1])
# write.csv(transposeddata, file="data2_01_PRAD_RNASeq2GeneNorm-20160128_transposed.csv")

transposeddata_PRAD <- read.csv("data2_01_PRAD_RNASeq2GeneNorm-20160128_transposed.csv", check.names=FALSE)
coldata_PRAD <- read.csv("data2_colData.csv", check.names=FALSE)

coldata_PRAD <- coldata_PRAD %>% as.data.frame %>%
  select(patientID, vital_status, patient.days_to_death, patient.days_to_last_followup, days_to_death, days_to_last_followup) %>%
  # Convert days to integer
  mutate(patient.days_to_death = as.integer(patient.days_to_death)) %>%
  mutate(patient.days_to_last_followup  = as.integer(patient.days_to_last_followup)) %>%
  mutate(days_to_death  = as.integer(days_to_death)) %>%
  mutate(days_to_last_followup  = as.integer(days_to_last_followup)) %>%
  # Find max time between all days (ignoring missings)
  rowwise %>%
  mutate(time = max(patient.days_to_death, patient.days_to_last_followup, days_to_death, days_to_last_followup, na.rm = TRUE)) %>%
  # Keep only survival variables and codes
  select(patientID, status = vital_status, time) %>% 
  # Discard individuals with survival time less or equal to 0
  filter(!is.na(time) & time > 0) %>% as.data.frame
print("---------------------------------------------------------------------------------------")

rownames(transposeddata_PRAD) <- transposeddata_PRAD[,1]
colnames(transposeddata_PRAD)[1] <- "patient_colname"
transposeddata_PRAD$patient_colname <- NULL

# Set index as the patientID
rownames(coldata_PRAD) <- coldata_PRAD$patientID

# Get matches between survival and assay data
transposeddata_PRAD <- transposeddata_PRAD[TCGAbarcode(rownames(transposeddata_PRAD)) %in% rownames(coldata_PRAD),]

censored_PRAD = subset(coldata_PRAD, coldata_PRAD$status == 0)
not_censored_PRAD = subset(coldata_PRAD, coldata_PRAD$status == 1)
# Splits the censored data randomly into 80% and 20%
censored20_PRAD = sort(sample(nrow(censored_PRAD), nrow(censored_PRAD)*.2))
ytest_censored_PRAD <- censored_PRAD[censored20_PRAD,]
ytrain_censored_PRAD <- censored_PRAD[-censored20_PRAD,]
# Splits the not_censored data randomly into 80% and 20%
not_censored80_PRAD = sort(sample(nrow(not_censored_PRAD), nrow(not_censored_PRAD)*.8))
ytrain_notcensored_PRAD <- not_censored_PRAD[not_censored80_PRAD,]
ytest_notcensored_PRAD <- not_censored_PRAD[-not_censored80_PRAD,]
# Merge the censored and not censored parts
ytrain_PRAD = rbind(ytrain_censored_PRAD, ytrain_notcensored_PRAD)
ytest_PRAD = rbind(ytest_censored_PRAD, ytest_notcensored_PRAD)

transposeddata000_PRAD <- transposeddata_PRAD
new_column_PRAD <- TCGAbarcode(rownames(transposeddata000_PRAD))
rownames(transposeddata000_PRAD) <- new_column_PRAD
xtrain_PRAD <- transposeddata000_PRAD[rownames(ytrain_PRAD),]

# Scaling xtrain
preprocess_train_PRAD <- preProcess(xtrain_PRAD, method = c("center", "scale"))
xtrain_PRAD <- predict(preprocess_train_PRAD, xtrain_PRAD)
variables_with_variance_zero_or_close_PRAD <- nearZeroVar(xtrain_PRAD)
xtrain_PRAD <- xtrain_PRAD[,-variables_with_variance_zero_or_close_PRAD]
xtrain_PRAD <- as.matrix(xtrain_PRAD)

xtest_PRAD <- transposeddata000_PRAD[rownames(ytest_PRAD),]

# Scaling xtest
xtest_PRAD <- predict(preprocess_train_PRAD, xtest_PRAD)
xtest_PRAD <- xtest_PRAD[,-variables_with_variance_zero_or_close_PRAD]
xtest_PRAD <- as.matrix(xtest_PRAD)

# Use all genes in the dataset divided in train and test sets
# xtrain
# xtest
ytrain_PRAD <- ytrain_PRAD %>% select(time, status)
ytest_PRAD <- ytest_PRAD %>% select(time, status)

end_time111 <- Sys.time()
total_time111 <- end_time111 - start_time111
print("TIME TO LOAD DATA")
print(total_time111)

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# TCox
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

start_time222 <- Sys.time()

print("Dissimilarity matrix calculation starting...")

print("TRAIN/TEST SETS")

print(dim(xtrain_BRCA))
print(dim(xtrain_PRAD))
print(dim(xtest_BRCA))
print(dim(xtest_PRAD))

xtrain_BRCA <- xtrain_BRCA[, which(colnames(xtrain_BRCA) %in% colnames(xtrain_PRAD))]
xtrain_PRAD <- xtrain_PRAD[, which(colnames(xtrain_PRAD) %in% colnames(xtrain_BRCA))]
xtest_BRCA <- xtest_BRCA[, which(colnames(xtest_BRCA) %in% colnames(xtest_PRAD))]
xtest_PRAD <- xtest_PRAD[, which(colnames(xtest_PRAD) %in% colnames(xtest_BRCA))]

print("-----------")

print(dim(xtrain_BRCA))
print(dim(xtrain_PRAD))
print(dim(xtest_BRCA))
print(dim(xtest_PRAD))

print("END TRAIN/TEST SETS")

x_BRCA_all <- rbind(xtrain_BRCA, xtest_BRCA)
x_PRAD_all <- rbind(xtrain_PRAD, xtest_PRAD)

print(dim(x_BRCA_all))
print(dim(x_PRAD_all))

x_BRCA_less <- x_BRCA_all[, which(colnames(x_BRCA_all) %in% colnames(x_PRAD_all))]
x_PRAD_less <- x_PRAD_all[, which(colnames(x_PRAD_all) %in% colnames(x_BRCA_all))]

print(dim(x_BRCA_less))
print(dim(x_PRAD_less))

## correlation matrices

suppressMessages(library("propagate"))
suppressMessages(library("lsa"))
suppressMessages(library("gmm"))

# My pc
# options(fftempdir = "~/Documents/msc-beatriz-r-correia/survival-methods/tcox/temporary/")
# Server
options(fftempdir = "~/msc-beatriz-correia/temporary/")

x_BRCA_cor <- bigcor(x_BRCA_less, y = NULL, fun = "cor", size = 2000, verbose=FALSE)
x_BRCA_cor <- as.data.frame(as.ffdf(x_BRCA_cor))

x_PRAD_cor <- bigcor(x_PRAD_less, y = NULL, fun = "cor", size = 2000, verbose=FALSE)
x_PRAD_cor <- as.data.frame(as.ffdf(x_PRAD_cor))

# angular distance
ang_weight <- vector()

for (i in 1:dim(x_BRCA_less)[2]){
  ang_weight[i] <- acos(cosine(x_BRCA_cor[,i],x_PRAD_cor[,i]))/pi
}

# To save an histogram
# My pc
# pdf('~/Documents/msc-beatriz-r-correia/survival-methods/tcox/results/hist_ang_wight.pdf')
# Server
pdf('~/msc-beatriz-correia/results/hist_ang_wights.pdf')

# normalized weights (between 0 and 1)
tumornormal_weights <- ang_weight / max(ang_weight)
hist(ang_weight)

# w
# hist(tumornormal_weights, main="normal weight") # COMMENTED FOR HISTOGRAM

# 1/w
tumornormal_weights6 <- 1 / tumornormal_weights
# hist(tumornormal_weights6, main="1 / normal weight") # COMMENTED FOR HISTOGRAM

print("DISSIMILARITY MATRIX CALCULATED ......................................................................")

end_time222 <- Sys.time()
total_time222 <- end_time222 - start_time222
print("TIME TO CALCULATE DISSIMILARITY MATRIX")
print(total_time222)

# Tumoral weights
if (weights.used == "w") {
  tumoral.weights <- tumornormal_weights # w
} else if (weights.used == "1/w") {
  tumoral.weights <- tumornormal_weights6 # 1/w
} else {
  print("!!!!! Unknown weights")
}

# Dataset to fit to
if (fit.to.dataset == "BRCA") { # breast
  xtrain_set <- xtrain_BRCA
  xtest_set <- xtest_BRCA
  ytrain_set <- ytrain_BRCA
  ytest_set <- ytest_BRCA
} else if (fit.to.dataset == "PRAD") { # prostate
  xtrain_set <- xtrain_PRAD
  xtest_set <- xtest_PRAD
  ytrain_set <- ytrain_PRAD
  ytest_set <- ytest_PRAD
} else {
  print("!!!!! Unknown dataset to fit to")
}

run.counter <- 1
all.outputs.info <- list()
start.counter <- 1
nr.of.null.runs <- 0
x <- 1
start_time <- Sys.time()

while (x <= nr.of.iterations.to.perform) {
  for (k in 1:length(weights.used)) {
    list.of.fits <- list()
    for (i in 1:length(alpha.values.list)) {
      tryCatch({
        print(run.counter)
        start_time2 <- Sys.time()
        output.info <- vector("list", length = 9)
        names(output.info) <- c("weights.used", "alpha", "nr.selected.genes", "list.of.genes", "pvalue.train", "pvalue.test", "iteration.nr", "cindex.train", "cindex.test")
        
        fit.name <- paste0("BRCA_", weights.used[k])
        fit.name <- paste0(fit.name, "_alpha")
        fit.name <- paste0(fit.name, alpha.values.list[i])
        pdf.name <- paste(fit.name, "pdf", sep=".")
        txt.name <- paste(fit.name, "txt", sep=".")
        
        output.info[1] <- weights.used[i]
        output.info[2] <- alpha.values.list[i]
        
        cv.fit_twii3 <- twiner_glmnet(xtrain_set, Surv(ytrain_set$time, ytrain_set$status), alpha.values.list[i], tumoral.weights)
        
        list.of.fits[[fit.name]] <- cv.fit_twii3
        
        coefs.v_twii3 <- twiner_coefs(list.of.fits[[fit.name]])
        coefs.v <- coefs.v_twii3
        
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
                                                     xtrain_set[, names(coefs.v)],
                                                     ytrain_set,
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
                                                     xtest_set[, names(coefs.v)],
                                                     ytest_set,
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
          c.index.train <- c.index.calculator(coefs.v, xtrain_set, ytrain_set)
          output.info[8] <- c.index.train
          c.index.test <- c.index.calculator(coefs.v, xtest_set, ytest_set)
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

print("FINISHED PROGRAM ...................................................................")

end_time333 <- Sys.time()
total_time333 <- end_time333 - start_time111
print("TIME TO EXECUTE EVERYTHING")
print(total_time333)
