library(curatedTCGAData)
library(TCGAutils)

extract_data <- function(data_code) {
  data <- curatedTCGAData(diseaseCode = data_code, assays = "RNASeq2GeneNorm", dry.run = FALSE)
  head(getSubtypeMap(data))
  print(data)
  td <- "/home/beatriz/Documents/msc-beatriz-correia/data"
  dataset_name <- paste(data_code, "original", sep="_")
  tempd <- file.path(td, dataset_name)
  if (!dir.exists(tempd))
    dir.create(tempd, mode = "0777")
  exportClass(data, dir = tempd, fmt = "csv", ext = ".csv")
  return(paste(data_code, "dataset successfully downloaded", sep=" "))
}

extract_primary_solid_tumor_data <- function(data_code) {
  data1 <- curatedTCGAData(diseaseCode = data_code, assays = "RNASeq2GeneNorm", dry.run = FALSE)
  data2 <- TCGAutils::splitAssays(data1, '01')
  head(getSubtypeMap(data2))
  print(data2)
  td <- "/home/beatriz/Documents/msc-beatriz-correia/data"
  dataset_name <- paste(data_code, "primary_solid_tumor", sep="_")
  tempd <- file.path(td, dataset_name)
  if (!dir.exists(tempd))
    dir.create(tempd, mode = "0777")
  exportClass(data2, dir = tempd, fmt = "csv", ext = ".csv")
  return(paste(data_code, "primary solid tumor dataset successfully downloaded", sep=" "))
}