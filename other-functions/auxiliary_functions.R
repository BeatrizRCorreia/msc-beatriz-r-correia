library(hash)
library(survcomp)

get_genes_selected <- function(nr.of.iterations, list.with.all.outputs) {
  dict <- hash()
  for (output in list.with.all.outputs) {
    for (gene in output[4]) {
      if (identical(gene, character(0))) {
        next
      }
      for (i in 1:length(gene)) {
        new.gene <- unlist(gene[i])
        if (gene[i] %in% names(dict)) {
          dict[[new.gene]] <- dict[[new.gene]] + 1
        }
        else {
          dict[[new.gene]] <- 1
        }
      }
    }
  }
  val <- unlist(as.list(dict))
  # genes.selected <- names(val[val >= nr.of.iterations/2])
  genes.selected <- sort(val[val >= nr.of.iterations/2], decreasing = TRUE)
  return(genes.selected)
}

calculate_standard_deviation <- function(string.with.element.name, average.number, total.n, list.with.all.outputs) {
  numerator.sum <- 0
  for (output in list.with.all.outputs) {
    if (unlist(output[3]) != 0) {
      if (string.with.element.name == "genes.selected") {
        numerator.sum <- numerator.sum + (unlist(output[3]) - average.number)^2
      }
      if (string.with.element.name == "pvalue.train") {
        numerator.sum <- numerator.sum + (unlist(output[5]) - average.number)^2
      }
      if (string.with.element.name == "pvalue.test") {
        numerator.sum <- numerator.sum + (unlist(output[6]) - average.number)^2
      }
      if (string.with.element.name == "cindex.train") {
        if (!is.na(unlist(output[8])))
          numerator.sum <- numerator.sum + (unlist(output[8]) - average.number)^2
      }
      if (string.with.element.name == "cindex.test") {
        if (!is.na(unlist(output[9])))
          numerator.sum <- numerator.sum + (unlist(output[9]) - average.number)^2
      }
    }
  }
  fraction <- numerator.sum/total.n
  standard.dev <- sqrt(fraction)
  return(standard.dev)
}

average_number_of_genes_selected <- function(nr.of.iterations, list.with.all.outputs) {
  sum <- 0
  for (output in list.with.all.outputs) {
    if (unlist(output[3]) != 0) {
      sum <- sum + unlist(output[3])
    }
  }
  average <- sum/nr.of.iterations
  return(average)
}

average_pvalue_train <- function(nr.of.iterations, list.with.all.outputs) {
  sum <- 0
  for (output in list.with.all.outputs) {
    if (!is.na(output[5])) {
      sum <- sum + unlist(output[5])
    }
  }
  average <- sum/nr.of.iterations
  return(average)
}

average_pvalue_train_significant <- function(list.with.all.outputs, info.wanted) {
  sum <- 0
  counter <- 0
  for (output in list.with.all.outputs) {
    if ((!is.na(output[5])) & (unlist(output[5]) < 0.05)) {
      sum <- sum + unlist(output[5])
      counter <- counter + 1
    }
  }
  if (counter == 0)
    average <- NA
  else
    average <- sum/counter
  if (info.wanted == "counter")
    return(counter)
  else if (info.wanted == "average")
    return(average)
}

average_pvalue_train_non_significant <- function(list.with.all.outputs, info.wanted) {
  sum <- 0
  counter <- 0
  for (output in list.with.all.outputs) {
    if ((!is.na(output[5])) & (unlist(output[5]) >= 0.05)) {
      sum <- sum + unlist(output[5])
      counter <- counter + 1
    }
  }
  if (counter == 0)
    average <- NA
  else
    average <- sum/counter
  if (info.wanted == "counter")
    return(counter)
  else if (info.wanted == "average")
    return(average)
}

average_pvalue_test <- function(nr.of.iterations, list.with.all.outputs) {
  sum <- 0
  for (output in list.with.all.outputs) {
    if (!is.na(output[6])) {
      sum <- sum + unlist(output[6])
    }
  }
  average <- sum/nr.of.iterations
  return(average)
}

average_pvalue_test_significant <- function(list.with.all.outputs, info.wanted) {
  sum <- 0
  counter <- 0
  for (output in list.with.all.outputs) {
    if ((!is.na(output[6])) & (unlist(output[6]) < 0.05)) {
      sum <- sum + unlist(output[6])
      counter <- counter + 1
    }
  }
  if (counter == 0)
    average <- NA
  else
    average <- sum/counter
  if (info.wanted == "counter")
    return(counter)
  else if (info.wanted == "average")
    return(average)
}

average_pvalue_test_non_significant <- function(list.with.all.outputs, info.wanted) {
  sum <- 0
  counter <- 0
  for (output in list.with.all.outputs) {
    if ((!is.na(output[6])) & (unlist(output[6]) >= 0.05)) {
      sum <- sum + unlist(output[6])
      counter <- counter + 1
    }
  }
  if (counter == 0)
    average <- NA
  else
    average <- sum/counter
  if (info.wanted == "counter")
    return(counter)
  else if (info.wanted == "average")
    return(average)
}

run_with_lowest_pvalue_train <- function(list.with.all.outputs) {
  smallest.run <- list(0, 0, 0, 0, 9999, 9999, 0, 0, 0)
  for (output in list.with.all.outputs) {
    if (!is.na(unlist(output[5]))) {
      if (unlist(output[5]) < smallest.run[5]) {
        smallest.run <- output
      }
    }
  }
  return(smallest.run)
}

run_with_lowest_pvalue_test <- function(list.with.all.outputs) {
  smallest.run <- list(0, 0, 0, 0, 9999, 9999, 0, 0, 0)
  for (output in list.with.all.outputs) {
    if (!is.na(unlist(output[6]))) {
      if (unlist(output[6]) < smallest.run[6]) {
        smallest.run <- output
      }
    }
  }
  return(smallest.run)
}

c.index.calculator <- function(coefs, xdata, ydata) {
  names.genes.selected <- names(coefs)
  xdata.genes.selected <- as.matrix(xdata[, c(names.genes.selected)])
  fitted.risk <- exp(as.vector(xdata.genes.selected %*% coefs))
  c.res <- concordance.index(fitted.risk, ydata[, "time"], ydata[, "status"] * 1, method = "noether")
  return(c.res$c.index)
}

c.index.average.train <- function(nr.of.iterations, list.with.all.outputs, info.wanted) {
  sum <- 0
  for (output in list.with.all.outputs) {
    if (!is.na(output[8])) {
      sum <- sum + unlist(output[8])
    }
    else if (is.na(output[8]) & (output[3] > 0)) {
      nr.of.iterations <- nr.of.iterations - 1
    }
  }
  if (nr.of.iterations == 0) {
    average <- NA
  }
  else {
    average <- sum/nr.of.iterations
  }
  if (info.wanted == "average") {
    return(average)
  }
  if (info.wanted == "iterations") {
    return(nr.of.iterations)
  }
}

c.index.average.test <- function(nr.of.iterations, list.with.all.outputs, info.wanted) {
  sum <- 0
  for (output in list.with.all.outputs) {
    if (!is.na(output[9])) {
      sum <- sum + unlist(output[9])
    }
    else if (is.na(output[9]) & (output[3] > 0)) {
      nr.of.iterations <- nr.of.iterations - 1
    }
  }
  if (nr.of.iterations == 0) {
    average <- NA
  }
  else {
    average <- sum/nr.of.iterations
  }
  if (info.wanted == "average") {
    return(average)
  }
  if (info.wanted == "iterations") {
    return(nr.of.iterations)
  }
}
