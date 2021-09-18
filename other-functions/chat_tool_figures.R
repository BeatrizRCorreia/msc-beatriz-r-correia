# library(glmSparseNet)
library(httr)
library(magrittr)
library(dplyr)
library(ggplot2)

#' Retrieve hallmarks of cancer count for genes
#'
#' @param genes gene names
#' @param metric see below
#' @param hierarchy see below
#' @param generate.plot flag to indicate if return object has a ggplot2 object
#' @param show.message flag to indicate if run.cache method shows messages
#'
#' @return data.frame with choosen metric and hierarchy
#' It also returns a vector with genes that do not have any
#' hallmarks.
#'
#' See http://chat.lionproject.net/api for more details on the
#' metric and hallmarks parameters
#'
#' To standardize the colors in the gradient you can use
#' scale_fill_gradientn(limits=c(0,1), colours=topo.colors(3)) to
#' limit between 0 and 1 for cprob and -1 and 1 for npmi
#'
#' @export
#'
#' @examples
#' \donttest{
#'     hallmarks(c('MOB1A', 'RFLNB', 'SPIC'))
#'     hallmarks(c('MOB1A', 'RFLNB', 'SPIC'), metric = 'cprob')
#' }
hallmarks <- function(
  genes, 
  metric = 'count', 
  hierarchy = 'full',
  generate.plot = TRUE, 
  show.message = FALSE
) {
  #
  valid.measures <- c('count', 'cprob', 'pmi', 'npmi')
  if (!metric %in% valid.measures) {
    stop('measure argument is not valid, it must be one of the following: ',
         paste(valid.measures, collapse = ', ')
    )
  }
  
  
  all.genes <- sort(unique(genes))
  
  #
  # necessary due to https://github.com/cambridgeltl/chat/issues/6
  if (metric == 'cprob') {
    temp.res        <- hallmarks(all.genes, metric = 'count',
                                 hierarchy = 'full', show.message = FALSE,
                                 generate.plot = FALSE)
    good.ix         <- Matrix::rowSums(temp.res$hallmarks) != 0
    all.genes       <- sort(unique(rownames(temp.res$hallmarks[good.ix,])))
    df.no.hallmarks <- temp.res$no.hallmakrs
    #
    message('There is a bug in the Hallmarks\' API that requires the ',
            'function to wait around 5 additional seconds to finish.\n',
            'Sorry. bug report: ',
            'https://github.com/cambridgeltl/chat/issues/6\n')
    Sys.sleep(5.5)
  } else {
    df.no.hallmarks <- NULL
  }
  
  # build base url for call
  baseUrl <- 'http://chat.lionproject.net/chartdata?measure=%s&hallmarks=%s'
  baseUrl <- sprintf(baseUrl, metric, hierarchy)
  # add genes
  call.url <- sprintf('%s&q=%s', baseUrl, paste(all.genes, collapse = '&q='))
  
  lines <- NULL
  conn <- NULL
  tryCatch({
    httr::set_config(httr::config(ssl_verifypeer = 0L))
    response <- httr::GET(call.url)
    csv_text <- content(response)
    # csv_table <- readr::read_delim(csv_text,"\t")
    # as.list(strsplit(csv_text, '\n')[[1]])
    lines <- as.list(strsplit(csv_text, '\n')[[1]])
    # lines <- strsplit(csv_text, split="\n")
    print(lines)
  }, error = function(err) {
    warning('Cannot call Hallmark API, please try again later.')
  })
  
  if (is.null(lines)) {
    return(list(hallmarks = NULL,
                no.hallmakrs = NULL,
                heatmap = NULL))
  }
  
  item_group <- cumsum(grepl(sprintf("^[A-Za-z0-9\\._,-]+\t%s", metric),
                             lines))
  all.items  <- list()
  col.names  <- c()
  
  for (ix in split(lines, item_group)) {
    item.id <- gsub(sprintf("\t%s", metric),"", ix[1])
    # prepare results
    item.val  <- list()
    my.names  <- c('gene.name')
    my.values <- c(item.id)
    for (line in ix[-1]) {
      if (line == '') {
        next
      }
      my.split  <- strsplit(line, '\t')[[1]]
      my.names  <- c(my.names, my.split[1])
      my.values <- c(my.values, as.numeric(my.split[2]))
      col.names <- c(col.names, my.split[[1]])
    }
    names(my.values) <- my.names
    all.items[[item.id]] <- my.values
  }
  
  col.names <- unique(col.names)
  df <- data.frame()
  for (ix in all.items) {
    # convert to numeric
    new.ix <- as.numeric(ix[names(ix) != 'gene.name'])
    # set previous names
    names(new.ix) <- names(ix)[names(ix) != 'gene.name']
    # create temporary data frame with controlled column names
    temp.df <- data.frame(t(new.ix[col.names]))
    rownames(temp.df) <- ix['gene.name']
    df <- rbind(df, temp.df)
  }
  
  df.scaled <- t(scale(t(df)))
  na.ix     <- which(apply(df.scaled, 1, function(col) {
    return(all(is.nan(col)))
  }))
  df.scaled <- df # use counts
  
  if (is.null(df.no.hallmarks)) {
    df.no.hallmarks <- data.frame(
      gene.name = sort(rownames(df.scaled)[na.ix]),
      stringsAsFactors = FALSE)$gene.name
  }
  
  df.scaled <- cbind(gene.name = rownames(df.scaled), df.scaled)
  
  #
  # Generate heatmap
  if (generate.plot) {
    df.scaled$gene.name <- rownames(df.scaled)
    
    g1 <- reshape2::melt(df.scaled, id.vars = c('gene.name')) %>%
      dplyr::filter(!!(as.name('value')) > 0) %>%
      ggplot2::ggplot(ggplot2::aes_string('gene.name', 'variable', fill = 'value')) +
      ggplot2::geom_raster() + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)) +
      ggplot2::ggtitle('Hallmarks heatmap',
                       subtitle = stringr::str_wrap(sprintf(
                         'Selected genes without hallmarks (%d): %s',
                         length(df.no.hallmarks),
                         paste(df.no.hallmarks, collapse = ', ')),
                         width = 50)) +
      ggplot2::xlab('External Gene Name') + ggplot2::ylab('') +
      ggplot2::scale_fill_gradientn(colours = rev(grDevices::topo.colors(2)))
  } else {
    g1 = NULL
  }
  
  df.scaled$gene.name <- NULL
  
  return(list(hallmarks = df.scaled, no.hallmakrs = df.no.hallmarks,
              heatmap = g1))
}
# Pdf for the plots
# pdf('/Users/beatrizcorreia/Documents/testes_comR/plots/plotsnovas_50_C_0.5.pdf')
pdf('/Users/beatrizcorreia/Documents/msc-beatriz-correia/results/BRCA/BRCA_hallmarks_tool/test7.pdf')

hallmarks(c('PGK1', 'AMOTL1', 'CD5', 'ATR', 'LOC100128977', 'CDK1', 'PPFIA3', 'SHCBP1'), metric='count')
