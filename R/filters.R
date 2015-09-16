#' Fold change filter
#'
#' Applies a simple fold change filter on the data. 
#' 
#' @param data - transposed gene data
#' @param theta - if > 1, theta is the target number of ranked variables;
#' otherwise, it is the top theta fraction of ranked variables.
#' @param class1.inxs - the row indices of class 1
#' @param class2.inxs - the row indices of class 2
#' 
#' @return filtered dataset
#' 
#' @examples
#' download.example.data()
#' webg.pathways = load.WebGestalt("Marra_0_wiki_protein_enrichment.tsv",'Wiki')
#' gene.ids = get.genes.wiki(webg.pathways)
#' Cattaneo.rna = load.gene.data("Cattaneo_array.csv",5)
#' prioritized.data = list.filter(Cattaneo.rna$transposed.data,gene.ids)
#' prioritized.filtered.data = fold.change.filter(prioritized.data,100,
#'                                                which(prioritized.data$labels=='adenoma'),
#'                                                which(prioritized.data$labels=='normal'))
#' 
#' @export
#' 
fold.change.filter <- function(data,theta,class1.inxs,class2.inxs) {
  fcs = apply(data[,-ncol(data)],2,function(col) {
    m1 = max(mean(col[class1.inxs]),min(col))
    m2 = max(mean(col[class2.inxs]),min(col))
    return(max(m1/m2,m2/m1)) } ) 
  fc.order = order(fcs)
  if (theta > 1) { # Assume this is the number of top features to select
    return(data[,c(fc.order[(length(fc.order)-theta+1):length(fc.order)],ncol(data))])
  } else {
    return(data[,-fc.order[1:round((ncol(data)-1)*theta)]])
  }  
}

#' Overall mean filter
#'
#' Applies the overall mean filter on the columns of the dataset.
#' 
#' @param data - transposed gene data
#' @param theta - if > 1, theta is the target number of ranked variables;
#' otherwise, it is the top theta fraction of ranked variables.
#' @return an array of inxs
#' 
#' @examples
#' download.example.data()
#' webg.pathways = load.WebGestalt("Marra_0_wiki_protein_enrichment.tsv",'Wiki')
#' gene.ids = get.genes.wiki(webg.pathways)
#' Cattaneo.rna = load.gene.data("Cattaneo_array.csv",5)
#' prioritized.data = list.filter(Cattaneo.rna$transposed.data,gene.ids)
#' prioritized.filtered.data = overall.mean.filter(prioritized.data,100)
#' 
#' @export
#' 
overall.mean.filter <- function(data,theta) {
  means = colMeans(data[,-ncol(data)]) # Compute the column means
  mean.order = order(means)
  if (theta > 1) { # Assume this is the number of top features to select
    return(data[,c(mean.order[(length(mean.order)-theta+1):length(mean.order)],ncol(data))])
  } else {
    return(data[,-mean.order[1:round((ncol(data)-1)*theta)]])
  }  
}

colVars <- function(x, na.rm=FALSE, dims=1, unbiased=TRUE, SumSquares=FALSE,
                    twopass=FALSE) {
  if (SumSquares) return(colSums(x^2, na.rm, dims))
  N <- colSums(!is.na(x), FALSE, dims)
  Nm1 <- if (unbiased) N-1 else N
  if (twopass) {x <- if (dims==length(dim(x))) x - mean(x, na.rm=na.rm) else
    sweep(x, (dims+1):length(dim(x)), colMeans(x,na.rm,dims))}
  (colSums(x^2, na.rm, dims) - colSums(x, na.rm, dims)^2/N) / Nm1
}

#' Overall variance filter
#'
#' Performs an overall variance filter using each column.
#' 
#' @param data - transposed gene data
#' @param theta - if > 1, theta is the target number of ranked variables;
#' otherwise, it is the top theta fraction of ranked variables.
#' @return an array of inxs
#' 
#' @examples
#' download.example.data()
#' webg.pathways = load.WebGestalt("Marra_0_wiki_protein_enrichment.tsv",'Wiki')
#' gene.ids = get.genes.wiki(webg.pathways)
#' Cattaneo.rna = load.gene.data("Cattaneo_array.csv",5)
#' prioritized.data = list.filter(Cattaneo.rna$transposed.data,gene.ids)
#' prioritized.filtered.data = overall.var.filter(prioritized.data,100)
#' 
#' @export
#' 
overall.var.filter <- function(data,theta) {
  vars = colVars(data[,-ncol(data)]) # Compute the column means
  var.order = order(vars)
  if (theta > 1) { # Assume this is the number of top features to select
    return(data[,c(var.order[(length(var.order)-theta+1):length(var.order)],ncol(data))])
  } else {
    return(data[,-var.order[1:round((ncol(data)-1)*theta)]])
  }  
}

#' Random filter
#'
#' Performs a random filter of the variables. For comparison purposes.
#' 
#' @param data - transposed gene data
#' @param theta - if > 1, theta is the target number of ranked variables;
#' otherwise, it is the top theta fraction of ranked variables.
#' @return an array of inxs
#' 
#' @examples
#' download.example.data()
#' webg.pathways = load.WebGestalt("Marra_0_wiki_protein_enrichment.tsv",'Wiki')
#' gene.ids = get.genes.wiki(webg.pathways)
#' Cattaneo.rna = load.gene.data("Cattaneo_array.csv",5)
#' prioritized.data = list.filter(Cattaneo.rna$transposed.data,gene.ids)
#' prioritized.filtered.data = random.filter(prioritized.data,100)
#' 
#' @export
#' 
random.filter <- function(data,theta) {
  if (theta > 1) { # Assume this is the number of top features to select
    return(data[,c(sample.int(ncol(data)-1)[(ncol(data)-1-theta+1):(ncol(data)-1)],ncol(data))])
  } else {
    return(data[,-sample.int(ncol(data)-1)[1:round((ncol(data)-1)*theta)]])
  }    
}


#' Returns a filtered dataset
#'
#' Helper function that takes an array of gene.ids and uses those to filter.
#' 
#' @param data - transposed gene data
#' @param gene.ids - result of get.genes.*(...)
#' @return an array of inxs
#' 
#' @examples
#' download.example.data()
#' webg.pathways = load.WebGestalt("Marra_0_wiki_protein_enrichment.tsv",'Wiki')
#' gene.ids = get.genes.wiki(webg.pathways)
#' Cattaneo.rna = load.gene.data("Cattaneo_array.csv",5)
#' filtered.data = list.filter(Cattaneo.rna$transposed.data,gene.ids)
#' 
#' @export
#' 
list.filter <- function(data,gene.ids) {
  gene.ixs = list.filter.inxs(data,gene.ids)
  return(data[,c(gene.ixs,ncol(data))])
}

#' Returns the gene.ixs for filtering
#'
#' A helper function that filters using a list of gene IDs and returns an
#' array of indices of the genes that passed the filter.
#' 
#' @param data - transposed gene data
#' @param gene.ids - result of get.genes.*(...)
#' @return an array of inxs
#' 
#' @examples
#' download.example.data()
#' webg.pathways = load.WebGestalt("Marra_0_wiki_protein_enrichment.tsv",'Wiki')
#' gene.ids = get.genes.wiki(webg.pathways)
#' Cattaneo.rna = load.gene.data("Cattaneo_array.csv",5)
#' gene.ixs = list.filter.inxs(Cattaneo.rna$transposed.data,gene.ids)
#'  
#' @export
#' 
list.filter.inxs <- function(data,gene.ids) {
  is.in <- function(gene.number) { return(gene.number %in% gene.ids) }
  gene.ixs = which(sapply(colnames(data), is.in))
  return(as.numeric(gene.ixs))
}

# Assumes that the order in the list is of decreasing importance
ordered.filter <- function(data,file,theta) {
  genes = read.csv(file)
  if (theta > 1) { # Assume this is the number of top features to select
    is.in <- function(gene.number) { return(gene.number %in% genes[1:theta,1]) }
  } else {
    is.in <- function(gene.number) { return(gene.number %in% genes[1:(round((1-theta)*ncol(data))),1]) }
  }    
  gene.ixs = which(sapply(colnames(data), is.in))
  return(data[,c(gene.ixs,ncol(data))])  
}
