#' Downloads example data files for package demonstration.
#'
#' @param dest The destination path
#' @param allow.cache Allow the files to be cached to prevent downloading multiple times.
#' @param base.url Source URL location
#' 
#' @examples
#' download.example.data()
#' 
#' @export
#' 
download.example.data <- function(dest=".",allow.cache=T,base.url="http://freyja.cs.cofc.edu/downloads/ComplementaryDomainPrioritization/") {
  files = c("Cattaneo_array.csv","Marra_0_kegg_protein_enrichment.tsv",
            "Marra_0_tf_protein_enrichment.tsv","Marra_0_wiki_protein_enrichment.tsv",
            "TCGA_protein.csv","TCGA_rna.csv","Uzozie_protein.csv","Marra_5percentFDR_tf_protein_enrichment.tsv",
            "Marra_5percentFDR_wiki_protein_enrichment.tsv","Marra_5percentFDR_kegg_protein_enrichment.tsv"
            )
  for (i in 1:length(files)) {
    dest.file = paste(dest,"/",files[i],sep="")
    if (allow.cache && file.exists(dest.file)) {
      print(paste("Using cached copy of",dest.file))
      next
    } else {
      print(paste("Downloading to",dest.file))
      download.file(paste(base.url,"/",files[i],sep=""),dest.file,method='wget')
    }
  }
}


#' Loading a correctly formatted gene file.
#'
#' Returns a list containing two data frames.
#' The second data frame is a formatted and transposed version 
#' of the original data with an additional labels column added.
#' For data format instructions see https://github.com/Anderson-Lab/ComplementaryDomainPrioritization/wiki.
#' 
#' @param file A string pointing to the location of the data file
#' @param start.data.inx An integer >= 2 that indicates where the actual data starts.
#' @return a list containing two data frames.
#' 
#' @examples
#' download.example.data()
#' TCGA.rna = load.gene.data("TCGA_rna.csv",3)
#' Cattaneo.rna = load.gene.data("Cattaneo_array.csv",5)
#' 
#' @export
#' 
load.gene.data <- function(file,start.data.inx) {
  gene.data = read.csv(file,skip=1)
  labels.df = read.csv(file,nrows=1,header=F)
  
  labels = c()
  for (i in start.data.inx:ncol(labels.df)) {
    labels[i-start.data.inx+1] = as.character(labels.df[[1,i]])
  }
  
  # Make sure we remove duplicates, keeping only one gene. Sort by overall mean
  gene.data$AvgValue = rowSums(gene.data[,start.data.inx:ncol(gene.data)])
  aa <- gene.data[order(gene.data$Entrez.Gene.Number, -abs(gene.data$AvgValue) ), ] #sort by id and reverse of abs(value)
  gene.data <- aa[ !duplicated(gene.data$Entrez.Gene.Number), -ncol(aa)]
  
  # Now go ahead and make a transposed and corrected version
  data = gene.data[,(start.data.inx-1):ncol(gene.data)]
  data <- as.data.frame(t(data))
  colnames(data) <- data[1, ]
  data <- data[-1, ]
  cmeans = colMeans(data)
  ixs = which(is.na(cmeans))
  if (length(ixs) > 0) {
    data = data[,-ixs]
  }
  data$labels <- factor(labels)
  return(list(data=gene.data,transposed.data=data))
}

#' Loading a correctly formatted protein file.
#'
#' Returns a list containing two data frames.
#' The second data frame is a formatted and transposed version 
#' of the original data with an additional labels column added.
#' For data format instructions see https://github.com/Anderson-Lab/ComplementaryDomainPrioritization/wiki.
#' 
#' @param file A string pointing to the location of the data file. 
#' @param start.data.inx An integer >= 2 that indicates where the actual data starts.
#' @return a list containing two data frames.
#' 
#' @examples
#' download.example.data()
#' TCGA.protein = load.protein.data("TCGA_protein.csv",3)
#' 
#' @export
#' 
load.protein.data <- function(file,start.data.inx) {
  protein.data = read.csv(file,skip=1)
  protein.labels.df = read.csv(file,nrows=1,header=F)
  protein.labels = c()

  for (i in start.data.inx:ncol(protein.labels.df)) {
    protein.labels[i-start.data.inx+1] = as.character(protein.labels.df[[1,i]])
  }
  
  # Now go ahead and make a transposed and corrected version of the protein data
  data.protein = protein.data[,(start.data.inx-1):ncol(protein.data)]
  data.protein <- as.data.frame(t(data.protein))
  colnames(data.protein) <- data.protein[1, ]
  data.protein <- data.protein[-1, ]
  cmeans = colMeans(data.protein)
  data.protein$labels <- factor(protein.labels)
  
  return(list(data=protein.data,transposed.data=data.protein))  
}

#' Loading the pathway data from WebGestalt output.
#'
#' Returns an array of pathways.
#' For data format instructions see https://github.com/Anderson-Lab/ComplementaryDomainPrioritization/wiki.
#' 
#' @param file A string pointing to the location of the data file
#' @param db Wiki, Kegg, or TF
#' @return a list containing the extracted pathways from WebGestalt output.
#' 
#' @examples
#' download.example.data()
#' pathways = load.WebGestalt("Marra_0_tf_protein_enrichment.tsv",'TF')
#' pathways = load.WebGestalt("Marra_0_wiki_protein_enrichment.tsv",'Wiki')
#' pathways = load.WebGestalt("Marra_0_kegg_protein_enrichment.tsv",'Kegg')
#' @export
#' 
load.WebGestalt <- function(file,db) {
  lines = readLines(file)
  if (db == 'Kegg') {
    pathway_lines = lines[grep('KEGG pathway\t',lines)]
  } else if (db == 'Wiki') {
    pathway_lines = lines[grep('Wikipathways pathway\t',lines)]
  } else if (db == 'TF') {
    pathway_lines = lines[grep('Transcrription Target\t',lines)]
  }
  pathways = as.data.frame(matrix(sapply(pathway_lines,function(line) { 
    fields = strsplit(line,'\t')[[1]]
    return(c(fields[2],fields[3])) },simplify=T),ncol=2,byrow=2))
  colnames(pathways) = c('Pathway','ID')
  return(pathways)
}