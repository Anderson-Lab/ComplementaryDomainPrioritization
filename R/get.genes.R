# Load the required libraries
library(GSEABase)	#Gene Set data structures
library(XML)
library(KEGGREST)
library(RCurl)

#' Downloading all genes in all WikiPathways pathways
#'
#' Returns the gene ids in WikiPathways for a given species
#' For more information on WikiPathways see http://webservice.wikipathways.org/ and
#' http://www.wikipathways.org/index.php/Help:WikiPathways_Webservice/API
#' 
#' @param url WikiPathways base URL for webservice calls.
#' @param species Species for pathway analysis
#' @return a dataframe
#' 
#' @examples
#' gene.ids = get.all.genes.wiki(webg.pathways)
#' 
#' @export
#' 
get.all.genes.wiki <- function(url="http://webservice.wikipathways.org",species="Homo sapiens") {
  # Get a list of all pathways from WikiPathways
  XML =getURL(paste(url,"/listPathways",sep=""))    
  # Parse the xml document
  doc = xmlParse(XML, asText=TRUE)
  
  # Extract the pathway info
  pathwayNodes = xmlElementsByTagName(xmlRoot(doc), "pathways", TRUE)
  pathways = lapply(pathwayNodes, function(n) {
    children = xmlChildren(n, addNames= TRUE)
    if(xmlValue(children$species) == species) {
      p = list()
      p[["id"]] = xmlValue(children$id)
      p[["name"]] = xmlValue(children$name)
      p[["species"]] = xmlValue(children$species)
      p[["url"]] = xmlValue(children$url)
      return(p)
    } else {
      return() # Skip non-species pathways
    }
  })
  
  # Remove NULL entries (non-human pathways)
  pathways = pathways[!sapply(pathways, is.null)]
  #new.pathways = pathways[sapply(pathways, function(p) {
  #  return(p$id %in% webg.pathways)
  #})]
  
  # A function that downloads a GeneSet for a pathway
  createGS = function(p) {
    print(p[["id"]])
    # Download the gene list (translated to Entrez ids)
    XML = getURL(paste(url,"/getXrefList?pwId=",p[['id']],"&code=L",sep=""))
    doc = xmlParse(XML, asText=TRUE)
    # Find the xref nodes with an xpath query
    resultNodes = xmlElementsByTagName(xmlRoot(doc), "xrefs", TRUE)
    # Extract the gene ids
    geneIds = sapply(resultNodes, xmlValue)
    if(length(geneIds) > 0) { # Skip empty lists
      # Create a GeneSet object
      geneSet = GeneSet(geneIds, geneIdType=EntrezIdentifier(), 
                        setName=paste(p[["id"]], " (", p[["name"]], ")", sep="")
      )
      print(geneSet)
      return(geneSet)
    }
  }
  geneSets = lapply(pathways, createGS) #Apply the createGS function on all pathways
  geneSets = geneSets[!sapply(geneSets, is.null)] #Remove empty sets
  all.gene.ids = c()
  for (i in 1:length(geneSets)) {
    all.gene.ids = c(all.gene.ids,as.vector(geneIds(geneSets[[i]])))
  }
  all.gene.ids = unique(all.gene.ids)
  return(as.numeric(all.gene.ids))
  
}


#' Downloading WikiPathways to access gene ids
#'
#' Returns the gene ids associated with select WikiPathways
#' For more information on WikiPathways see http://webservice.wikipathways.org/ and
#' http://www.wikipathways.org/index.php/Help:WikiPathways_Webservice/API
#' 
#' @param webg.pathways Dataframe where the second column contains the pathway ids
#' @param url WikiPathways base URL for webservice calls.
#' @param species Species for pathway analysis
#' @return a dataframe
#' 
#' @examples
#' download.example.data()
#' webg.pathways = load.WebGestalt("Marra_0_wiki_protein_enrichment.tsv",'Wiki')
#' gene.ids = get.genes.wiki(webg.pathways)
#' 
#' @export
#' 
get.genes.wiki <- function(webg.pathways,url="http://webservice.wikipathways.org",species="Homo sapiens") {
  # Get a list of all pathways from WikiPathways
  XML =getURL(paste(url,"/listPathways",sep=""))    
  # Parse the xml document
  doc = xmlParse(XML, asText=TRUE)
  
  # Extract the pathway info
  pathwayNodes = xmlElementsByTagName(xmlRoot(doc), "pathways", TRUE)
  pathways = lapply(pathwayNodes, function(n) {
    children = xmlChildren(n, addNames= TRUE)
    if(xmlValue(children$species) == species) {
      p = list()
      p[["id"]] = xmlValue(children$id)
      p[["name"]] = xmlValue(children$name)
      p[["species"]] = xmlValue(children$species)
      p[["url"]] = xmlValue(children$url)
      return(p)
    } else {
      return() # Skip non-species pathways
    }
  })
  
  #webg.pathways = read.table(paste(dir,'/',prefix,"wiki_filter.tsv",sep=""),sep="\t",header=F,colClasses=c(rep("factor",2)))
  webg.pathways = as.character(webg.pathways[,2])
  
  # Remove NULL entries (non-human pathways)
  pathways = pathways[!sapply(pathways, is.null)]
  
  new.pathways = pathways[sapply(pathways, function(p) {
    return(p$id %in% webg.pathways)
  })]
  
  # A function that downloads a GeneSet for a pathway
  createGS = function(p) {
    print(p[["id"]])
    # Download the gene list (translated to Entrez ids)
    XML = getURL(paste(url,"/getXrefList?pwId=",p[['id']],"&code=L",sep=""))
    doc = xmlParse(XML, asText=TRUE)
    # Find the xref nodes with an xpath query
    resultNodes = xmlElementsByTagName(xmlRoot(doc), "xrefs", TRUE)
    # Extract the gene ids
    geneIds = sapply(resultNodes, xmlValue)
    if(length(geneIds) > 0) { # Skip empty lists
      # Create a GeneSet object
      geneSet = GeneSet(geneIds, geneIdType=EntrezIdentifier(), 
                        setName=paste(p[["id"]], " (", p[["name"]], ")", sep="")
      )
      return(geneSet)
    }
  }
  geneSets = lapply(new.pathways, createGS) #Apply the createGS function on all pathways
  geneSets = geneSets[!sapply(geneSets, is.null)] #Remove empty sets
  all.gene.ids = c()
  for (i in 1:length(geneSets)) {
    all.gene.ids = c(all.gene.ids,as.vector(geneIds(geneSets[[i]])))
  }
  all.gene.ids = unique(all.gene.ids)
  return(as.numeric(all.gene.ids))
  
}

#' Downloading Kegg Pathways to find the gene ids
#'
#' Returns the gene ids
#' 
#' @param webg.pathways dataframe where the second column contains the pathway ids
#' @param sleep Number of seconds to sleep between Kegg requests
#' @return an array of gene ids
#' 
#' @examples
#' download.example.data()
#' webg.pathways = load.WebGestalt("Marra_0_kegg_protein_enrichment.tsv",'Kegg')
#' gene.ids = get.genes.kegg(webg.pathways)
#' 
#' @export
#' 
get.genes.kegg <- function(webg.pathways,sleep=10) {
  pathways = webg.pathways
  
  genes = c()
  #genes = read.csv('genes.csv')
  #genes = as.character(genes$x)
  z = 1
  #start = read.csv('currenti.csv')
  #start = start$x[1]+1
  pathways.with.no.genes = c()
  k = 1
  while (z <= nrow(pathways)) {
    submissions = c()
    for (j in z:nrow(pathways)) {
      submissions[j-z+1] = paste("path:hsa",pathways[j,2],sep="")
      if (length(submissions) == 10) # Kegg restrictions
        break
    }
    
    results = keggGet(submissions)
    for (i in 1:length(results)) {
      if (length(results) > 0 && sum(attributes(results[[i]])$name=="GENE") > 0) {
        if (length(results[[i]]$GENE) == 0) {
          pathways.with.no.genes[k] = pathways[i,2]
          k = k + 1
          next
        }
        for (j in seq(1,length(results[[i]]$GENE),2)) { # Skip gene name and just go with Entrez ID that is returned
          fields = strsplit(results[[i]]$GENE[j],";")
          genes = c(genes,fields[[1]][1])
        }
      } else {
        pathways.with.no.genes[k] = pathways[i,2]
        k = k + 1
      }
    }
    
    #    write.csv(genes,'temp.genes.csv')
    #    write.csv(z,'temp.z.csv')
    Sys.sleep(sleep)
    z = z + 10;
    #print(z)
    for (i in 1:length(submissions)) {
      print(submissions[i])
    }
  }
  
  genes = unique(genes)
  return(as.numeric(genes))

}

#' Downloading gene ids from all Kegg Pathways
#'
#' Writes the gene ids to a file
#' 
#' @param sleep Number of seconds to sleep between Kegg requests
#' @param species Species identifier for Kegg. Default is hsa.
#' @param start Index of pathway to start processing. Only necessary because of frequent Kegg timeouts.
#' @param file File to read gene ids from and 
#'  
#' @examples
#' gene.ids = get.all.genes.kegg(start=1)
#' 
#' @export
#' 
get.all.genes.kegg <- function(sleep=10,species="hsa",start=1,file="kegg.genes.csv") {
  pathways.kegg.list = keggList(paste("pathway/",species,sep=""))
  pathways = matrix(1:(2*length(names(pathways.kegg.list))),ncol=2)
  pathways[,2] = names(pathways.kegg.list)
  pathways[,1] = pathways.kegg.list
  
  if (file.exists(file)) {
    df = read.csv(file)
    genes = as.vector(df$x)
  } else {
    genes = c()
  }
  z = start
  pathways.with.no.genes = c()
  k = 1
  while (z <= nrow(pathways)) {
    submissions = c()
    for (j in z:nrow(pathways)) {
      submissions[j-z+1] = paste(pathways[j,2],sep="")
      if (length(submissions) == 10) # Kegg restrictions
        break
    }

    results = keggGet(submissions)
    for (i in 1:length(results)) {
      if (length(results) > 0 && sum(attributes(results[[i]])$name=="GENE") > 0) {
        if (length(results[[i]]$GENE) == 0) {
          pathways.with.no.genes[k] = pathways[i,2]
          k = k + 1
          next
        }
        for (j in seq(1,length(results[[i]]$GENE),2)) { # Skip gene name and just go with Entrez ID that is returned
          fields = strsplit(results[[i]]$GENE[j],";")
          genes = c(genes,fields[[1]][1])
        }
      } else {
        pathways.with.no.genes[k] = pathways[i,2]
        k = k + 1
      }
    }
    
    write.csv(as.numeric(unique(genes)),file=file)
    #Sys.sleep(sleep)
    z = z + 10;
    for (i in 1:length(submissions)) {
      print(submissions[i])
    }
    print(paste("Starting ",z))
  }
  
  genes = unique(genes)
  return(as.numeric(genes))
  
}

#' Filter for only those genes in the data.
#'
#' Returns a corrected gene ID vector
#' 
#' @param gene.ids A vector of Entrez Gene IDs
#' @param gene.data The gene data loaded with the helper functions.
#' @return a vector of Gene IDs
#' 
#' @examples
#' download.example.data()
#' Catteno.rna = load.gene.data("Catteno_array.csv",5)
#' webg.pathways = load.WebGestalt("Marra_0_kegg_protein_enrichment.tsv",'Kegg')
#' gene.ids = get.genes.kegg(webg.pathways)
#' genes.in.data = get.genes.in.data(gene.ids,Catteno.rna$data)
#' 
#' @export
#' 
get.genes.in.data <- function(gene.ids,gene.data) {
  is.in <- function(gene.number) { return(gene.number %in% gene.ids) }
  gene.ixs = which(sapply(gene.data$Entrez.Gene.Number, is.in))
    
  return(unique(gene.data$Entrez.Gene.Number[gene.ixs]))
}

#' Finds the genes ids from the Transcription Factor Targets database
#'
#' Returns ...
#' 
#' @param webg.pathways dataframe where the second column contains the pathway ids
#' @param db TF database to use
#' @return an array of gene ids
#' 
#' @examples
#' download.example.data()
#' webg.pathways = load.WebGestalt("Marra_0_tf_protein_enrichment.tsv",'TF')
#' gene.ids = get.genes.tf(webg.pathways)
#' 
#' @export
#' 
get.genes.tf <- function(webg.pathways,db=system.file("extdata", "c3.tft.v4.0.entrez.gmt", package = "GeneListPrioritization")) {
  lines = readLines(db) 
  gene.ids = c()
  for (i in 1:nrow(webg.pathways)) {
    key = gsub('hsa_','',webg.pathways[i,1])
    result = grep(key,lines,fixed=T)
    results = lines[result]
    for (line in results) {
      fields = strsplit(line,'\t')[[1]]
      gene.ids = c(gene.ids,fields[3:length(fields)])
    }
  }
  return(as.numeric(gene.ids))
}