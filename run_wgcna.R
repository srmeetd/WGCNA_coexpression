library(WGCNA)
library(igraph)
library(network)
library(networkD3)
library(dynamicTreeCut)
library(flashClust)
library(intergraph)
library(plyr)
library(dplyr)
library(optparse)

wgcna = function(x, output_root)
{
  data = read.table(x,header = T, sep = "\t")
  raw_counts = data[, sapply(data, function (x) any(grepl('raw_count|normalized_count', x)))]
  rownames(raw_counts) = make.names(data$Hybridization.REF, unique=T)
  df = raw_counts[-1,]
  RNAseq = df[apply(df,1,function(x) sum(x==0))<ncol(df)*0.8,]
  df2 <- mutate_all(RNAseq, function(x) as.numeric(as.character(x)))
  rownames(df2) = rownames(RNAseq)
  WGCNA_matrix = df2[order(apply(df2,1,mad), decreasing = TRUE)[1:5000],]
  WGCNA_matrix[grep ("TRIB3.57761",rownames(WGCNA_matrix)),]
  WGCNA_matrix = t(WGCNA_matrix)
  gene.class=colnames(WGCNA_matrix) 
  Data=as.data.frame(WGCNA_matrix)
  powers = c(c(1:10), seq(from = 12, to=20, by=2)) 
  sft = pickSoftThreshold(log(Data+0.1), powerVector = powers, verbose = 5,blockSize=8504) 
  softPower = 5
  adjacency = adjacency(log(Data+0.1), power = softPower,type = "signed")
  w = 1-adjacency
  geneTree = hclust(as.dist(w), method = 'average')
  modules = cutreeDynamic(dendro = geneTree, distM = w, deepSplit = 4, pamRespectsDendro = FALSE,
                          minClusterSize = 30)
  module.colours = labels2colors(modules)
  ENSEM = colnames(WGCNA_matrix)
  allLLIDs = ENSEM
  intModules = module.colours
  
  output_modules <- data.frame(gene_id=ENSEM, module=module.colours)
  write.table(output_modules, paste0(output_root, "_all_modules.tsv"), sep="\t", row.names = FALSE, quote=FALSE)
  trib1_module = module.colours[grep("TRIB3.57761", allLLIDs)]
  trib1_module_genes = output_modules[output_modules$module==trib1_module,]
  write.table(trib1_module_genes, paste0(output_root, "_trib3_module.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
}

options_list <- list(
  make_option(c("-i", "--input-file"), action="store", dest="input_file"),
  make_option(c("-o", "--ouput-root"), action="store", dest="output_root"))

opt <- parse_args(OptionParser(option_list=options_list))

wgcna(opt$input_file, opt$output_root)



