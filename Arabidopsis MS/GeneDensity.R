###from a data frame containing chromosome and gene positions calculate gene density
plottedGenes <- ggplot(genes) + geom_histogram(aes(x=pos),binwidth=100000) + facet_wrap(~chr,ncol=2)
allPlottingData <- ggplot_build(plottedGenes)
df <- allPlottingData[["data"]][[1]]

gene_density <- df[,c(2,8,12)] 
colnames(gene_density) <- c("GeneCount", "Chromosome", "Position")
