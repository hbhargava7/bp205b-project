BiocManager::install("monocle")
library(monocle)
library(tidyverse)
library(MASS)
library(reshape)
library(ggplot2)


expr_mtx <- read_csv("./pathtodata", col_names=FALSE)#expression matrix
expr_mtx <- t(as.matrix(dataset_expr_mtx))
sample_sheet <- read_csv("./pathtodata" ) #pheno data
row.names(expr_mtx) <- colnames(expr_mtx)
genes_ann <- as.data.frame(read_csv("./pathtodata")) #feature data
row.names(genes_ann) <- row.names(expr_mtx)

pd <- new("AnnotatedDataFrame", data = sample_sheet)
fd <- new("AnnotatedDataFrame", data = genes_ann)
CellDataSet <- newCellDataSet(as.matrix(expr_mtx),
                         phenoData = pd, featureData = fd, expressionFamily = negbinomial.size())

HCC_AA <- estimateSizeFactors(HCC_AA)
HCC_AA <- estimateDispersions(HCC_AA)

HCC_AA <- detectGenes(HCC_AA, min_expr = 0.1)
expressed_genes <- row.names(subset(fData(HCC_AA),
                                    num_cells_expressed >= 10))

# Log-transform each value in the expression matrix.
L <- log(exprs(HCC_AA[expressed_genes,]))

# Standardize each gene, so that they are all on the same scale,
# Then melt the data with plyr so we can plot it easily
melted_dens_df <- melt(Matrix::t(scale(Matrix::t(L))))

# Plot the distribution of the standardized gene expression values.
qplot(value, geom = "density", data = melted_dens_df) +
  stat_function(fun = dnorm, size = 0.5, color = 'red') +
  xlab("Standardized log(FPKM)") +
  ylab("Density")

#Proceed with Trajectory Analysis: technically parental vs LM2 are beginning vs. end of experimental selection
#can use as timepoints for trajectory analysis #might need to adapt for MDA dataset
HCC_AA <- detectGenes(HCC_AA, min_expr = 0.1)
fData(HCC_AA)$use_for_ordering <-
  fData(HCC_AA)$num_cells_expressed > 0.05 * ncol(HCC_AA)
plot_pc_variance_explained(HCC_AA, return_all = F)

HCC_AA <- reduceDimension(HCC_AA,
                          max_components = 2,
                          norm_method = 'log',
                          num_dim = 3,
                          reduction_method = 'tSNE', 
                          residualModelFormulaStr = "num_genes_expressed",
                          verbose = T)
HCC_AA <- clusterCells(HCC_AA, verbose=F)
#sanity check for clustering
plot_cell_clusters(HCC_AA, color_by = 'Cluster')
plot_cell_clusters(HCC_AA, color_by = 'SampleType')

#DEG
clustering_DEG_genes <-
  differentialGeneTest(HCC_AA[expressed_genes,],
                       fullModelFormulaStr = '~Cluster',
                       cores = 2)
#Selecting top 1000 DE genes for ordering
HCC_AA_ordering_genes <-
  row.names(clustering_DEG_genes)[order(clustering_DEG_genes$qval)][1:1000]
HCC_AA_ordering_genes <-setOrderingFilter(HCC_AA,ordering_genes = HCC_AA_ordering_genes)
HCC_AA <-reduceDimension(HCC_AA, method = 'DDRTree')
#order along pseudotime
HCC_AA <-orderCells(HCC_AA)
#write function to define initial state; otherwise monocle chooses arbitrarily
GM_state <- function(cds){
  if (length(unique(pData(cds)$State)) > 1){
    T0_counts <- table(pData(cds)$State, pData(cds)$SampleType)[,"parental"]
    return(as.numeric(names(T0_counts)[which(T0_counts == max(T0_counts))]))
  } else {
    return (1)
  }
}
HCC_AA <-orderCells(HCC_AA, root_state = GM_state(HCC_AA))
plot_cell_trajectory(HCC_AA, color_by = "SampleType")
plot_cell_trajectory(HCC_AA, color_by = "Pseudotime")
plot_cell_trajectory(HCC_AA, color_by = "SampleType") + facet_wrap("~State")

#Finding genes that change as a function of branching
#__ first critical branch along pseudotime between parental and LM2
BEAM_res <- BEAM(HCC_AA, branch_point = X, cores = 1)
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]

plot_genes_branched_heatmap(HCC_AA[row.names(subset(BEAM_res,
                                                    qval < 1e-4)),],
                            branch_point = X,
                            num_clusters = Y,
                            cores = 1,
                            use_gene_short_name = T,
                            show_rownames = T)
#gotta figure out cluster optimizations: what should it correspond to?

#repeat BEAM for each critical branch point
