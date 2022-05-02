library(corrplot)
library(RColorBrewer)

source("r_scripts/gene_expression/01_gene_pheatmap_non_normalised.R")

# draw correlation between the genes
gene_correlation <- as.matrix(cor(Challenge %>% select(Genes), use="pairwise.complete.obs"))


##Combining correlogram with the significance test
## Computing the p-value of correlations
## To compute the matrix of p-value, a custom R function is used 
# ... : further arguments to pass to the native R cor.test function
#tutorial: http://www.sthda.com/english/wiki/visualize-correlation-matrix-using-correlogram

# mat : is a matrix of data
cor.mtest <- function(mat, ...) {
    mat <- as.matrix(mat)
    n <- ncol(mat)
    p.mat<- matrix(NA, n, n)
    diag(p.mat) <- 0
    for (i in 1:(n - 1)) {
        for (j in (i + 1):n) {
            tmp <- cor.test(mat[, i], mat[, j], ...)
            p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
        }
    }
    colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
    p.mat
}

# matrix of the p-value of the correlatio
p.mat <- cor.mtest(gene_correlation)

jpeg("output_data/gene_expression/02_Corrplot_gene_lab_non_normalised.jpg", width = 1400, height = 1000)


corrplot(gene_correlation, 
         method = "circle",  #method of the plot, "color" would show colour gradient
         tl.col = "black", tl.srt=45, #colour of labels and rotation
         col = brewer.pal(n = 8, name ="RdYlBu"), #colour of matrix
         order="hclust") #hclust reordering
dev.off()


jpeg("output_data/gene_expression/03_Corrplot_gene_lab_significant_non_normalised.jpg", width = 1400, height = 1000)


corrplot(gene_correlation, 
         method = "circle",  #method of the plot, "color" would show colour gradien
         tl.col = "black", tl.srt=45, #colour of labels and rotation
         col = brewer.pal(n = 8, name ="RdYlBu"), #colour of matrix
         order="hclust", #hclust reordering
         p.mat = p.mat, sig.level = 0.01, insig = "blank") #Add significance level to the correlogram
        #remove the values that are insignificant
         
dev.off()

#switch off all dev devices
while (!is.null(dev.list()))  dev.off()         
