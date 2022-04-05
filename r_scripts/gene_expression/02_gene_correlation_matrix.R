library(corrplot)
library(RColorBrewer)

# draw correlation between the genes
gene_correlation <- as.matrix(cor(Challenge[, grepl("_N", colnames(Challenge))], use="pairwise.complete.obs"))

jpeg("output_data/Corrplot_gene_lab.jpg", width = 1400, height = 1000)

# ... : further arguments to pass to the native R cor.test function
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




head(p.mat[, 1:5])

corrplot(gene_correlation, 
         method = "circle",  #method of the plot, "color" would show colour gradient
         tl.col = "black", tl.srt=45, #colour of labels and rotation
         col = brewer.pal(n = 8, name ="RdYlBu"), #colour of matrix
         order="hclust" #hclust reordering ) 
         
         dev.off()
         