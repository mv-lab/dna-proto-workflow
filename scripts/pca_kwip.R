library(rgl)

print (snakemake@input[[1]])
print (snakemake@output[[1]])

y0<-read.delim(snakemake@input[[1]], header=T)
data<-as.matrix(y0[, -1])
labels<-colnames(data)
pca <- prcomp(data, scale=F)
print ('PCA done!')
plot3d(pca$x[,c(1,2,3)], size=10)
text3d(pca$x[,c(1,2,3)], text=labels)
rgl.postscript(snakemake@output[[1]],"pdf")
rgl.postscript(snakemake@output[[1]],"svg")

