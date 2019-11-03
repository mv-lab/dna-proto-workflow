
library(ggplot2)
library(vcfR)


VcfToTable <- function(vcf.file){
  vcf <- read.vcfR(vcf.file, verbose=F, convertNA=F)
  vcf.table <- data.frame(vcf@fix, vcf@gt, stringsAsFactors=F)
  return(vcf.table)
}


PlotVcf <- function(vcf.table, plot.type=c("sep","joint"), samples=NULL,
                    sites=NULL, scaled=T, allele.colors=NULL, site.color=NULL){
  # Re-order genovcf object by list of site IDs, if specified
  if (is.null(sites)==F) {
    vcf.table <- vcf.table[which(vcf.table$ID %in% sites),]
    vcf.table <- vcf.table[order(match(vcf.table$ID, sites)),]
    scaled <- F
  }
  # Re-order genovcf object by list of samples/individuals, if specified
  if (is.null(samples)==F) {
    site.info <- vcf.table[,c(1:9)]
    vcf.table <- vcf.table[,which(colnames(vcf.table) %in% samples)]
    vcf.table <- cbind(site.info, vcf.table[,order(match(colnames(vcf.table),
                                                         samples))])
  }
  # Get positions halfway between each variant and the variant before/after it
  vcf.table$POS <- as.numeric(vcf.table$POS)
  pos.chrom <- list()
  pos.start <- list()
  pos.stop <- list()
  for (i in 1:length(unique(vcf.table$CHROM))) {
    chrom <- vcf.table[which(vcf.table$CHROM==unique(vcf.table$CHROM)[i]),]
    pos.chrom[[i]] <- chrom$CHROM
    # If scaled is FALSE, plot sites in sequential order (not to scale)
    if (scaled==F) {
      pos.start[[i]] <- c(1:nrow(chrom))
      pos.stop[[i]] <- pos.start[[i]]+1
    }else{
      # If scaled is TRUE, plot sites to scale based on physical position
      pos1 <- NULL
      pos2 <- NULL
      for(j in 2:(nrow(chrom)-1)){
        here <- chrom[j, "POS"]
        up <- chrom[j+1, "POS"]
        down <- chrom[j-1, "POS"]
        pos1[j] <- here-((here-down)/2)
        pos2[j] <- here+((up-here)/2)
      }
      pos1[1] <- chrom[1, "POS"]
      pos2[1] <- pos1[2]
      pos1 <- c(pos1, pos2[length(pos2)])
      pos2 <- c(pos2, chrom[nrow(chrom), "POS"])
      pos.start[[i]] <- pos1
      pos.stop[[i]] <- pos2
    }
  }

  # Create genotype table for plotting alleles, where:
  #   sample column contains sample/individual names
  #   chr column contains site chromosomes
  #   start/end columns will be used as x-axis positions
  #   allele column has allele calls, which will be assigned distinct colors
  #   index column specifies y-axis positions
  #   marker column has site positions, which may be plotted as vertical lines
  pos.start <- unlist(pos.start)
  pos.stop <- unlist(pos.stop)
  pos.chrom <- unlist(pos.chrom)
  samples <- colnames(vcf.table)[10:ncol(vcf.table)]
  gen <- as.data.frame(matrix(nrow=(nrow(vcf.table)*length(samples)), ncol=7))
  colnames(gen) <- c("sample", "chr", "start", "end", "allele", "index",
                     "marker")
  gen$sample <- rep(samples, nrow(vcf.table))
  gen$sample <- gen$sample[order(match(gen$sample, samples))]
  gen$chr <- rep(pos.chrom, length(samples))
  gen$start <- rep(pos.start, length(samples))
  gen$end <- rep(pos.stop, length(samples))
  gen$index <- sort(rep(1:length(samples), nrow(vcf.table)))
  gen$marker <- rep(as.numeric(vcf.table$POS), length(samples))
  gt <- NULL
  for (i in 1:length(samples)) {
    gt <- c(gt, vcf.table[,which(colnames(vcf.table)==unique(gen$sample)[i])])
  }
  for (i in 1:length(gt)) {
    gen[i,"allele"] <- strsplit(gt[i][[1]][1],":")[[1]][1]
  }
  # Use genotype table to plot alleles across samples/individuals
  if (is.null(allele.colors)==T) {
    allele.colors <- rainbow(length(unique(gen$allele)))
  }
  genplot <- ggplot(gen) +
    geom_rect(aes(xmin=start, xmax=end, ymin=index, ymax=index+1, fill=allele))+
    scale_fill_manual(values=allele.colors) +
    theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
          panel.background=element_blank(),
          axis.line=element_line(colour="black")) +
    scale_x_continuous(expand=c(0,0), name="Position") +
    scale_y_continuous(breaks=unique(gen$index)+0.5, labels=unique(gen$sample),
                       expand=c(0,0), name="Individual")

  # If "joint" plot.type, generate one genotype plot with all chromosomes
  if (plot.type=="joint") {
    # If specified, plot vertical lines at site positions
    if (is.null(site.color)==F) {
      genplot <- genplot + geom_vline(xintercept=unique(gen$marker),
                                      colour=site.color)
    }
    genplot <- genplot + facet_grid(. ~ chr)
    return(genplot)
  # If "sep" plot.type, generate list of chromosome-specific genotype plots
    }else{
      genplot.list <- list()
    for (i in 1:length(unique(vcf.table$CHROM))) {
      genchr <- gen[which(gen$chr==unique(vcf.table$CHROM)[i]), ]
      genplot <- ggplot(genchr) +
        geom_rect(aes(xmin=start, xmax=end, ymin=index, ymax=index+1,
                      fill=allele)) +
        scale_fill_manual(values=allele.colors) +
        theme(panel.grid.major=element_blank(),
              panel.grid.minor=element_blank(),
              panel.background=element_blank(),
              axis.line=element_line(colour="black")) +
        scale_x_continuous(expand=c(0,0), name="Position") +
        scale_y_continuous(breaks=unique(genchr$index)+0.5, expand=c(0,0),
                           labels=unique(genchr$sample), name="Individual") +
        ggtitle(unique(vcf.table$CHROM)[i])
      #If specified, plot vertical lines at site positions
      if(is.null(site.color)==F){
        genplot <- genplot + geom_vline(xintercept=unique(gen$marker),
                                        colour=site.color)
      }
      genplot.list[[i]]=genplot
    }
      return(genplot.list)
  }
}



table <- VcfToTable(snakemake@input[[1]])
allele_plot <- PlotVcf(table)


pdf(snakemake@output[[1]],width=6,height=4,paper='special')
PlotVcf(table)
dev.off()

