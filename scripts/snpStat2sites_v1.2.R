#!/usr/bin/env Rscript
library(dplyr, quietly=T, warn.conflicts=F)
library(tidyr, quietly=T, warn.conflicts=F)

## This version is faster and more generic
## Relies on directory structure including "split" folder within dataset folder.
## Outputs to a single directory

convertHWE <- function(hweFile, maxHo=0.8, minHWEPval=0) {
    read.delim(hweFile) %>%
        mutate(He = 2*Freq*(1-Freq),
               Ho = He-(He * F)) %>%
        # filter(Ho <= maxHo, p.value >= minHWEPval) %>%
        select(-LRT, -Major, -Minor)
}

convertHo <- function(speciesFile) {
  read.delim(speciesFile) %>%
    mutate(He = 2*Freq*(1-Freq),
           Ho = He-(He * F)) %>%
    # filter(Ho <= maxHo, p.value >= minHWEPval) %>%
    select(Chromo,Position,Ho)
}

filterSnpStat <- function(snpstatFile, minMAD=2) {
    read.delim(snpstatFile, skip=1, col.names=c("Chromo", "Position",
               "majminor", "SB", "HWE", "baseQ", "mapQ", "edgeZ", "edgeZ2")) %>%
        separate(majminor, into=c("pMaj", "pMin", "mMaj", "mMin"),
                 sep=' ', convert=TRUE) %>%
        mutate(countMajor = pMaj + mMaj,
               countMinor = pMin + mMin) %>%
        select(Chromo, Position, countMajor, countMinor) %>%
        filter(countMinor >= minMAD,
               countMajor+countMinor < quantile(countMajor+countMinor,0.99))
}


main <- function(species_list, outbase, hweFile,snpstatFile,chr) {
  
    hwe = convertHWE(hweFile)
    snpstat = filterSnpStat(snpstatFile)
    merged = inner_join(hwe, snpstat, by=c("Chromo", "Position"))

    files=paste0(file.path(datadir,species,"split",paste0(chr,".hwe.gz")))
    for(i in 1:length(species)){
      speciesHo=convertHo(files[i])
      names(speciesHo)[3]=paste0("Ho_",species[i])
      merged = left_join(merged, speciesHo, by=c("Chromo", "Position"))
      write.table(merged, paste0(outbase, ".sitestats0.tsv"), sep="\t",
                  row.names=F, quote=F)
    }
    
    exclude_ho=(rowSums(merged[,grepl("Ho_",names(merged))]>maxHo, na.rm = T)/
      rowSums(!is.na(merged[,grepl("Ho_",names(merged))]), na.rm = T))>maxPopFail
    merged=merged[-which(exclude_ho),]

    write.table(merged[,c(1, 2)], paste0(outbase, ".sites"), quote=F, sep="\t",
                row.names=F, col.names=F)
    write.table(merged, paste0(outbase, ".sitestats.tsv"), sep="\t",
                row.names=F, quote = F)
}


args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
    cat("Error! usage is snpStat2sites.R <OUTPUT_DIR> <DATA_DIR> <SPECIES_LIST> <CHR>\n")
    q(save="no", status=1)
}

# args=c("output","data","data/species.list","Chr01")
outdir = args[1]
datadir = args[2]
species_list = args[3]
chr=args[4]

outbase=file.path(outdir,chr)
species=scan(species_list, what="char")
hweFile=file.path(datadir,paste0("everything/split/",chr,".hwe.gz"))
snpstatFile=file.path(datadir,paste0("everything/split/",chr,".snpStat.gz"))
maxHo=0.8
maxPopFail=0.8

main(species_list, outbase,hweFile,snpstatFile,chr)

