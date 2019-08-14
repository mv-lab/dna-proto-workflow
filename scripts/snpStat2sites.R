#!/usr/bin/env Rscript
library(dplyr, quietly=T, warn.conflicts=F)
library(tidyr, quietly=T, warn.conflicts=F)


filterHWE <- function(hweFile, maxHo=0.8, minHWEPval=0) {
    read.delim(hweFile) %>%
        mutate(He = 2*Freq*(1-Freq),
               Ho = He-(He * F)) %>%
        filter(Ho <= maxHo, p.value >= minHWEPval) %>%
        select(-LRT, -Major, -Minor)
}

filterSnpStat <- function(snpstatFile, minMAD=2) {
    read.delim(snpstatFile, skip=1, col.names=c("Chromo", "Position",
               "majminor", "SB", "HWE", "baseQ", "mapQ", "edgeZ", "edgeZ2")) %>%
        separate(majminor, into=c("pMaj", "pMin", "mMaj", "mMin"),
                 sep=' ', convert=TRUE) %>%
        mutate(countMajor = pMaj + mMaj,
               countMinor = pMin + mMin) %>%
        select(Chromo, Position, countMajor, countMinor) %>%
        filter(countMinor >= minMAD)
}


main <- function(basename, outbase) {
    hwe = filterHWE(paste0(basename, ".hwe.gz"))
    snpstat = filterSnpStat(paste0(basename, ".snpStat.gz"))
    merged = inner_join(hwe, snpstat, by=c("Chromo", "Position"))

    write.table(merged[,c(1, 2)], paste0(outbase, ".sites"), quote=F, sep="\t",
                row.names=F, col.names=F)
    write.table(merged, paste0(outbase, ".sitestats.tsv"), sep="\t",
                row.names=F)
}


args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
    cat("Error! usage is snpStat2sites.R <OUTPUT_BASENAME> <INPUT_BASENAME>\n")
    q(save="no", status=1)
}
outbase = args[1]
basename = args[2]
main(basename, outbase)
