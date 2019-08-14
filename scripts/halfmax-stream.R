#!/usr/bin/env Rscript
library(tidyverse)
library(furrr)
library(tictoc)
library(tibble)

options(echo=TRUE)
args <- commandArgs(trailingOnly = TRUE)
IN=args[1]
SAMPLES=args[2]
OUT=ifelse(length(args) >= 3, args[3], paste0(basename(IN), ".halfmax.tsv"))

# DEBUGGING
# IN="buxeales_ld_4kb-every-1kb.ld.tsv"
# SAMPLES="src/serlist/buxeales.txt"
# END DEBUG


ldModel = function(ld, n) {
	LD.data = ld$r2
	distance = ld$dist
	HW.st = c(C=0.1)
	fit = NULL
	try(
	    fit <- summary(nls( LD.data~ ((10+C*distance)/((2+C*distance)*(11+C*distance))) *(1+((3+C*distance)*(12+12*C*distance+(C*distance)^2))/(n*(2+C*distance)*(11+C*distance))),
			   start=c(C=0.1), control=nls.control(maxiter=100000, warnOnly=T))),
	    silent=T
	)
	if (is.null(fit)) {
		return(data.frame(rho=NA, halfmax=NA))
	}
	new.rho<-fit$parameters[1]

	dist.pred = (1:100000) / 10
	pred.r2 = ((10+new.rho*dist.pred)/((2+new.rho*dist.pred)*(11+new.rho*dist.pred)))*(1+((3+new.rho*dist.pred)*(12+12*new.rho*dist.pred+(new.rho*dist.pred)^2))/(n*(2+new.rho*dist.pred)*(11+new.rho*dist.pred)))
	halfmax = dist.pred[min(which(pred.r2 < max(pred.r2)/2))]
	data.frame(rho=new.rho, halfmax)
}

samps = read_tsv(SAMPLES, col_names=F)
nsamps = nrow(samps)

procChunk = function (chunk, ...) {
    chunk %>%
		nest(-chrom, -pos, .key="data") %>%
		mutate(data=map(data, ldModel, ...)) %>%
		unnest(data) %>%
		as.data.frame()
}

plan(multisession, gc=TRUE)
fh = gzfile(IN, open="r")
current = NULL
print(paste("Output to:" , OUT))
unlink(OUT) # as we append below
#halfmax = NULL
chunks=1e6
repeat {
	got = tryCatch(read.delim(fh, nrows=chunks, header=is.null(current)),
		       error=function (x) NULL)
	if (!is.null(current) && !is.null(got)) {
		colnames(got) = colnames(current)
	}
	current = rbind(current, got)
	tmp = current
	if (length(unique(current$SNP_A)) > 2) {
		tmp = current %>%
			filter(SNP_A != dplyr::last(SNP_A))
		current = current %>%
			filter(SNP_A == dplyr::last(SNP_A))
	}
    if (nrow(tmp) < 1) {
        next
    }
	tmp = tmp %>%
		transmute(chrom=CHR_A, pos=BP_A, window=floor(pos/5e4), dist=BP_B - BP_A,
				r2=R2, dp=DP) %>%
		nest(-window) %>%
		mutate(data=future_map(data, procChunk, n=nsamps)) %>%
		unnest(data) %>%
        select(-window) %>%
		as.data.frame()
    write_tsv(tmp, OUT, append=file.exists(OUT))
	#halfmax = rbind(halfmax, tmp)
	if (is.null(got) || nrow(current) == 0) {
		break
	}
    print(cbind(current[1,c("CHR_A", "BP_A")], NREC=nrow(current)))
}
#write_tsv(halfmax, OUT)
