#!/bin/bash -login

#generate the mafs and safs needed for downstream analysis
#files for use in NGStools analysis
#all bam files have been sorted with Samtools first

#PBS -N species_saf
#PBS -P xe2
#PBS -q normal
#PBS -l ncpus=2,walltime=150:00:00,mem=80GB
#PBS -l wd
#PBS -M randre20@une.edu.au
#PBS -m abe

#NOT#Running angsd to create genotype likelihood files (beagle + 
#NOT#print actual genotypes)

#Data (sorted and merged BAM files) stored /media/disk1/
#Reference stored /g/data1/xe2/references/eucalyptus/grandis - not bgz?
#Angsd stored ???

#NOT#setting current directory as PBS working directory 
#NOT#(default is /home/583/jxj583)
#NOT#cd $PBS_O_WORKDIR **NOTE: PBS -l wd above does the same thing

#setting variables
SPECIES=albens
CHR=Chr01
REF=/g/data1/xe2/references/eucalyptus/grandis/Egrandis_v2.fasta
ANC=/g/data1/xe2/references/eucalyptus/grandis/Egrandis_v2.fasta

BAMLIST=/g/data1/xe2/projects/euc_hybrid_adapt/workspace/olddata/step2/bamlists/$SPECIES.bamfile
SITES=/g/data1/xe2/projects/euc_hybrid_adapt/workspace/olddata/step2/sites/$CHR.sites
OUT=/g/data1/xe2/projects/euc_hybrid_adapt/workspace/olddata/step2/split/$SPECIES.$CHR

#NOT#create bam list of all bam files in directory
#NOT#ls *.bam > bam2.list

#doMajorMinor 1 = use genotype likelihood scores, doGL 2 = use GATK method
#doMaf 2 = fixed major, unknown minor; doGeno 36 = print genotypes as binary
#and print genotypes as direct alleles(32 does binary; 4 does full alleles);
#SNP pval increases speed for ngscovar
#rather than using all sites; doPost 1 means use HWE based prior
#-doGlf = generate beagle files
#-dosaf = generate saf files

#angsd -P 2 -bam bam2.list -out All_noreps_nobad_binary -ref $REF -anc $ANC -GL 2 -doMajorMinor 1 -doMaf 2 -doGeno 32 -dosaf 1 -skipTriallelic 1 -C 50 -baq 1 -minMapQ 30 -minQ 20 -SNP_pval 1e-3 -doPost 1
angsd -P 2 -bam $BAMLIST -out $OUT -ref $REF -anc $ANC -GL 2 -doMajorMinor 5 -doMaf 1 -doSaf 1 -C 50 -baq 1 -minMapQ 30 -minQ 20 -sites $SITES

