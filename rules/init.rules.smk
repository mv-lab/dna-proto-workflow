#######################################################################
#                            Init rules                               #
#######################################################################

from utils import snkmk
import os

configfile: "configs/toolconfig.yml"
configfile: "configs/runconfig.yml"
report: "../report/workflow.rst"

shell.prefix = "set -euo pipefail; "

wildcard_constraints:
    run="[^/]+",
    lib="[^/]+",
    aligner="[^/]+",
    sample="[^/]+",
    ref="[^/]+",
    type="[^/]+",

snkmk.create_fai()
snkmk.create_contigs_file()

#RUNLIB2SAMP, SAMP2RUNLIB = snkmk.make_runlib2samp("metadata/sample2runlib.csv")
#SAMPLESETS = snkmk.make_samplesets(s2rl_file="metadata/sample2runlib.csv", setfile_glob="metadata/samplesets/*.txt")
#VARCALL_REGIONS = snkmk.make_regions(config["refs"], window=config["varcall"]["chunksize"])

rule prepare_ref:
    input:
        "genomes_and_annotations/genomes/{ref}/genome.fa",
        name = lambda wc: config['init']['refs'][wc.ref],
        ref = lambda wc: config['refs'][name],
    output:
        "genomes_and_annotations/genomes/{ref}/genome.fa.fai"
    log:
        "log/init/create_fai.log"
    shell:
        "samtools faidx {input} 2> {log}"

rule contigs:
    input:
        "genomes_and_annotations/genomes/{ref}/genome.fa.fai",
        ref = lambda wc: config['init']['refs'][wc.ref],
    output:
        "metadata/contigs_of_interest.bed"
    log:
        "log/init/create_contigs.log"
    shell:
        "awk r'BEGIN {FS='\t\t'}; {print $1 FS '0' FS $2}' {input} > {output} 2> {log}"


rule align_prepare_reference:
    input:
        "genomes_and_annotations/genomes/Sorghum/genome.fa"
    output:
        expand("genomes_and_annotations/genomes/Sorghum/genome.fa.{ext}", ext=["amb", "ann", "bwt", "pac", "sa"])
    shell:
        "bwa index {input}"

rule init:
    input:
        #rules.contigs.output,
        #rules.prepare_ref.output,
