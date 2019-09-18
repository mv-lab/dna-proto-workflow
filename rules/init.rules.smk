#######################################################################
#                            Init rules                               #
#######################################################################

from utils import snkmk
import os

configfile: "configs/toolconfig.yml"
configfile: "configs/runconfig.yml"
report: "../report/workflow.rst"

snkmk.create_fai()
snkmk.create_contigs_file()

#RUNLIB2SAMP, SAMP2RUNLIB = snkmk.make_runlib2samp("metadata/sample2runlib.csv")
#SAMPLESETS = snkmk.make_samplesets(s2rl_file="metadata/sample2runlib.csv", setfile_glob="metadata/samplesets/*.txt")
#VARCALL_REGIONS = snkmk.make_regions(config["refs"], window=config["varcall"]["chunksize"])

rule prepare_ref:
    input:
        "rawdata/reference/genome.fa"
    output:
        "rawdata/reference/genome.fa.fai"
    log:
        "data/log/reference/fai.txt"
    shell:
        "samtools faidx {input} 2> {log}"

rule contigs:
    input:
        "rawdata/reference/genome.fa.fai"
    output:
        "metadata/contigs_of_interest.bed"
    log:
        "data/log/contigs/contigs.txt"
    shell:
        "awk r'BEGIN {FS='\t\t'}; {print $1 FS '0' FS $2}' {input} > {output} 2> {log}"

shell.prefix = "set -euo pipefail; "

wildcard_constraints:
    run="[^/]+",
    lib="[^/]+",
    aligner="[^/]+",
    sample="[^/]+",
    ref="[^/]+",
    type="[^/]+",


rule init:
    input:
        rules.contigs.output,
        rules.prepare_ref.output,
