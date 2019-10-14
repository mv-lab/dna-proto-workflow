#######################################################################
#                            Init rules                               #
#######################################################################

from utils import snkmk
import os
import pandas as pd

configfile: "config.yml"

shell.prefix = "set -euo pipefail;"

wildcard_constraints:
    run="[^/]+",
    lib="[^/]+",
    aligner="[^/]+",
    sample="[^/]+",
    ref="[^/]+",
    type="[^/]+",

snkmk.create_fai()
snkmk.create_contigs_file()

RUNLIB2SAMP, SAMP2RUNLIB = snkmk.make_runlib2samp(config['samples'])
SAMPLESETS = snkmk.make_samplesets(s2rl_file=config['samples'], setfile_glob="metadata/samplesets/*.txt")
VARCALL_REGIONS = snkmk.make_regions(config["refs"], window=config["varcall"]["chunksize"])

units = pd.read_csv(config['samples'], dtype=str).set_index(["run", "library"], drop=False)
units.index = units.index.set_levels([i.astype(str) for i in units.index.levels])  # enforce str in index
#print (units)


def get_il_fastq(wildcards):
    """Get fastq files of given metadata."""
    fastqs = units.loc[(wildcards.run, wildcards.lib), ["il_fq"]].dropna()
    return {"r1": fastqs.il_fq}


def get_fr_fastq(wildcards):
    """Get fastq files of given metadata."""
    fastqs = units.loc[(wildcards.run, wildcards.lib), ["fq1", "fq2"]].dropna()
    return {"r1": fastqs.fq1, "r2": fastqs.fq2}


rule align_prepare_reference:
    input:
        ref = lambda wc: config['refs'][wc.ref]
    output:
        "genomes_and_annotations/{ref}/{ref}.fa.amb",
        "genomes_and_annotations/{ref}/{ref}.fa.ann",
        "genomes_and_annotations/{ref}/{ref}.fa.bwt",
        "genomes_and_annotations/{ref}/{ref}.fa.pac",
        "genomes_and_annotations/{ref}/{ref}.fa.sa"
    shell:
        "bwa index -a bwtsw {input.ref}"
