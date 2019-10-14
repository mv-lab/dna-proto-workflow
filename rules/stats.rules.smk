#######################################################################
#                            Statistics                               #
#######################################################################


##### Target rules #####

rule stats:
    input:
        "output/plots/stats/quals.svg",
        "output/plots/stats/allele_filtered.tsv",

##### Actual rules #####

rule plot_quals:
    input:
        expand("output/variants/final/{caller}~{aligner}~{ref}~{sampleset}~filtered-{filter}.vcf.gz",
               caller=config["varcall"]["callers"],
               aligner=config["varcall"]["aligners"],
               ref=config["varcall"]["refs"],
               sampleset=config["varcall"]["samplesets"],
               filter=config["varcall"]["filters"]),
    output:
        "output/plots/stats/quals.svg"
    script:
        "../scripts/plot-quals.py"


rule tables_from_vcf:
    input:
        expand("output/variants/final/freebayes~{aligner}~{ref}~{sampleset}~filtered-{filter}.vcf.gz",
               caller=config["varcall"]["callers"],
               aligner=config["varcall"]["aligners"],
               ref=config["varcall"]["refs"],
               sampleset=config["varcall"]["samplesets"],
               filter=config["varcall"]["filters"]),
    output:
        "output/plots/stats/allele_filtered.tsv"
    shell:
        'printf "#Samples:con-all,D2,D2_F2_tt,D2_F2_TT\nCHROM\\tPOS\\tREF\\tALT\\tRO\\tAO\\tGT\\tGQ\\tSampleRO\\tSampleAO\n" > {output}'
