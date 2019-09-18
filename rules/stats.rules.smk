#######################################################################
#                            Statistics                               #
#######################################################################


##### Target rules #####

rule stats:
    input:
        "plots/quals.svg",

##### Actual rules #####

rule plot_quals:
    input:
        "data/variants/final/freebayes~bwa~genome~all_samples~filtered-strict.vcf.gz"
    output:
        "plots/quals.svg"
    script:
        "scripts/plot-quals.py"
