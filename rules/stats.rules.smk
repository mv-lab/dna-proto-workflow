
#######################################################################
#                            Statistics                               #
#######################################################################


##### Target rules #####

rule stats:
    input:
        "plots/quals.svg",
        "tables/calls.tsv.gz",  # rules.vcf_to_tsv.output,
        "plots/depths.svg",     # rules.plot_stats.output,
        "plots/allele-freqs.svg",


##### Actual rules #####

rule vcf_to_tsv:
    input:
        "data/variants/final/freebayes~bwa~genome~all_samples~filtered-strict.vcf.gz"
    output:
        report("tables/calls.tsv.gz", caption="../report/calls.rst", category="Calls")
    shell:
        "bcftools view --apply-filters PASS --output-type u {input} | "
        "rbt vcf-to-txt -g --fmt DP AD --info ANN | "
        "gzip > {output}"

rule plot_stats:
    input:
        "tables/calls.tsv.gz"
    output:
        depths=report("plots/depths.svg", caption="../report/depths.rst", category="Plots"),
        freqs=report("plots/allele-freqs.svg", caption="../report/freqs.rst", category="Plots")
    script:
        "../scripts/plot-depths.py"

rule plot_quals:
    input:
        "data/variants/final/freebayes~bwa~genome~all_samples~filtered-strict.vcf.gz"
    output:
        "plots/quals.svg"
    script:
        "scripts/plot-quals.py"
