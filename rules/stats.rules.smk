#######################################################################
#                            Statistics                               #
#######################################################################


##### Target rules #####

rule stats:
    input:
        "output/plots/stats/quals.svg",
        "output/plots/stats/allele_filtered.tsv",
        "output/plots/stats/allele_plot.pdf",
        #expand("output/variants/lbimpute/{caller}~{aligner}~{ref}~{sampleset}~filtered-{filter}.vcf", caller=config["varcall"]["callers"],aligner=config["varcall"]["aligners"], ref=config["varcall"]["refs"],sampleset=config["varcall"]["samplesets"],filter=config["varcall"]["filters"]),

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
        expand("output/variants/final/{caller}~{aligner}~{ref}~{sampleset}~filtered-{filter}.vcf.gz",
               caller=config["varcall"]["callers"],
               aligner=config["varcall"]["aligners"],
               ref=config["varcall"]["refs"],
               sampleset=config["varcall"]["samplesets"],
               filter=config["varcall"]["filters"]),
    output:
        "output/plots/stats/allele_filtered.tsv"
    shell:
        'printf "#Samples:con-all,D2,D2_F2_tt,D2_F2_TT\nCHROM\\tPOS\\tREF\\tALT\\tRO\\tAO\\tGT\\tGQ\\tSampleRO\\tSampleAO\n" > {output}'


rule lb_impute:
    input:
        "output/variants/final/{caller}~{aligner}~{ref}~{sampleset}~filtered-{filter}.vcf",
        #set = "output/alignments/sets/{aligner}~{ref}~{sampleset}.bam",
        #ref=lambda wc: config['refs'][wc.ref],
    output:
        "output/variants/lbimpute/{caller}~{aligner}~{ref}~{sampleset}~filtered-{filter}.vcf"
    log:
        "output/log/lbimpute/{caller}~{aligner}~{ref}~{sampleset}~filtered-{filter}.log"
    params:
        mem = config['LB-impute']['mem'],
        release = config['LB-impute']['release'],
        method = config['LB-impute']['method'],
        extra = config['LB-impute']['extra'],
	parents = config['LB-impute']['parents'],
    shell:
        "( java"
        "   -{params.mem}"
        "   -jar {params.release}"
        "   -method {params.method}"
        "   -f {input}"
        "   -parents {params.parents}"
        "   {params.extra}"
        "   -o {output}"
        ") >{log} 2>&1"


rule allele_plot:
    input:
        expand("output/variants/final/{caller}~{aligner}~{ref}~{sampleset}~filtered-{filter}.vcf",
               caller=config["varcall"]["callers"],
               aligner=config["varcall"]["aligners"],
               ref=config["varcall"]["refs"],
               sampleset=config["varcall"]["samplesets"],
               filter=config["varcall"]["filters"]),
    output:
        "output/plots/stats/allele_plot.pdf",
    script:
        "../scripts/allele_plot.R"
