#######################################################################
#                            Annotation                               #
#######################################################################

rule snpeff:
    input:
        expand("output/variants/final/freebayes~{aligner}~{ref}~{sampleset}~filtered-{filter}.vcf.gz",
               aligner=config["varcall"]["aligners"],
               ref=config["varcall"]["refs"],
               sampleset=config["varcall"]["samplesets"],
               filter = config['snpeff']['filter']),
    output:
        vcf="output/snpeff/annotated/all.vcf.gz",
        csvstats="output/snpeff/annotated/all.csv"
    log:
        "output/log/snpeff/snpeff.log"
    params:
        reference=config["snpeff"]['name'],
        extra="-Xmx6g -v"
    wrapper:
        "0.27.1/bio/snpeff"


rule prepare_database:
    input:
        name= config['snpeff']['name'],
    output:
        dir(expand("genomes_and_annotations/snpeffdata/{name}", name= config['snpeff']['name']))
    shell:
        "mkdir genomes_and_annotations/snpeffdata/{input.name}"
