
rule snpeff:
    input:
        "data/variants/final/freebayes~bwa~genome~all_samples~filtered-default.vcf.gz",
    output:
        vcf="annotated/all.vcf.gz",
        csvstats="snpeff/all.csv"
    log:
        "data/log/snpeff/snpeff.log"
    params:
        reference=config["snpeff"]['name'],
        extra="-Xmx6g"
    wrapper:
        "0.27.1/bio/snpeff"
