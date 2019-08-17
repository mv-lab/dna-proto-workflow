#######################################################################
#                           Variant Calling                           #
#######################################################################

# So I've rejigged the variant calling to use the megabam of all samples for
# each aligner/ref, as then we don't need 15 different massive subset bams.

rule freebayes:
    input:
        bam="data/alignments/sets/{aligner}~{ref}~all_samples.bam",  # use the megabam, see above
        bai="data/alignments/sets/{aligner}~{ref}~all_samples.bam.bai",
        sset="data/samplelists/{sampleset}.txt",
        ref=lambda wc: config['refs'][wc.ref],
    output:
        bcf="data/variants/raw_split/freebayes~{aligner}~{ref}~{sampleset}/{region}.bcf",
    log:
        "data/log/freebayes/{aligner}~{ref}~{sampleset}/{region}.log"
    benchmark:
        "data/log/freebayes/{aligner}~{ref}~{sampleset}/{region}.benchmark"
    priority: 1  # get them done earlier, normalisation is super quick
    params:
        theta=config["varcall"].get("theta_prior", 0.01),
        minmq=lambda wc: config["varcall"]["minmapq"].get(wc.aligner, 5),
        minbq=config["varcall"]["minbq"],
    shell:
        "(  freebayes"
        "   --theta {params.theta}"
        "   --samples {input.sset}"
        "   --ploidy 2"
        "   --use-best-n-alleles 3"
        "   --min-mapping-quality {params.minmq}"
        "   --min-base-quality {params.minbq}"
        "   --read-max-mismatch-fraction 0.1"
        "   --min-alternate-fraction 0"
        "   --min-alternate-count 2" # per sample
        "   --min-alternate-total 5" # across all samples
        "   --min-coverage 10" # across all samples
        "   --prob-contamination 1e-6"
        "   --use-mapping-quality"
        "   --strict-vcf"
        "   --region '{wildcards.region}'"
        "   --fasta-reference {input.ref}"
        "   {input.bam}"
        " | bcftools view"
        "   -O b  -o {output.bcf}"
        " ) >{log} 2>&1"


rule mpileup:
    input:
        bam="data/alignments/sets/{aligner}~{ref}~all_samples.bam",  # use the megabam, see above
        bai="data/alignments/sets/{aligner}~{ref}~all_samples.bam.bai",
        sset="data/samplelists/{sampleset}.txt",
        ref=lambda wc: config['refs'][wc.ref],
    output:
        bcf="data/variants/raw_split/mpileup~{aligner}~{ref}~{sampleset}/{region}.bcf",
    log:
        "data/log/mpileup/{aligner}~{ref}~{sampleset}/{region}.log"
    benchmark:
        "data/log/mpileup/{aligner}~{ref}~{sampleset}/{region}.benchmark"
    params:
        theta=config["varcall"].get("theta_prior", 0.01),
        minmq=lambda wc: config["varcall"]["minmapq"].get(wc.aligner, 5),
        minbq=config["varcall"]["minbq"],
    priority: 1  # get them done earlier, normalisation is super quick
    shell:
        "( bcftools mpileup"
        "   --adjust-MQ 50"
        "   --redo-BAQ"
        "   --max-depth 20000" # the default per file max (250x) is insane, i.e. <1x for most sets. new limit of 20000x  equates to a max. of 20x across all samples.
        "   --min-MQ {params.minmq}"
        "   --min-BQ {params.minbq}"
        "   --fasta-ref {input.ref}"
        "   --samples-file {input.sset}"
        "   --annotate FORMAT/DP,FORMAT/AD,FORMAT/SP,INFO/AD" #output extra tags
        "   --region '{wildcards.region}'"
        "   --output-type u" #uncompressed
        "   {input.bam}"
        " | bcftools call"
        "   --targets '{wildcards.region}'" # might not be needed
        "   --multiallelic-caller"
        "   --prior {params.theta}"
        "   -O b"
        "   -o {output.bcf}"
        " ) >{log} 2>&1"


rule bcfnorm:
    input:
        bcf="data/variants/raw_split/{caller}~{aligner}~{ref}~{sampleset}/{region}.bcf",
        ref=lambda wc: config['refs'][wc.ref],
    output:
        # Not a pipe! can't run multiple filters if a pipe
        bcf=temp("data/variants/norm_split/{caller}~{aligner}~{ref}~{sampleset}/{region}.bcf"),
    log:
        "data/log/bcfnormalise/{caller}~{aligner}~{ref}~{sampleset}/{region}.log"
    shell:
        "( bcftools norm"
        "   --fasta-ref {input.ref}"
        "   -O u"
        "   {input.bcf}"
        " | vt decompose_blocksub + -o -" # decompose MNP to multipe SNPs
        " | bcftools norm" # Split multi-alleics
        "   --fasta-ref {input.ref}"
        "   --do-not-normalize"
        "   --multiallelics -snps"
        "   -O u  -o {output.bcf}"
        " ) >{log} 2>&1"

rule bcffilter:
    input:
        bcf="data/variants/norm_split/{caller}~{aligner}~{ref}~{sampleset}/{region}.bcf",
        ref=lambda wc: config['refs'][wc.ref],
    output:
        # Not a pipe! can't run all regions separately if this is a pipe into merge
        bcf=temp("data/variants/filter_split/{caller}~{aligner}~{ref}~{sampleset}_filtered~{filter}/{region}.bcf"),
    log:
        "data/log/bcffilter/{caller}~{aligner}~{ref}~{sampleset}/{filter}/{region}.log"
    params:
        filtarg=lambda wc: config["varcall"]["filters"][wc.filter].replace('\n', ' ')
    shell:
        "( bcftools view"
        "   {params.filtarg}"
        "   -O u"
        "   {input.bcf}"
        " | bcftools norm" # We normalise here to re-join multi-allelic sites, after filtering with multi-allelics split
        "   --fasta-ref {input.ref}"
        "   --do-not-normalize"
        "   --multiallelics +snps" # Split multi-alleic sites
        "   -O b  -o {output.bcf}"
        " ) >{log} 2>&1"

localrules: bcfmerge_fofn
rule bcfmerge_fofn:
    input:
        bcf=lambda wc: expand("data/variants/filter_split/{caller}~{aligner}~{ref}~{sampleset}_filtered~{filter}/{region}.bcf",
                              caller=wc.caller, aligner=wc.aligner, ref=wc.ref, sampleset=wc.sampleset, filter=wc.filter,
                              region=sorted(VARCALL_REGIONS[wc.ref])),
    output:
        fofn=temp("data/variants/final/{caller}~{aligner}~{ref}~{sampleset}~filtered-{filter}.bcf.INPUT_FOFN"),
    run:
        with open(output[0], "w") as fh:
            for s in sorted(input):
                print(s, file=fh)

rule bcfmerge:
    input:
        bcf=lambda wc: expand("data/variants/filter_split/{caller}~{aligner}~{ref}~{sampleset}_filtered~{filter}/{region}.bcf",
                              caller=wc.caller, aligner=wc.aligner, ref=wc.ref, sampleset=wc.sampleset, filter=wc.filter,
                              region=sorted(VARCALL_REGIONS[wc.ref])),
        fofn="data/variants/final/{caller}~{aligner}~{ref}~{sampleset}~filtered-{filter}.bcf.INPUT_FOFN",
    output:
        bcf="data/variants/final/{caller}~{aligner}~{ref}~{sampleset}~filtered-{filter}.bcf",
    log:
        "data/log/mergebcf/{caller}~{aligner}~{ref}~{sampleset}_filtered~{filter}.log"
    threads: 4
    shell:
        "( bcftools concat"
        "   --threads {threads}"
        "   -O b"
        "   -o {output.bcf}"
        "   --file-list {input.fofn}"
        " ) >{log} 2>&1"


rule bcf2vcf:
    input:
        bcf="{path}.bcf",
    output:
        vcf="{path}.vcf.gz",
    log:
        "data/log/bcf2vcf/{path}.log"
    threads: 4
    shell:
        "( bcftools view"
        "   {input.bcf}"
        "   -O z"
        "   --threads {threads}"
        "   -o {output.vcf}"
        " ) >{log} 2>&1"

rule variantidx:
    input:
        "{path}"
    output:
        "{path}.csi"
    shell:
        "bcftools index -f {input}"

rule stats:
    input:
        "data/variants/{path}"
    output:
        "data/stats/variants/{path}.varstats"
    shell:
        "bcftools stats -s - -d 0,1000,2 --threads {threads} {input} >{output}"


def raw_variant_calls_input(wildcards):
    inputs = []
    for caller in config["varcall"]["callers"]:
        for aligner in config["varcall"]["aligners"]:
            for sampleset in config["varcall"]["samplesets"]:
                for ref in config["varcall"]["refs"]:
                    this_rawfiles = expand("data/variants/raw_split/{caller}~{aligner}~{ref}~{sampleset}/{region}.bcf",
                                           caller=caller, aligner=aligner, ref=ref, sampleset=sampleset, region=VARCALL_REGIONS[ref])
                    inputs.extend(this_rawfiles)
    return inputs


rule raw_variant_calls:
    input: raw_variant_calls_input

rule filtered_variants:
    input:
        expand("data/variants/final/{caller}~{aligner}~{ref}~{sampleset}~filtered-{filter}.{ext}",
               ext=["bcf", "bcf.csi", "vcf.gz", "vcf.gz.csi"],
               caller=config["varcall"]["callers"],
               aligner=config["varcall"]["aligners"],
               ref=config["varcall"]["refs"],
               sampleset=config["varcall"]["samplesets"],
               filter=config["varcall"]["filters"]),

rule varcall:
    input:
        rules.filtered_variants.input,