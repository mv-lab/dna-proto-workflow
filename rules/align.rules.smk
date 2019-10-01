
#######################################################################
#                       Alignment to Reference - Mapping              #
#######################################################################


##### Target rules #####

localrules: align_samples
rule align_samples:
    input:
        expand("data/alignments/samples/{aligner}/{ref}/{sample}.bam",
               ref=config["mapping"]["refs"],
               aligner=config["mapping"]["aligners"],
               sample=SAMP2RUNLIB),


localrules: align_stats
rule align_stats:
    input:
        expand("data/alignments/bamstats/sample/{aligner}~{ref}~{sample}.tsv",
               aligner=config["mapping"]["aligners"],
               ref=config["mapping"]["refs"],
               sample=SAMPLESETS["all_samples"]),
    output:
        expand("data/alnstats/everything_{type}.csv",
               type=["SN", "IS", "COV"])
    log: "log/align/bamstats/mergeallbamstats.log"
    shell:
        "python3 scripts/tidybamstat.py"
        "   -o data/alnstats/everything"  # prefix
        "   {input}"
        ">{log} 2>&1"


rule align:
   input:
        rules.align_samples.input,
        rules.align_stats.output,
        expand("data/alignments/sets/{aligner}~{ref}~all_samples.bam",
               ref=config["mapping"]["refs"],
               aligner=config["mapping"]["aligners"]),
        expand("data/alignments/bamlists/{aligner}~{ref}~all_samples.bamlist",
               ref=config["mapping"]["refs"],
               aligner=config["mapping"]["aligners"]),
        expand("data/alignments/sets/{aligner}~{ref}~all_samples.bam.bai",
               ref=config["mapping"]["refs"],
               aligner=config["mapping"]["aligners"]),
        expand("data/alignments/qualimap/samples/{aligner}~{ref}~{sample}/",
               aligner=config["mapping"]["aligners"],
               ref=config["mapping"]["refs"],
               sample=SAMPLESETS["all_samples"]),


##### Actual rules #####

rule ngmap:
    input:
        reads="data/reads/runs/{run}/{lib}.fastq.gz",
        ref=lambda wc: config['refs'][wc.ref],
    output:
        bam=temp("data/alignments/byrun.raw/ngm/{ref}/{run}/{lib}.bam"),
    log:
        "log/align/ngm/{ref}/{run}/{lib}.log"
    threads:
        8
    params:
        sample=lambda wc: RUNLIB2SAMP.get((wc.run, wc.lib), "{}~{}".format(wc.run, wc.lib)),
        sensitivity=config["mapping"]["ngm"]["sensitivity"],
    shell:
        "( ngm"
        "   -q {input.reads}"
        "   --paired --broken-pairs"
        "   -r {input.ref}"
        "   -t {threads}"
        "   --rg-id {wildcards.run}_{wildcards.lib}"
        "   --rg-sm {params.sample}"
        "   --sensitivity {params.sensitivity}" # this is the mean from a bunch of different runs
        "| samtools view -Suh - >{output.bam}"
        " ) >{log} 2>&1"

rule bwamem:
    input:
        expand("genomes_and_annotations/genomes/Sorghum/genome.fa.{ext}", ext=["amb", "ann", "bwt", "pac", "sa"]),
        reads="data/reads/runs/{run}/{lib}.fastq.gz",
        ref=lambda wc: config['refs'][wc.ref],
    output:
        bam=temp("data/alignments/byrun.raw/bwa/{ref}/{run}/{lib}.bam"),
    log: "log/align/bwa/{ref}/{run}/{lib}.log"
    threads:
        8
    params:
        sample=lambda wc: RUNLIB2SAMP.get((wc.run, wc.lib), "{}~{}".format(wc.run, wc.lib)),
    shell:
        "( bwa mem"
        "   -p" # paired input
        "   -t {threads}"
        "   -R '@RG\\tID:{wildcards.run}_{wildcards.lib}\\tSM:{params.sample}'"
        "   {input.ref}"
        "   {input.reads}"
        "| samtools view -Suh - >{output.bam}"
        " ) >{log} 2>&1"


rule bam_markdups_sort:
    input:
        bam="data/alignments/byrun.raw/{aligner}/{ref}/{run}/{lib}.bam",
        ref=lambda wc: config['refs'][wc.ref],
    output:
        bam=temp("data/alignments/byrun/{aligner}/{ref}/{run}/{lib}.bam"),
    threads: 4
    log: "log/align/markdup/{aligner}/{ref}/{run}/{lib}.log"
    shell:
        "( samtools fixmate "
        "   -m"
        "   -@ {threads}"
        "   --output-fmt bam,level=0"
        "   {input.bam}"
        "   -"
        " | samtools sort"
        "   -T ${{TMPDIR:-/tmp}}/{wildcards.run}_{wildcards.lib}_sort_$RANDOM"
        "   --output-fmt bam,level=0"
        "   -@ {threads}"
        "   -m 1g"
        "   -"
        " | samtools markdup"
        "   -T ${{TMPDIR:-/tmp}}/{wildcards.run}_{wildcards.lib}_markdup_$RANDOM"
        "   -s" # report stats
        "   -@ {threads}"
        "   --output-fmt bam,level=3"
        "   -"
        "   {output.bam}"
        " ) >{log} 2>&1"


rule mergebam_samp:
    input:
        lambda wc: ["data/alignments/byrun/{aln}/{ref}/{run}/{lib}.bam".format(
                            run=r, lib=l, aln=wc.aligner, ref=wc.ref)
	                for r, l in SAMP2RUNLIB[wc.sample]]
    output:
        bam="data/alignments/samples/{aligner}/{ref}/{sample}.bam",
    log:
        "log/align/mergesamplebam/{aligner}/{ref}/{sample}.log"
    threads: 8
    priority: 1 # so the temps get cleaned sooner
    shell:
        "( samtools merge"
        "   -@ {threads}"
        "   --output-fmt bam,level=4"
        "   {output.bam}"
        "   {input}"
        " ) >{log} 2>&1"


rule qualimap_samp:
    input:
        bam="data/alignments/samples/{aligner}/{ref}/{sample}.bam",
    output:
        "data/alignments/qualimap/samples/{aligner}~{ref}~{sample}/",
    log:
        "log/align/qualimap_sample/{aligner}~{ref}~{sample}.log"
    threads: 4
    shell:
        "( unset DISPLAY; qualimap bamqc"
        "   --java-mem-size=4G"
        "   -bam {input.bam}"
        "   -nr 10000"
        "   -nt {threads}"
        "   -outdir {output}"
        "   {input}"
        " ) >{log} 2>&1"


#localrules: bamlist
rule bamlist:
    input:
        lambda wc: expand("data/alignments/samples/{aligner}/{ref}/{sample}.bam",
                          aligner=wc.aligner, ref=wc.ref, sample=SAMPLESETS[wc.sampleset]),

    output:
        "data/alignments/bamlists/{aligner}~{ref}~{sampleset}.bamlist",
    run:
        with open(output[0], "w") as fh:
            for s in input:
                print(s, file=fh)


rule mergebam_set:
    input:
        lambda wc: expand("data/alignments/samples/{aligner}/{ref}/{sample}.bam",
                          aligner=wc.aligner, ref=wc.ref, sample=SAMPLESETS[wc.sampleset]),

    output:
        bam="data/alignments/sets/{aligner}~{ref}~{sampleset}.bam",
        bai="data/alignments/sets/{aligner}~{ref}~{sampleset}.bam.bai",
    log:
        "log/align/mergesetbam/{aligner}/{ref}/{sampleset}.log"
    threads: 4
    shell:
        "( samtools merge"
        "   --output-fmt bam,level=7"
        "   -@ {threads}"
        "   -"
        "   {input}"
        " | tee {output.bam}"
        " | samtools index - {output.bai}"  # indexing takes bloody ages, we may as well do this on the fly
        " ) >{log} 2>&1"


localrules: bamstat_samps
rule bamstat_samps:
    input:
        "data/alignments/samples/{aligner}/{ref}/{sample}.bam",
    output:
        "data/alignments/bamstats/sample/{aligner}~{ref}~{sample}.tsv",
    log:
        "log/align/bamstats_sample/{aligner}~{ref}~{sample}.tsv"
    shell:
        "(samtools stats -i 5000 -x {input} >{output}) >{log}"



# Utils
ruleorder: mergebam_set > bamidx # as we index on the fly for sets
rule bamidx:
    input:
        "{path}.bam"
    output:
        "{path}.bam.bai"
    log:
        "log/align/bamindex/{path}.log"
    shell:
        "samtools index {input}"

rule sam2bam:
    input:
        "{path}.sam"
    output:
        "{path}.bam"
    shell:
        "samtools view -Suh {input} >{output}"


### Align targets


allsets = set(
    list(config["mapping"]["samplesets"]) +
    list(config["varcall"]["samplesets"]) +
    list(config["sample_sets"])
)

rule align_librun:
    input:
        lambda wc: ["data/alignments/byrun/{aln}/{ref}/{run}/{lib}.bam".
                        format(run=r, lib=l, aln=a, ref=ref)
                        for r, l in RUNLIB2SAMP
                        for a in config["mapping"]["aligners"]
                        for ref in config["mapping"]["refs"]],


rule align_samplesets_all:
    # Currently there's no real need for this, as all variant calling uses the mega-bam
    # and ANGSD uses a list of filenames (which will be generated on the fly
    # for each ANGSD sample set)
    input:
        expand("data/alignments/sets/{aligner}~{ref}~{sampleset}.bam",
               ref=config["mapping"]["refs"],
               aligner=config["mapping"]["aligners"],
               sampleset=allsets),
        expand("data/alignments/bamlists/{aligner}~{ref}~all_samples.bamlist",
               ref=config["mapping"]["refs"],
               aligner=config["mapping"]["aligners"],
               sampleset=allsets),
