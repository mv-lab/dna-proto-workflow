#######################################################################
#                      De-novo Distance analysis                      #
#######################################################################


##### Target rules #####

rule kwip:
    input:
        expand("output/denovo/kwip/k{ksize}-s{sketchsize}/{set}.dist",
               ksize=config["denovodist"]["ksize"],
               sketchsize=config["denovodist"]["kwip_sketchsize"],
               set=config["denovodist"]["kwip_sets"]),

rule sourmash:
    input:
        expand("output/denovo/sourmash/k{ksize}-s{sketchsize}/{set}.dist",
               ksize=config["denovodist"]["ksize"],
               sketchsize=config["denovodist"]["sourmash_sketchsize"],
               set=config["denovodist"]["sourmash_sets"]),

rule mash:
    input:
        expand("output/denovo/mash/k{ksize}-s{sketchsize}/{set}.dist",
               ksize=config["denovodist"]["ksize"],
               sketchsize=config["denovodist"]["mash_sketchsize"],
               set=config["denovodist"]["mash_sets"]),

rule pca_mash:
    input:
        rules.mash.input, 
    output:
        "output/plots/denovo/mash/pca.pdf", 
    script:
        "../scripts/pca_mash.R"

rule pca_kwip:
    input:
        rules.kwip.input,
    output:
        "output/plots/denovo/kwip/pca.pdf",
    script:
        "../scripts/pca_kwip.R"



rule denovo:
    input:
        rules.kwip.input,
        rules.mash.input,
        rules.sourmash.input,
        rules.pca_mash.output,
	rules.pca_kwip.output,

##### Actual rules #####

rule mashsketch:
    input:
        lambda wc: expand("output/reads/samples/{sample}.fastq.gz",
                          sample=SAMPLESETS[wc.set]),
    output:
        temp("output/denovo/mash/k{ksize}-s{sketchsize}/{set}.msh"),
    log:
        "output/log/denovo/mash/sketch/k{ksize}-s{sketchsize}-{set}.log"
    threads: 27
    shell:
        " mash sketch"
        "   -k {wildcards.ksize}"
        "   -s {wildcards.sketchsize}"
        "   -p {threads}"
        "   -o {output}"
        "   {input}"
        " >{log} 2>&1"


rule mashdist:
    input:
        "output/denovo/mash/k{ksize}-s{sketchsize}/{set}.msh"
    output:
        dist="output/denovo/mash/k{ksize}-s{sketchsize}/{set}.dist",
    log:
        "output/log/denovo/mash/dist/k{ksize}-s{sketchsize}-{set}.log"
    threads: 27
    shell:
        "mash dist"
        "   -p {threads}"
        "   -t" # tabular format
        "   {input} {input}" # needs input twice
        " >{output}"
        " 2>{log}"

rule countsketch:
    input:
        "output/reads/samples/{sample}.fastq.gz",
    output:
        ct=temp("output/denovo/kwip/sketch/k{ksize}-s{sketchsize}/{sample}.ct.gz"),
        info="output/denovo/kwip/sketch/k{ksize}-s{sketchsize}/{sample}.ct.gz.info",
        tsv="output/denovo/kwip/sketch/k{ksize}-s{sketchsize}/{sample}.ct.gz.info.tsv",
    log:
        "output/log/denovo/kwip/sketch/k{ksize}-s{sketchsize}-{sample}.log"
    threads:
        3
    shell:
        "load-into-counting.py"
        "   -N 1"
        "   -x {wildcards.sketchsize}"
        "   -k {wildcards.ksize}"
        "   -b"
        "   -f"
        "   -s tsv"
        "   -T {threads}"
        "   {output.ct}"
        "   {input}"
        " >{log} 2>&1"

rule kwipdist:
    input:
        lambda wc: expand("output/denovo/kwip/sketch/k{ksize}-s{sketchsize}/{sample}.ct.gz",
                            ksize=wc.ksize, sketchsize=wc.sketchsize,
                            sample=SAMPLESETS[wc.set]),
    output:
        d="output/denovo/kwip/k{ksize}-s{sketchsize}/{set}.dist",
        k="output/denovo/kwip/k{ksize}-s{sketchsize}/{set}.kern",
    log:
        "output/log/denovo/kwip/dist/k{ksize}-s{sketchsize}-{set}.log"
    threads:
        4
    shell:
        "kwip"
        " -d {output.d}"
        " -k {output.k}"
        " -t {threads}"
        " {input}"
        " >{log} 2>&1"

rule unique_kmers:
    input:
        lambda wc: expand("output/reads/samples/{sample}.fastq.gz",
                          sample=SAMPLESETS[wc.set]),
    output:
        "output/readstats/unique-kmers/{set}.tsv",
    threads:
        27
    params:
        kmersize=config["denovodist"]["ksize"],
    log:
        "output/log/denovo/readstats/unique-kmers/{set}.log",
    shell:
        "( kdm-unique-kmers.py"
        "    -t {threads}"
        "    -k {params.kmersize}"
        "    {input}"
        "    >{output}"
        " ) 2>{log}"


rule sourmash_sketch:
    input:
        "output/reads/samples/{sample}.fastq.gz",
    output:
        temp("output/denovo/sourmash/sketch/k{ksize}-s{sketchsize}/{sample}.smh"),
    log:
        "output/log/denovo/sourmash/sketch/k{ksize}-s{sketchsize}-{sample}.log"
    shell:
        "( sourmash compute"
        "   --name '{wildcards.sample}'"
        "   -k {wildcards.ksize}"
        "   -n {wildcards.sketchsize}"
        "   -o {output}"
        "   {input}"
        ") >{log} 2>&1"

rule sourmash_dist:
    input:
        lambda wc: expand("output/denovo/sourmash/sketch/k{ksize}-s{sketchsize}/{sample}.smh",
                            ksize=wc.ksize, sketchsize=wc.sketchsize,
                            sample=SAMPLESETS[wc.set]),
    output:
        "output/denovo/sourmash/k{ksize}-s{sketchsize}/{set}.dist",
    log:
        "output/log/denovo/sourmash/dist/k{ksize}-s{sketchsize}-{set}.log"
    threads: 1
    shell:
        "(sourmash compare -k {wildcards.ksize} -o {output} {input} ) >{log} 2>&1"
