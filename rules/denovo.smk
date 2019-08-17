#######################################################################
#                      De-novo Distance analysis                      #
#######################################################################

rule mashsketch:
    input:
        lambda wc: expand("data/reads/samples/{sample}.fastq.gz",
                          sample=SAMPLESETS[wc.set]),
    output:
        temp("data/mash/k{ksize}-s{sketchsize}/{set}.msh"),
    log:
        "data/log/mash/sketch/k{ksize}-s{sketchsize}-{set}.log"
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
        "data/mash/k{ksize}-s{sketchsize}/{set}.msh"
    output:
        dist="data/mash/k{ksize}-s{sketchsize}/{set}.dist",
    log:
        "data/log/mash/dist/k{ksize}-s{sketchsize}-{set}.log"
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
        "data/reads/samples/{sample}.fastq.gz",
    output:
        ct=temp("data/kwip/sketch/k{ksize}-s{sketchsize}/{sample}.ct.gz"),
        info="data/kwip/sketch/k{ksize}-s{sketchsize}/{sample}.ct.gz.info",
        tsv="data/kwip/sketch/k{ksize}-s{sketchsize}/{sample}.ct.gz.info.tsv",
    log:
        "data/log/kwip/sketch/k{ksize}-s{sketchsize}-{sample}.log"
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
        lambda wc: expand("data/kwip/sketch/k{ksize}-s{sketchsize}/{sample}.ct.gz",
                            ksize=wc.ksize, sketchsize=wc.sketchsize,
                            sample=SAMPLESETS[wc.set]),
    output:
        d="data/kwip/k{ksize}-s{sketchsize}/{set}.dist",
        k="data/kwip/k{ksize}-s{sketchsize}/{set}.kern",
    log:
        "data/log/kwip/dist/k{ksize}-s{sketchsize}-{set}.log"
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
        lambda wc: expand("data/reads/samples/{sample}.fastq.gz",
                          sample=SAMPLESETS[wc.set]),
    output:
        "data/readstats/unique-kmers/{set}.tsv",
    threads:
        27
    params:
        kmersize=config["denovodist"]["ksize"],
    log:
        "data/log/readstats/unique-kmers/{set}.log",
    shell:
        "( kdm-unique-kmers.py"
        "    -t {threads}"
        "    -k {params.kmersize}"
        "    {input}"
        "    >{output}"
        " ) 2>{log}"


rule sourmash_sketch:
    input:
        "data/reads/samples/{sample}.fastq.gz",
    output:
        temp("data/sourmash/sketch/k{ksize}-s{sketchsize}/{sample}.smh"),
    log:
        "data/log/sourmash/sketch/k{ksize}-s{sketchsize}-{sample}.log"
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
        lambda wc: expand("data/sourmash/sketch/k{ksize}-s{sketchsize}/{sample}.smh",
                            ksize=wc.ksize, sketchsize=wc.sketchsize,
                            sample=SAMPLESETS[wc.set]),
    output:
        "data/sourmash/k{ksize}-s{sketchsize}/{set}.dist",
    log:
        "data/log/sourmash/dist/k{ksize}-s{sketchsize}-{set}.log"
    threads: 1
    shell:
        "(sourmash compare -k {wildcards.ksize} -o {output} {input} ) >{log} 2>&1"


rule kwip:
    input:
        expand("data/kwip/k{ksize}-s{sketchsize}/{set}.dist",
               ksize=config["denovodist"]["ksize"],
               sketchsize=config["denovodist"]["kwip_sketchsize"],
               set=config["denovodist"]["kwip_sets"]),

rule sourmash:
    input:
        expand("data/sourmash/k{ksize}-s{sketchsize}/{set}.dist",
               ksize=config["denovodist"]["ksize"],
               sketchsize=config["denovodist"]["sourmash_sketchsize"],
               set=config["denovodist"]["sourmash_sets"]),

rule mash:
    input:
        expand("data/mash/k{ksize}-s{sketchsize}/{set}.dist",
               ksize=config["denovodist"]["ksize"],
               sketchsize=config["denovodist"]["mash_sketchsize"],
               set=config["denovodist"]["mash_sets"]),

rule denovo:
    input:
        rules.kwip.input,
        rules.mash.input,
        rules.sourmash.input,
