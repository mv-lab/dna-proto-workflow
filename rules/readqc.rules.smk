#######################################################################
#                            Read-level QC                            #
#######################################################################

##### Target rules #####

rule qc_runlib:
    input:
        ["data/reads/runs/{run}/{lib}.fastq.gz".format(run=run, lib=lib)
		for run, lib in RUNLIB2SAMP],

rule read_stats:
    input:
        "data/stats/reads/readnum_librun.tsv",
        "data/stats/reads/readnum_samples.tsv",

rule qc_samples:
    input:
        expand("data/reads/samples/{sample}.fastq.gz", sample=SAMP2RUNLIB),

rule readqc:
    input:
        rules.qc_runlib.input,
        rules.read_stats.input,
        rules.qc_samples.input,


##### Actual rules #####


ruleorder: qcreads_il > qcreads
rule qcreads:
    input:
        unpack(get_fr_fastq)
        #r1="rawdata/runs/{run}/{lib}_R1.fastq.gz",
        #r2="rawdata/runs/{run}/{lib}_R2.fastq.gz",
    output:
        reads="data/reads/runs/{run}/{lib}.fastq.gz",
    log:
        log="log/adapterremoval/{run}/{lib}.log",
        settings="data/stats/adapterremoval/{run}/{lib}.txt",
    threads:
        7
    params:
        adp1=lambda wc: config["qc"].get(wc.run, config["qc"]["_DEFAULT_"])["adapter1"],
        adp2=lambda wc: config["qc"].get(wc.run, config["qc"]["_DEFAULT_"])["adapter2"],
        minqual=lambda wc: config["qc"].get(wc.run, config["qc"]["_DEFAULT_"])["minqual"],
    shell:
        "( AdapterRemoval"
        "   --file1 {input.r1}"
        "   --file2 {input.r2}"
        "   --adapter1 {params.adp1}"
        "   --adapter2 {params.adp2}"
        "   --combined-output"
        "   --interleaved-output"
        "   --trimns"
        "   --trimqualities"
        "   --trimwindows 10"
        "   --minquality {params.minqual}"
        "   --threads 2"
        "   --settings {log.settings}"
        "   --output1 /dev/stdout"
        " | seqhax pairs"
        "   -l 20"
        "   -b >(pigz -p 5 >{output.reads})"
        "   /dev/stdin"
        ") >{log.log} 2>&1"


#localrules: qcreads
rule qcreads_il:
    input:
        unpack(get_il_fastq),
    output:
        reads="data/reads/runs/{run}/{lib}.fastq.gz",
    log:
        log="log/adapterremoval/{run}/{lib}.log",
        settings="data/stats/adapterremoval/{run}/{lib}.txt",
    threads:
        7
    params:
        adp1=lambda wc: config["qc"].get(wc.run, config["qc"]["_DEFAULT_"])["adapter1"],
        adp2=lambda wc: config["qc"].get(wc.run, config["qc"]["_DEFAULT_"])["adapter2"],
        minqual=lambda wc: config["qc"].get(wc.run, config["qc"]["_DEFAULT_"])["minqual"],
    shell:
        "( AdapterRemoval"
        "   --file1 {input}"
        "   --adapter1 {params.adp1}"
        "   --adapter2 {params.adp2}"
        "   --combined-output"
        "   --interleaved"
        "   --interleaved-output"
        "   --trimns"
        "   --trimqualities"
        "   --trimwindows 10"
        "   --minquality {params.minqual}"
        "   --threads 2"
        "   --settings {log.settings}"
        "   --output1 /dev/stdout"
        " | seqhax pairs"
        "   -l 20"
        "   -b >(pigz -p 5 >{output.reads})"
        "   /dev/stdin"
        ") >{log.log} 2>&1"


localrules: samplefastq
rule samplefastq:
    input:
        lambda wc: ["data/reads/runs/{run}/{lib}.fastq.gz".format(run=r, lib=l) for r, l in SAMP2RUNLIB[wc.sample]],
    output:
        "data/reads/samples/{sample}.fastq.gz"
    log:
        "log/samplefastq/{sample}.log"
    threads: 1
    shell:
        "cat {input} > {output}"

rule read_count_librun:
    input:
        ["data/reads/runs/{run}/{lib}.fastq.gz".format(run=run, lib=lib)
		for run, lib in RUNLIB2SAMP],
    output:
        "data/stats/reads/readnum_librun.tsv",
    threads:
        28
    log:
        "log/readstats/seqhax-stats-librun.log",
    shell:
        "( seqhax stats"
        "    -t {threads}"
        "    {input}"
        "    >{output}"
        " ) 2>{log}"


rule read_count_sample:
    input:
    	expand("data/reads/samples/{sample}.fastq.gz", sample=SAMP2RUNLIB),
    output:
        "data/stats/reads/readnum_samples.tsv",
    threads:
        27
    log:
        "log/readstats/seqhax-stats-sample.log",
    shell:
        "( seqhax stats"
        "    -t {threads}"
        "    {input}"
        "    >{output}"
        " ) 2>{log}"
