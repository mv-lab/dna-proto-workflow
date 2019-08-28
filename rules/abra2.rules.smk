rule abra2:
    input:
        bam="data/alignments/sets/bwa~genome~all_samples.bam",
        #ref= lambda wc: config['refs'][wc.ref],
    output:
        "data/abra/bwa~genome~all_samples.bam",
    log:
        log = "data/log/bwa~genome.log",
    threads:
        3
    params:
        region = config['abra2']['regions'],
        abra_temp = config['abra2']['temp'],
        abra_release = config['abra2']['release'],
        mem= config['abra2']['memory'],
    shell:
        "( java"
        "   -{params.mem}"
        "   -jar {params.abra_release}"
        "   --in {input.bam}"
        "   --out {output}"
        "   --ref rawdata/reference/genome.fa"
        "   --threads {threads}"
        "   --targets {params.region}"
        "   --tmpdir {params.abra_temp}"
        ") >{log.log} 2>&1"
