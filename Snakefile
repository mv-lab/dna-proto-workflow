
include: "rules/init.rules.smk"

##### Modules #####

include: "rules/denovo.rules.smk"
include: "rules/readqc.rules.smk"
include: "rules/align.rules.smk"
include: "rules/varcall.rules.smk"
include: "rules/stats.rules.smk"
include: "rules/abra2.rules.smk"

##### Target rules #####


rule all:
    input:
        #rules.denovo.input,
        rules.reads.input,
        rules.align.input,
        #rules.abra2.output,
        rules.varcall.input,
        #rules.stats.input,


