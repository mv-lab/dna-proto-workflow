
include: "rules/init.rules.smk"

##### Modules #####

include: "rules/denovo.rules.smk"
include: "rules/readqc.rules.smk"
include: "rules/align.rules.smk"
include: "rules/varcall.rules.smk"
include: "rules/stats.rules.smk"
include: "rules/snpeff.rules.smk"

##### Target rules #####

rule all:
    input:
        rules.denovo.input,
        rules.readqc.input,
        rules.align.input,
        rules.varcall.input,
        rules.stats.input,
