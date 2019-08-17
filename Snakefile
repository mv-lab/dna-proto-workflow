
include: "rules/common.smk"

##### Modules #####

include: "rules/denovo.smk"
include: "rules/readqc.smk"
include: "rules/mapping.smk"
include: "rules/calling.smk"

##### Target rules #####

rule all:
    input:
        rules.denovo.input,
        rules.reads.input,
        rules.align.input,
        rules.varcall.input,
