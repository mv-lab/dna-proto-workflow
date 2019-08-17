from utils import snkmk

configfile: "config.yml"

RUNLIB2SAMP, SAMP2RUNLIB = snkmk.make_runlib2samp("metadata/sample2runlib.csv")

SAMPLESETS = snkmk.make_samplesets(s2rl_file="metadata/sample2runlib.csv",
                                   setfile_glob="metadata/samplesets/*.txt")

VARCALL_REGIONS = snkmk.make_regions(config["refs"], window=config["varcall"]["chunksize"])

shell.prefix = "set -euo pipefail; "

wildcard_constraints:
    run="[^/]+",
    lib="[^/]+",
    aligner="[^/]+",
    sample="[^/]+",
    ref="[^/]+",
    type="[^/]+",
