#######################################################################
#                            Init rules                               #
#######################################################################

from utils import snkmk
import os

configfile: "configs/toolconfig.yml"
configfile: "configs/runconfig.yml"
report: "../report/workflow.rst"


# Check if the reference is ready
is_fai = os.path.exists('rawdata/reference/genome.fa.fai')
print ('Reference is ready? ', is_fai)
if not is_fai:
    print ('>> Preparing reference...')
    os.system('samtools faidx rawdata/reference/genome.fa')
    print ('Done!')


# Prepare contigs of interest (bed)
print ('Prepating contigs of interest ...')
os.system(r"awk r'BEGIN {FS='\t'}; {print $1 FS '0' FS $2}' rawdata/reference/genome.fa.fai > metadata/contigs_of_interest.bed")


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
