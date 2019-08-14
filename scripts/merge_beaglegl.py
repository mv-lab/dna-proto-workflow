import gzip
from tempfile import mktemp
from subprocess import Popen, PIPE, STDOUT, check_output
import os

infiles = list(sorted(snakemake.input))
output = snakemake.output.beagle
nsnps_file = snakemake.output.nsnps
threads = snakemake.threads
wcfile = mktemp()

if output.endswith(".gz"):
    outcmd = "pigz -p {}".format(threads)
else:
    outcmd = "cat -"

with Popen("{} >{}".format(outcmd, output), shell=True, stdin=PIPE) as outproc:
    out = outproc.stdin
    # Grab header
    hdr = check_output("zcat {} | head -n 1".format(infiles[0]), shell=True)
    out.write(hdr)
    # Grab all but header for remainder of files
    for infile in infiles:
        Popen("zcat {} | tail -n +2 | tee >(wc -l >>{})".format(infile, wcfile),
              executable="/bin/bash", shell=True, stdout=out).wait()

with open(nsnps_file, "w") as ofh:
    nsnps = sum([int(l.strip()) for l in open(wcfile)])
    print(nsnps, file=ofh)
os.unlink(wcfile)
