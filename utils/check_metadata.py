import snkmk
from check_config import *

print ('\nMetadata ...\n')

config = readconfig('config.yml')

RUNLIB2SAMP, SAMP2RUNLIB = snkmk.make_runlib2samp("metadata/sample2runlib.csv")

print ('RUNLIB2SAMP\n',RUNLIB2SAMP)
print ('\nSAMP2RUNLIB\n',SAMP2RUNLIB)

SAMPLESETS = snkmk.make_samplesets(s2rl_file="metadata/sample2runlib.csv",
                                   setfile_glob="metadata/samplesets/*.txt")
print ('\nSAMPLESETS\n',SAMPLESETS)

VARCALL_REGIONS = snkmk.make_regions(config["refs"], window=config["varcall"]["chunksize"])
print ('\nVARCALL_REGIONS\n',VARCALL_REGIONS,'\n')

