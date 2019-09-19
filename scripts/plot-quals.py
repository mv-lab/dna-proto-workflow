import matplotlib
import matplotlib.pyplot as plt
from pysam import VariantFile
import numpy as np

###matplotlib.use("Agg")

print ("$ Script: plot-quals")
print ("INPUT: ", snakemake.input)
quals = [record.qual for record in VariantFile(snakemake.input[0])]
plt.hist(quals)
plt.savefig(snakemake.output[0])
