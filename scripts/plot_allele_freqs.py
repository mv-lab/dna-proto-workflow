import os
import matplotlib
from pathlib import Path
import pandas as pd
from IPython.display import display
import numpy as np
import seaborn as sns
import re
from collections import OrderedDict
from matplotlib import pyplot as plt
sns.set()
pd.options.display.max_columns = None # default is 20
pd.options.display.max_rows = 60 # default is 60

### Set input ###

# Define input table paths
tsv_file = Path('freebayes_D2.filtered.tsv')
raw_output_pdf = Path('freebayes_allele_freq_D2.pdf')
window_output_pdf = Path('freebayes_allele_freq_window_D2.pdf')


# Order as: control, mutant, dwarf/mutant F2 pool, WT F2 pool
samples = ['con-all', 'D2', 'D2_F2_tt', 'D2_F2_TT']
sample_F2_WT = 'D2_F2_TT'
sample_F2_mutant = 'D2_F2_tt'

# These lengths are from the Sbicolor_454_v3.0.1.fa genome
# All chromosomes that should be included in the plots.
chromosome_lengths = OrderedDict([('Chr01', 80884392), 
                      ('Chr02', 77742459),
                      ('Chr03', 74386277),
                      ('Chr04', 68658214),
                      ('Chr05', 71854669),
                      ('Chr06', 61277060),
                      ('Chr07', 65505356),
                      ('Chr08', 62686529),
                      ('Chr09', 59416394),
                      ('Chr10', 61233695)])

# Set window six
windowsize = 5000000
stepsize = 1000000


# This function generates a tidy/longform table, which makes is easier to plot with seaborn.
# Will have one row per sample per variant, sample name in "sample", and rows for control/mutant observations and
# frequency of mutant alleles
#

def getTidyTable(raw_table, samples):

    # collects new rows
    row_accumulator = []

    def splitListToRows(row):
        '''Split one row into long form, one new row per sample'''
        if row['CHROM'] not in chromosome_lengths:
            # Only keep markers on the included chromosomes
            return

        split_gt = row['GT'].split(',')
        if '.' in split_gt:
            # Remove rows where not all samples are genotyped
            return
        split_gq = row['GQ'].split(',')
        split_ro = row['SampleRO'].split(',')
        split_ao = row['SampleAO'].split(',')

        control_ref = split_gt[0] == '0/0'
        new_rows = []

        for i in range(4):
            new_row = row.to_dict()

            # Remove rows unnecessary for plotting
            new_row.pop('REF')
            new_row.pop('ALT')
            new_row.pop('RO')
            new_row.pop('AO')

            new_row['GT'] = split_gt[i]
            new_row['GQ'] = int(split_gq[i])
            new_row['SampleRO'] = split_ro[i]
            new_row['SampleAO'] = split_ao[i]
            new_row['sample'] = samples[i]

            # Get observations of WT/mutant alleles
            if control_ref:
                new_row['controlO'] = int(new_row['SampleRO'])
                new_row['mutantO'] = int(new_row['SampleAO'])
            else:
                new_row['controlO'] = int(new_row['SampleAO'])
                new_row['mutantO'] = int(new_row['SampleRO'])

            if (new_row['mutantO'] + new_row['controlO']) == 0:
                # Remove rows where there are no observations for a sample
                break
            else:
                # Get mutant allele frequency
                new_row['mutant_freq'] = new_row['mutantO'] / (new_row['mutantO'] + new_row['controlO'])
            new_rows.append(new_row)

        if len(new_rows) == 4:
            # only keep rows if there is one for each sample
            row_accumulator.extend(new_rows)

    raw_table.apply(splitListToRows, axis=1)
    table = pd.DataFrame(row_accumulator)
    table = table[['CHROM', 'POS', 'sample', 'GT', 'GQ',
                   'SampleRO', 'SampleAO', 'controlO', 'mutantO', 'mutant_freq']]
#     table = table[table['CHROM'].str.startswith('Chr')]
    return table



# Create a table with averaged mutant frequencies for the F2 pools with a sliding window

def get_window_table(table, sample_tt, sample_TT, windowsize, stepsize):

    # positions to include in window before/after the current pos, i.e. half the window size
    w = int(windowsize/2)

    rows = []
    for chrom in chromosome_lengths:

        # get F2 pool genotypes for current chromosomes
        chrom_table = table[(table['CHROM'] == chrom)
                               & ((table['sample'] == sample_TT)
                                  | (table['sample'] == sample_tt))].reset_index(drop=True)

        for i in range(1, chromosome_lengths[chrom]+1, stepsize):
            # Loop over windows

            wstart = max(1, i-w)
            wend = min(chromosome_lengths[chrom], i+w)

            # Get table for window
            wtable = chrom_table[(chrom_table['POS'] >= wstart)
                                 & (chrom_table['POS'] <= wend)]

            freqs_tt = wtable[wtable['sample'] == sample_tt]['mutant_freq']
            freqs_TT = wtable[wtable['sample'] == sample_TT]['mutant_freq']
            varcount = int(wtable.shape[0] / 2) # because of two samples per marker
            avg_TT = np.mean(freqs_TT)
            avg_tt = np.mean(freqs_tt)
            std_TT = np.std(freqs_TT)
            std_tt = np.std(freqs_tt)

            wsize = wend - wstart
            row_tt = {'CHROM': chrom,
                      'POS': i,
                      'sample': sample_tt,
                      'avg_mutant_freq': avg_tt,
                      'std': std_tt,
                      'varcount': varcount,
                      'window_size': wsize}
            rows.append(row_tt)
            row_TT = {'CHROM': chrom,
                      'POS': i,
                      'sample': sample_TT,
                      'avg_mutant_freq': avg_TT,
                      'std': std_TT,
                      'varcount': varcount,
                      'window_size': wsize}
            rows.append(row_TT)
    window_table = pd.DataFrame(rows)
    window_table = window_table[['CHROM', 'POS', 'sample', 'avg_mutant_freq', 'std', 'varcount', 'window_size']]
    return window_table



#----------------------  LOAD INPUT ---------------------------

# Load TSV table
raw_table = pd.read_csv(tsv_file, sep='\t', na_values=['.'], comment='#')

# Remove leading commas
for col in ['GT', 'GQ', 'SampleRO', 'SampleAO']:
    raw_table.loc[:,col] = raw_table[col].str[1:]

table = getTidyTable(raw_table, samples)


#---------------------- Plot raw frequencies ------------------

# Plot raw mutant allele frequencies of F2 pools
plot_table = table[(table['sample'] == sample_F2_mutant) | (table['sample'] == sample_F2_WT)]
plot = sns.relplot(data=plot_table, x='POS', y='mutant_freq', row='CHROM', style='sample',
                   hue='sample', aspect=7.0, height=4., kind='line', markers=True, dashes=False, linewidth=0.5)

plot.savefig(raw_output_pdf)



### Plot sliding-window averaged frequencies ###

windowsize=5000000
stepsize=1000000

wtable = get_window_table(table, sample_F2_mutant, sample_F2_WT, windowsize, stepsize)

# Plot the per-window average allele frequencies
# Create one subplot per chromosome
wgrid = sns.FacetGrid(wtable, row='CHROM',aspect=7.0, height=4., hue='sample',legend_out=True)

# Plot lines
wgrid = wgrid.map(sns.lineplot, 'POS', 'avg_mutant_freq', linewidth=0.5).add_legend()

# Overlay with scatterplot with densities per chromosome
for i, chrom in enumerate(chromosome_lengths):
    plotdata = wtable[wtable['CHROM']==chrom]
    mincount = min(plotdata['varcount'])
    maxcount = max(plotdata['varcount'])
    sizes=(np.log2(mincount+1)*100,np.log2(maxcount+1)*100)
    sns.scatterplot(data=wtable[wtable['CHROM']==chrom],
                    x='POS', y='avg_mutant_freq', size='varcount', hue='sample', legend=False,
                    ax=wgrid.axes[i,0], sizes=sizes)

wgrid.savefig(window_output_pdf)
