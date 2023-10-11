# -*- coding: utf-8 -*-
"""
Created on Fri Oct  6 13:38:18 2023

@author: 2424069M
"""

from datetime import datetime
from pathlib import Path
from glob import glob
import subprocess
import re

from pydna.amplify import pcr
from Bio import AlignIO, SeqIO
from Bio import Align
from Bio.Seq import Seq
from pyteomics import fasta

from matplotlib import colors
from itertools import combinations
import matplotlib.pyplot as plt
import numpy as np
import math

# import warnings
# warnings.filterwarnings("ignore")

# Directories definition
root_dir = Path(__file__).resolve().parent.parent

# Variables definition
dAmpl16S = {}                                                                    # Dictionary with amplicon 16S
now = datetime.now()
time = now.strftime("%m-%d-%Y_%H;%M;%S")
aligner = Align.PairwiseAligner()

# Primers
primers = open('primers\main_PrimersList.txt', 'r')
fwd = primers.readline()
rev = primers.readline()
primers.close()

# Feature read
featureRead = open('feature/feature.txt', 'r')
get_product_str = featureRead.readline()
get_product = [get_product_str]
featureRead.close()

# Max alignment score (local)
local_count = fwd.count('A') + fwd.count('T') + fwd.count('C') + fwd.count('G') # Only count ATCG nucleotides

# Fasta file creation
fastafilename = 'out_fasta/fa (' + str(time) + ').txt'
fastafilealignmentname = 'out_fasta/fa_alignment (' + str(time) + ').fa'
fastafilealignmentsortname = 'out_fasta/fa_alignment_sort (' + str(time) + ').fa'
fastafile = open(fastafilename, 'a')

# Unique amplicon
unique_amplicons = 0

# Branch gb
print('>> Getting amplicons from (.gb):')
for file in glob('in_gb\*.gb'):
    dAmpl16S[file] = []
    iAmpl = 0
    id_strain = str(re.findall(r'\(.*?\)', str(file)))
    for gb_record in SeqIO.parse(open(file, 'r'), 'genbank'):
        # Selection of 16S ribosomal RNAs
        for iF in gb_record.features:
            i_product = iF.qualifiers.get('product')
            if i_product == get_product:
                iAmpl += 1
                i_start = iF.location.nofuzzy_start
                i_end = iF.location.nofuzzy_end
                i_seq = gb_record.seq[i_start:i_end]
                i_seq_len = len(i_seq)
                
                # Get (unique) amplicons
                t = pcr((fwd, rev), i_seq)
                c_t = aligner.align(t.seq[0:len(fwd)-1], fwd)
                # Flip amplicon (if necessary)
                if c_t.score < local_count:
                    t_saved = t.seq.reverse_complement()
                else:
                    t_saved = t.seq
                
                # Unique amplicons (if necessary)
                # Lorem ipsum...
                dAmpl16S[file].append(t_saved_unique)
                
                # Write amplicons to FASTA file
                fastafile.write('>' + id_strain[3:-3] + '_' + str(iAmpl) + '\n')
                fastafile.write(str(t_saved) + '\n')
        
        if unique_amplicons == 0:
            gb_description = ' > %s (%s), %i copies of %s' % (gb_record.description, gb_record.name, iAmpl, get_product_str)
        else:
            gb_description = ' > %s (%s), %i unique copies of %s' % (gb_record.description, gb_record.name, iAmpl, get_product_str)
        print(gb_description)
        
# Branch FASTA
# Lorem ipsum...

fastafile.close()

# Alignment (Next-generation MUSCLE v5)
muscle_exe = 'muscle\muscle5.1.win64.exe'
input_sequences = fastafilename
output_alignment = fastafilealignmentname
print('\n>> Alignment in progress...')
subprocess.run([muscle_exe, '-align', input_sequences, '-output', output_alignment], shell=True, capture_output=True, text=True)
# Sorting FASTA file
with fasta.read(fastafilealignmentname) as f:
    fasta.write(sorted(f), fastafilealignmentname)
print('>> Alignment done!')

# Alignment analysis (#mismatches)
print('\n>> Plotting in progress...')
file_misAnalysis = fastafilealignmentname

# for file in glob('out_fasta\*.fa'):
a = AlignIO.read(file_misAnalysis, 'fasta')
filesave = file_misAnalysis.replace('out_fasta/', 'out_misAnalysis/')
filesave = filesave.replace('fa_alignment', 'mismatchTable')
filesave = filesave.replace('.fa', '.txt')

# Alignment dimensions
n_seq = a.__len__()                                                         # Num of sequences (rows)
alig_len = a.get_alignment_length()                                         # Max length of alignment (columns)
comb = np.array(list(combinations(np.arange(1, n_seq+1), 2)))
num_comb = math.comb(n_seq, 2)
indx1 = comb[:,0]
indx2 = comb[:,1]

# Exporting sequences
mId = 'IDs'
mseq = np.zeros((1, alig_len))
for record in a:
    iseq = list(record.seq)
    mseq = np.vstack([mseq, iseq])
    mId += ' ' + record.id

mseq = mseq[1:n_seq+1, :]
mId = mId.split()[1:n_seq+1]

# Matrix of mismatches
mMismatch = np.ones((n_seq, n_seq)) * (-1)
for iM in range(0, num_comb):
    iIndx1 = indx1[iM] - 1
    iIndx2 = indx2[iM] - 1        
    mseq1 = mseq[iIndx1]
    mseq2 = mseq[iIndx2]
    mismatch_arr = np.compare_chararrays(mseq1, mseq2, "!=", rstrip = True)
    iMismatch = mismatch_arr.sum()
    mMismatch[iIndx2, iIndx1] = iMismatch

# Save results
np.savetxt(filesave, mMismatch, fmt='%.0f')

# Plotting
nrows = n_seq
ncols = n_seq
nTicks = np.arange(0.5, n_seq+0.5)
Z = mMismatch
x = np.arange(ncols + 1)
y = np.flip(np.arange(nrows + 1))
cmap = colors.LinearSegmentedColormap.from_list('', ['blue', 'yellow'])
cmap.set_under('black')
mX, mY = np.meshgrid(x, y)

fig, ax = plt.subplots()
plt.plot(mX, mY, c='k', linewidth = '0.5')
plt.plot(np.transpose(mX), np.transpose(mY), c='k', linewidth = '0.5')
colormesh = ax.pcolormesh(x, y, Z, cmap = cmap, vmin = 0, vmax = 10)        # Z.min() \\ Z.max()
cbar = fig.colorbar(colormesh, label = '# mismatches')
plt.xticks(nTicks, mId, rotation=90)
plt.yticks(nTicks, np.flip(mId))
plt.rcParams.update({'font.size': 45})

# Plot size
fig.set_figheight(n_seq)
fig.set_figwidth(n_seq+10)

# Save plot
figurefilename = 'out_misAnalysis/alignmismatches (' + str(time) + ').png'
fig.savefig(figurefilename)

print('>> Plotting done!')