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

import warnings
warnings.filterwarnings("ignore")

# Directories definition
root_dir = Path(__file__).resolve().parent.parent

# Variables definition
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

# Unique amplicons
unique_amplicons = 1                                                            # Unique amplicons [0/1]
if unique_amplicons == 1: dAmpl16S = {}                                         # Dictionary with amplicon 16S

# Fasta file creation
if unique_amplicons == 1:
    fastafilename = 'out_fasta/fa_unique (' + str(time) + ').txt'
    fastafilealignmentname = 'out_fasta/fa_unique_alignment (' + str(time) + ').fa'
    fastafilealignmentsortname = 'out_fasta/fa_unique_alignment_sort (' + str(time) + ').fa'
else:
    fastafilename = 'out_fasta/fa (' + str(time) + ').txt'
    fastafilealignmentname = 'out_fasta/fa_alignment (' + str(time) + ').fa'
    fastafilealignmentsortname = 'out_fasta/fa_alignment_sort (' + str(time) + ').fa'

fastafile = open(fastafilename, 'a')

# Branch gb
print('>> Getting amplicons from (.gb):')
for file in glob('in_gb\*.gb'):
    if unique_amplicons == 1: dAmpl16S[file] = []
    iAmpl = 0
    if unique_amplicons == 1: iAmpl_unique = 0 
    id_strain = str(re.findall(r'\(.*?\)', str(file)))
    for gb_record in SeqIO.parse(open(file, 'r'), 'genbank'):
        # Selection of feature's region (no ITS)
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
                
                # Write (unique) amplicons to FASTA file
                if unique_amplicons == 1:
                    if not dAmpl16S[file].count(t_saved):
                        iAmpl_unique += 1
                        dAmpl16S[file].append(t_saved)
                        fastafile.write('>' + id_strain[3:-3] + '_' + str(iAmpl) + '\n')
                        fastafile.write(str(t_saved) + '\n')
                else:   
                    fastafile.write('>' + id_strain[3:-3] + '_' + str(iAmpl) + '\n')
                    fastafile.write(str(t_saved) + '\n')
        
        # Printing results
        if unique_amplicons == 0:
            if iAmpl == 1:
                copystr = 'copy'
            else:
                copystr = 'copies'
            gb_description = ' > %s (%s), %i %s of %s' % (gb_record.description, gb_record.name, iAmpl, copystr, get_product_str)
        else:
            if iAmpl_unique == 1:
                copystr = 'copy'
            else:
                copystr = 'copies'
            gb_description = ' > %s (%s), %i unique %s of %s' % (gb_record.description, gb_record.name, iAmpl_unique, copystr, get_product_str)
        print(gb_description)
        
# Branch FASTA
print('\n>> Getting amplicons from (.fa):')
for file in glob('in_fasta\*.fa'):
    if unique_amplicons == 1: dAmpl16S[file] = []
    iAmpl = 0
    if unique_amplicons == 1: iAmpl_unique = 0 
    id_strain = str(re.findall(r'\(.*?\)', str(file)))
    f = open(file, 'r')
    lines = f.readlines() 
    hre = re.compile('>(\S+)')

    # Detect if 'file' is a full sequence or fwd/rev read and if a sequence is from fwd/rev
    frre = re.compile('__')
    outfr = frre.search(file)
    
    if outfr:
        # FWD and REV read
        # Lorem ipsum...
        print('')
    else:
        # Full sequence read
        for line in lines:
            outh = hre.search(line)
            if outh:     
                iAmpl += 1
                id_strain_saved = '>' + id_strain[3:-3] + '_' + str(iAmpl) + '\n'
            else:
                # Get (unique) amplicons
                i_seq = Seq(line)
                t = pcr((fwd, rev), i_seq)
                c_t = aligner.align(t.seq[0:len(fwd)-1], fwd)
                # Flip amplicon (if necessary)
                if c_t.score < local_count:
                    t_saved = t.seq.reverse_complement()
                else:
                    t_saved = t.seq
                    
                # Write (unique) amplicons to FASTA file
                if unique_amplicons == 1:
                    if not dAmpl16S[file].count(t_saved):
                        iAmpl_unique += 1
                        dAmpl16S[file].append(t_saved)
                        fastafile.write(id_strain_saved)
                        fastafile.write(str(t_saved) + '\n')
                else:
                    fastafile.write(id_strain_saved)
                    fastafile.write(str(t_saved) + '\n')
        f.close()
        
        # Printing results
        fa_name = file.replace('in_fasta\\', '')
        fa_name = fa_name.replace('.fa', '')
        if unique_amplicons == 0:
            if iAmpl == 1:
                copystr = 'copy'
            else:
                copystr = 'copies'
            fa_description = ' > %s, %i %s of %s' % (fa_name, iAmpl, copystr, get_product_str)
        else:
            if iAmpl_unique == 1:
                copystr = 'copy'
            else:
                copystr = 'copies'
            fa_description = ' > %s, %i unique %s of %s' % (fa_name, iAmpl_unique, copystr, get_product_str)
        print(fa_description)

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
filesave = file_misAnalysis.replace('out_fasta/', 'out_misAnalysis/')
if unique_amplicons == 1:
    filesave = filesave.replace('fa_unique_alignment', 'mismatchTable_unique')
    figurefilename = 'out_misAnalysis/alignmismatches_unique (' + str(time) + ').png'
else:
    filesave = filesave.replace('fa_alignment', 'mismatchTable')
    figurefilename = 'out_misAnalysis/alignmismatches (' + str(time) + ').png'
filesave = filesave.replace('.fa', '.txt')

# Alignment dimensions
a = AlignIO.read(file_misAnalysis, 'fasta')
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
fig.savefig(figurefilename)

print('>> Plotting done!')