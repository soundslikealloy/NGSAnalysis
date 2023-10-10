# -*- coding: utf-8 -*-
"""
Created on Fri Oct  6 13:38:18 2023

@author: 2424069M
"""

from datetime import datetime
import re

from glob import glob
from pydna.amplify import pcr
from Bio import SeqIO
from Bio import Align
from Bio.Seq import Seq, complement
from Bio.Align.Applications import MuscleCommandline
from pyteomics import fasta
# from StringIO import StringIO
# from Bio import AlignIO
import subprocess

from pathlib import Path

# Directories definition
root_dir = Path(__file__).resolve().parent.parent

# variables definition
get_product = ['16S ribosomal RNA']
dAmpl16S = {}                                                                    # Dictionary with amplicon 16S
now = datetime.now()
time = now.strftime("%m-%d-%Y_%H;%M;%S")
aligner = Align.PairwiseAligner()

# Primers
primers = open('primers\main_PrimersList.txt', 'r')
fwd = primers.readline()
rev = primers.readline()
primers.close()

# Max alignment score (local)
local_count = fwd.count('A') + fwd.count('T') + fwd.count('C') + fwd.count('G') # Only count ATCG nucleotides

fastafilename = 'fasta_out/fa (' + str(time) + ').txt'
fastafilealignmentname = 'fasta_out/fa_alignment (' + str(time) + ').fa'
fastafilealignmentsortname = 'fasta_out/fa_alignment_sort (' + str(time) + ').fa'

fastafile = open(fastafilename, 'a')

# Branch gb
print('>> Getting amplicons from (.gb):')
for file in glob('gb_in\*.gb'):
    dAmpl16S[file] = []
    iAmpl = 0
    id_strain = str(re.findall(r'\(.*?\)', str(file)))
    for gb_record in SeqIO.parse(open(file, 'r'), 'genbank'):
        gb_description = ' > %s (%s), %i features' % (gb_record.description, gb_record.name, len(gb_record.features))
        print(gb_description)
        
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
                    # t_saved = complement(t.seq)
                    t_saved = t.seq.reverse_complement()
                else:
                    t_saved = t.seq
                # unique amplicons here (if necessary)
                
                # Write amplicons to FASTA file
                fastafile.write('>' + id_strain[3:-3] + '_' + str(iAmpl) + '\n')
                fastafile.write(str(t_saved) + '\n')    # t_saved = t.seq
                dAmpl16S[file].append(t_saved)            # t_saved = t.seq

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

# Alignment analysis (#missmatches)
# Lorem ipsum...