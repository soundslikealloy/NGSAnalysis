# -*- coding: utf-8 -*-
"""
Created on Fri Oct  6 13:38:18 2023

@author: 2424069M
"""

from glob import glob
from pydna.amplify import pcr
from Bio import SeqIO
from Bio.Seq import Seq, complement

# variables definition
get_product = ['16S ribosomal RNA']
fwd = 'CCTAYGGGRBGCASCAG'
rev = 'GGACTACNNGGGTATCTAAT'

dict16S = dict.fromkeys(glob('*.gb'), [])
print(dict16S)

for file in glob('*.gb'):
    for gb_record in SeqIO.parse(open(file,'r'), 'genbank'):
        print('>> %s (%s), %i features' % (gb_record.description, gb_record.name, len(gb_record.features)))
        print()
        
        # Selection of 16S ribosomal RNAs
        for iF in gb_record.features:
            i_product = iF.qualifiers.get('product')
            if i_product == get_product:
                i_start = iF.location.nofuzzy_start
                i_end = iF.location.nofuzzy_end
                i_seq = gb_record.seq[i_start:i_end]
                i_seq_len = len(i_seq)
                
                # get unique amplicons
                t = pcr((fwd, rev), i_seq)
                if not any([dict16S[file].count(t.seq), dict16S[file].count(complement(t.seq))]):
                    dict16S[file].append(t.seq)
                    print(t.seq)                                                # DEBUGGING
                    print()                                                     # DEBUGGING