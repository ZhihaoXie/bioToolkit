#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Version:   V1.0.0
# Author:    Zhihao Xie  \(#^o^)/
# Date:      2019/8/8 16:53
# CopyRight: Copyright (c) Zhihao Xie, All rights reserved.

import os
import sys
import re
from Bio import SeqIO
import gzip

# <---------- main ----------> #
if len(sys.argv) <= 1:
    print("Note: calculate fasta length, include length of exclude N")
    print(f"Usage: python {sys.argv[0]} fasta1 [fasta2 ..] > fasta_length.out")
    sys.exit(1)

for fasta in sys.argv[1:]:
    #fasta = sys.argv[1]
    if fasta.endswith('.gz'):
        fhandle = gzip.open(fasta, 'rt')
    else:
        fhandle = open(fasta)

    for rec in SeqIO.parse(fhandle, 'fasta'):
        rec_id = rec.id
        myl = len(rec)
        tmpseq = str(rec.seq)
        nNumber = tmpseq.count('N')
        nNumber += tmpseq.count('n')
        mylnoN = myl - nNumber
        print(f"{rec_id}\t{myl}\t{mylnoN}")
        
    fhandle.close()
