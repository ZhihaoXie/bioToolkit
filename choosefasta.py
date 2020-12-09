#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# FileName:  truncate_seq.py
# Author:    Zhihao Xie  \(#^o^)/
# Date:      2019/5/6 9:39
# Version:   v1.0.0
# CopyRight: Copyright @Zhihao Xie, All rights reserved.

import sys
import re
import os
from Bio import SeqIO
#from Bio.Seq import Seq
#from Bio.SeqRecord import SeqRecord
#from Bio.Alphabet import IUPAC


if __name__ == '__main__':
    if len(sys.argv) < 3:
        print("Usage:\n\tpython3 {} <list> <fasta> [> <output>]".format(sys.argv[0]))
        sys.exit()

    ref_seq = os.path.abspath(sys.argv[2])
    listFile = os.path.abspath(sys.argv[1])
    pos_hash = {}
    # position list
    with open(listFile, 'r') as f:
        for line in f:
            if re.search(r'^\s*$', line):
                continue
            else:
                arrs = line.strip().split("\t")
                if len(arrs) == 0:
                    continue
                if len(arrs) > 1:
                    pos_hash.setdefault(arrs[0],[]).append(arrs[1:])
                else:
                    pos_hash.setdefault(arrs[0], [])
    
    total_s_number = len(pos_hash.keys())
    tmp_number = 0
    # get subseq
    for seq_record in SeqIO.parse(ref_seq, "fasta"):
        seq_id = seq_record.id
        seq_desc = seq_record.description.replace(seq_record.id, "").strip()
        seq_seq = seq_record.seq
        seq_length = len(seq_seq)
        if tmp_number >= total_s_number:
            break
        if seq_id in pos_hash:
            tmp_number += 1
            if len(pos_hash[seq_id]) == 0:
                print(">{} {}\n{}".format(seq_id, seq_desc, seq_seq))  # 整个序列
                continue
            else:
                for tmp_list in pos_hash[seq_id]:
                    if len(tmp_list) == 3:
                        start = int(tmp_list[0])
                        end = int(tmp_list[1])
                        sub_seq_name = tmp_list[2]
                        strand = 1
                        if start > end:
                            strand = -1
                            start, end = end, start
                        start = start - 1
                        start = 0 if start < 0 else start
                        end = seq_length if end > seq_length else end
                        sub_seq = seq_seq[start:end]
                        if strand == -1:
                            sub_seq = sub_seq.reverse_complement()
                        #print("{}\t{}\t{}\t{}\t{}".format(seq_id, start+1, end, sub_seq_name, sub_seq))
                        print(">{}\n{}".format(sub_seq_name, sub_seq))

