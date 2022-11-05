#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
'''
@Time    :   2019/11/28 16:45:02
@Author  :   Xie Zhihao 
@Version :   1.0.0
@Contact :   xiezhihao1122@outlook.com
@License :   (C)Copyright 2019, BSD License
@Desc    :   Split a multiFASTA file by number or chunks or file size
'''

import sys, os, re
import argparse
import gzip
from Bio import SeqIO

def is_gzip(myfile):
    '''determine a file is a gzip file?
    '''
    gzip_yes = False
    try:
        tmp_hd = open(myfile, 'rb')
    except:
        sys.stderr.write("Warning: can't read %s\n" % (myfile))
    else:
        tmp_b2b = tmp_hd.read(2)
        if tmp_b2b == b'\x1f\x8b':
            gzip_yes = True
        else:
            gzip_yes = False
        tmp_hd.close()
    return gzip_yes

def num_with_si_suffix(num):
    unit_to_value = {
            'k': 2**10,    # kilo
            'M': 2**20,    # mega
            'G': 2**30,    # giga
            'T': 2**40,   # tera
            }

    if num[-1].isdigit():
        value=int(num)
    else:
        value=int(num[0:-1])
        suffix=num[-1]
        if suffix not in unit_to_value.keys():
            msg = "unknown suffix '{0}'".format(suffix)
            raise argparse.ArgumentTypeError(msg)
        value *= unit_to_value[suffix]

    return value   

def count_entries(your_file):
    '''count total number of fasta file.
    return count
    '''
    count = 0
    if is_gzip(your_file):
        stream = gzip.open(your_file, 'rt')
    else:
        stream = open(your_file, 'r')
    for rec in SeqIO.parse(stream, 'fasta'):
        count += 1
    stream.close()
    return count

def split_by_chunks(fasta_file, chunk_size, output_dir):
    (fasta_file_basename, fasta_file_extension)  = os.path.splitext(os.path.basename(fasta_file))
    output_file_prefix = os.path.join(output_dir, fasta_file_basename)
    gz_yes = is_gzip(fasta_file)
    try:
        if gz_yes:
            f_hd = gzip.open(fasta_file, 'rt')
        else:
            f_hd = open(fasta_file, 'r')
        
        cur_chunk = 1
        entries_in_chunk = 0
        chunk_file = "{0}_{1}{2}".format(output_file_prefix, cur_chunk, fasta_file_extension)
        if gz_yes:
            chunk_stream = gzip.open(chunk_file, 'wt')
        else:
            chunk_stream = open(chunk_file, 'wt')
        for rec in SeqIO.parse(f_hd, 'fasta'):
            if (entries_in_chunk == chunk_size):
                # 关闭文件, 新打开文件
                cur_chunk = cur_chunk + 1
                entries_in_chunk = 0
                chunk_stream.close()
                chunk_file = "{0}_{1}{2}".format(output_file_prefix, cur_chunk, fasta_file_extension)
                if gz_yes:
                    chunk_stream = gzip.open(chunk_file, 'wt')
                else:
                    chunk_stream = open(chunk_file, 'wt')

            SeqIO.write(rec, chunk_stream, 'fasta')
            entries_in_chunk = entries_in_chunk + 1
        chunk_stream.close()
        f_hd.close()
    except IOError:
        sys.exit("Error cannot open {0}".format(fasta_file))

def split_by_size(fasta_file, max_file_size, output_dir):
    (fasta_file_basename, fasta_file_extension)  = os.path.splitext(os.path.basename(fasta_file))
    output_file_prefix = os.path.join(output_dir, fasta_file_basename)
    gz_yes = is_gzip(fasta_file)
    try:
        if gz_yes:
            f_hd = gzip.open(fasta_file, 'rt')
        else:
            f_hd = open(fasta_file, 'r')

        cur_chunk = 1
        file_size = 0
        chunk_file = "{0}_{1}{2}".format(output_file_prefix, cur_chunk, fasta_file_extension)
        if gz_yes:
            chunk_stream = gzip.open(chunk_file, 'wt')
        else:
            chunk_stream = open(chunk_file, 'wt')
        for rec in SeqIO.parse(f_hd, 'fasta'):
            if (file_size >= max_file_size):
                cur_chunk = cur_chunk + 1
                chunk_stream.close()
                del(chunk_file)
                file_size = 0
                chunk_file = "{0}_{1}{2}".format(output_file_prefix, cur_chunk,fasta_file_extension)
                if gz_yes:
                    chunk_stream = gzip.open(chunk_file, 'wt')
                else:
                    chunk_stream = open(chunk_file, 'wt')

            output = "{1}{0}{2}{0}".format(os.linesep, rec.id, rec.seq)
            SeqIO.write(rec, chunk_stream, 'fasta')
            file_size += len(output)
        chunk_stream.close()
        f_hd.close()
    except IOError:
        sys.exit("Error cannot open {0}".format(fasta_file))

def get_parameters():
    parser = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--input-file', dest='fasta_file', required=True,
            default=argparse.SUPPRESS, type=str,
            help='Input FASTA file to split, support for gzip compression')
    parser.add_argument('-n', '--num', dest='number', type=int,
            help="Number of sequences of every subfile")
    parser.add_argument('-c', '--chunks', dest='num_chunks', type=int,
            help='Number of chunks to split')
    parser.add_argument('-m', '--max_file_size', dest='max_file_size',
            type=num_with_si_suffix, help='Instead of a precise number,'
            ' indicate an expected file size for each chunk in bytes'
            '(k, M, G, T suffixes are accepted).')
    parser.add_argument('-o', '--output-dir', dest='output_dir', type=str,
            default=os.curdir + os.sep, help='Output directory')
    return parser.parse_args(), parser


# <---------- main ----------> #
if __name__ == "__main__":
    args, parser = get_parameters()
    if not args.number and not args.num_chunks and not args.max_file_size:
        print("Error: Please indicate the number of subfile or the number of chunks or the max file size to split.")
        sys.exit(parser.print_help())
    
    fasta_file = os.path.abspath(args.fasta_file)
    if not os.path.isfile(fasta_file):
        print("Error: {0} don't exist.".format(fasta_file))
        sys.exit(1)
    out_dir = os.path.abspath(args.output_dir)
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)
    
    if args.number and args.number > 0:
        print("Dividing {0} in subfiles which has {1} entries".format(fasta_file, args.number), flush=True)
        split_by_chunks(fasta_file, args.number, out_dir)
    elif args.num_chunks and args.num_chunks > 1:
        print("Start reading {0}".format(fasta_file), flush=True)
        num_entries = count_entries(fasta_file)
        print("{0} has {1} FASTA entries".format(fasta_file, num_entries), flush=True)
        chunk_size = int((num_entries + args.num_chunks -1)/(args.num_chunks))
        print("Dividing {0} in {1} chunks of {2} entries".format(fasta_file, args.num_chunks, chunk_size), flush=True)
        split_by_chunks(fasta_file, chunk_size, out_dir)
    elif args.max_file_size:
        print("Start create chunks", flush=True)
        split_by_size(fasta_file, args.max_file_size, out_dir)
    
    print("Done.")

