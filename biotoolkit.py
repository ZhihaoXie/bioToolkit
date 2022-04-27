#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
'''
@Filename:    biotoolkit.py
@Version:     V0.0.0.1
@Date:        2021/12/02 14:09:51
@Description: biotoolkit
'''

import sys
import os
import re
import argparse
import gzip
import pyfastx
from collections import OrderedDict
from Bio import SeqIO  # required biopython >=1.78
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def gbk2fa(args):
    '''genbank convert to fasta format
    '''
    gbk, fa = args.gbk, args.fasta
    SeqIO.convert(gbk, 'genbank', fa, 'fasta')


def gbk2gff(args):
    '''genbank convert to gff3
    '''
    from BCBio import GFF
    gbk, gff3 = args.gbk, args.gff3
    with open(gbk, 'r') as f1, open(gff3, 'w') as f2:
        GFF.write(SeqIO.parse(f1, "genbank"), f2)


def fq2fa(args):
    '''fastq convert to fasta format
    '''
    fq, fa = args.fastq, args.fasta
    if fq.endswith('.gz') and fa.endswith('.gz'):
        try:
            fq_handle = gzip.open(fq, 'rt')
            fa_handle = gzip.open(fa, 'wt')
        except Exception as e:
            sys.stderr.write("Warning: %s\n" (str(e)))
        else:
            for rec in SeqIO.parse(fq_handle, 'fastq'):
                SeqIO.write(rec, fa_handle, 'fasta')
        finally:
            fq_handle.close()
            fa_handle.close()
    elif fq.endswith('.gz') and not fa.endswith('.gz'):
        try:
            fq_handle = gzip.open(fq, 'rt')
            fa_handle = open(fa, 'w')
        except Exception as e:
            sys.stderr.write("Warning: %s\n" (str(e)))
        else:
            for rec in SeqIO.parse(fq_handle, 'fastq'):
                SeqIO.write(rec, fa_handle, 'fasta')
        finally:
            fq_handle.close()
            fa_handle.close()
    elif not fq.endswith('.gz') and fa.endswith('.gz'):
        try:
            fq_handle = open(fq, 'r')
            fa_handle = gzip.open(fa, 'wt')
        except Exception as e:
            sys.stderr.write("Warning: %s\n" (str(e)))
        else:
            for rec in SeqIO.parse(fq_handle, 'fastq'):
                SeqIO.write(rec, fa_handle, 'fasta')
        finally:
            fq_handle.close()
            fa_handle.close()
    else:
        SeqIO.convert(fq, 'fastq', fa, 'fasta')


def fa2bed(args):
    '''fasta file convert to bed file
    '''
    fa, bed_file = args.fasta, args.bed
    if fa.endswith('.gz'):
        fhandle = gzip.open(fa, 'rt')
    else:
        fhandle = open(fa, 'r')

    with open(bed_file, 'w') as f:
        for rec in SeqIO.parse(fhandle, 'fasta'):
            rec_id = rec.id
            rec_len = len(rec)
            f.write(f"{rec_id}\t0\t{rec_len}\n")

    fhandle.close()


def lengthStats(args):
    '''length stats of fasta file
    '''
    fa = args.fasta
    try:
        if fa.endswith('.gz'):
            fhandle = gzip.open(fa, 'rt')
        else:
            fhandle = open(fa, 'r')
    except Exception as e:
        sys.stderr.write(str(e))
    else:
        for rec in SeqIO.parse(fhandle, 'fasta'):
            rec_id = rec.id
            rec_len = len(rec)
            tmpseq = str(rec.seq)
            nNumber = tmpseq.count('N')
            nNumber += tmpseq.count('n')
            mylnoN = rec_len - nNumber
            print(f"{rec_id}\t{rec_len}\t{mylnoN}")
    finally:
        fhandle.close()


def stats(args):
    '''stats of fasta file
    '''
    fa = args.fasta
    num_seqs, sum_len, min_len, avg_len, max_len, total_gc, gap_num, gap_len, N50_num, N50 = 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
    try:
        if fa.endswith('.gz'):
            fhandle = gzip.open(fa, 'rt')
        else:
            fhandle = open(fa, 'r')
    except Exception as e:
        sys.stderr.write(str(e))
    else:
        len_list = []
        for rec in SeqIO.parse(fhandle, 'fasta'):
            rec_id = rec.id
            rec_len = len(rec.seq)
            len_list.append(rec_len)
            num_seqs += 1
            sum_len += rec_len
            if min_len == 0:  # init min_len
                min_len = rec_len
            elif min_len > rec_len:
                min_len = rec_len
            if max_len < rec_len:
                max_len = rec_len
            g_num = rec.seq.count('G')
            c_num = rec.seq.count('C')
            s_num = rec.seq.count('S')
            total_gc += g_num + c_num + s_num
            n_list = re.findall('N+', str(rec.seq), re.I)
            gap_num += len(n_list)
            gap_len += sum(list(map(lambda x: len(x), n_list)))

        avg_len = round(sum_len / num_seqs, 2)
        gc_ratio = round(total_gc / sum_len, 2)
        len_list.sort(reverse=True)
        temp_n50 = 0
        for i, v in enumerate(len_list):
            temp_n50 += v
            if temp_n50 >= 0.5 * sum_len:
                N50 = v
                N50_num = i + 1
                break

        print("#file\tnum_seqs\tsum_len\tmin_len\tavg_len\tmax_len\tgc\tgap_num\tgap_len\tN50_num\tN50")
        print(
            f"{os.path.basename(fa)}\t{num_seqs}\t{sum_len}\t{min_len}\t{avg_len}\t{max_len}\t{gc_ratio}\t{gap_num}\t{gap_len}\t{N50_num}\t{N50}"
        )
    finally:
        fhandle.close()


def cutHead10Fastq(args):
    '''剪去FASTQ序列前10bp
    '''
    fq, new_fq, head = args.fastq, args.out, args.trim_num
    if not new_fq.endswith(".gz"):
        sys.stderr.write(f"Warning: {new_fq} not endswith .gz, will add .gz for it.\n")
        new_fq += ".gz"
    with gzip.open(new_fq, 'wt') as fout:
        for name, seq, qual in pyfastx.Fastq(fq, build_index=False):
            new_seq = seq[head:]
            new_qual = qual[head:]
            fout.write(f"@{name}\n{new_seq}\n+\n{new_qual}\n")


def getGeneFromGFF3(args):
    '''get gene from gff3 and genome
    '''
    gff3, genomeFasta, cds_fna = args.gff3, args.genome, args.cds
    gene_dict = {}
    with open(gff3, "r") as fhd:
        for line in fhd:
            if re.search(r'^##|^#|^\s*$', line):
                continue
            cols = line.rstrip("\n").split("\t")
            seq_id, feature_type, start_site, end_site, strand = cols[0], cols[2], int(cols[3]), int(cols[4]), cols[6]
            # only CDS
            if feature_type == "CDS":
                if start_site > end_site:
                    start_site, end_site = end_site, start_site

                cds_infos = cols[8].split(';')
                locus_tag, cds_name, cds_nameID, cds_id = '', '', '', ''
                tmp_list = []
                for i in cds_infos:
                    if re.search(r'locus_tag=', i):
                        locus_tag = re.findall('locus_tag=(.+)', i)[0]
                    if re.search(r'ID=', i):
                        cds_nameID = re.findall('ID=(.+)', i)[0].replace('cds-', '')
                    if re.search(r'Name=', i):
                        cds_name = re.findall('Name=(.+)', i)[0]
                    if re.search(r'gene=', i):
                        cds_genesymbol = i
                        tmp_list.append(cds_genesymbol)
                    if re.search(r'product=', i):
                        cds_product = i
                        tmp_list.append(cds_product)
                if locus_tag:
                    cds_id += locus_tag + ":"
                if cds_name:
                    cds_id += cds_name
                elif cds_nameID:
                    cds_id += cds_nameID
                cds_id = cds_id.strip().strip(':')
                cds_desc = " ".join(tmp_list)  # as cds description
                # 要注意，存在cds名称一样但基因位置不一的CDS，因此要保存所有重复的CDS
                gene_dict.setdefault(seq_id, OrderedDict()).setdefault(cds_id, []).append(
                    (start_site, end_site, strand, cds_desc))

    with open(cds_fna, 'w') as fout:
        for seq_record in SeqIO.parse(genomeFasta, "fasta"):
            seq_seq = seq_record.seq
            seq_id = seq_record.id
            if seq_id in gene_dict:
                for cds_id, tmp_list in gene_dict[seq_id].items():
                    if len(tmp_list) == 1:
                        for start_site, end_site, strand, cds_desc in tmp_list:
                            start_site = start_site - 1
                            start_site = 0 if start_site < 0 else start_site
                            end_site = len(seq_seq) if end_site > len(seq_seq) else end_site
                            sub_seq = seq_seq[start_site:end_site]  # sub_seq is Seq type
                            if strand == "-" or strand == "-1":
                                sub_seq = sub_seq.reverse_complement()
                            sub_seq_id = f"{seq_id}:{cds_id}"
                            sub_record = SeqRecord(sub_seq, id=sub_seq_id, description=cds_desc)
                            SeqIO.write(sub_record, fout, 'fasta')
                    else:
                        for start_site, end_site, strand, cds_desc in tmp_list:
                            start_site = start_site - 1
                            start_site = 0 if start_site < 0 else start_site
                            end_site = len(seq_seq) if end_site > len(seq_seq) else end_site
                            sub_seq = seq_seq[start_site:end_site]  # sub_seq is Seq type
                            if strand == "-" or strand == "-1":
                                sub_seq = sub_seq.reverse_complement()
                            sub_seq_id = f"{seq_id}:{cds_id}:site={start_site+1}-{end_site}"
                            sub_record = SeqRecord(sub_seq, id=sub_seq_id, description=cds_desc)
                            SeqIO.write(sub_record, fout, 'fasta')


def gbkGetGeneRegionByName(args):
    '''get sequence region between two genes
    '''
    gbk, select_genes, region_fa = args.gbk, args.gene_names, args.region
    two_genes = select_genes.strip().lower().split(",")  # gene lower
    if len(two_genes) < 2:
        raise ValueError("gene string error. please provide two genes name, separated by commas")

    had_found = OrderedDict()  # 标注基因是否查找到
    gene_sites = {}
    gbk_dist = {}
    try:
        if gbk.endswith('.gz'):
            gbk_hd = gzip.open(gbk, 'rt')
        else:
            gbk_hd = open(gbk, 'r')
    except Exception as e:
        sys.stderr.write("Error, open gbk file failed: ", str(e))
    else:
        for gb_rec in SeqIO.parse(gbk_hd, 'genbank'):
            # genbank information
            gb_id, gb_name = gb_rec.id, gb_rec.name
            gb_desc = gb_rec.description
            gb_seq = gb_rec.seq
            had_found.setdefault(gb_id, set())
            for feature in gb_rec.features:
                if feature.type == "CDS":
                    gene_locus = feature.qualifiers.get("locus_tag", [])
                    gene_name = feature.qualifiers.get("gene", [])
                    if len(gene_name) == 0:
                        continue
                    gene_name = gene_name[0]
                    gene_locus = gene_locus[0]
                    if gene_name.lower() in two_genes:
                        gene_location = feature.location
                        gene_sites.setdefault(gb_id, {}).setdefault(gene_name, (gene_locus, gene_location))
                        gbk_dist.setdefault(gb_id, (gb_name, gb_desc, gb_seq))
                        had_found[gb_id].add(gene_name)
    finally:
        gbk_hd.close()

    with open(region_fa, 'w') as fout:
        for gb_id, genes in had_found.items():
            if len(genes) < 2:
                sys.stderr.write(f"Warning: found {genes} of {gb_id} in {gbk}\n")
                continue
            else:
                min_site, max_site = 0, 0
                for gg in gene_sites[gb_id]:
                    _, gs = gene_sites[gb_id][gg]
                    if int(gs.end) > max_site:
                        max_site = int(gs.end)
                    if min_site == 0:
                        min_site = int(gs.start)
                    elif int(gs.start) < min_site:
                        min_site = int(gs.start)
                sub_seq = gbk_dist[gb_id][2][min_site:max_site]
                select_genes = select_genes.replace(',', '-')
                fout.write(f">{select_genes}|location:{min_site+1}-{max_site}|ACCESSION:{gb_id}\n{sub_seq}\n")


def chooseseq(args):
    '''choose sequences by sequence's id
    '''
    id_list, seq_file, seq_type = args.ids, args.seq, args.seq_type
    out_file = args.out
    pos_hash = {}
    # read id file
    with open(id_list, 'r') as f:
        for line in f:
            if re.search(r'^\s*$', line):
                continue
            else:
                cols = line.strip().split("\t")
                if len(cols) == 0:
                    continue
                if len(cols) == 1:
                    pos_hash.setdefault(cols[0], [])
                elif len(cols) == 4:
                    pos_hash.setdefault(cols[0], []).append(cols[1:])
                else:
                    sys.stderr.write(
                        f"Warning: {line} was non-standard line format, must be: 'id start end new_id' or 'id', tab-separated.")
                    continue

    total_s_number = len(pos_hash.keys())
    tmp_number = 0
    with open(out_file, 'w') as fout:
        if seq_type == "fasta":
            for name, seq in pyfastx.Fasta(seq_file, build_index=False):
                if tmp_number >= total_s_number:
                    break
                if name in pos_hash:
                    tmp_number += 1
                    if len(pos_hash[name]) == 0:
                        new_seq = SeqRecord(Seq(seq), id=name, description="")
                        SeqIO.write(new_seq, fout, 'fasta')
                        continue
                    else:
                        for temp_list in pos_hash[name]:
                            start = int(temp_list[0])
                            end = int(temp_list[1])
                            sub_seq_name = temp_list[2]
                            strand = 1
                            if start > end:
                                strand = -1
                                start, end = end, start
                            start = start - 1
                            start = 0 if start < 0 else start
                            end = len(seq) if end > len(seq) else end
                            sub_seq = seq[start:end]
                            if strand == -1:
                                new_sub_seq = Seq(sub_seq)
                                sub_seq = str(new_sub_seq.reverse_complement())
                            new_seq = SeqRecord(Seq(sub_seq), id=sub_seq_name, description="")
                            SeqIO.write(new_seq, fout, 'fasta')
        elif seq_type == 'fastq':
            for name, seq, qual in pyfastx.Fastq(seq_file, build_index=False):
                if tmp_number >= total_s_number:
                    break
                if name in pos_hash:
                    tmp_number += 1
                    if len(pos_hash[name]) == 0:
                        fout.write(f"@{name}\n{seq}\n+\n{qual}\n")
                        continue
                    else:
                        for temp_list in pos_hash[name]:
                            start = int(temp_list[0])
                            end = int(temp_list[1])
                            sub_seq_name = temp_list[2]
                            strand = 1
                            if start > end:
                                strand = -1
                                start, end = end, start
                            start = start - 1
                            start = 0 if start < 0 else start
                            end = len(seq) if end > len(seq) else end
                            sub_seq = seq[start:end]
                            sub_qual = qual[start:end]
                            if strand == -1:
                                new_sub_seq = Seq(sub_seq)
                                sub_seq = str(new_sub_seq.reverse_complement())
                            fout.write(f"@{sub_seq_name}\n{sub_seq}\n+\n{sub_qual}\n")
        else:
            sys.stderr.write("Error: Unknown format for sequence. only support fasta/fastq of DNA.")
            sys.exit()


def search(args):
    '''search the location of sub-seq in the genome. print to stdout
    '''
    seq_file, pattern_match = args.seq, args.pattern
    pattern_match = pattern_match.upper()
    pattern_matchNew = ""
    for name, seq in pyfastx.Fasta(seq_file, build_index=False):
        temp_list = [m for m in re.finditer(pattern_match, str(seq).upper())]
        if len(temp_list) == 0:
            temp_seq = Seq(pattern_match)
            temp_seq = temp_seq.reverse_complement()  # reverse pattern
            pattern_matchNew = str(temp_seq)
            temp_list = [m for m in re.finditer(pattern_matchNew, str(seq).upper())]
        if len(temp_list) == 0:
            sys.stdout.write(f"{name}\t{pattern_match}\tNot-found!\n")
            continue
        if len(pattern_match) == 1:
            temp_list = [str(m.start() + 1) for m in temp_list]
            sys.stdout.write(f"{name}\t{pattern_match}\t" + ",".join(temp_list) + "\n")
        else:
            temp_list = [f"{m.start()+1}-{m.end()}" for m in temp_list]
            if pattern_matchNew:
                # 模式被反向互补了
                sys.stdout.write(f"{name}\t{pattern_match}(rc:{pattern_matchNew})\t" + ",".join(temp_list) + "\n")
            else:
                sys.stdout.write(f"{name}\t{pattern_match}\t" + ",".join(temp_list) + "\n")


def getGeneFromGBK(args):
    '''extract protein or nucleotide seq of gene from genbank file
    '''
    gbk, outfile, molecular_type = args.gbk, args.out, args.molecular
    with open(outfile, 'w') as fout:
        for gb_record in SeqIO.parse(gbk, 'genbank'):
            gb_id = gb_record.id
            for ele in gb_record.features:
                if ele.type == "CDS":
                    cds_seq = ele.extract(gb_record.seq)  # extract func return a Seq object
                    locus_tag = ele.qualifiers['locus_tag'][0]
                    gene_symbol = ele.qualifiers.get("gene", [])
                    if len(gene_symbol) > 0:
                        gene_symbol = gene_symbol[0]
                    else:
                        gene_symbol = ""
                    product = ele.qualifiers.get("product", [])
                    if len(product) > 0:
                        product = product[0]
                    else:
                        product = ""
                    cds_description = f"{gene_symbol} {product}".strip()
                    if molecular_type in ['DNA', 'dna']:
                        cds_id = f"{gb_id}_cds_{locus_tag}"
                        SeqIO.write(SeqRecord(cds_seq, id=cds_id, description=cds_description), fout, 'fasta')
                    elif molecular_type in ['protein', 'pro', 'p']:
                        protein_id = ele.qualifiers.get("protein_id", [])
                        if len(protein_id) > 0:
                            protein_id = protein_id[0]
                        else:
                            protein_id = ""
                        cds_id = f"{gb_id}_cds_{locus_tag}_{protein_id}".strip("_")
                        # 优先提取"translation"字符串
                        cds_protein = ele.qualifiers.get("translation", [])
                        if len(cds_protein) > 0:
                            cds_protein = cds_protein[0]
                        else:
                            transl_table = ele.qualifiers.get("transl_table", [])
                            if len(transl_table) == 0:
                                transl_table = 1  # default translate table is 1
                            else:
                                transl_table = int(transl_table[0])
                            if len(cds_seq) % 3 != 0:
                                sys.stderr.write(
                                    f"Warning: {locus_tag}'s sequence length {len(cds_seq)} is not a multiple of three, maybe pseudogene!\n"
                                )
                                continue
                            cds_protein = ele.extract(gb_record.seq).translate(table=transl_table, to_stop=True)
                            cds_protein = str(cds_protein)
                        SeqIO.write(SeqRecord(Seq(cds_protein), id=cds_id, description=cds_description), fout, 'fasta')


def geneStats(args):
    """gene stats of genbank
    """
    from Bio.SeqUtils import GC, GC_skew

    def gene_stats_gbk(gbk_file, stats_out):
        genome_name = os.path.basename(gbk_file)
        genome_len, gene_len, gene_num, cds_len, cds_num, tRNA_num, rRNA_num = 0, 0, 0, 0, 0, 0, 0
        gene_ratio, gene_density = 0, 0
        genome_gc, genome_gcskew = 0, 0
        rec_num = 0
        for rec in SeqIO.parse(gbk_file, "genbank"):  # every fragment
            rec_num += 1
            genome_len += len(rec)
            genome_gc += GC(rec.seq)
            genome_gcskew += GC_skew(rec.seq, window=len(rec))[0]
            for seq_feature in rec.features:
                if seq_feature.type == "gene":
                    gene_num += 1
                    gene_seq = seq_feature.extract(rec.seq)
                    gene_len += len(gene_seq)
                if seq_feature.type == "CDS":
                    cds_num += 1
                    cds_seq = seq_feature.extract(rec.seq)
                    cds_len += len(cds_seq)
                if seq_feature.type == "tRNA":
                    tRNA_num += 1
                if seq_feature.type == "rRNA":
                    rRNA_num += 1

        if gene_num == 0:
            gene_num = cds_num + tRNA_num + rRNA_num
            gene_len = cds_len
        gene_ratio = round(gene_len / genome_len * 100, 2)
        gene_density = round(gene_num / genome_len * 1000, 2)  # 每1Kb的基因密度
        genome_gc = round(genome_gc / rec_num, 2)
        genome_gcskew = round(genome_gcskew / rec_num, 2)
        # output header
        if os.path.isfile(stats_out) and os.path.getsize(stats_out) >= 107:
            with open(stats_out, 'a', encoding='utf-8') as fout:  # 追加
                print(
                    f"{genome_name}\t{genome_len}\t{gene_len}\t{gene_num}\t{cds_len}\t{cds_num}\t{tRNA_num}\t{rRNA_num}\t{gene_ratio}\t{gene_density}\t{genome_gc}\t{genome_gcskew}",
                    file=fout)
        else:
            with open(stats_out, 'w', encoding='utf-8') as fout:
                print(
                    "#Genome\tGenome_len\tGene_len\tGene_num\tCDS_len\tCDS_num\ttRNA_num\trRNA_num\tGene_ratio\tGene_density\tGC\tGC_skew",
                    file=fout)
                print(
                    f"{genome_name}\t{genome_len}\t{gene_len}\t{gene_num}\t{cds_len}\t{cds_num}\t{tRNA_num}\t{rRNA_num}\t{gene_ratio}\t{gene_density}\t{genome_gc}\t{genome_gcskew}",
                    file=fout)

    gbk, gbk_list = "", ""
    if args.genbank:
        gbk = os.path.abspath(args.genbank)
    if args.list:
        gbk_list = os.path.abspath(args.list)
    output_file = os.path.abspath(args.output)
    if os.path.isfile(output_file) and os.path.getsize(output_file) > 0:
        os.remove(output_file)
    if os.path.isfile(gbk) and os.path.getsize(gbk) > 0:
        gene_stats_gbk(gbk, output_file)
    elif os.path.isfile(gbk_list) and os.path.getsize(gbk_list) > 0:
        with open(gbk_list, 'r') as f:
            for line in f:
                if re.search(r'^\s*$|^#', line):
                    continue
                temp_gbk = line.strip()
                gene_stats_gbk(temp_gbk, output_file)
    else:
        sys.stderr.write("Please provide a genbank file or a list of genbank.\n")
        sys.exit()


def main():
    parser = argparse.ArgumentParser(prog='biotoolkit')
    parser.add_argument('-v', '--version', action='version', version=__doc__)
    subparsers = parser.add_subparsers(title='subcommands', help='Desired action to perform')
    # 添加子命令
    parser_gbk2fa = subparsers.add_parser('gbk2fa', help="genbank convert to fasta format. gzip format is not supported.")
    parser_gbk2fa.add_argument('-g', dest='gbk', required=True, help="input of genbank format")
    parser_gbk2fa.add_argument('-f', dest='fasta', required=True, help="output of fasta format")
    parser_gbk2fa.set_defaults(func=gbk2fa)

    parser_gbk2gff = subparsers.add_parser('gbk2gff', help="genbank convert to gff3")
    parser_gbk2gff.add_argument('-g', dest='gbk', required=True, help="input of genbank format")
    parser_gbk2gff.add_argument('-f', dest='gff3', required=True, help="output of gff3 format")
    parser_gbk2gff.set_defaults(func=gbk2gff)

    parser_fq2fa = subparsers.add_parser('fq2fa', help="fastq convert to fasta format. gzip format is supported.")
    parser_fq2fa.add_argument('-q', dest='fastq', required=True, help="input of fastq format")
    parser_fq2fa.add_argument('-f', dest='fasta', required=True, help="output of fasta format")
    parser_fq2fa.set_defaults(func=fq2fa)

    parser_fa2bed = subparsers.add_parser('fa2bed', help="fasta file convert to bed file")
    parser_fa2bed.add_argument('-f', dest='fasta', required=True, help="input of fasta format")
    parser_fa2bed.add_argument('-b', dest='bed', required=True, help="output of bed format")
    parser_fa2bed.set_defaults(func=fa2bed)

    parser_lengthStats = subparsers.add_parser('lengthStats', help="length stats of fasta file. output to stdout")
    parser_lengthStats.add_argument('-f', dest='fasta', required=True, help='input of fasta format')
    parser_lengthStats.set_defaults(func=lengthStats)

    parser_stats = subparsers.add_parser('stats', help="stats of fasta file. output to stdout")
    parser_stats.add_argument('-f', dest='fasta', required=True, help="input of fasta format")
    parser_stats.set_defaults(func=stats)

    parser_cutHead10Fastq = subparsers.add_parser('cutHead10Fastq', help="trim 10 bp from head ends for fastq")
    parser_cutHead10Fastq.add_argument('-q', dest='fastq', required=True, help="input of fastq")
    parser_cutHead10Fastq.add_argument('-o', dest='out', required=True, help="output for fastq, which had trim")
    parser_cutHead10Fastq.add_argument('-n', dest='trim_num', default=10, type=int, help="trim number from head ends")
    parser_cutHead10Fastq.set_defaults(func=cutHead10Fastq)

    parser_getgenegff3 = subparsers.add_parser('getGeneFromGFF3', help="get gene by CDS from gff3 and genome")
    parser_getgenegff3.add_argument('-f', dest='gff3', required=True, help="gff3 file")
    parser_getgenegff3.add_argument('-g', dest='genome', required=True, help='genome, fasta format')
    parser_getgenegff3.add_argument('-c', dest='cds', required=True, help='gene fasta as output')
    parser_getgenegff3.set_defaults(func=getGeneFromGFF3)

    parser_getGeneFromGBK = subparsers.add_parser(
        'getGeneFromGBK', help="extract protein or nucleotide seq of gene from genbank file. gzip format is not supported.")
    parser_getGeneFromGBK.add_argument('-g', dest='gbk', required=True, help="input of genbank file")
    parser_getGeneFromGBK.add_argument('-o', dest="out", required=True, help="ouput of gene")
    parser_getGeneFromGBK.add_argument('-m',
                                       dest="molecular",
                                       default='dna',
                                       help="molecular type, dna or protein, default is %(default)s")
    parser_getGeneFromGBK.set_defaults(func=getGeneFromGBK)

    parser_gbkGetGeneRegionByName = subparsers.add_parser('gbkGetGeneRegionByName',
                                                          aliases=['geneRegion'],
                                                          help="get sequence region between two genes")
    parser_gbkGetGeneRegionByName.add_argument('-g', dest='gbk', required=True, help="input of genbank")
    parser_gbkGetGeneRegionByName.add_argument('-n', dest='gene_names', required=True, help="two genes name, separated by commas")
    parser_gbkGetGeneRegionByName.add_argument('-o', dest='region', required=True, help="output for selected gene region")
    parser_gbkGetGeneRegionByName.set_defaults(func=gbkGetGeneRegionByName)

    parser_chooseseq = subparsers.add_parser('chooseseq', help="choose sequences by sequence's id")
    parser_chooseseq.add_argument('-i',
                                  dest="ids",
                                  required=True,
                                  help="select id file, one-line must be: 'id start end new_id' or 'id', tab-separated")
    parser_chooseseq.add_argument('-s',
                                  dest='seq',
                                  required=True,
                                  help="input of fasta or fastq of DNA. only support SE reads for fastq")
    parser_chooseseq.add_argument('-o', dest="out", required=True, help="ouput, gzip format is not supported")
    parser_chooseseq.add_argument('-f',
                                  dest="seq_type",
                                  default='fasta',
                                  help="seq format, support fasta|fastq, default is fasta")
    parser_chooseseq.set_defaults(func=chooseseq)

    parser_search = subparsers.add_parser('search', help="search the location of sub-seq in the genome. print to stdout")
    parser_search.add_argument('-s', dest="seq", required=True, help="input of fasta format")
    parser_search.add_argument('-p', dest='pattern', required=True, help="pattern for search")
    parser_search.set_defaults(func=search)

    parser_geneStats = subparsers.add_parser('geneStats', help="gene stats of genbank")
    parser_geneStats.add_argument('-g', dest='genbank', help="genbank file for stats")
    parser_geneStats.add_argument('-l', dest='list', help="a list of genbank files, per line a genbank file")
    parser_geneStats.add_argument('-o', dest="output", help="output of gene stats")
    parser_geneStats.set_defaults(func=geneStats)

    # print help
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit()

    # parse args and run subcmd
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
