# biotoolkit instruction

tag： bioinformatics

---

Some scripts for bioinformatics.

# Dependencies

- perl (>= v5.10)
  - bioperl
- python3 (>= v3.6)
  - biopython >=1.78
  - pyfastx >=0.8
  - bcbio-gff >=0.6.6
  - loguru

 Note: Before using these scripts, please install Perl and Python 3. 

# Usage

## biotoolkit.py

Usage: `python biotoolkit.py -h`

```
usage: biotoolkit [-h] [-v]
                  {gbk2fa,gbk2gff,fq2fa,fa2bed,lengthStats,stats,cutHead10Fastq,getGeneFromGFF3,getGeneFromGBK,gbkGetGeneRegionByName,geneRegion,chooseseq,search,geneStats,translate}
                  ...

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit

subcommands:
  {gbk2fa,gbk2gff,fq2fa,fa2bed,lengthStats,stats,cutHead10Fastq,getGeneFromGFF3,getGeneFromGBK,gbkGetGeneRegionByName,geneRegion,chooseseq,search,geneStats,translate}
                        Desired action to perform
    gbk2fa              genbank convert to fasta format. gzip format is not
                        supported.
    gbk2gff             genbank convert to gff3
    fq2fa               fastq convert to fasta format. gzip format is
                        supported.
    fa2bed              fasta file convert to bed file
    lengthStats         length stats of fasta file. output to stdout
    stats               stats of fasta file, include total length, average
                        length, gc, N50. output to stdout
    cutHead10Fastq      trim 10 bp from head ends for fastq
    getGeneFromGFF3     get gene by CDS from gff3 and genome
    getGeneFromGBK      extract protein or nucleotide seq of gene from genbank
                        file. gzip format is not supported.
    gbkGetGeneRegionByName (geneRegion)
                        get sequence region between two genes
    chooseseq           choose sequences by sequence's id
    search              search the location of sub-seq in the genome. print to
                        stdout
    geneStats           gene stats of genbank
    translate           translate cds to protein
```

The subcommand can obtain parameter help information through -h.

Example:

```
$ python biotoolkit.py translate -h
usage: biotoolkit translate [-h] [-i CDS] [-o PROTEIN] [-c CODING_TABLE]

optional arguments:
  -h, --help       show this help message and exit
  -i CDS           cds sequences, fasta format
  -o PROTEIN       protein sequences as output, fasta format
  -c CODING_TABLE  codon table, default is 1
  
$ python biotoolkit.py translate -i ./sequence.fasta -o p.fa -c 11
```

## Other scripts

### `scripts/get_metadata_cif.py`

批量分析mmCIF文件获取蛋白质元数据并输出表格

Usage：

```
# 分析单个文件
python get_metadata_cif.py 1ake.cif.gz -o results.csv
  
# 分析多个文件
python get_metadata_cif.py file1.cif file2.cif.gz file3.cif -o summary.xlsx

# 分析整个目录（递归）
python get_metadata_cif.py /path/to/pdb/structures -o pdb_catalog.csv

# 分析当前目录（不递归）
python get_metadata_cif.py . --no-recursive -o local.csv
```

### download_kegg_picture.pl

When you use ko number to path mapping for KEGG, you may want to download the path map. This script will help you.
Usage: ` perl download_kegg_picture.pl -i mapid_file -u url -o out_dir`
- mapid_ File is a file containing ko number, such as ko00710
- the URL looks like this:< http://www.genome.jp/kegg-bin/show_pathway?144541224825059 >
out_ Dir is the directory of image output results.

expand: [KEGG](http://www.genome.jp/kegg/)

### calc_SNP_Coregene.pl

Count the number of SNP of core genes of different species and reference species after multi-sequence alignment of multi species. It is also applicable to statistics of SNP number after sequences alignment of other polynucleic acids. It is recommended to use muscle software for multi-sequences alignment. The alignment result is in the fasta format after alignment, with ". mus" as the suffix.

Fasta format seq ID is "specifics_id | gene_id", such as "E.coli | gene1". 

The list file contains only species id, and each ID is a line, and sequence ID of reference must be in the first line.

Usage： `perl calc_SNP_Coregene.pl <list> <in_dir> > output`

### split_fasta.py

Split a multiFASTA file by number or chunks or file size

```
python split_fasta.py -h
arguments:
  -h, --help            show this help message and exit
  -i FASTA_FILE, --input-file FASTA_FILE
                        Input FASTA file to split, support for gzip
                        compression
  -n NUMBER, --num NUMBER
                        Number of sequences of every subfile (default: None)
  -c NUM_CHUNKS, --chunks NUM_CHUNKS
                        Number of chunks to split (default: None)
  -m MAX_FILE_SIZE, --max_file_size MAX_FILE_SIZE
                        Instead of a precise number, indicate an expected file
                        size for each chunk in bytes(k, M, G, T suffixes are
                        accepted). (default: None)
  -o OUTPUT_DIR, --output-dir OUTPUT_DIR
                        Output directory (default: .)
```

### name2taxid.py

Query taxonomy id by scientific name from names.dmp

```
usage: name2taxid.py [-h] [-n NAMEDMP] [-q QUERY] [-qf QUERYFILE]

scientific name to taxonmy id

optional arguments:
  -h, --help     show this help message and exit
  -n NAMEDMP     names.dmp of taxonomy as query database
  -q QUERY       query sci name
  -qf QUERYFILE  query sci name file, one line as a name, -q and -qf conflict
```

### run_multitask.py

Execute shell scripts in parallel

```
Usage: python run_multitask.py <bash_cmd.sh> [task_number]
```

### multiprocessing_run_cmd.py

Execute shell scripts in parallel, each line is a command.

```
Usage: python multiprocessing_run_cmd.py cmd.sh [threads_num]
```
