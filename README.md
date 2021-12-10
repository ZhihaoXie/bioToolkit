# bio_scripts instruction

tag： bioinformatics

---

Some scripts for bioinformatics.

# 安装依赖

- perl (>= v5.10)
  - bioperl
- python3 (>= v3.6)
  - biopython >=1.78
  - pyfastx >=0.8
  - bcbio-gff >=0.6.6
  - loguru

备注，在使用这些脚本之前，请先安装Perl、Python3。

# 使用说明

## biotoolkit.py

查看说明： `python biotoolkit -h`

```
usage: biotoolkit [-h] [-v]
                  {gbk2fa,gbk2gff,fq2fa,fa2bed,lengthStats,stats,cutHead10Fastq,getGeneFromGFF3,getGeneFromGBK,gbkGetGeneRegionByName,geneRegion,chooseseq,search}
                  ...

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit

subcommands:
  {gbk2fa,gbk2gff,fq2fa,fa2bed,lengthStats,stats,cutHead10Fastq,getGeneFromGFF3,getGeneFromGBK,gbkGetGeneRegionByName,geneRegion,chooseseq,search}
                        Desired action to perform
    gbk2fa              genbank convert to fasta format. gzip format is not
                        supported.
    gbk2gff             genbank convert to gff3
    fq2fa               fastq convert to fasta format. gzip format is
                        supported.
    fa2bed              fasta file convert to bed file
    lengthStats         length stats of fasta file. output to stdout
    stats               stats of fasta file. output to stdout
    cutHead10Fastq      trim 10 bp from head ends for fastq
    getGeneFromGFF3     get gene by CDS from gff3 and genome
    getGeneFromGBK      extract protein or nucleotide seq of gene from genbank
                        file. gzip format is not supported.
    gbkGetGeneRegionByName (geneRegion)
                        get sequence region between two genes
    chooseseq           choose sequences by sequence's id
    search              search the location of sub-seq in the genome. print to
                        stdout
```


## gbk-summary.pl

gbk-summary.pl 用于统计GBK文件中基因个数、基因平均长度、rRNA和tRNA数量等。

```
用法： perl gbk-summary.pl <gbk.file>  [> out.file]
```


## download_kegg_picture.pl

当你用 ko number 做完 pathway mapping，你可能想要下载 pathway 图，那么这个脚本将会帮助你。

用法： `perl download_kegg_picture.pl -i mapid_file -u url -o out_dir`

mapid_file 是包含 ko number 的文件，如 ko00710

URL 类似这样： http://www.genome.jp/kegg-bin/show_pathway?144541224825059

out_dir 则是图片输出结果的目录。

拓展： [KEGG](http://www.genome.jp/kegg/)


## calc_SNP_Coregene.pl

统计多物种多序列比对后不同物种与参考物种的核心基因的SNP数目。也适用于其他多核酸序列比对后统计SNP数目。多序列比对建议使用软件muscle，比对结果是比对后的fasta格式，以".mus"为后缀。

fasta seq ID 格式为"species_id|gene_id"，如："E.coli|gene1"。list file只包含species_id，且每一个species_id为一行，参考序列的species_id须在第一行。

用法： `perl calc_SNP_Coregene.pl <list> <in_dir> > output`


## split_fasta.py

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


## name2taxid.py

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


## run_multitask.py

Execute shell scripts in parallel

```
Usage: python run_multitask.py <bash_cmd.sh> [task_number]
```

## multiprocessing_run_cmd.py

Execute shell scripts in parallel, each line is a command.

```
Usage: python multiprocessing_run_cmd.py cmd.sh [threads_num]
```

