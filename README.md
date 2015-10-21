# bio_scripts instruction

标签： git_instruction

---

some small tools for bioinformatics

在使用这些脚本工具前，请先安装Perl和bioperl模块。

## choosefasta.pl
用法： `perl choosefasta.pl <list> <fasta_file> [> <output>]` 

说明：

1. list file
list file是以tab符分割的csv文件或文本文件；第一列是源序列的ID，第二列是起始位置，第三列是终止位置，第四列是取出的子序列的ID（自定义）。如果只有第一列则是取出整个序列。如：
scaff1  1   100 scaff_1
scaff2  25  500 scaff_2
scaff3
2. fasta_file
fasta seq文件是标准的fasta序列格式.

[fasta wiki](https://en.wikipedia.org/wiki/FASTA_format)


## choosefastq.pl
用法： `perl choosefastq.pl <list> <fastq> [> out]` 

说明：
list file是所要提取fastq序列的ID（@可有可无）。
fastq 文件是标准的fastq序列格式.

[fastq wiki](https://en.wikipedia.org/wiki/FASTQ_format)


## gbk-summary.pl
gbk-summary.pl 用于统计GBK文件中基因个数、基因平均长度、rRNA和tRNA数量等。

用法： `perl gbk-summary.pl <gbk.file>  [> out.file]`


## get_geneseq_fromGBK.pl
get_geneseq_fromGBK.pl从GBK文件中提取基因的核苷酸序列和蛋白序列。

用法： `perl get_geneseq_fromGBK.pl genbank_file prefix`


## download_kegg_picture.pl
当你用 ko number 做完 pathway mapping，你可能想要下载 pathway 图，那么这个脚本将会帮助你。

用法： `perl download_kegg_picture.pl -i mapid_file -u url -o out_dir`

mapid_file 是包含 ko nukber 的文件，如 ko00710

URL 可能是这样： http://www.genome.jp/kegg-bin/show_pathway?144541224825059

out_dir 则是图片输出结果的目录。

拓展： [KEGG](http://www.genome.jp/kegg/)
