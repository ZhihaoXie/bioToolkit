#!/usr/bin/perl -w

########################################################################################
# Usage: summary the number of gene,length of gene,etc for every sequence in gbk file. #
# Copyright: Zhihao Xie      \\(#^0^)/                                                 #
########################################################################################

use Bio::SeqIO;
use Bio::SeqFeature::Generic;
use File::Basename;
my $Usage = "\n Usage: perl $0 <gbk.file>  [> out.file]\n Copyright: Zhihao Xie              \\(#^0^)/\n";
if (@ARGV < 1) {
    print STDERR "$Usage\n";
    exit;
}
my $gbk_file = shift;
my $fgbk = basename $gbk_file;

print "Id\tnumber_of_genes\tlength_of_genes\tavg_length_genes\tCDS_ratio\tGC_of_genes\tnumber_of_rrna\tnumber_of_trna\n";

my ($all_gene,$all_length,$all_ave,$all_rrna,$all_trna,$all_gc,$gc);
$all_rrna=0;$all_trna=0;
my $sall_length;

my $seq_obj=Bio::SeqIO->new(-file=>"< $gbk_file",-format=>'genbank');
while (my $seq=$seq_obj->next_seq()) 
{
    my $num_cds=0;my $length_cds=0;my $num_rrna=0;my $num_trna=0;my $agv_cds=0;my $gc_num=0;my $gc_cds=0;
    my $id = $seq->id;
    my $seq_length = $seq->length;
    $sall_length += $seq_length;
    my @features = $seq->all_SeqFeatures;
    for my $f (@features) {
        if ($f->primary_tag eq 'CDS') {
            $num_cds += 1;
            my $l_cds = $f->length;
            $length_cds += $l_cds;
            $all_gene += 1;
            $all_length += $l_cds;
            my $start = $f->start;
            my $end = $f->end;
            my $str = $seq->subseq($start,$end);
            my $number = $str =~ tr/GCgc/GCgc/;
            $gc_num += $number;
            $all_gc += $number;
        }
        if ($f->primary_tag eq 'rRNA') {
            $num_rrna += 1;
            $all_rrna += 1;
        }
        if ($f->primary_tag eq 'tRNA') {
            $num_trna += 1;
            $all_trna += 1;
        }
    }
    if ($num_cds > 0) {
    $agv_cds = $length_cds/$num_cds;
    }
    if ($length_cds > 0) {
    $gc_cds = $gc_num/$length_cds*100; 
    }
    my $ratio = $length_cds/$seq_length*100;
    print "$id\t$num_cds\t$length_cds\t";
    printf "%0.2f\t%0.2f%%\t%0.2f%%", $agv_cds,$ratio,$gc_cds;
    print "\t$num_rrna\t$num_trna\n";
}
if ($all_gene > 0) {
    $all_ave = $all_length/$all_gene;
    $gc = $all_gc/$all_length*100;
}
my $all_ratio = $all_length/$sall_length*100;
print "#SUM:\n";
print "$fgbk\t$all_gene\t$all_length\t";
printf "%0.2f\t%0.2f%%\t%0.2f%%", $all_ave,$all_ratio,$gc;
print "\t$all_rrna\t$all_trna\n";
