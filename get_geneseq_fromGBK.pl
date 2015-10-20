#!/usr/bin/perl -w 
# 
# Extract genes' protein and nucleotide seq from genbank file.
# Copyright: Zhihao Xie      2014.2.18
#
use Bio::SeqIO;
use Bio::SeqFeature::Generic;

my $Usage = 'This script will extract genes\' protein and nucleotide seq from genbank file!';
if (@ARGV < 2) {
	print STDERR "\t$Usage\n\tperl $0 genbank_file prefix\n";
	exit;
}

my $file = shift;
my $prefix = shift;
my $obj = Bio::SeqIO->new(-file=>"< $file",-format=>'genbank');
open OUT1,"> $prefix.faa";
open OUT2,"> $prefix.ffn";

while (my $seq = $obj->next_seq()) 
{
	my @feat = $seq->get_SeqFeatures('CDS');
	for $f (@feat) {
		unless ($f->has_tag('locus_tag')) {
			$i_n++;
			my $id = $f->seq_id;
			my $start = $f->start;
			my $end = $f->end;
			$f->add_tag_value('locus_tag',"$id-$start-$end-$i_n");
		}
		my @tag = $f->get_tag_values('locus_tag');
		my $start = $f->start;
		my $end = $f->end;
		my @pro;
		my $str = $seq->subseq($start,$end);
		if ($f->strand == -1) {
			my $m = Bio::Seq->new(-seq=>$str,-id=>'1');
			$m = $m->revcom;
			$str = $m->seq;
		}

		if ($f->has_tag('translation')) {
			@pro = $f->get_tag_values('translation');
		}
        my $protein=$pro[0];
		print OUT2 ">$tag[0]\n$str\n";
		if ($protein) {
			print OUT1 ">$tag[0]\n$protein\n";
		}
	}
}
close OUT1;
close OUT2;
