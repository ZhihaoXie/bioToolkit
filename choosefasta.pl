#!/usr/bin/perl -w
#
# perl choosefasta.pl <list> <fasta> [> <output>]
# Copyright: Zhihao Xie
#
use Bio::SeqIO;

if(@ARGV < 2) {
    print STDERR "perl choosefasta.pl <list> <fasta> [> <output>]\n";
    exit;
}

my $list = shift;
my $input = shift;
my $fasta = Bio::SeqIO->new(-file=>$input,-format=>'fasta');
while(my $seq = $fasta->next_seq){
	my $id = $seq->id;
	$hash{$id} = $seq;
}
open FILE,$list or die "Can not open $list, please check input!'\n\n";
while(<FILE>){
	chomp;
	my @a = split /\t/,$_;
	unless($hash{$a[0]}){
		print STDERR "The sequence name $a[0] not exite in the fasta file~\n";
		next;
	}
	my $seq = $hash{$a[0]};
	my $id;
	my $str;
	unless($a[1]){
		$str = $seq->seq;
		$id = $seq->id;
		print ">$id\n";
		print "$str\n";
		next;
	}
	if($a[2]>=$a[1]){
		if($a[2] > $seq->length){
			$a[2]=$seq->length;
		}
		$str = $seq->subseq($a[1],$a[2]);
	}
	else{
		if($a[1] > $seq->length){
			$a[1]=$seq->length;
		}
		$str = $seq->subseq($a[2],$a[1]);
		my $n = Bio::Seq->new(-seq=>$str,-id=>'1');
		$n = $n->revcom;
		$str = $n->seq;
	}
	$id = $a[3];
	print ">$id\n";
	print "$str\n";
}
