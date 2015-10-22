#!/usr/bin/perl -w

if (@ARGV<2) {
    print STDERR "Usage:\n";
    print STDERR "\tperl $0 list in_dir > output\n\n";
    print STDERR "\tsum snp of core gene to reference.\n\tlist file contain seq name, and reference name is in first line.\n";
    print STDERR "\tcopyright by Zhihao Xie\t 2015.5.15\n";
    exit;
}
use Bio::AlignIO;
use Statistics::Descriptive;
use Data::Dumper;

my $list = shift;
my $in_dir = shift;
my @file_list = glob("$in_dir/*.mus");

my $ref_name;
my @another_spec;

open LIST,"$list";
$ref_name = <LIST>;
chomp($ref_name);
while (<LIST>) {
    chomp;
    push @another_spec, $_;
}
close LIST;

print "Ref\t";                 # print output head line
print join "\t", @another_spec;
print "\tMin\tMax\tMean\n";

for my $file (@file_list) {
    my $alignin = Bio::AlignIO->new(-file=>"$file", -format=>'fasta');
    my $aln = $alignin->next_aln;
    my $new_aln = $aln->remove_columns(['match','gaps']);
    #print Dumper \$new_aln;
    my $ref_gene;
    my %hash;
    my %sum_snp;
    foreach my $seq ($new_aln->each_seq()) {
        (my $str = $seq->display_id) =~ s/\|.*//;
        if ($str =~ /$ref_name/) {
          $ref_gene = $seq->display_id;
        }
        if ($seq->seq) {
            @{$hash{$str}} = split //,$seq->seq;
        } else {
            $sum_snp{$str} = 0;
        }
    }

    foreach my $i (@another_spec) {
        if (defined $hash{$i}) {
            $sum_snp{$i} = 0;
            for (my $j = 0; $j<=$#{$hash{$i}}; $j++) {
                if ($hash{$i}[$j] ne $hash{$ref_name}[$j]) {
                    $sum_snp{$i}+=1;
                }
            }
        }
    }
    
    my @sum;
    for my $i (@another_spec) {
        push @sum, $sum_snp{$i}; 
    }
    my $stat = Statistics::Descriptive::Full->new();
    $stat->add_data(\@sum);
    my $mean = $stat->mean();
    my $min=$stat->min();
    my $max=$stat->max();
    print $ref_gene,"\t";
    print join "\t", @sum;
    print "\t",$min,"\t",$max,"\t",$mean,"\n";
}
