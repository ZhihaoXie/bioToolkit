#!/usr/bin/perl
#
# Author:    Zhihao Xie  \(#^o^)/
#
use File::Basename;
use Data::Dumper;

#summary the status of an assembly

if (@ARGV < 1) {
    die(qq/
    # summary the status of an assembly
    Usage:
        perl $0 <scaffold_seq_list> [large_length_filter]\n
    /);
}
my $scaffold = shift;
my $large = shift;
unless ($large) {
    $large = 1;
}

my %stats;
#my @all_gname;
print "Sample\tTotal_number\tTotal_len(bp)\tAVG_len(bp)\tLongest_len(bp)\tShortest_len(bp)\tN50\tN90\tN_len(bp)\tGC(%)\n";

open TABLE, "< $scaffold" or die $!;
while (<TABLE>) {
    next if /^#|^\s*$/;
    chomp;
    my @list = split /\t/, $_;
    $list[0] =~ s/^\s*|\s*$//g;
    $list[1] =~ s/^\s*|\s*$//g;
    #push @all_gname, $list[0];

    my (@names, %seq);
    open (SCAF, "< $list[1]") or die $!;
    while (<SCAF>) {
        if($_ =~ /^>(\S+)/) {
            $name = $1;
            $names[@names] = $name;
        } else {
            chomp $_;
            $seq{$name} .= $_;
        }
    }
    close SCAF;

    @names = sort {length($seq{$b}) <=> length($seq{$a})} @names;

    my $total_number = 0;
    my $total_length = 0;
    my $largest_length;
    my $shortest_length;
    my $N50_length;
    my $N90_length;
    my ($GC_total, $GC_percent, $N_length, $N_len_percent, $gaps);

    ### largest scaffold and largest length of scaffold ###
    my $largest_scaf = $names[0];
    $largest_length = length($seq{$names[0]});
    if ($largest_length < $large) {
        $largest_length = 0;
    }

    ### scaffold summary ###
    for (my $i=0; $i<@names; $i++) {
        if (length($seq{$names[$i]}) >= $large) {
            $total_number++;
            $total_length += length($seq{$names[$i]});
            if (defined $shortest_length && length($seq{$names[$i]}) < $shortest_length) {
                $shortest_length = length($seq{$names[$i]});
            } elsif (!defined $shortest_length) {
                $shortest_length = length($seq{$names[$i]});
            }
            my $tmp_N_num = $seq{$names[$i]}=~tr/Nn/Nn/;
            $N_length += $tmp_N_num;
            my $tmp_gap = $seq{$names[$i]}=~s/(N+)/$1/ig;
            $gaps += $tmp_gap;
            my $GC_num = $seq{$names[$i]}=~tr/GCgc/GCgc/;
            $GC_total += $GC_num;
        }
    }
    $GC_percent = sprintf("%.2f", $GC_total/$total_length*100);
    my $mean_length = sprintf("%.2f", $total_length/$total_number);
    $N_len_percent = sprintf("%.3f", $N_length/$total_length*100);

    my $length;
    my ($N50_counts, $N90_counts);
    for (my $i=0; $i<@names; $i++) {
        $length += length($seq{$names[$i]});
        if ($length >= 0.5*$total_length && $N50_length eq "") {
            $N50_length = length($seq{$names[$i]});
            $N50_counts = $i + 1;
        }
        if ($length >= 0.9*$total_length && $N90_length eq "") {
            $N90_length = length($seq{$names[$i]});
            $N90_counts = $i + 1;
        }
    }

    #print "SampleID\tTotal number(#)\tTotal len.(bp)\tAverage len.(bp)\tMax len.(bp)\tN50 len.(bp)\tN90 len.(bp)\tGC(%)\n";
    print "$list[0]\t$total_number\t$total_length\t$mean_length\t$largest_length\t$shortest_length\t$N50_length\t$N90_length\t$N_length\t$GC_percent\n";

    #print "Sample\tTotal_num\tTotal_length(bp)\tMean_length(bp)\tLongest_scaf_name\tLongest(bp)\tShortest(bp)\tN50_length(bp)\tN50_counts\tN90_length(bp)\tN90_counts\tGC(%)\tN_length\tN_percent(%)\tGaps\tContigs_counts\tLongest(bp)\tN50_length(bp)\tN50_counts\tN90_length(bp)\tN90_counts\n";
    #print "$basename\t$total_number\t$total_length\t$mean_length\t$largest_scaf\t$largest_length\t$shortest_length\t$N50_length\t$N50_counts\t$N90_length\t$N90_counts\t$GC_percent\t$N_length\t$N_len_percent\t$gaps\t";

}

