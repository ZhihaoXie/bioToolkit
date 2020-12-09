#!/usr/bin/perl -w
#
# perl choosefastq <list> <fastq> [> out]
# copyright: Zhihao Xie
#

if (@ARGV < 2) {
    print STDERR "perl $0 <ID_list> <fastq> [> out]\nAuthor: xzh \t \@2014.4.22\t  \\(#^0^)/\n";
    exit;
}

my ($list_f,$fastq_f) = @ARGV;

my %choose;
open LIST,"< $list_f";
while (<LIST>) {
    chomp;
    if (/^@/) {
        $choose{$_} = 1;
    } else {
        $choose{"@".$_} = 1;
    }
}
close LIST;

open FASTQ,"< $fastq_f";
while (<FASTQ>) {
    if (/^@/) {
        my @arr = split;
        if ($choose{$arr[0]}) {
            my $position = tell(FASTQ);
            print $_;
            my $line1 = <FASTQ>; print $line1;
            my $line2 = <FASTQ>; print $line2;
            my $line3 = <FASTQ>; print $line3;
            seek(FASTQ, $position, 0);
        }
    }
}
close FASTQ;
