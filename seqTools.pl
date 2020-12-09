#!/usr/bin/perl 
 
use strict;
use warnings;
 
use Getopt::Std;
use Bio::SeqIO;
 
&main;
exit;
 
sub main {
    &usage if ( @ARGV < 1 );
    my $command  = shift(@ARGV);
    my %function = (
        format  => \&format,
        revcom  => \&revcom,
        length  => \&length,
        content => \&content,
        search  => \&search,
    );
    die("Unknown command \"$command\"\n") if ( !defined( $function{$command} ) );
    &{$function{$command}};
}
 
sub format {
    my %opts = ( w => 60, l => 0 );
    getopts( 'i:o:w:l:', \%opts );
    die(
        qq/
Usage:    seqTools.pl format [options]
 
Options:  -i STR   Input file in FASTA format.
          -o STR   Output file in FASTA format.
          -w INT   The line width for FASTA output. [$opts{w}]
          -l INT   Filters contigs shorter than [int] nt. [don't filter]
\n/
    ) unless ( exists $opts{i} && exists $opts{o} );
 
    my ( $fi, $fo, $width ) = ( $opts{i}, $opts{o}, $opts{w} );
    my $in  = Bio::SeqIO->new( -file => "$fi",  -format => 'fasta' );
    my $out = Bio::SeqIO->new( -file => ">$fo", -format => 'fasta' );
    while ( my $seq = $in->next_seq() ) {
        if ($seq->length >= $opts{l}) {
            $out->width($width);
            $out->write_seq($seq);
        }
    }
}
 
sub revcom {
    my %opts = ( w => 60 );
    getopts( 'i:o:w:', \%opts );
    die(
        qq/
Usage:    seqTools.pl revcom [options]
 
Options:  -i STR   Input file in FASTA format.
          -o STR   Output file in FASTA format.
          -w INT   The line width for FASTA output. [$opts{w}]
\n/ ) unless ( exists $opts{i} && exists $opts{o} );
 
    my ( $fi, $fo, $width ) = ( $opts{i}, $opts{o}, $opts{w} );
    my $in  = Bio::SeqIO->new( -file => "$fi",  -format => 'fasta' );
    my $out = Bio::SeqIO->new( -file => ">$fo", -format => 'fasta' );
    while ( my $seq = $in->next_seq() ) {
        my $rc_seq = $seq->revcom;
        $out->width($width);
        $out->write_seq($rc_seq);
    }
}
 
sub length {
    my %opts;
    getopts( 'i:', \%opts );
    die(
        qq/
Usage:    seqTools.pl length [options]
 
Options:  -i STR   Input file in FASTA format.
\n/
    ) unless ( exists $opts{i} );
 
    my $fi = $opts{i};
    my $in = Bio::SeqIO->new( -file => "$fi", -format => 'fasta' );
    while ( my $seq = $in->next_seq() ) {
        print $seq->id, "\t", $seq->length, "\n";
    }
}
 
sub content {
    my %opts;
    getopts( 'i:p:', \%opts );
    die(
        qq/
Usage:    seqTools.pl content [options]
 
Options:  -i STR   Input file in FASTA format.
\n/
    ) unless ( exists $opts{i} );
 
    my $fi = $opts{i};
    my $in = Bio::SeqIO->new( -file => "$fi", -format => 'fasta' );
    while ( my $seq = $in->next_seq() ) {
        print $seq->id, "\t";
        my $length   = $seq->length;
        my $sequence = $seq->seq;
        my $num_a    = $sequence =~ tr/A//;
        my $num_c    = $sequence =~ tr/C//;
        my $num_g    = $sequence =~ tr/G//;
        my $num_t    = $sequence =~ tr/T//;
        my $num_o    = $length - $num_a - $num_c - $num_g - $num_t;
        my $gc       = sprintf( "%.2f", ( $num_c + $num_g ) / $length * 100 );
        print "Length=$length,GC=$gc,A=$num_a,C=$num_c,G=$num_g,T=$num_t,Others=$num_o\n";
    }
}
 
sub search {
    my %opts;
    getopts( 'i:p:', \%opts );
    die(
        qq/
Usage:    seqTools.pl search [options]
 
Options:  -i STR   Input file in FASTA format.
          -p STR   Pattern you want to search.
 
Notes:    1. You can use RegExp supported by Perl.
          2. It uses case-sensitive mode.
\n/
    ) unless ( exists $opts{i} && $opts{p} );
 
    my ( $fi, $pattern ) = ( $opts{i}, $opts{p} );
    my $in = Bio::SeqIO->new( -file => "$fi", -format => 'fasta' );
    while ( my $seq = $in->next_seq() ) {
        print $seq->id, "\t";
        my @starts;
        my $sequence = $seq->seq;
        while ( $sequence =~ m/($pattern)/gi ) {
            my $start = pos($sequence) - CORE::length($1) + 1;
            pos($sequence) = $start;
            push @starts, $start;
        }
        my $num = @starts;
        my $str = join ",", @starts;
        print "$num\t$str\n";
    }
}
 
sub usage {
    die(
        qq/
Usage:   seqTools.pl <command> [<arguments>]\n
Command: format   Format the FASTA record(s).
         revcom   Reverse and complement the FASTA record(s).
         length   Get the length of FASTA record(s).
         content  Calculate the GC content of FASTA record(s).
         search   Find subseq in the FASTA record(s).
 
Author:  Zhihao Xie, 
Version: 20160824
Notes:  1. A similar webtool developed by lh3: http:\/\/lh3lh3.users.sourceforge.net\/fasta.shtml
        2. When you search, you can use RegExp supported by Perl.
        3. The search function use case-sensitive mode.
\n/
    );
}
