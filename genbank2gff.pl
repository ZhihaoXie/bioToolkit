#!/usr/bin/perl -w

# Author:    Zhihao Xie  \(#^o^)/
# Date:      2018.03.20
# Version:   v1.0.0
# CopyRight: Copyright ©Zhihao Xie, All rights reserved.

use strict;
use Bio::SeqIO;
use Bio::Tools::GFF;

if (@ARGV < 2) {
    print "Usage: perl $0 <genbank> <gff_output>\n";
    exit 1;
}
my ($INFILE, $OUTFILE) = @ARGV;

my $seqi = Bio::SeqIO->new(
    -file   => $INFILE,
    -format => 'genbank',
);

my $gff = Bio::Tools::GFF->new(
    -gff_version    => 3,
    -file           => ">$OUTFILE",
);

while (my $seq = $seqi->next_seq) {
    for my $feature ($seq->get_SeqFeatures()) {
        if ($feature->has_tag('score')) {
            (my $score) = $feature->get_tag_values('score');
            $feature->remove_tag('score');
            $feature->score($score);
        }
        my $gene_id;
        my $method = $feature->primary_tag;
        if ($method =~ /(gene|RNA|CDS|exon)/) {
            if ($feature->has_tag('locus_tag')) {
                ($gene_id) = $feature->get_tag_values('locus_tag');
            } elsif ($feature->has_tag('old_locus_tag')) {
                ($gene_id) = $feature->get_tag_values('old_locus_tag');
            } elsif ($feature->has_tag('gene')) {
                ($gene_id) = $feature->get_tag_values('gene');
            }
            if (! $feature->has_tag('ID')) {
                $feature->add_tag_value("ID", $gene_id);
            }
            if (! $feature->has_tag('locus_tag')) {
                $feature->add_tag_value("locus_tag", $gene_id);
            }
        }
        
        $gff->write_feature($feature);
    }
}
$seqi->close();
$gff->close();

open OUTFILE2, ">> $OUTFILE";
print OUTFILE2 "##FASTA\n";
close OUTFILE2;

my $out_obj = Bio::SeqIO->new(-file=>">> $OUTFILE", -format=>'fasta');
$seqi = Bio::SeqIO->new(
    -file   => $INFILE,
    -format => 'genbank',
);
while (my $seq = $seqi->next_seq()) {
    $out_obj->write_seq($seq);
}
