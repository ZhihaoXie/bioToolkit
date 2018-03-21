#!/usr/bin/perl
# copyright: Zhihao Xie
#
use strict; 
use warnings;
use LWP;
use LWP::Simple;
use HTML::TreeBuilder;
use Getopt::Long;

my $url = "";
my ($i, $j, @kegg_map);
my $fmid;
my $url_pre;
my $prefix;
my $out_dir = './';
my $help;
my $verbose;
GetOptions(
    "i:s"=>\$fmid,
    "u|url:s"=>\$url_pre,
    "o|out:s"=>\$out_dir,
    "h|help!"=>\$help,
    "v!"=>\$verbose,
);

my $usage = "
    Usage: perl $0 -i mapid_file -u url -o out_dir
    Options:
        -i        pathway mapping id (default ko number)
        -u|url    pathway mapping result url
        -o|out    output dir (default current directory)
        -h|help   Detailed help
        -v        Detailed help and version\n\n";
if ($help || $verbose) {
    print STDERR $usage;
    exit;
}
unless ($fmid && $url_pre) {
    print STDERR $usage;
    exit;
}

unless ($url_pre =~ /\/$/) {
    $url_pre = $url_pre . '/';
}
unless ($url_pre =~ /^http:\/\//) {
    $url_pre = 'http://' . $url_pre;
}
if ($url_pre =~ /(http:\/\/[^\/]+)\//) {
    $prefix = $1;
}
#print "$fmid\t$url_pre\t$prefix\n";

#input kegg map_id list
open(HEHE,"<$fmid") or die $!;
$i = -1;
while(<HEHE>){
    $_ =~ s/[\r\n\s]+//g;
    $i ++;
    $kegg_map[$i] = $_;
}
close HEHE;

for $j (0 .. $#kegg_map) {
    $url = $url_pre.$kegg_map[$j].'.args';
    my $html = get($url);
    my $root = HTML::TreeBuilder->new_from_content($html);
    my %images;
    foreach my $node ($root->find_by_tag_name('img')) {
        $images{ $node->attr('src') }++;
    }
    foreach my $pic (sort keys %images) {
        if($pic =~ /$kegg_map[$j]/){
            my $agent     =   LWP::UserAgent->new();
            my $request   =   HTTP::Request->new(GET=>$prefix.$pic);
            print "$pic\n";
            my $response  =   $agent->request($request);
            open(HEHE, ">$out_dir/$kegg_map[$j].png") or die $!;
            binmode HEHE;
            print HEHE $response->content();
            close HEHE;
        }
    }
}
