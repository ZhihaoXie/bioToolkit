#!/bin/bash

# ftp site
FTP_SERVER="ftp://ftp.ncbi.nlm.nih.gov/blast/db"

if [ $# -lt 1 ];then
    1>&2 echo "download nt/nr for blast from ncbi blast db"
    1>&2 echo "sh $0 <out_dir> [nr|nt]"
    exit 1
fi

aria2c=`which aria2c`
if [ ! -n $aria2c ];then
    echo "aria2c not in \$PATH. install it or check"
    exit 1
fi

out_dir=`readlink -f $1`
library_name="$2"
if [ ! -d $out_dir ];then
    mkdir -p $out_dir
fi

case $library_name in
    "nr")
        cd $out_dir
        1>&2 echo -n "Downloading $library_name database from FTP..."
        wget -q --no-remove-listing --spider ${FTP_SERVER}/
        awk '{ print $NF }' .listing | perl -ple 'tr/\r//d' | grep 'nr' |grep '\.tar\.gz' > ${library_name}.list
        #cat ${library_name}.list | xargs -n1 -I{} wget -q -c -t 0 $FTP_SERVER/{}
        cat ${library_name}.list | xargs -n1 -I{} aria2c -c -x 10 -u 100M --max-download-limit 100M -t 120 -m 0 --retry-wait 1 -o {} $FTP_SERVER/{}
        rm -f .listing
        1>&2 echo " done."
        ;;
    "nt")
        cd $out_dir
        1>&2 echo -n "Downloading $library_name database from FTP..."
        wget -q --no-remove-listing --spider ${FTP_SERVER}/
        awk '{ print $NF }' .listing | perl -ple 'tr/\r//d' | grep 'nt' |grep '\.tar\.gz' > ${library_name}.list
        #cat ${library_name}.list | xargs -n1 -I{} wget -q -c -t 0 $FTP_SERVER/{}
        cat ${library_name}.list | xargs -n1 -I{} aria2c -c -x 10 -u 100M --max-download-limit 100M -t 120 -m 0 --retry-wait 1 -o {} $FTP_SERVER/{}
        rm -f .listing
        1>&2 echo " done."
        ;;
    *)
        1>&2 echo "Unsupported library. Valid options are: nr nt"
        exit 1
        ;;
esac

