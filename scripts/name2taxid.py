#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
'''
@Filename:    name2taxid.py
@Version:     V0.0.0.1
@Date:        2021/10/12 09:47:00
@Description: scientific name to taxonmy id
'''

import sys
import os
import re
import argparse


def name2taxidQ(name_dmp, q_name):
    col_delimiter = '\t|\t'
    row_delimiter = '\t|\n'
    with open(name_dmp, 'r') as names_file:
        for line in names_file:
            if re.search(r'^\s*$|^#', line):
                continue
            elif len(line) == 0:
                break
            line = line.rstrip(row_delimiter)
            values = line.split(col_delimiter)
            tax_id, name_txt, _, name_type = values[:4]
            if name_txt == q_name:
                return q_name, tax_id
            else:
                continue


def name2taxidF(name_dmp, q_file):
    col_delimiter = '\t|\t'
    row_delimiter = '\t|\n'
    scientific_names = {}
    common_names = {}
    # Reading taxonomy
    with open(name_dmp, 'r') as names_file:
        for line in names_file:
            if re.search(r'^\s*$|^#', line):
                continue
            elif len(line) == 0:
                break
            line = line.rstrip(row_delimiter)
            values = line.split(col_delimiter)
            tax_id, name_txt, _, name_type = values[:4]
            if name_type == 'scientific name':
                scientific_names[name_txt] = tax_id
            elif name_type == 'common name':
                common_names[name_txt] = tax_id

    if os.path.isfile(q_file):
        with open(q_file, 'r') as f:
            for line in f:
                if re.search(r'^\s*$|^#', line):
                    continue
                elif len(line) == 0:
                    break
                sci_name = line.strip().replace("_", " ")
                if sci_name in scientific_names:
                    yield sci_name, scientific_names[sci_name]
                elif sci_name in common_names:
                    yield sci_name, common_names[sci_name]


def main():
    parser = argparse.ArgumentParser(description="scientific name to taxonmy id")
    parser.add_argument('-n', dest='namedmp', help="names.dmp of taxonomy as query database")
    parser.add_argument('-q', dest='query', help="query sci name")
    parser.add_argument('-qf', dest="queryfile", help="query sci name file, one line as a name, -q and -qf conflict")
    args = parser.parse_args()
    if len(sys.argv) <= 1:
        parser.print_help()
        sys.exit()
    namedmp = args.namedmp
    query = args.query
    query_file = args.queryfile
    if query:
        query = query.replace("_", " ")
        sci_name, sci_taxid = name2taxidQ(namedmp, query)
        print(f"{sci_name}\t{sci_taxid}")
    elif query_file:
        for sci_name, sci_taxid in name2taxidF(namedmp, query_file):
            print(f"{sci_name}\t{sci_taxid}")


if __name__ == "__main__":
    main()
