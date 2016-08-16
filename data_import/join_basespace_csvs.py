'''
roots around in Illumina App dirs for csv results and joins various samples (provided by SampleSheet) into table

Example invocation
python3 join_basespace_csvs.py  --infile  ~/dummy  --outfile  test10.csv  --AppResults_dir  '/shared/soniat/misc/basespace_FT/mount/Projects/MiSeq M02127 2016 Run206/AppResults'  --SampleSheet  /shared/soniat/misc/basespace_FT/SampleSheet.csv
'''
__author__ = "Sonia Timberlake"

import pandas as pd
import os 
import sys
import argparse 

from MiSeq_utils import *

def parse_args(argv):
    parser = argparse.ArgumentParser(usage=__doc__,
                                     description='''
                                     get sample info from Sample Sheet,
                                     join csvs from basespace App dir 
                                     ''')
    parser.add_argument('--infile', required=True,
                        type=argparse.FileType('r'), default=sys.stdin)
    parser.add_argument('--outfile', required=True,
                        type=argparse.FileType('w'), default=sys.stdout)
    parser.add_argument('--AppResults_dir', default='./',
                        type=str )
    parser.add_argument('--SampleSheet', required=True,
                        type=str)
    parser.add_argument('--test', action='store_true')
    parser.add_argument('--samp_id_suffix', default=' (2)' )
    return parser.parse_args(argv)

def main(argv):
    #if not argv: 0
    opts = parse_args(argv[1:])
    samp_sheet_df = parse_SampleSheet(opts.SampleSheet)
    samp_ids = list(samp_sheet_df.Sample_ID)

    if opts.test:
        all = join_sample_csvs(opts.AppResults_dir, samp_id_suffix=opts.samp_id_suffix)
    else:
        all = join_sample_csvs(opts.AppResults_dir, samp_ids=samp_ids, samp_id_suffix=opts.samp_id_suffix)
    all.to_csv(opts.outfile)
    return all

if __name__ == '__main__':
    main(sys.argv)
