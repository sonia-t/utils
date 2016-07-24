'''
roots around in Illumina App dirs for csv results and joins various samples (provided by SampleSheet) into table

Example invocation
%run -d join_basespace_csvs.py  --infile  ~/dummy  --outfile  test10.csv  --AppResults_dir  '/shared/soniat/misc/basespace_FT/mount/Projects/MiSeq M02127 2016 Run206/AppResults'  --SampleSheet  /shared/soniat/misc/basespace_FT/SampleSheet.csv
'''
__author__ = "Sonia Timberlake"

import pandas as pd
import os 
import sys
import argparse 



def parse_SampleSheet(SampleSheet_file):
    ''' rudimentary Miseq SampleSheet.csv --> df parser ''' 
    with open(SampleSheet_file,'r') as f: lines=f.readlines()
    start_line = [i for i in range(len(lines)) if 'Data' in lines[i] ]
    assert(len(start_line)==1)
    table = [line.strip().split(',') for line in lines[1+start_line[0]:]]
    df = pd.DataFrame( data=table[1:], columns = table[0])
    return df 
            

def __old():
    indir = '/shared/soniat/misc/basespace_FT/mount/Projects/MiSeq M02127 2016 Run206/AppResults/'
    dfs={}
    for samp_id, bcnum in dict(zip(['UC-23 B','UC-23 A','UC-22 C'], [73,85,97])).items():
        infile = '{}/{} (2)/Files/{}_S{}.summary.csv'.format(indir, samp_id, samp_id.replace(' ','-'), bcnum)
        dfs[samp_id] = pd.read_table(infile, sep=',')
        dfs[samp_id].set_index( list(dfs[samp_id].columns[:-2]) , inplace=True)
        dfs[samp_id].drop( "%_hits", axis=1, inplace=True)
    all = pd.concat(dfs, axis=1)
    return all

def join_sample_csvs(AppResults_dir, samp_ids=['UC-23 B','UC-23 A','UC-22 C']):
    ''' find and read csvs from basespace AppResults_dir dir '''
    dfs = {}
    for samp_id in samp_ids:
        sample_path = '{}/{} (2)/Files/'.format(AppResults_dir,samp_id)
        try:
            file_paths = [sample_path + fn for fn in os.listdir(sample_path)
                          if (fn.startswith(samp_id.replace(' ','-'))
                              and fn.endswith('csv')) ]
        except:
            print('WARNING: no csv not found in {}'.format(samp_id ,sample_path ), file=sys.stderr)
            continue
        else:
            print('reading csv {} ...'.format(file_paths[0]), file=sys.stderr)        
            dfs[samp_id] = pd.read_table(file_paths[0], sep=',')
            dfs[samp_id].set_index( list(dfs[samp_id].columns[:-2]) , inplace=True)
            dfs[samp_id].drop( "%_hits", axis=1, inplace=True)
    try:
        all = pd.concat(dfs, axis=1)
        all.columns = all.columns.get_level_values(0)
    except:
        all = pd.DataFrame()
    return all

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
    return parser.parse_args(argv)

def main(argv):
    #if not argv: 0
    opts = parse_args(argv[1:])
    if opts.test:
        all = join_sample_csvs(opts.AppResults_dir)
        all.to_csv(opts.outfile)
        return
    samp_sheet_df = parse_SampleSheet(opts.SampleSheet)
    all = join_sample_csvs(opts.AppResults_dir, samp_ids = list(samp_sheet_df.Sample_ID)[:3])
    all.to_csv(opts.outfile)
    return all

if __name__ == '__main__':
    main(sys.argv)
