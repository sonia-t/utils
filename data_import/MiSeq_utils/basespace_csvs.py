'''
roots around in Illumina App dirs for csv results and joins various samples (provided by SampleSheet) into table
'''
__author__ = "Sonia Timberlake"

import pandas as pd
import os 
import sys



def parse_SampleSheet(SampleSheet_file):
    ''' rudimentary Miseq SampleSheet.csv --> df parser ''' 
    with open(SampleSheet_file,'r') as f: lines=f.readlines()
    start_line = [i for i in range(len(lines)) if 'Data' in lines[i] ]
    assert(len(start_line)==1)
    table = [line.strip().split(',') for line in lines[1+start_line[0]:]]
    df = pd.DataFrame( data=table[1:], columns = table[0])
    df.Sample_ID = df.Sample_ID.apply(lambda x: x.strip())
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

def join_sample_csvs(AppResults_dir,
                     samp_ids=['08-120','UC-23 B','UC-23 A','UC-22 C'], 
                     samp_id_suffix=' (2)' ):
    ''' find and read csvs from basespace AppResults_dir dir '''
    dfs = {}
    for samp_id in samp_ids:
        sample_path = '{}/{}{}/Files/'.format(AppResults_dir, samp_id, samp_id_suffix)
        try:
            file_path = [sample_path + fn for fn in os.listdir(sample_path)
                          if fn.endswith('csv') ][0]
            # (fn.startswith(samp_id.replace(' ','-'))
            print('reading csv {} ...'.format(file_path), file=sys.stderr)        
            dfs[samp_id] = pd.read_table(file_path, sep=',')
            dfs[samp_id].set_index( list(dfs[samp_id].columns[:-2]) , inplace=True)
            dfs[samp_id].drop( "%_hits", axis=1, inplace=True)
        except:
            print('WARNING: no csv for samp_id {} found in {}'.format(samp_id ,sample_path ), file=sys.stderr)
    try:
        all = pd.concat(dfs, axis=1)
        all.columns = all.columns.get_level_values(0)
    except:
        all = pd.DataFrame()
    return all

