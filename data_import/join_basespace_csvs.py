import pandas as pd
import os 


indir = '/shared/soniat/misc/basespace_FT/mount/Projects/MiSeq M02127 2016 Run206/AppResults/'

dfs={}


for samp_id, bcnum in dict(zip(['UC-23 B','UC-23 A','UC-22 C'], [73,85,97])).items():
    infile = '{}/{} (2)/Files/{}_S{}.summary.csv'.format(indir, samp_id, samp_id.replace(' ','-'), bcnum)
    dfs[samp_id] = pd.read_table(infile, sep=',')
    dfs[samp_id].set_index( list(dfs[samp_id].columns[:-2]) , inplace=True)
    dfs[samp_id].drop( "%_hits", axis=1, inplace=True)


all = pd.concat(dfs, axis=1)


for samp_id in ['UC-23 B','UC-23 A','UC-22 C']:
    relevant_path = '{}/{} (2)/Files/'.format(indir,samp_id)
    included_extenstions = ['csv']
    file_names = [fn for fn in os.listdir(relevant_path)
                  if  (fn.startswith(samp_id.replace(' ','-'))
                  and fn.endswith('csv')) ]
    print(file_names)
    assert (len(file_names)==1)
    infile = file_names[0]


