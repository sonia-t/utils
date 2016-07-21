import pandas as pd


dfs={}



for samp_id, bcnum in dict(zip(['UC-23 B','UC-23 A','UC-22 C'], [73,85,97])).items():
    infile = '{} (2)/Files/{}_S{}.summary.csv'.format(samp_id, samp_id.replace(' ','-'), bcnum)
    dfs[samp_id] = pd.read_table(infile, sep=',')
    dfs[samp_id].set_index( list(dfs[samp_id].columns[:-2]) , inplace=True)
    dfs[samp_id].drop( "%_hits", axis=1, inplace=True)


all = pd.concat(dfs, axis=1)


