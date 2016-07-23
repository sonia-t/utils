__author__ = "Sonia Timberlake"

import pandas as pd
import os 


indir = '/shared/soniat/misc/basespace_FT/mount/Projects/MiSeq M02127 2016 Run206/AppResults/'

dfs={}

def parse_SampleSheet(sampsheetfile):
    ''' rudimentary Miseq SampleSheet.csv --> df parser ''' 
    with open(sampesheetfile,'r') as f: lines=readlines()
    start_line = [i for i in range(len(lines)) if 'Data' in lines[i] ]
    table = [line.strip().split(',') for line in lines[1+start_lines:]]
    df = pd.DataFrame( data=table[1:], columns = table[0])
    return df 
            

def __old():
    for samp_id, bcnum in dict(zip(['UC-23 B','UC-23 A','UC-22 C'], [73,85,97])).items():
        infile = '{}/{} (2)/Files/{}_S{}.summary.csv'.format(indir, samp_id, samp_id.replace(' ','-'), bcnum)
        dfs[samp_id] = pd.read_table(infile, sep=',')
        dfs[samp_id].set_index( list(dfs[samp_id].columns[:-2]) , inplace=True)
        dfs[samp_id].drop( "%_hits", axis=1, inplace=True)


    all = pd.concat(dfs, axis=1)

def read_three_samples(samp_ids=['UC-23 B','UC-23 A','UC-22 C']):
    for samp_id in samp_ids:
        relevant_path = '{}/{} (2)/Files/'.format(indir,samp_id)
        included_extenstions = ['csv']
        file_names = [fn for fn in os.listdir(relevant_path)
                      if  (fn.startswith(samp_id.replace(' ','-'))
                      and fn.endswith('csv')) ]
        print(file_names)
        assert (len(file_names)==1)
        infile = file_names[0]
        dfs[samp_id] = pd.read_table(infile, sep=',')
        dfs[samp_id].set_index( list(dfs[samp_id].columns[:-2]) , inplace=True)
        dfs[samp_id].drop( "%_hits", axis=1, inplace=True)

    all = pd.concat(dfs, axis=1)
    return all

def parse_args(argv):
    parser = argparse.ArgumentParser(usage=__doc__,
                                     description='''
                                     get sample info from Sample Sheet,
                                     join csvs from basespace App dir 
                                     ''')
    parser.add_argument('--in-clone-stats', 
                        type=argparse.FileType('r'), default=sys.stdin)
    parser.add_argument('--in-premode',
                        type=argparse.FileType('r')


def main(argv=None):
    if not argv == 0: sys.exit()
    opts = parse_args(argv)
    df = parse_SampleSheet(sampsheetfile)

if __name__ == '__main__':
    main(sys.argv[1:])
