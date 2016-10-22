'''
select fastq reads - probably better to use samtools faidx or pyfaidx but coulnd't get it to work for complexer- pattern matching
'''
__author__ = "Sonia Timberlake"

import os
import sys
import argparse
from os.path import basename, dirname, join
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# Globals ? 


# Functions 

def select_records_by_substring(inh, suffix, substring_list, outh):
    ''' probably reinventing the wheel :( '''
    count = 0
    for record in SeqIO.parse(inh, suffix):
        if any(x in record.description for x in substring_list):
            SeqIO.write(record, outh, suffix)
            count += 1
    return count

def parse_args(argv):
    parser = argparse.ArgumentParser(usage=__doc__,
                                     description='''
                                     count reads in a fastq/fasta / usearch file 
                                     ''')
    parser.add_argument('--infile', 
                        type=argparse.FileType('r'), default=sys.stdin)
    parser.add_argument('--outfile', 
                        type=argparse.FileType('w'), default=sys.stdout)
    
    parser.add_argument('--pattern',
                        default=sys.stdin , help='supports substrings but not regex')
    parser.add_argument('--file_format', choices=["tab","fasta","fastq"],
                        default='fasta' , help='tab or fast[aq]')
    parser.add_argument('--test', action='store_true')
    return parser.parse_args(argv)



def main(argv=sys.argv):
    opts = parse_args(argv[1:])

    prefix, suffix = basename(opts.infile.name).split('.', 1)

    if suffix not in ('fastq','fasta','fa','fq','fst', 'fna', 'faa'):
        raise Exception("Encountered file {} in rule counts which doesn't end "
                        "with fastq, fasta, fa, fq, faa, fna, fst.".format(infile.name))
    suffix = 'fastq' if suffix in ('fastq', 'fq') else 'fasta'


    with open(opts.pattern, "r") as inh:
        substring_list = [line.rstrip('\n') for line in inh]

    count = select_records_by_substring(opts.infile, suffix, substring_list, opts.outfile)

    print("{} records match between {}, {} .....  output to {}".format(
        count, opts.pattern, opts.infile.name, opts.outfile.name))
    return (count, opts.infile.name, opts.outfile.name)

if __name__ == '__main__':
    main()
