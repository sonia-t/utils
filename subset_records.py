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
    count_matches_found =0
    all_strings_matched = set()
    for record in SeqIO.parse(inh, suffix):
        strings_matched = [x for x in substring_list if x in record.description]
        if any(strings_matched):
            SeqIO.write(record, outh, suffix)
            count_matches_found += 1
            all_strings_matched |=  set(strings_matched)
    return (count_matches_found, all_strings_matched)

def select_records_by_exact(inh, suffix, string_list, outh):
    ''' probably reinventing the wheel :( '''
    count_matches_found =0
    all_strings_matched = set()
    for record in SeqIO.parse(inh, suffix):
        strings_matched = [x for x in string_list if x==record.description]
        if any(strings_matched):
            SeqIO.write(record, outh, suffix)
            count_matches_found += 1
            all_strings_matched |=  set(strings_matched)
    return (count_matches_found, all_strings_matched)

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
                        default=sys.stdin,
                        help='supports substrings but not regex')
    parser.add_argument('--file_format', choices=["tab","fasta","fastq"],
                        default='fasta' , help='tab or fast[aq]')
    parser.add_argument('--exact', action='store_true')
    parser.add_argument('--test', action='store_true')
    return parser.parse_args(argv)



def main(argv=sys.argv):
    opts = parse_args(argv[1:])

    prefix, suffix = os.path.splitext(opts.infile.name) # basename(opts.infile.name).split('.', 1)
    suffix =suffix[1:]
    if suffix not in ('fastq','fasta','fa','fq','fst', 'fna', 'faa'):
        raise Exception("Encountered file {} in rule counts with suffix {} not in "
                        "set: (fastq, fasta, fa, fq, faa, fna, fst) ".format(opts.infile.name, suffix))
    suffix = 'fastq' if suffix in ('fastq', 'fq') else 'fasta'


    with open(opts.pattern, "r") as inh:
        substring_list = [line.rstrip('\n') for line in inh]

    select_records = select_records_by_exact if opts.exact else select_records_by_substring
    count_matches_found, all_strings_matched= select_records(opts.infile, suffix, substring_list, opts.outfile)

    print("matching {}, {} .....  output to {}".format(
        opts.pattern, opts.infile.name, opts.outfile.name))
    
    print("searched for {} substrings \n{} substrings match {} records".format(
        len(set(substring_list)),len(all_strings_matched), count_matches_found))
    print( 'substrings not found: ' + ",".join(set(substring_list) - set(all_strings_matched) ))
    
    return count_matches_found, all_strings_matched, opts.pattern, opts.infile.name, opts.outfile.name

if __name__ == '__main__':
    main()
