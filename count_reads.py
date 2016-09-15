'''
count reads or records in usearch-formatted fastq and tab files 
'''
__author__ = "Sonia Timberlake"

import os
import sys
import argparse
from os.path import basename, dirname, join
from Bio import SeqIO
from collections import OrderedDict

# Globals ? 
default_delimiter={'fields': ';','key-value':'='}

# Functions 
def parseAnnotation(record, fields=None, delimiter=default_delimiter):
    """
    Extracts annotations from a FASTA/FASTQ sequence description

    Arguments:
      record : Description string to extract annotations from
      fields : List of fields to subset the return dictionary to;
               if None return all fields
      delimiter : a tuple of delimiters for (fields, values, value lists)

    Returns:
      OrderedDict : An OrderedDict of field/value pairs
    """
    annotation = record.split(delimiter[0])
    field_dict = OrderedDict([('ID', annotation.pop(0))])
    for ann in annotation:
        vals = ann.split(delimiter[1])
        field_dict[vals[0].upper()] = vals[1]

    # Subset return dictionary to requested fields
    if fields is not None:
        if not isinstance(fields, list):  fields = [fields]
        for f in set(field_dict).difference(fields):  del field_dict[f]

    return field_dict

def parse_usearch_header(header, 
                         delims={'fields':';','key-value':'='}):
    '''
    Extracts annotations from a FASTA/FASTQ sequence description
    Default delimiters are the usearch/qiime format
    '''
    fields = header.split(delims['fields'])
    field_dict = OrderedDict([('ID', fields.pop(0))])
    for field in fields:
        pair = field.split(delims['key-value'])
        if len(pair)>1: # in ase not every field is a key-value pair 
            field_dict[pair[0].lower()] = pair[1]
        else:
            field_dict[pair[0].lower()] = pair[0]
    return field_dict

def count_tab(infile, count_unit='reads', derep_field=None):
    ''' when wc -l won't do because the reads are dereplicated '''
    if count_unit == 'reads' and derep_field:
        count = sum( int(parse_usearch_header(line)[derep_field]) for line in infile)
    else:
        count = sum( 1 for line in infile)
    # TODO add in ID dereplication ? 
    return count, count_unit

def count_fastq(infile, count_unit='reads'):
    prefix, suffix = basename(infile.name).split('.', 1)
    if len(suffix)==2: suffix = {'fq':'fastq', 'fa':'fasta'}[suffix] # lazy exception handling
    if not suffix in ('fastq','fasta'):
        raise Exception("Encountered file {} in rule counts which doesn't end "
                        "with fastq, fasta, fa, fq.".format(infile.name))
    if 'derep' in prefix and count_unit=='reads': # have to use the size field # TODO: maybe should look for it by default ? 
        derep_field = 'size'
        count = sum( int(parse_usearch_header(record.description)[derep_field])
                                 for record in SeqIO.parse(infile, suffix))
    else:
        count = sum(1 for record in SeqIO.parse(infile, suffix))
    return count, count_unit

def parse_args(argv):
    parser = argparse.ArgumentParser(usage=__doc__,
                                     description='''
                                     count reads in a fastq/fasta / usearch file 
                                     ''')
    parser.add_argument('--infile', 
                        type=argparse.FileType('r'), default=sys.stdin)
    parser.add_argument('--outfile', 
                        type=argparse.FileType('w'), default=sys.stdout)
    
    parser.add_argument('--count_unit', choices=["reads","records"],
                        default='reads' , help='in dereplicated fasta, reads!=records')
    parser.add_argument('--file_format', choices=["tab","fasta","fastq"],
                        default='fasta' , help='tab or fast[aq]')
    parser.add_argument('--test', action='store_true')
    return parser.parse_args(argv)



def main(argv=sys.argv):
    opts = parse_args(argv[1:])

    prefix, suffix = basename(opts.infile.name).split('.', 1)

    if opts.file_format=='tab':
        count, count_unit =count_tab(opts.infile, opts.count_unit, derep_field='size')
    else:
        count, count_unit = count_fastq(opts.infile, opts.count_unit)
    print("{}\t{}\t{}".format(count, count_unit, opts.infile.name), file=opts.outfile)
    return (count, opts.infile.name)

if __name__ == '__main__':
    main()
