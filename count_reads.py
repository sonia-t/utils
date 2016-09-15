'''
fastq utils - remember not to reinvent the wheel! 
'''
__author__ = "Sonia Timberlake"

import os
import sys
import argparse
from os.path import basename, dirname, join

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
                         delims={'fields':' ','key-value':'='}):
    '''
    Extracts annotations from a FASTA/FASTQ sequence description
    Default delimiters are the usearch/qiime format
    '''
    fields = header.split(delims['fields'])
    field_dict = OrderedDict([('ID', fields.pop(0))])
    for field in fields:
        pairs = field.split(delims['key-value'])
        field_dict[pairs[0].lower()] = pairs[1]
    return field_dict

def count_fastq(filename, units='reads'):
    prefix, suffix = basename(filename).split('.', 1)
    if len(suffix)==2: suffix = {'fq':'fastq', 'fa':'fasta'}[suffix] # lazy exception handling
    if not assert(suffix in ('fastq','fasta')):
        raise Exception("Encountered file {} in rule counts which doesn't end "
                        "with fastq, fasta, fa, fq.".format(filename))
    if 'derep' in prefix and units=='reads': # have to use the size field # TODO: maybe should look for it by default ? 
        header_field = 'size'
        count = sum( int(parse_usearch_header(record.description)[header_field])
                                 for record in SeqIO.parse(open(str(fname)), suffix) )
    else:
        count = sum(1 for record in SeqIO.parse(open(str(fname)), suffix))
    return 

def parse_args(argv):
    parser = argparse.ArgumentParser(usage=__doc__,
                                     description='''
                                     count reads in a fastq/fasta / usearch file 
                                     ''')
    parser.add_argument('--infile', required=True,
                        type=argparse.FileType('r'), default=sys.stdin)
    parser.add_argument('--outfile', required=True,
                        type=argparse.FileType('w'), default=sys.stdout)
    parser.add_argument('--count_units', choices=["reads","records"],
                        default='reads' , help='in dereplicated fasta, reads!=records')
    parser.add_argument('--test', action='store_true')
    return parser.parse_args(argv)



def main(argv):
    #if not argv: 0
    opts = parse_args(argv[1:])
    count = count_fastq(opts.infile, opts.count_units)
    return (count, opts.infile.name)

if __name__ == '__main__':
    main(sys.argv)
