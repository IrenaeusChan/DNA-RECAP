from __future__ import division
from clint.textui import indent
from tqdm import tqdm
import pysam
import multiprocessing as mp
import pandas as pd

import dnarecap.utils.logger as log
import dnarecap.utils.tools as tools
import dnarecap.process.scarmapper as scarmapper

def trim_cigar_ends(cigar, read_reference_start, read_reference_end, read_query_start, read_query_end):
    """
    Trim leading and trailing matches from the cigar string.
    """
    # Trim leading matches
    while cigar and cigar[0][0] == 0:
        match_length = cigar[0][1]
        read_reference_start += match_length
        read_query_start += match_length
        del cigar[0]

    # Trim trailing matches
    while cigar and cigar[-1][0] == 0:
        match_length = cigar[-1][1]
        read_reference_end -= match_length
        read_query_end -= match_length
        del cigar[-1]

    return cigar, read_reference_start, read_reference_end, read_query_start, read_query_end

def read_to_indelinfo(read, cigar, regionSeq, regionStart, read_reference_start, read_query_start, verbose=False):
    # setup indelinfo dictionary
    indelinfo = {'source':'',
                'read':read.query_name,
                'chrom':'',
                'pos':None,
                'whichread': 'r1' if read.is_read1 else 'r2',
                'strands': '',
                'ref':'',
                'alt':'',
                'type':'',
                'mh':'',
                'mh_length': 0,
                'lft_tmplt':'',
                'rt_tmplt':'',
                'lft_tmplt_pos':'',
                'rt_tmplt_pos':'',
                'classification':''}

    # If there's only one cigar operation left, it's an insertion or deletion
    if len(cigar) == 1:
        # theres one event, so record the coordinates
        indelinfo['chrom'] = read.reference_name
        indelinfo['strands'] = "++" if read.is_reverse is False else "--"

        # if insertion
        if cigar[0][0] == 1:           
            indelinfo['pos'] = read_reference_start-1
            indelinfo['ref'] = regionSeq[read_reference_start-regionStart-1:read_reference_start-regionStart] #fasta.fetch(chr,read_reference_start-1-1,read_reference_start-1) # 0 based and position before the isertion
            indelinfo['alt'] = read.query_sequence[read_query_start-1:read_query_start+cigar[0][1]]
            indelinfo['type'] = 'INDEL'
        
        # if deletion
        elif cigar[0][0] == 2 or cigar[0][0] == 3:
            indelinfo['pos'] = read_reference_start-1
            indelinfo['ref'] = regionSeq[read_reference_start-regionStart-1:read_reference_start-regionStart + cigar[0][1]] #fasta.fetch(chr,read_reference_start-1,read_reference_start+cigar[0][1])
            indelinfo['alt'] = indelinfo['ref'][0]
            indelinfo['type'] = 'INDEL'

    # if its a complex indel
    else:
        indelinfo['chrom'] = read.reference_name
        indelinfo['pos'] = read_reference_start-1
        indelinfo['strands'] = "++" if read.is_reverse is False else "--"
        indelinfo['type'] = 'INDEL'

        # iterate through cigar and get indels
        varlen = sum([ x[1] for x in cigar ])
        indelinfo['ref'] = regionSeq[read_reference_start-regionStart-1:read_reference_start-regionStart+sum([ x[1] for x in cigar ])] #fasta.fetch(chr,read_reference_start-1,read_reference_start+sum([ x[1] for x in cigar ]))
        indelinfo['alt'] = read.query_sequence[read_query_start+1:read_query_start+sum([x[1] if x[0] == 1 or x[0] == 0 else 0 for x in cigar])]

    return indelinfo