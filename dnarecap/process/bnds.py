from __future__ import division
from clint.textui import indent
from tqdm import tqdm
import pysam, sys
import multiprocessing as mp
import pandas as pd
from enum import Enum

import dnarecap.utils.logger as log
import dnarecap.utils.tools as tools
import dnarecap.process.scarmapper as scarmapper

class StrandCombo(Enum):
    FORWARD_FORWARD = "++"
    REVERSE_REVERSE = "--"
    REVERSE_FORWARD = "-+"
    FORWARD_REVERSE = "+-"

class ComparisonType(Enum):
    START_TO_START = "start_to_start"
    START_TO_END = "start_to_end"
    END_TO_START = "end_to_start"
    END_TO_END = "end_to_end"

# Mapping of strand combinations to comparison types
STRAND_TO_COMPARISON = {
    ("leftSoftClip", StrandCombo.FORWARD_FORWARD): ComparisonType.START_TO_END,
    ("leftSoftClip", StrandCombo.REVERSE_REVERSE): ComparisonType.END_TO_START,
    ("leftSoftClip", StrandCombo.REVERSE_FORWARD): ComparisonType.END_TO_START,
    ("leftSoftClip", StrandCombo.FORWARD_REVERSE): ComparisonType.START_TO_END,

    ("rightSoftClip", StrandCombo.FORWARD_FORWARD): ComparisonType.END_TO_START,
    ("rightSoftClip", StrandCombo.REVERSE_REVERSE): ComparisonType.START_TO_END,
    ("rightSoftClip", StrandCombo.REVERSE_FORWARD): ComparisonType.START_TO_END,
    ("rightSoftClip", StrandCombo.FORWARD_REVERSE): ComparisonType.END_TO_START
}

def find_largest_overlap(seq1, seq2, comparison_type, min_match=2):
    """Find largest overlapping segment between sequences based on comparison type"""
    
    def get_sequence_parts(seq, direction):
        if direction == "start":
            return [seq[:i] for i in range(min_match, len(seq) + 1)]
        else:  # end
            return [seq[-i:] for i in range(min_match, len(seq) + 1)]
    
    # Get parts to compare based on comparison type
    if comparison_type == ComparisonType.START_TO_START:
        parts1 = get_sequence_parts(seq1, "start")
        parts2 = get_sequence_parts(seq2, "start")
    elif comparison_type == ComparisonType.START_TO_END:
        parts1 = get_sequence_parts(seq1, "start")
        parts2 = get_sequence_parts(seq2, "end")
    elif comparison_type == ComparisonType.END_TO_START:
        parts1 = get_sequence_parts(seq1, "end")
        parts2 = get_sequence_parts(seq2, "start")
    elif comparison_type == ComparisonType.END_TO_END:
        parts1 = get_sequence_parts(seq1, "end")
        parts2 = get_sequence_parts(seq2, "end")

    # Find largest matching overlap
    largest_overlap = ""
    for part1 in parts1:
        for part2 in parts2:
            if part1 == part2 and len(part1) > len(largest_overlap):
                largest_overlap = part1
                
    return largest_overlap

def print_sequence_alignment(rSeq, pSeq, sSeq, mh="", afterSeq="", comparison_type=None):
    """Print aligned sequences with colored matching regions"""
    
    def print_with_padding(seq, padding=0):
            print(" " * padding + seq, file=sys.stderr)
            
    print("\n=============== Sequence Alignment ===============", file=sys.stderr)
    print_with_padding(f"{rSeq}")
    if comparison_type == ComparisonType.START_TO_END:
        if afterSeq:
            print_with_padding(pSeq, len(sSeq))
            print_with_padding(f"{sSeq} {afterSeq}")
        else:
            print_with_padding(pSeq, len(sSeq)-len(mh))
            print_with_padding(sSeq, 0)
    elif comparison_type == ComparisonType.END_TO_START:
        if afterSeq:
            print_with_padding(pSeq, 0)
            print_with_padding(f"{afterSeq} {sSeq}", len(pSeq)-20)
        else:
            print_with_padding(pSeq, 0)
            print_with_padding(sSeq, len(pSeq)-len(mh))
            
    print("=" * 50, file=sys.stderr)

# Adapted from Dr. Dave Spencer
def read_to_bndinfo(read, read_strand, start, end, read_reference_start, read_reference_end, leftSoftClip, rightSoftClip, ref_fasta, maxNM, svDistanceThreshold, minSecMapQual, verbose=False):
    # setup indelinfo dictionary
    bndinfo = {'source':'',
                'read':read.query_name,
                'chrom':'',
                'pos':None,
                'chrom2':None,
                'pos2':None,
                'whichread': 'r1' if read.is_read1 else 'r2',
                'strands': '',
                'type':'',
                'mh':'',
                'mh_length': 0,
                'classification':''}
    
    sChr, sPos, sStrand, sCigar, sMq, sNm = read.get_tag('SA').split(';')[0].split(',') if read.has_tag('SA') else [None,None,None,None,None,None]

    sCigarTuples = tools.make_cigar_tuples(sCigar)
    sPos = int(sPos)
    sEnd = sPos + sum([x[1] for x in sCigarTuples if x[0] == 0 or x[0] == 2 or x[0] == 3]) - 1
    sLeftSoftClip = sCigarTuples[0][1] if sCigarTuples[0][0] >= 4 else 0
    sRightSoftClip = sCigarTuples[-1][1] if sCigarTuples[-1][0] >= 4 else 0
    lengthSA = sum([x[1] for x in sCigarTuples if x[0] == 0 or x[0] == 1])

    pSeq = read.query_alignment_sequence    # Removes Soft Clip and Hard Clip

    # Apply filters to secondary alignment
    if (int(sMq) >= minSecMapQual and int(sNm) <= maxNM and int(sMq)>0 and
        (sChr != read.reference_name or abs(int(sPos) - read_reference_start) >= svDistanceThreshold) or sStrand != read_strand or 
        (rightSoftClip > 0 and sEnd < read_reference_start) or (leftSoftClip > 0 and sPos > read_reference_end)):

        bndinfo['chrom'] = read.reference_name
        bndinfo['chrom2'] = sChr
        bndinfo['whichread'] = 'r1' if read.is_read1 else 'r2' 
        # Strands matters for directionality, if it's ++, it means both sStrand and read_strand are on the same strand, Rev or Fwd
        # If it is +-, it means they are on opposite strands so if the read_strand is on the Fwd strand, then sStrand is Rev
        if sStrand == read_strand and sStrand == '+':
            bndinfo['strands'] = '++'
        elif sStrand == read_strand and sStrand == '-':
            bndinfo['strands'] = '--'
        elif sStrand != read_strand and sStrand == '+':
            bndinfo['strands'] = '-+'
        elif sStrand != read_strand and sStrand == '-':
            bndinfo['strands'] = '+-'
        else:
            log.logit(f"ERROR: No proper strand orientation found: {read.query_name}: {str(read)}", color="red")
            return None
        bndinfo['type'] = 'BND'

        if bndinfo['strands'] == '++' or bndinfo['strands'] == '--':
            # if --->...> ...>---> or <---<... <...<---
            if rightSoftClip > leftSoftClip and sLeftSoftClip > sRightSoftClip:
                bndinfo['pos'] = read_reference_end
                bndinfo['pos2'] = sPos

            # if ...>---> --->...> or <...<--- <---<...
            elif leftSoftClip > rightSoftClip and sRightSoftClip > sLeftSoftClip:
                bndinfo['pos'] = read_reference_start - 1 
                bndinfo['pos2'] = sEnd + 1

            # This typically happens when there are S/H clips on both ends of the SA
            # e.g. 16S118M16S or 20H110M20H
            # Or the secondary alignment soft clip doesn't match the primary alignment soft clip
            # This read could potentially be misaligned and therefore we skip it
            else:
                log.logit(f"ERROR: No proper soft clip orientation found: {read.query_name}: {str(read)}", color="red")
                return None

        elif bndinfo['strands'] == '+-' or bndinfo['strands'] == '-+':
            # if --->...> <---<... or <---<... --->...>
            if rightSoftClip > leftSoftClip and sRightSoftClip > sLeftSoftClip:
                bndinfo['pos'] = read_reference_end
                bndinfo['pos2'] = sEnd

            # if ...>---> <...<--- or <...<--- ...>--->
            elif leftSoftClip > rightSoftClip and sLeftSoftClip > sRightSoftClip:
                bndinfo['pos'] = read_reference_start - 1
                bndinfo['pos2'] = sPos

            # See error message above for explanation
            else:
                log.logit(f"ERROR: No proper soft clip orientation found: {read.query_name}: {str(read)}", color="red")
                return None
        
        # See error message above for explanation
        else:
            log.logit(f"ERROR: No proper soft clip orientation found: {read.query_name}: {str(read)}", color="red")
            return None

        if verbose:
            log.logit(f"Processing: {read.query_name} --> {bndinfo['chrom']}:{bndinfo['pos']}:{bndinfo['chrom2']}:{bndinfo['pos2']}:{bndinfo['whichread']}:{bndinfo['strands']}", color ='cyan')
        # Since all sequences are always 5' to 3' from BAM/CRAM
        # We need to care about the strands FWD(+) and REV(-)

        # NOTE: query_sequence is the FULL 150 bp
        # NOTE: cigar_dict = {'M': 0, 'I': 1, 'D': 2, 'N': 3, 'S': 4, 'H': 5, 'P': 6, '=': 7, 'X': 8}

        # If the sequence in question alignes to FWD strand, and the SA strand is also "+"
        if read_strand == sStrand:
            # If both the read and SA are on the FWD/REV strand, then we can just take the sequence as is
            if leftSoftClip > rightSoftClip:
                overlap = leftSoftClip - lengthSA
                sSeq = read.query_sequence[sLeftSoftClip:leftSoftClip-overlap+sLeftSoftClip]
                # The +sLeftSoftClip is because sometimes the SA will have S/H clips on BOTH ends, we already account for the single clip via the IF statement
                # But when we calculate the sequence we need to account that the sequence may be shorter due to an additional S/H clip
            elif rightSoftClip > leftSoftClip:   # Right soft clips should be the end of the matched sequence until the end of the soft clip
                overlap = lengthSA - rightSoftClip
                sSeq = read.query_sequence[sum([x[1] for x in read.cigar if x[0] == 0 or x[0] == 1])-overlap-sRightSoftClip:rightSoftClip-sRightSoftClip+sum([x[1] for x in read.cigar if x[0] == 0 or x[0] == 1])]
        elif read_strand != sStrand:
            # If the read is on the REV strand, but the SA is on the FWD strand, we need to reverse complement it                    
            if leftSoftClip > rightSoftClip:
                overlap = leftSoftClip - lengthSA
                sSeq = tools.revcomp(read.query_sequence[sRightSoftClip:leftSoftClip-overlap+sRightSoftClip])
                # The +sRightSoftClip is because sometimes the SA will have S/H clips on BOTH ends, we already account for the single clip via the IF statement
                # But when we calculate the sequence we need to account that the sequence may be shorter due to an additional S/H clip e.g. 61H74M16H
                # The IF statement accounts for the 61H clip, but doesn't account for the 16H clip
            elif rightSoftClip > leftSoftClip:
                overlap = lengthSA - rightSoftClip
                sSeq = tools.revcomp(read.query_sequence[sum([x[1] for x in read.cigar if x[0] == 0 or x[0] == 1])-overlap-sLeftSoftClip:rightSoftClip-sLeftSoftClip+sum([x[1] for x in read.cigar if x[0] == 0 or x[0] == 1])])
                # This is the same on the right side, we need to account for the additional S/H clip on the left side e.g. 16H88M47H
                # The IF statement accounts for the 47H, but doesn't account for the 16H clip

        # #This was mainly used for testing purposes - These were reads with possible fusions
        # for check_read in bam.fetch(sChr, sPos, sEnd, multiple_iterators = True):
        #     if check_read.query_name == "consensus_read_CCCCCCATC_7_144982641_7_144982766_1" or \
        #         check_read.query_name == "consensus_read_ACCCTGCAG_4_93917621_4_93917755_1" or \
        #             check_read.query_name == "consensus_read_CCTAGACTA_12_91010217_12_91010350_1" or \
        #                 check_read.query_name == "consensus_read_GTCCTACTG_13_72645505_13_72645379_2" or \
        #                     check_read.query_name == "consensus_read_ATTACCTGG_13_72645397_13_72645539_1" or \
        #                         check_read.query_name == "consensus_read_TCACCTAAC_22_71977465_22_71977579_1":
        #         # These have singular mismatches between the read and the reference, won't really matter in the final script
        #         # But the chr13 and chr14 have Ns in their sequence... this MIGHT matter? We can account for when the N is a wildcard
        #         continue
        #     if check_read.query_name == read.query_name and check_read.reference_start+1 == sPos and check_read.cigarstring == sCigar:
        #         # if the read is on the same strand as the SA, then we can just take the sequence as is
        #         if check_read.query_alignment_sequence != sSeq:
        #             print(f"Error: {sSeq} != {check_read.query_alignment_sequence}", file=sys.stderr)
        #         break

        # e.g. chr1:13646284:chr1:15540216:r2:-+
        # AAAATTAACTGGGCATGGTGGTGTGTGCCTGTAATCCCAGCTACAGGTGCACGCCACCACGCCCAGCTAATTGTTGTATTTTTAGTAGAGACGGGGTTTCACCACGTTGGCCAGGATGGTCTCGATGTCTTGACCTCATGATCCACCCGCC
        #                                         CTACAGGTGCACGCCACCACGCCCAGCTAATTGTTGTATTTTTAGTAGAGACGGGGTTTCACCACGTTGGCCAGGATGGTCTCGATGTCTTGACCTCATGATCCACCCGCC
        # AAAATTAACTGGGCATGGTGGTGTGTGCCTGTAATCCCAGCTAC

        rSeq = read.query_sequence
        if read_strand == "-":
            pSeq = tools.revcomp(pSeq)
            rSeq = tools.revcomp(rSeq)
        if sStrand == "-":
            sSeq = tools.revcomp(sSeq)

        #if verbose:
        #    print(f"\t\tSequence: {rSeq}", file=sys.stderr)
        #    print(f"\t\tPrimary Sequence: {pSeq}", file=sys.stderr)
        #    print(f"\t\tSecondary Sequence: {sSeq}", file=sys.stderr)

        # This is the simple one to find overlap
        if leftSoftClip > rightSoftClip:
            comparison_type = STRAND_TO_COMPARISON[("leftSoftClip", StrandCombo(bndinfo['strands']))]
            #if verbose: print(f"\t\tleftSoftClip, {comparison_type}", file=sys.stderr)
        elif rightSoftClip > leftSoftClip:
            comparison_type = STRAND_TO_COMPARISON[("rightSoftClip", StrandCombo(bndinfo['strands']))]
            #if verbose: print(f"\t\trightSoftClip, {comparison_type}", file=sys.stderr)
        
        # Identify the possible Microhomology
        bndinfo['mh'] = find_largest_overlap(pSeq, sSeq, comparison_type)

        if bndinfo['strands'] == "--":
            # Because we are looking at the reverse strand, we need to reverse complement the MH
            bndinfo['mh'] = tools.revcomp(bndinfo['mh'])
        
        # If we cannot find any Microhomology, we will check the reference
        afterSecondaryAlignmentSeq, aSeq = "", ""
        if bndinfo['mh'] == "":
            if comparison_type == ComparisonType.START_TO_END:
                if sStrand == "+":
                    #if verbose: print(f"afterSecondaryAlignmentSeq ({sChr}:{sEnd}-{sEnd+(len(pSeq)/2)}): {afterSecondaryAlignmentSeq}", file=sys.stderr)
                    afterSecondaryAlignmentSeq = ref_fasta.fetch(sChr, sEnd, sEnd+(len(pSeq)/2))
                if sStrand == "-":
                    afterSecondaryAlignmentSeq = tools.revcomp(ref_fasta.fetch(sChr, sPos-(len(pSeq)/2), sPos))
                    #if verbose: print(f"afterSecondaryAlignmentSeq ({sChr}:{sPos-(len(pSeq)/2)}-{sPos}): {afterSecondaryAlignmentSeq}", file=sys.stderr)
                aSeq = afterSecondaryAlignmentSeq[0:20]
                bndinfo['mh'] = find_largest_overlap(pSeq, afterSecondaryAlignmentSeq, ComparisonType.START_TO_START)
            elif comparison_type == ComparisonType.END_TO_START:
                if sStrand == "+":
                    #if verbose:  print(f"{sChr}:{sPos-(len(pSeq)/2)}-{sPos}", file=sys.stderr)
                    afterSecondaryAlignmentSeq = ref_fasta.fetch(sChr, sPos-(len(pSeq)/2), sPos)
                    #if verbose: print(f"afterSecondaryAlignmentSeq ({sChr}:{sPos-(len(pSeq)/2)}-{sPos}): {afterSecondaryAlignmentSeq}", file=sys.stderr)
                if sStrand == "-":
                    #if verbose:  print(f"{sChr}:{sEnd}-{sEnd+(len(pSeq)/2)}", file=sys.stderr)
                    afterSecondaryAlignmentSeq = tools.revcomp(ref_fasta.fetch(sChr, sEnd, sEnd+(len(pSeq)/2)))
                    #if verbose: print(f"afterSecondaryAlignmentSeq ({sChr}:{sEnd}-{sEnd+(len(pSeq)/2)}): {afterSecondaryAlignmentSeq}", file=sys.stderr)
                aSeq = afterSecondaryAlignmentSeq[-20:]
                bndinfo['mh'] = find_largest_overlap(pSeq, afterSecondaryAlignmentSeq, ComparisonType.END_TO_END)

        #if verbose:
        #    print(f"\t\tMicrohomology: {indelinfo['mh']}", file=sys.stderr)
        #    # We need a way to visually make sure our MH are correct...
        #    print_sequence_alignment(
        #                rSeq=rSeq,
        #                pSeq=pSeq,
        #                sSeq=sSeq,
        #                mh=indelinfo['mh'],
        #                comparison_type=comparison_type,
        #                afterSeq=aSeq
        #            )
            
        # We add this here because we need to normalize things, but this is not FULLY representative of the strands...
        #indelinfo['strands'] = '++'if sStrand == read_strand else '+-'
        bndinfo['mh_length'] = len(bndinfo['mh']) if bndinfo['mh'] != '' else 0
    
    # These are the reads that have a supplementary alignment but do not pass the filters
    elif read_reference_start < end and read_reference_end > start:
       bndinfo['chrom'] = read.reference_name
       bndinfo['pos'] = start
       bndinfo['ref'] = '.'
       bndinfo['alt'] = '.'
       bndinfo['type'] = 'REF'

    # So there are reads that have SA tags, do not pass the filters, and do not overlap the region of interest
    # This is because our window may capture reads that have SAs elsewhere
    # We will return None for these
    else:
        return None

    return bndinfo  
