import pandas as pd
import pyranges as pr
import re
import pysam
import gzip
from tqdm import tqdm
from clint.textui import indent
import dnarecap.utils.logger as log
import biotite.sequence.align as align
import biotite.sequence as seq

def process_input_file(input_file, window, chromosome, verbose=False):
    """
    Load and validate the input file.
    """
    input_df = pd.read_csv(input_file)
    input_df['End'] = input_df['pam_site']
    input_df['Start'] = input_df['pam_site'] - 1
    input_df['guide_name'] = input_df['guide_name'].astype(str)
    input_df['info'] = input_df.apply(lambda row: f"{row['guide_name']},{row['target_sequence']},{row['PAM']},{row['Chromosome']},{row['pam_site']},{row['strand']}", axis=1)
    input_df['ontarget'] = input_df['on_target'] # This can be used because our multiguide RNA can have multiple off-target regions

    # make pyranges object
    input_pr = pr.PyRanges(input_df[['Chromosome','Start','End','guide_name','pam_site','info','ontarget']])

    # Cluster the intervals
    input_pr = input_pr.cluster(slack=window)
    merged_pr = input_pr.merge(by='Cluster', strand=False, slack=window)
    merged_df = merged_pr.df.join(input_pr.df.groupby('Cluster')['guide_name'].agg(list).reset_index().set_index('Cluster'),on='Cluster',how='left')
    merged_df = merged_df.join(input_pr.df.groupby('Cluster')['pam_site'].agg(list).reset_index().set_index('Cluster'),on='Cluster',how='left')
    merged_df = merged_df.join(input_pr.df.groupby('Cluster')['info'].agg(list).reset_index().set_index('Cluster'),on='Cluster',how='left')
    merged_df = merged_df.join(input_pr.df.groupby('Cluster')['ontarget'].agg(list).reset_index().set_index('Cluster'),on='Cluster',how='left')

    if chromosome is not None:
        log.logit(f"Processing ONLY chromosome {chromosome}")
        merged_df = merged_df[merged_df['Chromosome'] == chromosome]
    return merged_df

# make cigar tuples from a cigar string
def make_cigar_tuples(cigar):
    cigar_dict = {'M': 0, 'I': 1, 'D': 2, 'N': 3, 'S': 4, 'H': 5, 'P': 6, '=': 7, 'X': 8}
    # use regex to parse cigar string into tuples
    cigar_tuples = re.findall(r'(\d+)([MIDNSHP=X])', cigar)
    cigar_tuples = [(cigar_dict[operation],int(length)) for length, operation in cigar_tuples]
    return cigar_tuples

def cigar_summary(cigar):
    # use regex to parse cigar string into tuples
    cigar_tuples = re.findall(r'(\d+)([MIDNSHP=X])', cigar)
    cigar_sum = { 'M': 0, 'I': 0, 'D': 0, 'N': 0, 'S': 0, 'H': 0, 'P': 0, '=': 0, 'X': 0 }
    for length, operation in cigar_tuples:
        cigar_sum[operation] += int(length)

    return cigar_sum

def revcomp(seq):
    tab = str.maketrans('ACGTacgtRYMKrymkVBHDvbhd', 'TGCAtgcaYRKMyrkmBVDHbvdh') # maketrans <- maps the reverse complement
    return seq.translate(tab)[::-1] # translate(x)[::-1] <- works backward through the string, effectively reversing the string

# complement function
def comp(sequence):
    complement = {
        'A': 'T', 'T': 'A', 
        'C': 'G', 'G': 'C',
        'a': 't', 't': 'a', 
        'c': 'g', 'g': 'c',
        'N': 'N', 'n': 'n'
    }
    return ''.join(complement.get(base, base) for base in sequence[::-1])

def process_soft_clipping(read, cigar, leftSoftClip, rightSoftClip, read_strand, read_reference_start, read_reference_end, regionStart, regionSeq, window, maxNM, svDistanceThreshold, verbose):
    """
    Adjust alignment cigars for soft clips and supplementary alignments.
    """
    minSecMapQual = 10
    minSoftClipLength = 5

    # Unused parameters
    saAlignmentTolerance = 5
    deletionDistanceThreshold = 1_000

    # Get left soft clip length
    if cigar[0][0] == 4: leftSoftClip = cigar[0][1]
    # Get right soft clip length
    if cigar[-1][0] == 4: rightSoftClip = cigar[-1][1]

    # This adjusts alignment cigars for soft clips and supplementary alignments. 
    # It extends the cigar if the softclips map nearby and therefore the read has a deletion.

    # Process left soft clip, if present
    if leftSoftClip > rightSoftClip:
        # First check for local supplementary alignments 
        sChr, sPos, sStrand, sCigar, sMq, sNm = read.get_tag('SA').split(';')[0].split(',') if read.has_tag('SA') else [None,None,None,None,None,None]

        if (sChr != None and int(sMq)>=minSecMapQual and
            int(sNm)<maxNM and sChr == read.reference_name and
            sStrand == read_strand and int(sPos) < read_reference_start and
            abs(int(sPos) - read_reference_start) < svDistanceThreshold):
            
            sPos = int(sPos)
            sCigarTuples = make_cigar_tuples(sCigar)

            if sCigarTuples[-1][0] >= 4 and sCigarTuples[0][0] == 0 and sCigarTuples[0][1] == cigar[0][1]:
                cigar = sCigarTuples[:-1] + [(2,read_reference_start-sPos-1)] + cigar[1:]

            else:
                # If there is a supplementary alignment, then realign the entire read to the reference sequence
                cigar = make_cigar_tuples(align.write_alignment_to_cigar(align.align_optimal(seq.NucleotideSequence(regionSeq[sPos-regionStart:read_reference_end-regionStart+1]), #seq.NucleotideSequence(fasta.fetch(chr,sPos-1,read_reference_end)),
                                                                            seq.NucleotideSequence(read.query_sequence),
                                                                            matrix=align.SubstitutionMatrix.std_nucleotide_matrix(),
                                                                            gap_penalty=(-10,-1),local=False,max_number=1)[0]))
            read_reference_start = sPos

        else:
            # If the left soft clip is long enough, see if it supports a deletion within window
            # this performs an exact match and so isnt ideal, but should be good enough for now
            if leftSoftClip >= minSoftClipLength:
                clippedSeq = read.query_sequence[:leftSoftClip]
                refSeq = regionSeq[read_reference_start - regionStart - window + 1:read_reference_start - regionStart +  1] #fasta.fetch(chr,read_reference_start - window,read_reference_start)
                # Find the position of the left-most occurance of clippedSeq in refSeq
                leftClipPos = refSeq.rfind(clippedSeq)
                # Irenaeus Edit: leftClipPos < len(refSeq) - len(clippedSeq) to avoid cases where the clippedSeq is found at the very end and there is no detectable deletion
                if leftClipPos != -1 and leftClipPos < len(refSeq) - len(clippedSeq):
                    delSeq = refSeq[leftClipPos + len(clippedSeq)-1:]
                    cigar = [(0,len(clippedSeq))] + [(2,len(delSeq))] + cigar[1:]
                    read_reference_start = read_reference_end - sum([x[1] for x in cigar if x[0] == 0 or x[0] == 2 or x[0] == 3]) + 1

    # Process right soft clip, if present
    if rightSoftClip > leftSoftClip:
        # first check for local supplementary alignments 
        sChr, sPos, sStrand, sCigar, sMq, sNm = read.get_tag('SA').split(';')[0].split(',') if read.has_tag('SA') else [None,None,None,None,None,None]

        if (sChr != None):
            sPos = int(sPos)
            sCigarTuples = make_cigar_tuples(sCigar)

            if (int(sMq)>=minSecMapQual and sCigarTuples[0][0] >= 4 and sCigarTuples[-1][0] == 0 and
                int(sNm)<maxNM and sChr == read.reference_name and 
                sStrand == read_strand and int(sPos) > read_reference_end and 
                int(sPos) - read_reference_end < svDistanceThreshold):
            
                if sCigarTuples[0][0] >= 4 and sCigarTuples[-1][0] == 0 and sCigarTuples[-1][1] == cigar[-1][1]:
                    cigar = cigar[:-1] + [(2,sPos-read_reference_end)] + sCigarTuples[1:]

                else:
                    cigar = make_cigar_tuples(align.write_alignment_to_cigar(align.align_optimal(seq.NucleotideSequence(regionSeq[read_reference_start-regionStart:sPos-regionStart+sum([x[1] for x in sCigarTuples if x[0] == 0 or x[0] == 2 or x[0] == 3])]), #fasta.fetch(chr,read_reference_start-1,sPos-1+sum([x[1] for x in sCigarTuples if x[0] == 0 or x[0] == 2 or x[0] == 3]))),
                                                                                seq.NucleotideSequence(read.query_sequence),
                                                                                matrix=align.SubstitutionMatrix.std_nucleotide_matrix(),
                                                                                gap_penalty=(-10,-1),local=False,max_number=1)[0]))

            read_reference_end = read_reference_start + sum([x[1] for x in cigar if x[0] == 0 or x[0] == 2 or x[0] == 3])

        else:
            # If the left soft clip is long enough, see if it supports a deletion within deletionDistanceThreshold
            # Basically, if there is a deletion that is being mislabeled as a soft clip, this will find it
            # this performs an exact match and so isnt ideal, but should be good enough for now
            if rightSoftClip >= minSoftClipLength:                    
                clippedSeq = read.query_sequence[-rightSoftClip:]
                refSeq = regionSeq[read_reference_end-regionStart:read_reference_end-regionStart + window] #fasta.fetch(chr,read_reference_end,read_reference_end + window)
                # Find the position of the left-most occurance of clippedSeq in refSeq
                rightClipPos = refSeq.find(clippedSeq)
                # Irenaeus Edit: rightClipPos > 0 to avoid cases where the clippedSeq is found at position 0 and there is no detectable deletion
                if rightClipPos != -1 and rightClipPos > 0:
                    delSeq = refSeq[0:rightClipPos]
                    cigar = cigar[:-1] + [(2,len(delSeq))] + [(0,len(clippedSeq))]
                    read_reference_end = read_reference_start + sum([x[1] for x in cigar if x[0] == 0 or x[0] == 2 or x[0] == 3])

    return cigar, leftSoftClip, rightSoftClip, read_reference_start, read_reference_end

# Position-based filtering
def get_read_variant_pairs(reads, variants):
    pairs = []
    for read in reads:
        for idx, variant in variants.iterrows():
            # Check position overlap
            if (read.reference_start <= variant.pos + len(variant.ref) + 50 and 
                read.reference_end >= variant.pos - 50):
                pairs.append((read, idx))
    return pairs

def quick_cigar_check(read, variant_ref, variant_alt):
    """Quick check if read CIGAR could support this variant without alignment"""
    expected_del_size = len(variant_ref) - len(variant_alt)
    expected_ins_size = len(variant_alt) - len(variant_ref)
    
    if expected_del_size > 0:  # Deletion variant
        # Check if read has deletion of approximately the right size
        cigar_deletions = sum(length for op, length in read.cigartuples if op == 2)
        return abs(cigar_deletions - expected_del_size) <= 2  # Allow some tolerance
    
    elif expected_ins_size > 0:  # Insertion variant
        # Check if read has insertion of approximately the right size
        cigar_insertions = sum(length for op, length in read.cigartuples if op == 1)
        return abs(cigar_insertions - expected_ins_size) <= 2
    
    return True  # For other cases, allow through

def batch_align_reads(df, read_variant_pairs, alignment_matrix, gap_penalty):
    """Process alignments in batches to reduce object creation overhead"""
    results = []
    
    # Pre-convert sequences to reduce repeated conversions
    sequence_cache = {}
    
    for read, variant_idx in tqdm(read_variant_pairs, desc="Aligning reads to variants", unit="pairs"):
        variant = df.loc[variant_idx]
        
        # Cache sequence objects, if the sequence already exists, reuse it
        if read.query_name not in sequence_cache:
            sequence_cache[read.query_name] = seq.NucleotideSequence(read.query_sequence)
        
        ref_key = f"{variant_idx}_ref"
        alt_key = f"{variant_idx}_alt"
        
        if ref_key not in sequence_cache:
            sequence_cache[ref_key] = seq.NucleotideSequence(variant.refseq)
            sequence_cache[alt_key] = seq.NucleotideSequence(variant.altseq)
        
        read_seq = sequence_cache[read.query_name]
        ref_seq = sequence_cache[ref_key]
        alt_seq = sequence_cache[alt_key]
        
        # Quick pre-screen
        if not quick_cigar_check(read, variant.ref, variant.alt):
            continue
            
        # Perform alignments
        ref_align = align.align_optimal(ref_seq, read_seq, 
                                      matrix=alignment_matrix, gap_penalty=gap_penalty,
                                      local=True, terminal_penalty=False, max_number=1)[0]
        alt_align = align.align_optimal(alt_seq, read_seq,
                                      matrix=alignment_matrix, gap_penalty=gap_penalty,
                                      local=True, terminal_penalty=False, max_number=1)[0]
        
        results.append((variant_idx, read, ref_align, alt_align))
    
    return results

# Function to count the number of reads that support an indel or BND event
def add_normal_counts(df, control_bam_file, fasta, chrom, start, end, window, verbose, handicap=5, flank=150):
    """
    Adding control read counts supporting each INDEL/BND from control BAM file.
    """
    df = df.copy()

    # Early exit for empty input
    if len(df) == 0:
        return df
    
    # Filter out REF variants
    df = df[df['type'] != 'REF'].copy()
    if len(df) == 0:
        return df

    log.logit(f"Getting INDEL/BNDs counts from Control BAM file for interval {chrom}:{start}-{end}")

    control = pysam.AlignmentFile(control_bam_file, "rc", reference_filename=fasta)
    ref_fasta = pysam.FastaFile(fasta)

    # Make REF and ALT sequences for each INDEL/BND
    df['refseq'] = ''
    df['altseq'] = ''
    
    for it, row in df.iterrows():
        #df.at[it,'refseq'] = ref_fasta.fetch(row['chrom'],row['pos']-1-flank,row['pos']+len(row['ref'])+flank)
        refseq = ref_fasta.fetch(row['chrom'],row['pos']-1-flank,row['pos']+len(row['ref'])+flank)
        df.at[it,'refseq'] = refseq
        
        if row['type'] in ['INDEL','DEL','INS']:
            #df.at[it,'altseq'] = ref_fasta.fetch(row['chrom'],row['pos']-1-flank,row['pos']-1) + row['alt'] + ref_fasta.fetch(row['chrom'],row['pos']+len(row['ref'])-1,row['pos']+len(row['ref'])+flank)
            altseq = ref_fasta.fetch(row['chrom'],row['pos']-1-flank,row['pos']-1) + row['alt'] + ref_fasta.fetch(row['chrom'],row['pos']+len(row['ref'])-1,row['pos']+len(row['ref'])+flank)
            df.at[it,'altseq'] = altseq

        elif row['type'] == 'BND':
            # Adjust pos2 if its too close to the edge of the reference
            pos2 = row['pos2']
            if pos2-1-flank < 0:
                pos2 = flank + 1
            elif pos2 + flank > ref_fasta.get_reference_length(row['chrom2']):
                pos2 = ref_fasta.get_reference_length(row['chrom2']) - flank

            #df.at[it,'altseq'] = ref_fasta.fetch(row['chrom'],row['pos']-1-flank,row['pos']-1) + (ref_fasta.fetch(row['chrom2'],pos2-1,pos2+flank) if row['strands']=='++' else revcomp(ref_fasta.fetch(row['chrom2'],pos2-1-flank,pos2)))
            altseq = ref_fasta.fetch(row['chrom'],row['pos']-1-flank,row['pos']-1) + (ref_fasta.fetch(row['chrom2'],pos2-1,pos2+flank) if row['strands']=='++' else revcomp(ref_fasta.fetch(row['chrom2'],pos2-1-flank,pos2)))
            df.at[it,'altseq'] = altseq

    # PRE-FILTER: Group variants by type for faster processing
    deletion_variants = df[(df['type'] == 'INDEL') & (df['ref'].str.len() > df['alt'].str.len())].copy()
    insertion_variants = df[(df['type'] == 'INDEL') & (df['alt'].str.len() > df['ref'].str.len())].copy()
    bnd_variants = df[df['type'] == 'BND'].copy()

    # Early exit if no variants to check
    if len(deletion_variants) + len(insertion_variants) + len(bnd_variants) == 0:
        df['control_alt_counts'] = 0
        df['control_total_counts'] = 0
        return df

    all_reads = list(control.fetch(contig=chrom, start=start-window, end=end, multiple_iterators = True))

    # Multi-stage filtering
    quality_reads = [read for read in all_reads if not (
        read.is_mapped is False or read.is_duplicate is True or 
        read.is_secondary is True or read.is_supplementary is True or 
        read.mapping_quality == 0)]
    
    # Type-specific read filtering
    deletion_reads = [r for r in quality_reads if 'D' in r.cigarstring] if len(deletion_variants) > 0 else []
    insertion_reads = [r for r in quality_reads if 'I' in r.cigarstring] if len(insertion_variants) > 0 else []
    bnd_reads = [r for r in quality_reads if r.has_tag('SA')] if len(bnd_variants) > 0 else []

    log.logit(f"Filtered to {len(deletion_reads)} deletion reads, "f"{len(insertion_reads)} insertion reads, {len(bnd_reads)} BND reads")

    # Checks to see if the control read and variant could possibly match based on position
    # This is because if the control read is not overlapping the variant position +/- 50bp, it can't possibly support the variant
    deletion_pairs = get_read_variant_pairs(deletion_reads, deletion_variants)
    insertion_pairs = get_read_variant_pairs(insertion_reads, insertion_variants)
    bnd_pairs = get_read_variant_pairs(bnd_reads, bnd_variants)

    total_pairs = len(deletion_pairs) + len(insertion_pairs) + len(bnd_pairs)
    log.logit(f"Reduced to {total_pairs} read-variant pairs for alignment")

    df['control_alt_counts'] = 0
    df['control_total_counts'] = len(set(r.query_name for r in quality_reads))

    # Pre-compute alignment matrix and gap penalty (for all alignments)
    alignment_matrix = align.SubstitutionMatrix.std_nucleotide_matrix()
    gap_penalty = (-10,-1)

    # Process each type separately
    all_pairs = deletion_pairs + insertion_pairs + bnd_pairs

    # Batch process all read-variant pairs
    if len(all_pairs) > 0:
        results = batch_align_reads(df, all_pairs, alignment_matrix, gap_penalty)

        # Process alignment results
        for variant_idx, read, ref_align, alt_align in results:
            if alt_align.score > 0 and alt_align.score - handicap > ref_align.score:
                ref_align_cigar_sum = cigar_summary(align.write_alignment_to_cigar(ref_align))
                
                variant = df.loc[variant_idx]
                expected_del = len(variant.ref) - len(variant.alt)
                expected_ins = len(variant.alt) - len(variant.ref)
                
                if ((expected_del > 0 and ref_align_cigar_sum['D'] == expected_del) or 
                    (expected_ins > 0 and ref_align_cigar_sum['I'] == expected_ins)):
                    df.at[variant_idx, 'control_alt_counts'] += 1

    control.close()
    ref_fasta.close()

    return df.copy()

def summarise_repaired_reads(breaks, summary_outfile, input):
    log.logit(f"Summarising INDELs/BNDs from {input['guide_name']} to {summary_outfile}")
    total_reads, repaired_reads, control_total_reads, control_repaired_reads = 0, 0, 0, 0
    repaired_fraction, control_repaired_fraction = 0, 0
    indel_keys, bnd_keys = '.', '.'

    indel_count, bnd_count = 0, 0
    num_ins_nhej, num_ins_tins, num_ins_nmh_ej, num_ins_other = 0, 0, 0, 0
    num_del_nhej, num_del_tmej, num_del_nmh_ej, num_del_other = 0, 0, 0, 0
    num_indel_nhej, num_indel_tmej, num_indel_nmh_ej, num_indel_complex, num_indel_other = 0, 0, 0, 0, 0
    num_bnd_tmej, num_bnd_nmh_ej, num_bnd_other = 0, 0, 0

    indel_nhej_fraction, indel_tmej_fraction, indel_nmh_ej_fraction, indel_complex_fraction, indel_other_fraction = 0, 0, 0, 0, 0
    indel_ins_nhej_fraction, indel_ins_tins_fraction, indel_ins_nmh_ej_fraction, indel_ins_other_fraction = 0, 0, 0, 0
    indel_del_nhej_fraction, indel_del_tmej_fraction, indel_del_nmh_ej_fraction, indel_del_other_fraction = 0, 0, 0, 0
    bnd_tmej_fraction, bnd_nmh_ej_fraction, bnd_other_fraction = 0, 0, 0

    refs = breaks[breaks['type']=='REF'].copy()
    indels = breaks[breaks['type']=='INDEL'].copy()
    bnds = breaks[breaks['type']=='BND'].copy()
    all = pd.concat([refs, indels, bnds], ignore_index=True)

    # Handle empty DataFrame or missing columns
    if len(all) == 0:
        log.logit("WARNING: No data found for summarization", color = 'yellow')
        # Set all values to 0 and write empty summary
        total_reads = 0
        repaired_reads = 0
        control_repaired_reads = 0
        control_total_reads = 0
    else:
        total_reads = sum(all['counts'])
        repaired_reads = sum(all[all['type']!='REF']['counts'])
        # Handle NaN values safely
        control_alt_mean = all['control_alt_counts'].mean() if 'control_alt_counts' in all.columns else 0
        control_total_mean = all['control_total_counts'].mean() if 'control_total_counts' in all.columns else 0
        
        control_repaired_reads = int(control_alt_mean) if not pd.isna(control_alt_mean) else 0
        control_total_reads = int(control_total_mean) if not pd.isna(control_total_mean) else 0

    # If we don't have any INDELs, we can skip this
    if len(indels) > 0:
        indel_count = len(indels[indels['type']!='REF'])

        # DNA damage results
        num_indel_nhej = sum(indels['classification'] == 'NHEJ')
        num_indel_tmej = sum(indels['classification'] == 'TMEJ')
        num_indel_nmh_ej = sum(indels['classification'] == 'Non-MH EJ')
        num_indel_complex = sum(indels['classification'] == 'Complex')
        num_indel_other = sum(indels['classification'] == 'Other')

        # Insertions
        num_ins_nhej = sum((indels['indel_type'] == "Ins") & (indels['classification'] == 'NHEJ'))
        num_ins_tins = sum((indels['indel_type'] == "Ins") & (indels['classification'] == 'TINS'))
        num_ins_nmh_ej = sum((indels['indel_type'] == "Ins") & (indels['classification'] == 'Non-MH EJ'))
        num_ins_other = sum((indels['indel_type'] == "Ins") & (indels['classification'] == 'Other'))

        # Deletions
        num_del_nhej = sum((indels['indel_type'] == "Del") & (indels['classification'] == 'NHEJ'))
        num_del_tmej = sum((indels['indel_type'] == "Del") & (indels['classification'] == 'TMEJ'))
        num_del_nmh_ej = sum((indels['indel_type'] == "Del") & (indels['classification'] == 'Non-MH EJ'))
        num_del_other = sum((indels['indel_type'] == "Del") & (indels['classification'] == 'Other'))

        # Percentage Classification
        indel_nhej_fraction = round(num_indel_nhej/total_reads,4) if total_reads > 0 else 0
        indel_tmej_fraction = round(num_indel_tmej/total_reads,4) if total_reads > 0 else 0
        indel_nmh_ej_fraction = round(num_indel_nmh_ej/total_reads,4) if total_reads > 0 else 0
        indel_complex_fraction = round(num_indel_complex/total_reads,4) if total_reads > 0 else 0
        indel_other_fraction = round(num_indel_other/total_reads,4) if total_reads > 0 else 0

        indel_ins_nhej_fraction = round(num_ins_nhej/total_reads,4) if total_reads > 0 else 0
        indel_ins_tins_fraction = round(num_ins_tins/total_reads,4) if total_reads > 0 else 0
        indel_ins_nmh_ej_fraction = round(num_ins_nmh_ej/total_reads,4) if total_reads > 0 else 0
        indel_ins_other_fraction = round(num_ins_other/total_reads,4) if total_reads > 0 else 0

        indel_del_nhej_fraction = round(num_del_nhej/total_reads,4) if total_reads > 0 else 0
        indel_del_tmej_fraction = round(num_del_tmej/total_reads,4) if total_reads > 0 else 0
        indel_del_nmh_ej_fraction = round(num_del_nmh_ej/total_reads,4) if total_reads > 0 else 0
        indel_del_other_fraction = round(num_del_other/total_reads,4) if total_reads > 0 else 0

    if len(bnds) > 0:
        bnd_count = len(bnds[bnds['type']!='REF'])

        # DNA damage results
        num_bnd_tmej = sum(bnds['classification'] == 'TMEJ')
        num_bnd_nmh_ej = sum(bnds['classification'] == 'Non-MH EJ')
        num_bnd_other = sum(bnds['classification'] == 'Other')

        bnd_tmej_fraction = round(num_bnd_tmej/len(bnds),4) if len(bnds) > 0 else 0
        bnd_nmh_ej_fraction = round(num_bnd_nmh_ej/len(bnds),4) if len(bnds) > 0 else 0
        bnd_other_fraction = round(num_bnd_other/len(bnds),4) if len(bnds) > 0 else 0
    
    repaired_fraction = round(repaired_reads/total_reads,4) if total_reads > 0 else 0
    control_repaired_fraction = round(control_repaired_reads/control_total_reads,4) if control_total_reads > 0 else 0
    ontarget = ';'.join(input['ontarget']) if isinstance(input['ontarget'], list) else input['ontarget']
    offtargetsites = ';'.join(input['info']) if isinstance(input['info'], list) else input['info']
    
    indel_keys = ';'.join(indels['key'].tolist()) if len(indels) > 0 else '.'
    bnd_keys = ';'.join(bnds['key'].tolist()) if len(bnds) > 0 else '.'

    # Create the One Line Summary
    out_string=f"{input['Chromosome']}\t{input['Start']}\t{input['End']}\t{';'.join([ str(x) for x in input['pam_site'] ])}\t{total_reads}\t{repaired_reads}\t{repaired_fraction}\t{control_total_reads}\t{control_repaired_reads}\t{control_repaired_fraction}\t{indel_count}\t{indel_keys}\t{bnd_count}\t{bnd_keys}\t{offtargetsites}\t{ontarget}\t\
        {num_indel_nhej}\t{num_indel_tmej}\t{num_indel_nmh_ej}\t{num_indel_other}\t{num_indel_complex}\t{indel_nhej_fraction}\t{indel_tmej_fraction}\t{indel_nmh_ej_fraction}\t{indel_other_fraction}\t{indel_complex_fraction}\t\
        {num_ins_nhej}\t{num_ins_tins}\t{num_ins_nmh_ej}\t{num_ins_other}\t{indel_ins_nhej_fraction}\t{indel_ins_tins_fraction}\t{indel_ins_nmh_ej_fraction}\t{indel_ins_other_fraction}\t\
        {num_del_nhej}\t{num_del_tmej}\t{num_del_nmh_ej}\t{num_del_other}\t{indel_del_nhej_fraction}\t{indel_del_tmej_fraction}\t{indel_del_nmh_ej_fraction}\t{indel_del_other_fraction}\t\
        {num_bnd_tmej}\t{num_bnd_nmh_ej}\t{num_bnd_other}\t{bnd_tmej_fraction}\t{bnd_nmh_ej_fraction}\t{bnd_other_fraction}"
    open(summary_outfile, 'a').write(out_string + "\n")

def process_time(total_elapsed):
    """
    Format elapsed time into a human-readable string.
    """ 
    if total_elapsed < 60:
        time_str = f"{total_elapsed:.2f} seconds"
    elif total_elapsed < 3600:
        minutes = int(total_elapsed // 60)
        seconds = total_elapsed % 60
        time_str = f"{minutes}m {seconds:.1f}s"
    else:
        hours = int(total_elapsed // 3600)
        minutes = int((total_elapsed % 3600) // 60)
        seconds = total_elapsed % 60
        time_str = f"{hours}h {minutes}m {seconds:.1f}s"
    return time_str

def annotate_vcf_with_scars(vcf_file, df_with_scars, chromosome, outfile):
    """
    Annotate VCF file with scar information from the DataFrame.
    """
    # Create lookup dictionary from scarmapper results
    scar_lookup = {}
    for _, row in df_with_scars.iterrows():
        key = f"{row['chrom']}:{row['pos']}:{row['ref']}:{row['alt']}"
        scar_lookup[key] = {
            'classification': row['classification'],
            'mh': '.' if (row['mh'] == '' or row['mh'] is None) else row['mh'],
            'lft_tmplt': '.' if (row['lft_tmplt'] == '' or row['lft_tmplt'] is None) else row['lft_tmplt'],
            'rt_tmplt': '.' if (row['rt_tmplt'] == '' or row['rt_tmplt'] is None) else row['rt_tmplt']
        }

    input_vcf = pysam.VariantFile(vcf_file, 'r')
    input_vcf.header.info.add("DNARECAP_CLASS", number='A', type='String', description="DNA repair classification (NHEJ, TMEJ, etc.)")
    input_vcf.header.info.add("DNARECAP_MH", number='A', type='String', description="Microhomology sequence")
    input_vcf.header.info.add("DNARECAP_LFT_TMPLT", number='A', type='String', description="Left templated insertion")
    input_vcf.header.info.add("DNARECAP_RT_TMPLT", number='A', type='String', description="Right templated insertion")
    
    # Create output filename
    if outfile.endswith('.tsv'):
        output_vcf = outfile.replace('.tsv', '.annotated.vcf.gz')
    elif outfile.endswith('.vcf'):
        output_vcf = outfile + '.gz'
    else:
        output_vcf = outfile + '.annotated.vcf.gz'
    vcf_out = pysam.VariantFile(output_vcf, 'wz', header=input_vcf.header)

    for variant in input_vcf:
        classifications = []
        mh_seqs = []
        lft_tmplts = []
        rt_tmplts = []

        for alt in variant.alts:
            key = f"{variant.chrom}:{variant.pos}:{variant.ref}:{alt}"
            if key in scar_lookup:
                classifications.append(scar_lookup[key]['classification'])
                mh_seqs.append(scar_lookup[key]['mh'])
                lft_tmplts.append(scar_lookup[key]['lft_tmplt'])
                rt_tmplts.append(scar_lookup[key]['rt_tmplt'])
            else:
                # No scar data found (likely SNP)
                classifications.append('.')
                mh_seqs.append('.')
                lft_tmplts.append('.')
                rt_tmplts.append('.')
        
        # Add annotations
        if any(c != '.' for c in classifications):
            variant.info['DNARECAP_CLASS'] = classifications
            if "TMEJ" in classifications:
                variant.info['DNARECAP_MH'] = mh_seqs
            elif "TINS" in classifications:
                if any(lt != '' for lt in lft_tmplts):
                    variant.info['DNARECAP_LFT_TMPLT'] = lft_tmplts
                if any(rt != '' for rt in rt_tmplts):
                    variant.info['DNARECAP_RT_TMPLT'] = rt_tmplts
        
        vcf_out.write(variant)
    input_vcf.close()
    vcf_out.close()

    log.logit(f"Annotated VCF written to {output_vcf}")