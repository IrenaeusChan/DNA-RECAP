from __future__ import division
from clint.textui import indent
from tqdm import tqdm
import pysam
import os, sys, gzip
import multiprocessing as mp
import pandas as pd

import dnarecap.utils.logger as log
import dnarecap.utils.tools as tools
import dnarecap.process.scarmapper as scarmapper
import dnarecap.process.bnds as bnds
import dnarecap.process.indels as indels

# Create different keys based on variant type
def create_variant_key(row):
    if row['type'] == 'BND':
        return f"{row['chrom']}:{row['pos']}:{row['chrom2']}:{row['pos2']}:{row['strands']}:{row['counts']}:{row['control_alt_counts']}:{row['Distance']}:{row['mh_length']}:{row['classification']}"
    elif row['type'] == 'REF':
        return f"{row['chrom']}:{row['pos']}:REF:{row['counts']}:{row['control_alt_counts']}"
    else:  # INDEL types
        return f"{row['chrom']}:{row['pos']}:{row['ref']}:{row['alt']}:{row['counts']}:{row['control_alt_counts']}:{row['Distance']}:{row['mh_length']}:{row['classification']}"

# Adapted from Dr. Dave Spencer
def process_regions(bam_file, fasta, source, chrom, start, end, pam_site, window, distance, minreads, maxcontrol, verbose=False):
    """
    Process a BAM file to identify INDELs/BNDs within a specified genomic interval.
    """
    log.logit(f"Processing interval {chrom}:{start-window}-{end+window}")
    # Open BAM files
    bam = pysam.AlignmentFile(bam_file, "rc", reference_filename=fasta)
    # Open fasta file
    ref_fasta = pysam.FastaFile(fasta)

    # Define constants
    svDistanceThreshold = 100_000
    maxNM = 5
    minSecMapQual = 10

    # Initialize empty DataFrame to store read alignments
    indels_df = pd.DataFrame(columns=['source', 'read', 'chrom', 'pos', 'ref', 'alt', 'strands', 'type'])
    bnds_df = pd.DataFrame(columns=['source','read','chrom','pos','chrom2','pos2','strands','type'])

    # Grab reference sequence for region +/- svDistanceThreshold
    regionStart = max(start - svDistanceThreshold, 1)
    regionEnd = min(end + svDistanceThreshold, ref_fasta.get_reference_length(chrom))
    regionSeq = ref_fasta.fetch(chrom,regionStart-1,regionEnd)

    # Get reads that align within window of start and end
    all_reads = list(bam.fetch(chrom, start-window, end+window))
    # skip if not primary alignment or a duplicate or alignment doesnt overlap start,end
    valid_reads = [read for read in all_reads if not (
        read.is_mapped is False or read.is_duplicate is True or 
        read.is_secondary is True or read.is_supplementary is True or 
        read.mapping_quality == 0)]
    
    with indent(4, quote=" >"):
        for read in tqdm(valid_reads):#, desc=f"Processing read in {read.query_name}", disable=not verbose, leave=False):
            # Get cigar information
            cigar = read.cigartuples # get cigar info
            mate_cigar = tools.make_cigar_tuples(read.get_tag('MC')) if read.has_tag('MC') else None

            # Quality filter: Get mismatches from NM tag and subtract deleted/skipped bases from cigar and skip if too many mismatches
            if read.has_tag('NM') and int(read.get_tag('NM')) > maxNM:
                if (int(read.get_tag('NM')) - sum([x[1] for x in cigar if x[0] == 2 or x[0]==3])) > maxNM:
                    continue
            
            # Define read defaults
            read_strand = "+" if read.is_reverse is False else "-"
            leftSoftClip = 0
            rightSoftClip = 0
            # This is the start of the current read of interest
            read_reference_start = read.reference_start + 1
            read_reference_end = read.reference_end
            # This is the start of the mate read of interest
            read_next_reference_start = read.next_reference_start + 1
            read_next_reference_end = None if mate_cigar is None else read_next_reference_start + sum([x[1] for x in mate_cigar if x[0] == 0 or x[0] == 2 or x[0] == 3])
            read_query_start = 0
            read_query_end = len(read.query_sequence)
            
            cigar, leftSoftClip, rightSoftClip, read_reference_start, read_reference_end = tools.process_soft_clipping(read, cigar, leftSoftClip, rightSoftClip, read_strand, read_reference_start, read_reference_end, regionStart, regionSeq, window, maxNM, svDistanceThreshold, verbose)

            # Adjusted read has to be within distance of the start/end positions to consider
            if read_reference_start - leftSoftClip >= end + distance or read_reference_end + rightSoftClip <= start - distance:
                continue

            # Remove leading or trailing soft clips
            if cigar[0][0] >= 4:
                read_query_start += cigar[0][1]
                del cigar[0]
            
            if cigar[-1][0] >= 4:
                read_query_end -= cigar[-1][1]
                del cigar[-1]

            # There are multiple cigar operations bookened by matches, then process an indel                
            if len(cigar) > 1:
                cigar, read_reference_start, read_reference_end, read_query_start, read_query_end = indels.trim_cigar_ends(cigar, read_reference_start, read_reference_end, read_query_start, read_query_end)
                # Initialize indelinfo and bndinfo
                indelinfo = indels.read_to_indelinfo(read, cigar, regionSeq, regionStart, read_reference_start, read_query_start, verbose)
                indels_df = pd.concat([indels_df, pd.DataFrame([indelinfo])], ignore_index=True)
            # If the read has a supplementary alignment and soft clipping, process as BND
            elif read.has_tag('SA') is True and (leftSoftClip > 0 or rightSoftClip > 0):
                bndinfo = bnds.read_to_bndinfo(read, read_strand, start, end, read_reference_start, read_reference_end, leftSoftClip, rightSoftClip, ref_fasta, maxNM, svDistanceThreshold, minSecMapQual, verbose)
                if bndinfo is not None:
                    bnds_df = pd.concat([bnds_df, pd.DataFrame([bndinfo])], ignore_index=True)
                else:
                    log.logit(f"Read {read.query_name} has SA tag but did not meet criteria for BND, skipping")
            # If no indels but read overlaps region, mark as REF, these are reads without indication of "damage"
            elif read_reference_start < end and read_reference_end > start:
                indelinfo = {'source':'',
                    'read': read.query_name,
                    'chrom': chrom,
                    'pos': start,
                    'whichread': 'r1' if read.is_read1 else 'r2',
                    'strands': '',
                    'ref':'.',
                    'alt':'.',
                    'type':'REF',
                    'mh':'',
                    'mh_length': 0,
                    'lft_tmplt':'',
                    'rt_tmplt':'',
                    'classification':''}
                indels_df = pd.concat([indels_df, pd.DataFrame([indelinfo])], ignore_index=True)

            else:
                continue

    log.logit(f"Summarizing INDELs from reads")
    if verbose: log.logit(f"Total reads processed in region: {len(indels_df)}")
    # Getting unique indels and their counts
    # Sorting and ONLY keeping the first occurrence of each read
    indels_df = indels_df.sort_values(by=['read','chrom','pos','strands','ref','alt','type','mh'],key=lambda col: col != '',ascending=False).groupby('read').first().reset_index()
    # Count how many times each UNIQUE INDEL occurs
    indels_df = indels_df.groupby(['read','chrom','pos','strands','ref','alt','type','mh'],dropna=False).size().reset_index(name='counts')
    # Sorting by chrom, pos for easier viewing
    indels_df = indels_df.sort_values(by=['chrom','pos']).reset_index()

    # Ensure 'pos' is of integer type
    indels_df['pos'] = indels_df['pos'].astype(pd.Int64Dtype())
    # Adding distance from cut site
    indels_df['Positions'] = [pam_site] * len(indels_df)
    indels_df['Distance'] = indels_df.apply(lambda r: min([abs(r['pos'] - x) for x in r['Positions']] + [abs(r['pos'] + len(r['ref']) - 1 - x) for x in r['Positions']] ) if r['pos'] is not pd.NA else pd.NA,axis=1)
    indels_df['source'] = [source] * len(indels_df)

    # For consistency with BNDs, add chrom2 and pos2 columns with empty values
    indels_df['chrom2'] = ''
    indels_df['pos2'] = pd.NA
    
    if len(bnds_df) > 0:
        log.logit(f"Summarizing BNDs from reads")
        # Getting unique bnds and their counts
        # Sorting and ONLY keeping the first occurrence of each read
        bnds_df = bnds_df.sort_values(by=['read','chrom','pos','chrom2','pos2','strands','type','mh'],key=lambda col: col != '',ascending=False).groupby('read').first().reset_index()
        # Count how many times each UNIQUE INDEL occurs
        bnds_df = bnds_df.groupby(['read','chrom','pos','chrom2','pos2','strands','type','mh'],dropna=False).size().reset_index(name='counts')
        # Sorting by chrom, pos for easier viewing
        bnds_df = bnds_df.sort_values(by=['chrom','pos']).reset_index()
    
        # Ensure 'pos' is of integer type
        bnds_df['pos'] = bnds_df['pos'].astype(pd.Int64Dtype())
        bnds_df['pos2'] = bnds_df['pos2'].astype(pd.Int64Dtype())
        # Adding distance from cut site
        bnds_df['Positions'] = [pam_site] * len(bnds_df)
        bnds_df['Distance'] = bnds_df.apply(lambda r: min([abs(r['pos'] - x) for x in r['Positions']] + [abs(r['pos'] - x) for x in r['Positions']] ) if r['pos'] is not pd.NA else pd.NA,axis=1)
        bnds_df['source'] = [source] * len(bnds_df)

        # bndinfo needs to add ref, alt, lft_tmplt, rt_tmplt for consistency
        bnds_df['ref'] = '-'
        bnds_df['alt'] = '-'
        bnds_df['lft_tmplt'] = ''
        bnds_df['rt_tmplt'] = ''
        bnds_df['lft_tmplt_pos'] = ''
        bnds_df['rt_tmplt_pos'] = ''

    # Close BAM and fasta files
    bam.close()
    ref_fasta.close()

    # Combine indels and bnds
    if len(bnds_df) == 0:
        combined = indels_df
    else:
        combined = pd.concat([indels_df, bnds_df], ignore_index=True, sort=False)

    log.logit(f"Found {len(combined[combined['type'] == 'REF'])} REF calls")
    log.logit(f"Found {len(combined[combined['type'] == 'INDEL'])} INDELs")
    log.logit(f"Found {len(combined[combined['type'] == 'BND'])} BNDs")
    return combined

def process_breaks(bam_file, control_bam_file, merged_df, fasta, window, distance, minreads, maxcontrol, threads, verbose, outfile):
    results = []
    if threads == 1:
        # Single-threaded processing
        for index, row in merged_df.iterrows():
            guide_name = ';'.join(row['guide_name']) if isinstance(row['guide_name'], list) else row['guide_name']
            result = process_regions(bam_file, fasta, guide_name, row['Chromosome'], row['Start'], row['End'], row['pam_site'], window, distance, minreads, maxcontrol, verbose)
            results.append(result)
    else:
        # Set up multiprocessing pool, if threads > 1
        # Ideally, each process will handle one region from merged_df
        with mp.Pool(processes=threads) as pool:
            # Combine guide names in list comprehension too
            results = pool.starmap(process_regions, [
                (bam_file, fasta, 
                 ';'.join(row['guide_name']) if isinstance(row['guide_name'], list) else row['guide_name'], 
                 row['Chromosome'], row['Start'], row['End'], row['pam_site'], 
                 window, distance, minreads, maxcontrol, verbose) 
                for index, row in merged_df.iterrows()
            ])

    # Combine results from all processes
    breaks = pd.concat(results, ignore_index=True)
    
    breaks['control_alt_counts'] = 0
    breaks['control_total_counts'] = 0
    # Now we need to get control counts for each INDEL/BND, let's split this by source to speed it up
    results = []
    if threads == 1:
        for index, row in merged_df.iterrows():
            # Combine guide names here too for filtering
            guide_name = ';'.join(row['guide_name']) if isinstance(row['guide_name'], list) else row['guide_name']
            result = tools.add_normal_counts(breaks[breaks['source'] == guide_name], control_bam_file, fasta,
                row['Chromosome'], row['Start'], row['End'], 
                window, verbose)
            results.append(result)
    else:
        with mp.Pool(processes=threads) as pool:
            # Combine guide names in multiprocessing too
            results = pool.starmap(tools.add_normal_counts, [
                (breaks[breaks['source'] == (';'.join(row['guide_name']) if isinstance(row['guide_name'], list) else row['guide_name'])], 
                control_bam_file, fasta, 
                row['Chromosome'], row['Start'], row['End'],
                window, verbose) 
                for index, row in merged_df.iterrows()
            ])
    breaks = pd.concat(results, ignore_index=True)

    log.logit(f"Applying filters: minreads={minreads}, maxcontrol={maxcontrol}, distance={distance}")
    # Apply filters step by step (preserve REF calls with ref=='.')
    breaks = breaks[(breaks['counts'] >= minreads) | (breaks['ref']=='.')]
    breaks = breaks[(breaks['control_alt_counts'] <= maxcontrol) | (breaks['ref']=='.')]
    breaks = breaks[(breaks['Distance'] <= distance) | (breaks['ref']=='.')]
    breaks = breaks[~breaks['alt'].str.contains('N')]

    breaks = scarmapper.run_scarmapper(breaks, fasta)
    
    # Add MH Length
    breaks['mh_length'] = breaks['mh'].apply(lambda x: len(x) if x != '' else 0)
    
    # Create Key
    breaks['key'] = breaks.apply(create_variant_key, axis=1)

    log.logit(f"Found {len(breaks[(breaks['type'] == 'REF')])} REF calls after filtering")
    log.logit(f"Found {len(breaks[(breaks['type'] == 'INDEL')])} INDELs after filtering")
    log.logit(f"Found {len(breaks[(breaks['type'] == 'BND')])} BNDs after filtering")

    breaks_out = breaks[['read','source','chrom','pos','ref','alt','chrom2','pos2','strands','type','indel_type','counts','control_total_counts','control_alt_counts','Distance','classification','mh','mh_length','lft_tmplt','rt_tmplt','lft_tmplt_pos','rt_tmplt_pos','key']]

    log.logit(f"Writing final INDELs/BNDs to {outfile}")
    breaks_out.to_csv(outfile, sep='\t', index=False)

    return breaks

def process_vcf(vcf_file, sample_name, has_control, control_sample_name, fasta, chromosome, minreads, maxcontrol, ignore_filters, threads, verbose, outfile):
    # Read in VCF file
    input_vcf = pysam.VariantFile(vcf_file)
    # Check if VCF is indexed
    vcf_path = os.path.abspath(vcf_file)
    if not os.path.exists(vcf_path + ".tbi"):
        log.logit(f"WARNING: VCF file {vcf_file} is not indexed. We will index the VCF file.", color = "yellow")
        pysam.tabix_index(vcf_file, force=True, seq_col=0, start_col=1, end_col=1, preset="vcf", keep_original=True)
    if chromosome is not None:
        variants = input_vcf.fetch(chromosome)
    else:
        variants = input_vcf.fetch()              # This will have all the variants

    num_samples = len(input_vcf.header.samples)

    # There are 3 possibilities:
    # 1. Sample VCF - If it is just sample, we can either filter for PASS (if --ignore_filters is not set) and skip checking if the variant is in the control sample column
    # 2. Normal + Sample VCF - We can filter for PASS (if --ignore_filters is not set) and we WON'T check if the variant is in the control sample column because this is normal, so it is possible for somatic variants to be in the normal sample
    # 3. Control + Sample VCF - We can filter for PASS (if --ignore_filters is not set) and we will check if the variant is in the control sample column, but we will only use the control sample column for control counts if the --has_control flag is set, otherwise we will assume all samples are case samples and there are no control counts

    # There are 2 results, either:
    # 1) The VCF has a control sample column, in which case we need to check that the --has_control flag is set and then we will use the control sample column for control counts
    # 2) The VCF does not have a control sample column, in which case we need to check that the --has_control flag is not set and then we will assume all samples are case samples and there are no control counts
    if has_control and num_samples < 2:
        log.logit(f"ERROR: VCF file {vcf_file} does not have a control sample column, but --has_control flag is set. Please check your VCF file and try again.", color = "red")
        sys.exit(1)

    sample_names = list(input_vcf.header.samples)

    # If --has_control flag is set but --control_sample_name is not provided, we will use the first sample in the VCF file as the control sample
    if has_control and control_sample_name is None:
        log.logit(f"WARNING: --has_control flag is set but --control_sample_name is not provided. We will use the first sample in the VCF file as the control sample.", color = "yellow")
        control_sample_index = 0
    # If --has_control flag is set and --control_sample_name is provided, we need to find the index of the control sample column
    if has_control and control_sample_name is not None:
        if control_sample_name not in input_vcf.header.samples:
            log.logit(f"ERROR: Control sample name {control_sample_name} not found in VCF file {vcf_file}. Please check your VCF file and try again.", color = "red")
            sys.exit(1)
        control_sample_index = sample_names.index(control_sample_name)
    # If --has_control flag is not set, we will assume all samples are case samples and there are no control counts, so we will set control_sample_index to None
    if not has_control:
        control_sample_index = None

    if sample_name not in input_vcf.header.samples:
        log.logit(f"ERROR: Sample name {sample_name} not found in VCF file {vcf_file}. Please check your VCF file and try again.", color = "red")
        sys.exit(1)
    else:
        sample_index = sample_names.index(sample_name)

    variants_df = pd.DataFrame(columns=['variant', 'source', 'chrom', 'pos', 'ref', 'alt', 'type', 'indel_type', 'sample_alt_counts', 'sample_total_counts', 'control_alt_counts', 'control_total_counts', 'classification', 'mh', 'mh_length', 'lft_tmplt', 'rt_tmplt', 'lft_tmplt_pos', 'rt_tmplt_pos'])

    for variant in variants:
        key = f"{variant.chrom}:{variant.pos}:{variant.ref}:{variant.alts[0]}"
        # If --ignore_filters flag is set, ignore filters and process all variants. If not set, only process variants that PASS filters (i.e. FILTER column is 'PASS')
        # Ideally this should always be FALSE because we only want to process variants that PASS filters, but we are adding this flag in case there are some cases where users want to process all variants regardless of filters
        if ignore_filters is False and 'PASS' not in variant.filter.keys():
            continue
        # Note: For multiallelics, either the user processes them beforehand, or we will just take the first alt allele
        for alt in variant.alts:
            if len(variant.ref) == 1 and len(alt) == 1:
                var_type = 'SNP'
                continue
            else:
                var_type = 'INDEL'

            # Get control counts if control_sample_index is not None
            if control_sample_index is not None:
                control_sample = variant.samples[control_sample_index]
                control_alt_counts = control_sample.get('AD')[1] if control_sample.get('AD') is not None else 0
                control_total_counts = sum(control_sample.get('AD')) if control_sample.get('AD') is not None else 0
            else:
                control_alt_counts = 0
                control_total_counts = 0

            if sample_index is not None:
                sample = variant.samples[sample_index]
                sample_alt_counts = sample.get('AD')[1] if sample.get('AD') is not None else 0
                sample_total_counts = sum(sample.get('AD')) if sample.get('AD') is not None else 0
            else:
                sample_alt_counts = 0
                sample_total_counts = 0

            if minreads is not None and sample_alt_counts < minreads:
                continue
            if maxcontrol is not None and control_alt_counts > maxcontrol:
                continue

            variant_data = {
                'variant': f"{variant.chrom}:{variant.pos}:{variant.ref}:{alt}",
                'source': 'VCF',
                'chrom': variant.chrom,
                'pos': variant.pos,
                'ref': variant.ref,
                'alt': alt,
                'type': var_type,
                'indel_type': '',
                'sample_alt_counts': sample_alt_counts,
                'sample_total_counts': sample_total_counts,
                'control_alt_counts': control_alt_counts,
                'control_total_counts': control_total_counts,
                'classification': '',
                'mh': '',
                'mh_length': 0,
                'lft_tmplt': '',
                'rt_tmplt': '',
                'lft_tmplt_pos': '',
                'rt_tmplt_pos': ''
            }

            variants_df = pd.concat([variants_df, pd.DataFrame([variant_data])], ignore_index=True)
    input_vcf.close()
    
    df_with_scars = scarmapper.run_scarmapper(variants_df, fasta)
    df_with_scars['mh_length'] = df_with_scars['mh'].apply(lambda x: len(x) if x != '' else 0)
    df_with_scars['key'] = df_with_scars.apply(lambda row: f"{row['variant']}:{row['classification']}", axis=1)
    df_with_scars.to_csv(outfile, sep='\t', index=False)
    
    tools.annotate_vcf_with_scars(vcf_file, df_with_scars, chromosome, outfile)

    return df_with_scars