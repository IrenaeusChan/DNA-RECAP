import pysam
import dnarecap.utils.logger as log
import dnarecap.utils.tools as tools
from dnarecap.utils.variant import Variant

def indel_scars(variant, fasta):
    """
    Classify microhomology and templated insertions for a given INDEL variant using ScarMapper logic.
    1. For deletions, identify microhomology at the junction.
    2. For insertions longer than 4 bp, identify templated insertions from flanking sequences.
    3. Return microhomology sequence and templated insertion sequences (if any).
    4. Adapted from Dennis Simpson ScarMapper: https://pubmed.ncbi.nlm.nih.gov/33963863/ found here: https://github.com/Gaorav-Gupta-Lab/ScarMapper/
    """
    refseq = pysam.FastaFile(fasta)
    try:
        chrom_length = refseq.get_reference_length(variant.chrom)
    except KeyError:
        log.logit(f"Warning: Chromosome {variant.chrom} not found in reference")
        refseq.close()
        return "", "", "", None, None
    # Classify INDEL Type and get length
    if len(variant.ref) > 1:
        indel_type = 'Del'
        indel_length = len(variant.ref) - 1
    else:
        indel_type = 'Ins'
        indel_length = len(variant.alt)

    # Bounds checking for target sequence fetching
    left_start = max(0, variant.pos - indel_length)
    left_end = variant.pos
    right_start = variant.pos + indel_length
    right_end = min(chrom_length, variant.pos + indel_length * 2)
    
    # Get the target sequence with bounds checking
    if left_start >= left_end or right_start >= right_end:
        log.logit(f"Warning: Invalid coordinates for variant {variant.chrom}:{variant.pos}")
        refseq.close()
        return "", "", "", None, None

    # Get the target sequence
    right_target = refseq.fetch(variant.chrom, variant.pos+indel_length, variant.pos+indel_length*2).upper()
    left_target = refseq.fetch(variant.chrom, variant.pos-indel_length, variant.pos).upper()

    microhomology = ""
    mh_count = 0
    lft_template = None
    rt_template = None
    lft_pos_template = None
    rt_pos_template = None

    # Check for microhomology in the deletion.
    if indel_type == "Del":
        # Right target is the sequence after the deletion
        # Left target is the sequence before the deletion
        # We need to check for microhomology on both sides of the deletion
        # e.g. chr1:26934776:TCAGCCTC:T       The deletion is CAGCCTC = variant.ref[1:]
        # Right target is e.g. CAGTAGCTGGGAC
        # For each nucleotide in the right target, check if it matches the reference sequence
        for fnucleotide, rnucleotide in zip(right_target, variant.ref[1:]):
            if fnucleotide == rnucleotide:
                microhomology += fnucleotide
                mh_count += 1
            else:
                break
        if not microhomology:
            for fnuc, rnuc in zip(reversed(left_target), reversed(variant.ref[1:])):
                if fnuc == rnuc:
                    microhomology += fnuc
                    mh_count += 1
                else:
                    break
    # Adapted from https://github.com/Gaorav-Gupta-Lab/ScarMapper/
    elif indel_type == "Ins" and indel_length > 4:
        insertion = variant.alt[1:]                        # This is the insertion
        lft_query1 = tools.revcomp(insertion[:5])          # Reverse complement the first 5 bp
        lft_query2 = insertion[-5:]                        # Grab the last 5 bp
        rt_query1 = insertion[:5]                          # Grab the first 5 bp
        rt_query2 = tools.revcomp(insertion[-5:])          # Reverse complement of the last 5 bp
        
        search_space = 50 #(15+len(insertion))    # <-- This search space is currently subjective
        lower_limit = max(0, variant.pos - search_space)  # Prevent negative coordinates
        upper_limit = min(chrom_length, variant.pos + search_space)  # Prevent exceeding chromosome
        
        left_not_found = True
        right_not_found = True
        possible_left_not_found = True
        possible_right_not_found = True
        lft_template = ""
        rt_template = ""

        # Set starting positions and search for left template
        lft_position = max(0, variant.pos-len(insertion))
        rt_position = variant.pos

        # If it's the left side... it should ONLY be the reverse complement

        # E.g. chr4	49098604	T	('TTTCCA',)
        # TTCCATCCCATTCCATTCCATTCCATTCCATTCCGTTCCGTTCCGTTCC T TTCCATTCCATTCCATTCCATTGCATTCCATTCCATTCCATTCCATTCTA
        # target_segment = GTTCC
        # lft_query1 = TGGAA           lft_query2 = TTCCA
        #  rt_query1 = TTCCA            rt_query2 = TGGAA

        # We let lft_query2 and rt_query2 because the mutation could occur on either strand? 
        # Double check with John about this...
        # it could only happen from 5'-->3'... by allowing query2, we assume that the mutation could happen on either strand.

        while left_not_found and rt_position > lower_limit and lft_position >= 0 and lft_position < rt_position:
            target_segment = refseq.fetch(variant.chrom, lft_position, rt_position).upper()
            if lft_query1 == target_segment:
                lft_template = target_segment
                left_not_found = False
                pos_template = lft_position

            lft_position -= 1
            rt_position -= 1

        # Reset starting positions and search for right template
        lft_position = variant.pos
        rt_position = min(chrom_length, variant.pos + len(insertion))

        while right_not_found and lft_position < upper_limit and rt_position <= chrom_length and lft_position < rt_position:
            target_segment = refseq.fetch(variant.chrom, lft_position, rt_position).upper()
            
            if rt_query1 == target_segment:
                rt_template = target_segment
                right_not_found = False
                pos_template = rt_position

            lft_position += 1
            rt_position += 1

        # Since we don't know if 3'-->5' can happen
        # Keep track just in case we need to check for it.
        lft_position = max(0, variant.pos-len(insertion))
        rt_position = variant.pos

        while possible_left_not_found and rt_position > lower_limit and lft_position >= 0 and lft_position < rt_position:
            target_segment = refseq.fetch(variant.chrom, lft_position, rt_position).upper()
            if lft_query2 == target_segment:
                lft_pos_template = tools.revcomp(target_segment)
                possible_left_not_found = False
                pos_template = lft_position

            lft_position -= 1
            rt_position -= 1

        lft_position = variant.pos
        rt_position = min(chrom_length, variant.pos + len(insertion))

        while possible_right_not_found and lft_position < upper_limit and rt_position <= chrom_length and lft_position < rt_position:
            target_segment = refseq.fetch(variant.chrom, lft_position, rt_position).upper()
            if rt_query2 == target_segment:
                rt_pos_template = tools.revcomp(target_segment)
                possible_right_not_found = False
                pos_template = rt_position

            lft_position += 1
            rt_position += 1

    refseq.close()
    return microhomology, lft_template, rt_template, lft_pos_template, rt_pos_template

def run_scarmapper(df, fasta):
    """
    Run ScarMapper to identify microhomology and templated insertions for each INDEL.
    """
    for it, row in df.iterrows():
        if row['type'] == 'REF':
            df.at[it, 'indel_type'] = 'Ref'
            continue
        
        if row['type'] == 'INDEL':
            indel_string = f"{row['chrom']}:{row['pos']}:{row['ref']}:{row['alt']}"
            microhomology, lft_templated_ins, rt_templated_ins, lft_pos_templated_ins, rt_pos_templated_ins = indel_scars(Variant.createVariant(indel_string), fasta)

            df.at[it,'mh'] = microhomology
            df.at[it,'lft_tmplt'] = lft_templated_ins
            df.at[it,'rt_tmplt'] = rt_templated_ins
            df.at[it,'lft_tmplt_pos'] = lft_pos_templated_ins
            df.at[it,'rt_tmplt_pos'] = rt_pos_templated_ins

            # Classify the INDEL
            if len(row['ref']) > 1 and len(row['alt']) > 1:
                indel_type = "Complex"
                indel_size = 100
            elif len(row['ref']) > 1 and len(row['alt']) == 1:
                indel_type = "Del"
                indel_size = len(row['ref']) - 1
            elif len(row['ref']) == 1 and len(row['alt']) > 1:
                indel_type = "Ins"
                indel_size = len(row['alt']) - 1

            # Classify the INDEL
            df.at[it, 'indel_type'] = indel_type
            df.at[it, 'classification'] = classify_indel(indel_type,indel_size,microhomology,lft_templated_ins,rt_templated_ins)
        
        elif row['type'] == 'BND':
            indel_string = f"{row['chrom']}:{row['pos']}:{row['chrom2']}:{row['pos2']}:{row['strands']}"
            indel_type = "Fusion"
            microhomology = row['mh'] if 'mh' in row else ''

            # Classify the INDEL
            df.at[it, 'indel_type'] = indel_type
            df.at[it, 'classification'] = classify_indel(indel_type,0,microhomology,"","")

        elif row['type'] == 'SNV':
            df.at[it, 'mh'] = ''
            df.at[it, 'lft_tmplt'] = ''
            df.at[it, 'rt_tmplt'] = ''
            df.at[it, 'lft_tmplt_pos'] = None
            df.at[it, 'rt_tmplt_pos'] = None
            df.at[it, 'indel_type'] = ''
            df.at[it, 'classification'] = ''

    return df.copy()

def classify_indel(indel_type, indel_size, microhomology, lft_template, rt_template):
    if indel_type == "Del":
        if indel_size >= 1 and indel_size < 4:
            return "NHEJ"
        elif indel_size >= 4 and len(microhomology) >= 2:
            return "TMEJ"
        elif indel_size >= 4 and len(microhomology) <= 1:
            return "Non-MH EJ"
        else:
            return "Other"
    elif indel_type == "Ins":
        if indel_size >= 1 and indel_size < 4:
            return "NHEJ"
        # NOTE: I am switching this to >= 4 because DEL is at >= 4, and I want to keep consistent
        elif indel_size >= 4 and (lft_template or rt_template):
            return "TINS"
        # NOTE: This will also be switched to >= 4 to be consistent with DEL
        elif indel_size >= 4:
            return "Non-MH EJ"
        else:
            return "Other"
    elif indel_type == "Fusion":
        if len(microhomology) >= 2:
            return "TMEJ"
        elif len(microhomology) <= 1:
            return "Non-MH EJ"
        else:
            return "Other"
    elif indel_type == "Complex":
        return "Complex"
