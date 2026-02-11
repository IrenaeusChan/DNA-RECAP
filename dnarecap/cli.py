"""
DNA-RECAP: DNA Repair Event Classification and Analysis Pipeline
A tool for analyzing and classifying DNA repair outcomes from CIBAR-seq data.
Can also be used to add repair event classification to VCF files.
"""

import os, sys, signal, time
import click
from clint.textui import puts, colored, indent
import dnarecap.utils.logger as log

from dnarecap.version import __version__

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])

@click.group(context_settings=CONTEXT_SETTINGS)
@click.version_option(version=__version__)
def cli():
    '''Description.'''
    # to make this script/module behave nicely with unix pipes
    # http://newbebweb.blogspot.com/2012/02/python-head-ioerror-errno-32-broken.html
    signal.signal(signal.SIGPIPE, signal.SIG_DFL)

@cli.group('process', short_help="Tools to process CIBAR-seq data")
def process():
    signal.signal(signal.SIGPIPE, signal.SIG_DFL)

@process.command('breaks', short_help="Process INDELs/BNDs from CIBAR-seq data")
@click.option('--bam_file', '-b', type=click.Path(exists=True), required=True, help="Experimental BAM file")
@click.option('--control_bam_file', '-cb', type=click.Path(exists=True), required=True, help="Control BAM file")
@click.option('--input_file', '-i', type=click.Path(exists=True), required=True, help="Input file with off-target coordinates")
# These are the options for the indels command
@click.option('--fasta', '-f', type=click.Path(exists=True), required=True, help="Reference fasta file")
@click.option('--window', '-w', type=int, default=500, show_default=True, help="Distance between off-target sites for merging intervals.")
@click.option('--distance', '-d', type=int, default=25, show_default=True, help="Window size around off-target site to identify mutations")
@click.option('--minreads', '-m', type=int, default=1, show_default=True, help="Minimum supporting reads to report an indel/bnd event.")
@click.option('--chromosome', '-c', type=str, default=None, help="Chromosome to process.")
@click.option('--maxcontrol', '-x', type=int, default=0, show_default=True, help="Maximum supporting reads to to report an indel/bnd event.")
@click.option('--threads', '-t', type=int, default=1, show_default=True, help="Number of threads/cores to use for multiprocessing.")
@click.option('--verbose', '-v', is_flag=True, show_default=True, default=False, help="Print verbose output")
@click.option('--outfile', '-o', type=click.Path(), required=True, default=None, help="Output file")
def process_indels(bam_file, control_bam_file, input_file, fasta, window, distance, minreads, chromosome, maxcontrol, threads, verbose, outfile):
    """
    Process INDELs and BNDs from CIBAR-seq data
    """
    from dnarecap.process import process
    import dnarecap.utils.tools as tools
    start_time = time.time()
    puts(colored.green("Processing INDELs/BNDs from CIBAR-seq data..."))
    merged_df = tools.process_input_file(input_file, window, chromosome, verbose)
    break_results = process.process_breaks(bam_file, control_bam_file, merged_df, fasta, window, distance, minreads, maxcontrol, threads, verbose, outfile)
    outfile_summary = outfile.replace('.tsv', '.summary.tsv')
    open(outfile_summary, 'w').write("chrom\tstart\tend\tpam_positions\ttotal_reads\trepaired_reads\trepaired_fraction\tcontrol_reads\tcontrol_repaired_reads\tcontrol_repaired_fraction\tindel_count\tindel_info\tbnd_count\tbnd_info\ttarget_info\tis_target\t\
                                     num_indel_nhej\tnum_indel_tmej\tnum_indel_nmh_ej\tnum_indel_other\tnum_indel_complex\tindel_nhej_fraction\tindel_tmej_fraction\tindel_nmh_ej_fraction\tindel_other_fraction\tindel_complex_fraction\t\
                                     num_ins_nhej\tnum_ins_tins\tnum_ins_nmh_ej\tnum_ins_other\tindel_ins_nhej_fraction\tindel_ins_tins_fraction\tindel_ins_nmh_ej_fraction\tindel_ins_other_fraction\t\
                                     num_del_nhej\tnum_del_tmej\tnum_del_nmh_ej\tnum_del_other\tindel_del_nhej_fraction\tindel_del_tmej_fraction\tindel_del_nmh_ej_fraction\tindel_del_other_fraction\t\
                                     num_bnd_tmej\tnum_bnd_nmh_ej\tnum_bnd_other\tbnd_tmej_fraction\tbnd_nmh_ej_fraction\tbnd_other_fraction\n")
    for index, row in merged_df.iterrows():
        source = ';'.join(row['guide_name']) if isinstance(row['guide_name'], list) else row['guide_name']
        tools.summarise_repaired_reads(break_results[(break_results['source'] == source)], outfile.replace('.tsv', '.summary.tsv'), row)
    puts(colored.green("Finished processing INDELs/BNDs."))
    total_elapsed = tools.process_time(time.time() - start_time)
    puts(colored.green(f"Total time elapsed: {total_elapsed}."))

@process.command('vcf', short_help="Adds Repair Event Classification to VCF files")
@click.option('--vcf_file', '-i', type=click.Path(exists=True), required=True, help="Input VCF file")
@click.option('--sample_name', '-s', type=str, required=False, help="Sample name to add to the output VCF file. If not provided, will use the sample name from the input VCF file.")
@click.option('--has_control', '-n', is_flag=True, show_default=True, default=False, help="Whether the VCF file has a control sample. If not provided, will assume the VCF file does not have a control sample.")
@click.option('--control_sample_name', '-cs', type=str, required=False, help="Control sample name to use for control counts if --has_control flag is set. If not provided, but --has_control flag is set, will use the first sample in the VCF file as the control sample.")
@click.option('--fasta', '-f', type=click.Path(exists=True), required=True, help="Reference fasta file")
@click.option('--chromosome', '-c', type=str, default=None, help="Chromosome to process.")
@click.option('--minreads', '-m', type=int, default=1, show_default=True, help="Minimum supporting reads to report an indel/bnd event.")
@click.option('--maxcontrol', '-x', type=int, default=0, show_default=True, help="Maximum supporting reads to to report an indel/bnd event.")
@click.option('--ignore_filters', '-p', type=bool, is_flag=True, show_default=True, default=False, help="By default, only variants that pass all filters (i.e. FILTER column is 'PASS') will be processed. Set this flag to ignore filters and process all variants.")
@click.option('--threads', '-t', type=int, default=1, show_default=True, help="Number of threads/cores to use for multiprocessing.")
@click.option('--verbose', '-v', is_flag=True, show_default=True, default=False, help="Print verbose output")
@click.option('--outfile', '-o', type=click.Path(), required=True, default=None, help="Output TSV/VCF file with added Repair Event Classification")
def process_vcf(vcf_file, sample_name, has_control, control_sample_name, fasta, chromosome, minreads, maxcontrol, ignore_filters, threads, verbose, outfile):
    """
    Process VCF files to add Repair Event Classification
    """
    from dnarecap.process import process
    import dnarecap.utils.tools as tools
    start_time = time.time()
    puts(colored.green("Processing VCF file to add Repair Event Classification..."))
    vcf_df = process.process_vcf(vcf_file, sample_name, has_control, control_sample_name, fasta, chromosome, minreads, maxcontrol, ignore_filters, threads, verbose, outfile)
    
    puts(colored.green("Finished processing VCF file."))
    total_elapsed = tools.process_time(time.time() - start_time)
    puts(colored.green(f"Total time elapsed: {total_elapsed}."))