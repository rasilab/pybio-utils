#! /usr/bin/env python3

import argparse
import logging

import numpy as np
import pandas as pd

import pbio.utils.bed_utils as bed_utils
import pbio.utils.gtf_utils as gtf_utils
import pbio.misc.logging_utils as logging_utils
import pbio.misc.parallel as parallel
import pbio.misc.utils as utils

logger = logging.getLogger(__name__)

default_chr_name_file = None
default_exon_feature = "exon"
default_cds_feature = "CDS"

default_num_cpus = 3
default_num_groups = 500

def get_transcript_id(gtf_entry, attr_name, offset, is_quote):
    attributes_str = gtf_entry['attributes']
    start = attributes_str.find(attr_name) + offset

    # make sure we actually found the attribute
    if start == offset - 1:
        return None
    
    end = attributes_str.find(";", start)
    # some GFF3 attributes may end without a semi-colon!
    # we assume that otherwise all attributes will end with a semi-colon
    if end == -1:
        end = len(attributes_str)
    transcript_id = attributes_str[start+is_quote:end-is_quote]
    return transcript_id
    
def get_transcript_ids(gtf_entries, attr_name, offset, is_quote):
    ret = parallel.apply_df_simple(gtf_entries, get_transcript_id, attr_name, offset, is_quote)
    return ret

def get_bed12_entry(gtf_entries):
    starts = np.array(gtf_entries['start'])
    start = min(starts)
    
    rel_starts = starts - start
    rel_starts_str = ",".join(str(s) for s in rel_starts)
    
    # we now subtract 1 from the start because BED is base-0
    start -= 1
    
    # we do not subtract 1 from the end, though, because BED
    # is open on the "end"
    
    lengths_str = ",".join(gtf_entries['length'])
    
    ret = {    
        'seqname': gtf_entries['seqname'].iloc[0],
        'start': start,
        'end': max(gtf_entries['end']),
        'id': gtf_entries['transcript_id'].iloc[0],
        'score': 0,
        'strand': gtf_entries['strand'].iloc[0],
        'thick_start': -1,
        'thick_end': -1,
        'color': 0,
        'num_exons': len(gtf_entries),
        'exon_lengths': lengths_str,
        'exon_genomic_relative_starts': rel_starts_str
    }
    
    return ret

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="This script converts a GTF file to a BED12 file. In particular, "
        "it creates bed entries based on the exon features and transcript_id field. "
        "It uses the CDS regions to determine the \"thick_start\" and \"thick_end\" "
        "features of the BED12 file.")

    parser.add_argument('gtf', help="The GTF file")
    parser.add_argument('out', help="The (output) BED12 file")

    parser.add_argument('--chr-name-file', help="If this file is given, then the "
        "bed entries will be sorted according to the order of seqnames in this "
        "file. Presumably, this is the chrName.txt file from STAR.", 
        default=default_chr_name_file)

    parser.add_argument('--exon-feature', help="The name of features which are "
        "treated as exons", default=default_exon_feature)
    parser.add_argument('--cds-feature', help="The name of features which are "
        "treated as CDSs", default=default_cds_feature)

    parser.add_argument('-p', '--num-cpus', help="The number of CPUs to use",
        type=int, default=default_num_cpus)
    parser.add_argument('-g', '--num-groups', help="The number of groups to split "
        "into for parallelization", type=int, default=default_num_groups)

    logging_utils.add_logging_options(parser)
    args = parser.parse_args()
    logging_utils.update_logging(args)

    msg = "Reading annotations file"
    logger.info(msg)

    attr_name = "transcript_id "
    offset = len(attr_name)
    is_quote = 1
    use_gff3_specs = args.gtf.endswith('gff')
    if use_gff3_specs:
        msg = """The program will use GFF3 specifications to extract features. In addition, 
        the stop codon will be removed from the CDSs."""
        logger.info(msg)
        # CDS and exons have a "Parent" attribute whose value matches the "ID" attribute of the
        # corresponding transcript feature, equivalent to "transcript_id" for our purpose.
        attr_name = "Parent="
        offset = len(attr_name)
        is_quote = 0

    # gtf or gff, we use the same fields
    gtf = gtf_utils.read_gtf(args.gtf)

    msg = "Extracting exon and CDS features"
    logger.info(msg)

    m_exons = gtf['feature'] == args.exon_feature
    m_cds = gtf['feature'] == args.cds_feature
        
    exons = gtf[m_exons].copy()
    cds_df = gtf[m_cds].copy()


    msg = "Extracting CDS transcript ids"
    logger.info(msg)

    cds_transcript_ids = parallel.apply_parallel_split(cds_df, args.num_cpus, 
            get_transcript_ids, attr_name, offset, is_quote, progress_bar=True, num_groups=args.num_groups)
    cds_transcript_ids = utils.flatten_lists(cds_transcript_ids)
    cds_df['transcript_id'] = cds_transcript_ids


    msg = "Calculating CDS genomic start and end positions"
    logger.info(msg)

    cds_groups = cds_df.groupby('transcript_id')

    # we subtract 1 from start because gtf is 1-based
    cds_min_starts = cds_groups['start'].min()
    cds_start_df = pd.DataFrame()
    cds_start_df['id'] = cds_min_starts.index
    cds_start_df['cds_start'] = cds_min_starts.values - 1

    # we do not subtract 1 from end because bed is "open" on the end
    cds_max_end = cds_groups['end'].max()
    cds_end_df = pd.DataFrame()
    cds_end_df['id'] = cds_max_end.index
    cds_end_df['cds_end'] = cds_max_end.values

    # Under GFF3 specifications, START and STOP codons are included in the CDS,
    # so we need to exclude it from the coding region.
    if use_gff3_specs:
        strand = cds_groups['strand'].unique()
        strand_df = pd.DataFrame()
        strand_df['id'] = strand.index
        strand_df['strand'] = strand.values
        cds_start_df = cds_start_df.merge(strand_df, on='id', how='left')
        cds_end_df = cds_end_df.merge(strand_df, on='id', how='left')

        m_positive = cds_end_df['strand'].apply(lambda x: bool(x == '+'))
        cds_end_df.loc[m_positive, 'cds_end'] = cds_end_df.loc[m_positive, 'cds_end'] - 3
        cds_end_df.drop('strand', axis=1, inplace=True)
        m_negative = cds_start_df['strand'].apply(lambda x: bool(x == '-'))
        cds_start_df.loc[m_negative, 'cds_start'] = cds_start_df.loc[m_negative, 'cds_start'] + 3
        cds_start_df.drop('strand', axis=1, inplace=True)

    msg = "Extracting exon transcript ids"
    logger.info(msg)

    exon_transcript_ids = parallel.apply_parallel_split(exons, args.num_cpus, 
            get_transcript_ids, attr_name, offset, is_quote, progress_bar=True, num_groups=args.num_groups)
    exon_transcript_ids = utils.flatten_lists(exon_transcript_ids)
    exons['transcript_id'] = exon_transcript_ids

    exons['length'] = exons['end'] - exons['start'] + 1
    exons['length'] = exons['length'].astype(str)

    # store these for sorting later
    transcript_ids = np.array(exons['transcript_id'])

    msg = "Combining exons into BED12 entries"
    logger.info(msg)

    exons = exons.sort_values('start')
    exon_groups = exons.groupby('transcript_id')

    bed12_df = parallel.apply_parallel_groups(exon_groups, args.num_cpus, 
            get_bed12_entry, progress_bar=True)
    bed12_df = pd.DataFrame(bed12_df)


    msg = "Joining BED12 entries to CDS information"
    logger.info(msg)

    bed12_df = bed12_df.merge(cds_start_df, on='id', how='left')
    bed12_df = bed12_df.merge(cds_end_df, on='id', how='left')

    bed12_df = bed12_df.fillna(-1)

    bed12_df['thick_start'] = bed12_df['cds_start'].astype(int)
    bed12_df['thick_end'] = bed12_df['cds_end'].astype(int)

    msg = "Sorting BED12 entries"
    logger.info(msg)

    # We will break ties among transcripts by the order they appear 
    # in the GTF file. This is the same way star breaks ties.
    bed12_df = bed_utils.sort(bed12_df, seqname_order=args.chr_name_file, 
                    transcript_ids=transcript_ids)

    msg = "Writing BED12 to disk"
    logger.info(msg)

    bed_utils.write_bed(bed12_df[bed_utils.bed12_field_names], args.out)

if __name__ == '__main__':
    main()
