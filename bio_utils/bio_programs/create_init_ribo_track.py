#! /usr/bin/env python3

""" Create a bedGraph (bigWig) file from a riboseq BAM file
to visualise initiating ribosomes (p-site offset shifted).
"""

import os
import pandas as pd

import argparse
import logging
import misc.logging_utils as logging_utils
logger = logging.getLogger(__name__)

import yaml
import json

import misc.utils as utils
import misc.shell_utils as shell_utils

import riboutils.ribo_utils as ribo_utils
import riboutils.ribo_filenames as filenames

display_mode = 'full'
colour_positive = '128,128,0'
colour_negative = '218,165,32'

def _convert(bed, bw, args):

    in_files = [bed, args.chrSizes]
    out_files = [bw]

    try:
        msg = "Calling bedGraphToBigWig"
        logger.info(msg)

        cmd = "bedGraphToBigWig {} {} {}".format(bed, args.chrSizes, bw)
        shell_utils.call_if_not_exists(cmd, out_files, in_files=in_files,
                                       overwrite=args.overwrite, call=True)
    finally:
        if not args.keep_bed:
            try:
                os.remove(bed)
                msg = "Removing: {}".format(bed)
                logger.info(msg)
            except OSError:
                msg = "Could not remove: {}. Skipping.".format(bed)
                logger.info(msg)


def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="""Create 2 bedGraph files e.g. to use on a genome browser, to visualise
        the riboseq read counts at each genomic position. Periodic reads are mapped to 
        their 5' ends and shifted according to their p-site offset, one file for each strand. 
        This script does NOT determine the periodic lengths and offsets, thus 
        the 'periodic-offsets' file must be available. Unless specified, a file is created
        for every sample in the configuration file.""")

    parser.add_argument('config', help="The (yaml) config file.")
    parser.add_argument('outloc', help="The output directory. Temporary files are also"
                                    "written to this location if converting to BigWig.")

    parser.add_argument('--use-pretty-names', help="If this flag is given, then will use the names"
        "in 'riboseq_sample_name_map' if present.", action='store_true')
    parser.add_argument('--isoform-strategy', help="""See b-tea.cl_utils.
        Currently we only use 'merged'. By default, we use the 'unique' mapped reads,
        and no other option is implemented (e.g. stranded)""", default=None)
    parser.add_argument('--input-list', help="""A space-delimited list of sample names, 
        each quoted separately. They must match the keys in the configuration file. 
        Only these will be converted.""", nargs='*', type=str)
    # TO DO: add option to sum profiles of all samples of one condition, as in Rp-Bp.
    parser.add_argument('--add-chr', help="If this flag is present then 'chr' will be pre-pended"
        "to sequence names. This is done before any other changes to sequence names, so this"
        "must be taken into account if using [--chr-dict]", action='store_true')
    parser.add_argument('-d', '--chr-dict', help="""A dictionary mapping of 
        sequence names found in the data to the sequence names, as in 'chrom.sizes'. 
        The format is as follows: '{"key":"value"}'.""", type=json.loads)
    parser.add_argument('--overwrite', help="If this flag is present, then existing files"
        "will be overwritten.", action='store_true')

    subparser = parser.add_subparsers(help="""Optional arguments if converting to bigWig. 
        The program 'bedGraphToBigWig' must be available on the user's path, and can be"
        "downloaded from 'http://hgdownload.soe.ucsc.edu/admin/exe/'.""", dest='to_bigWig')
    parser_to_bigWig = subparser.add_parser('to-bigWig', help="""If given, will convert
        to bigWig.""")
    parser_to_bigWig.add_argument('chrSizes', help="""The 'chrom.sizes' file 
        for the UCSC database.""", type=str)
    parser_to_bigWig.add_argument('-k', '--keep-bed', help="If this flag is given, then "
        "the bedGraph file will not be deleted", action='store_true')

    logging_utils.add_logging_options(parser)
    args = parser.parse_args()
    logging_utils.update_logging(args)

    config = yaml.load(open(args.config))
    note = config.get('note', None)
    is_unique = not ('keep_riboseq_multimappers' in config)

    required_keys = [
        'riboseq_data',
        'riboseq_samples',
    ]
    utils.check_keys_exist(config, required_keys)

    if args.input_list:
        sample_names = {name: [name] for name in args.input_list}
    else:
        sample_names = config['riboseq_samples']

    if args.use_pretty_names:
        sample_name_map = ribo_utils.get_sample_name_map(config)
    else:
        sample_name_map = {name: [name] for name in config['riboseq_samples'].keys()}

    # check output path
    if os.path.exists(args.outloc):
        args.outloc = os.path.join(args.outloc, '')
    else:
        msg = "Invalid output path or wrong permission: {}. Quitting.".format(args.outloc)
        raise OSError(msg)

    for name in sorted(sample_names.keys()):

        msg = "Processing sample: {}".format(name)
        logger.info(msg)

        # get riboseq bam file
        bam = filenames.get_riboseq_bam(
            config['riboseq_data'],
            name,
            is_unique=is_unique,
            isoform_strategy=args.isoform_strategy,
            note=note
        )

        # get the lengths and offsets
        lengths, offsets = ribo_utils.get_periodic_lengths_and_offsets(
            config,
            name,
            isoform_strategy=args.isoform_strategy,
            is_unique=is_unique
        )

        if len(lengths) == 0:
            msg = "No periodic read lengths and offsets were found!"
            logger.critical(msg)
            return

        # now we don't modify the ribo scripts, since we assume that riboseq is always fr
        msg = "Finding P-sites"
        logger.info(msg)

        lengths = [int(l) for l in lengths]
        offsets = [int(l) for l in offsets]
        p_sites = ribo_utils.get_p_sites(bam, lengths, offsets)

        if args.add_chr:
            p_sites['seqname'] = 'chr' + p_sites['seqname'].astype(str)
        if args.chr_dict:
            for chrom_old, chrom_new in args.chr_dict.items():
                seqname_m = p_sites['seqname'] == str(chrom_old)
                p_sites.loc[seqname_m, 'seqname'] = str(chrom_new)

        # create the profiles separately for each strand
        p_sites = p_sites.groupby(p_sites.columns.tolist()).size().reset_index().rename(columns={0: 'count'})

        m_positive = p_sites['strand'] == '+'
        p_sites_positive = p_sites[m_positive]

        m_negative = p_sites['strand'] == '-'
        p_sites_negative = p_sites[m_negative]

        p_sites_positive = p_sites_positive[['seqname', 'start', 'end', 'count']]
        p_sites_negative = p_sites_negative[['seqname', 'start', 'end', 'count']]

        # write bedGraph to disk
        # currently header is not "customisable", add options later if necessary

        note_str = filenames.get_note_string(note)
        unique_str = filenames.get_unique_string(is_unique)
        iso_str = filenames.get_isoform_strategy_string(args.isoform_strategy)

        header = ('track type=bedGraph '
                  'name="{}" '
                  'description="{}" '
                  'visibility={} '
                  'color={} '
                  'autoScale=off '
                  'alwaysZero=on '
                  'graphType=bar').format(sample_name_map[name],
                                          sample_name_map[name],
                                          display_mode,
                                          colour_positive)

        fn = ''.join([name, note_str, iso_str, unique_str, '.sense.bedGraph'])
        bg_sense_output = os.path.join(str(args.outloc), fn)
        file_exists = (os.path.exists(bg_sense_output) and os.path.getsize(bg_sense_output) > 0)
        if args.overwrite or not file_exists:
            with open(bg_sense_output, 'w') as f:
                f.write("{}\n".format(header))
            p_sites_positive.to_csv(bg_sense_output, index=False, header=False, sep='\t', mode='a')
        else:
            msg = "The bedGraph file {} already exists. Skipping".format(bg_sense_output)
            logger.warning(msg)

        header = ('track type=bedGraph '
                  'name="{}" '
                  'description="{}" '
                  'visibility={} '
                  'color={} '
                  'autoScale=off '
                  'alwaysZero=on '
                  'graphType=bar').format(sample_name_map[name],
                                          sample_name_map[name],
                                          display_mode,
                                          colour_negative)

        fn = ''.join([name, note_str, iso_str, unique_str, '.anti-sense.bedGraph'])
        bg_antisense_output = os.path.join(str(args.outloc), fn)
        file_exists = (os.path.exists(bg_antisense_output) and os.path.getsize(bg_antisense_output) > 0)
        if args.overwrite or not file_exists:
            with open(bg_antisense_output, 'w') as f:
                f.write("{}\n".format(header))
            p_sites_negative.to_csv(bg_antisense_output, index=False, header=False, sep='\t', mode='a')
        else:
            msg = "The bedGraph file {} already exists. Skipping".format(bg_antisense_output)
            logger.warning(msg)

        if args.to_bigWig:

            programs = ['bedGraphToBigWig']
            shell_utils.check_programs_exist(programs)

            fn = ''.join([name, note_str, iso_str, unique_str, '.sense.bw'])
            output = os.path.join(str(args.outloc), fn)
            _convert(bg_sense_output, output, args)

            fn = ''.join([name, note_str, iso_str, unique_str, '.anti-sense.bw'])
            output = os.path.join(str(args.outloc), fn)
            _convert(bg_antisense_output, output, args)



if __name__ == '__main__':

    main()
