#! /usr/bin/env python3

import sys
import os

import argparse
import logging
import misc.logging_utils as logging_utils
logger = logging.getLogger(__name__)

import yaml
import json
import csv

import bio_utils.bed_utils as bed_utils
import misc.utils as utils
import misc.shell_utils as shell_utils
import misc.pandas_utils as pandas_utils

import riboutils.ribo_utils as ribo_utils
import riboutils.ribo_filenames as filenames


def _get_bed(input_filename, fields_to_keep, args):
    ''' Get bed12 file, sort and adjust features.
        This function contains hard coded fields if
        using color, etc.
    '''

    bed_df = bed_utils.get_bed_df(input_filename)

    # Adjust chrom field.
    if args.add_chr:
        bed_df['seqname'] = 'chr' + bed_df['seqname'].astype(str)
    if args.chr_dict:
        for chrom_old, chrom_new in args.chr_dict.items():
            seqname_m = bed_df['seqname'] == str(chrom_old)
            bed_df.loc[seqname_m, 'seqname'] = str(chrom_new)

    # Sort on the chrom field, and then on the chromStart field.
    bed_df.sort_values(['seqname', 'start'], ascending=[True, True], inplace=True)

    if args.use_color:
        # Adjust scores where canonical (and novel if present) are darker,
        # and where anything 'in between' is paler, and suspect, etc. are very pale.
        canonical_m = bed_df['orf_type'] == 'canonical'
        novel_m = bed_df['orf_type'] == 'novel'

        within_m = bed_df['orf_type'].str.contains('within')
        noncoding_m = bed_df['orf_type'].str.contains('noncoding')
        suspect_m = bed_df['orf_type'].str.contains('suspect')
        pale_m = within_m | noncoding_m | suspect_m

        # extended, truncated, overlap and all primes
        all_others_m = ~novel_m & ~canonical_m & ~pale_m

        bed_df.loc[canonical_m, 'score'] = 1000
        bed_df.loc[novel_m, 'score'] = 800
        bed_df.loc[all_others_m, 'score'] = 600
        bed_df.loc[pale_m, 'score'] = 400

        # Add color
        # all canonical but not novel are blue, the novels are all red,
        # the 5' not novel green, the others black
        canonical_m = bed_df['orf_type'].str.contains('canonical')
        novel_m = bed_df['orf_type'].str.contains('novel')
        five_m = bed_df['orf_type'].str.contains('five')
        three_m = bed_df['orf_type'].str.contains('three')

        canonical_not_novel_m = canonical_m & ~novel_m
        five_not_novel = five_m & ~novel_m
        three_not_novel = three_m & ~novel_m
        novel_m = novel_m & ~pale_m

        bed_df['color'] = '0,0,0'
        bed_df.loc[canonical_not_novel_m, 'color'] = '0,0,255'
        bed_df.loc[novel_m, 'color'] = '255,0,0'
        bed_df.loc[five_not_novel, 'color'] = '0,255,0'
        bed_df.loc[three_not_novel, 'color'] = '255,255,0'

    # remove unused fields
    bed_df = bed_df[fields_to_keep]

    # Writes bed file to output directory.
    _, name = os.path.split(input_filename)
    filen, ext = os.path.splitext(name)
    if ext == '.gz':
        filen, ext = os.path.splitext(filen)
    output_filename = args.out + filen
    pandas_utils.write_df(bed_df, str(output_filename + '.tmp.bed'), index=False, sep='\t',
    header=False, do_not_compress=True, quoting=csv.QUOTE_NONE)

    return output_filename

def _convert(bed, bb, use_config_fields, args):
    in_files = [bed, args.chrSizes]
    out_files = [bb]
    if use_config_fields:
        cmd = "bedToBigBed -as={} -type={} -extraIndex={} {} {} {}".format(use_config_fields['as_file'],
            use_config_fields['bed_type'], "name", bed, args.chrSizes, bb)
        in_files.append(use_config_fields['as_file'])
    else:
        cmd = "bedToBigBed {} {} {}".format(bed, args.chrSizes, bb)
    shell_utils.call_if_not_exists(cmd, out_files, in_files=in_files,
        overwrite=args.overwrite, call=True)
    try:
        os.remove(bed)
        msg = "Removing: {}".format(bed)
        logger.info(msg)
    except OSError:
        msg = "Could not remove: {}".format(bed)
        logger.info(msg)


def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="This script converts a set of bed files to bigBed files, e.g. to use"
        "on a genome browser, first by removing unused fields, and then by calling the executable"
        "program 'bedToBigBed'. This program must be available on the user's path, and can be"
        "downloaded from 'http://hgdownload.soe.ucsc.edu/admin/exe/'. Unless specified,"
        "all samples and conditions in the configuration file are used. Input files must either be"
        ".bed or .bed.gz. ** Some features are still experimental and not for general use**.")

    parser.add_argument('config', help="The (yaml) config file.")
    parser.add_argument('out', help="The output directory. All bed files are temporarily re-written "
        "to this location.")
    parser.add_argument('chrSizes', help="The 'chrom.sizes' file for the UCSC database.")

    parser.add_argument('--add-chr', help="If this flag is present then 'chr' will be pre-pended"
        "to sequence names. This is done before any other changes to sequence names, so this"
        "must be taken into account if giving a dictionary mapping", action='store_true')

    parser.add_argument('-d', '--chr-dict', help="""A dictionary mapping of sequence names found
        in the data to the sequence names, as in 'chrom.sizes'. The format is as follows:
        '{"key":"value"}' """, type=json.loads)

    parser.add_argument('--configure-fields', help="A file with comma-separated items (one per line)"
        "corresponding to fields that will be included in the bigBed file. The field names must"
        "correspond to the ones used the bed file. Each field name must be followed by"
        "'type','standard field name', 'description', as needed to generate the AutoSql format (.as)"
        "file describing these fields. Standard fields must be separated from any extra fields"
        "by an empty line. See e.g.3 here: https://genome.ucsc.edu/goldenpath/help/bigBed.html."
        "One extra index will be created on the name field by default. If multiple bed files are"
        "passed in, these will be used for all input files.", required='--use-color' in sys.argv)

    parser.add_argument('--use-color', help="If this flag is present then custom score and"
        "color (field 9) fields are added. These are currently not configurable, and presumably used"
        "only with the ORF predictions from the Rp-Bp pipeline. The [--custom-field] option is"
        "required if using colour, even if no extra fields are given.", action='store_true')

    parser.add_argument('--input-list', help="A space-delimited list of input files, sample names"
        "or conditions, each quoted separately. They must either be all files, or sample/condition"
        "names, and this must be specified with the [--input-type]. Only these will be converted.",
        nargs='*', required='--input-type' in sys.argv, type=str)

    parser.add_argument('--input-type', help="The 'type' of [--input-list], either f (files), "
        "s (samples) or c (conditions).", required='--input-list' in sys.argv, type=str,
        choices=['f', 's', 'c'])

    parser.add_argument('--overwrite', help="If this flag is present, then existing files" 
        "will be overwritten.", action='store_true')

    logging_utils.add_logging_options(parser)
    args = parser.parse_args()
    logging_utils.update_logging(args)

    config = yaml.load(open(args.config))

    required_keys = [
        'riboseq_data',
        'riboseq_samples',
        'riboseq_biological_replicates'
    ]
    utils.check_keys_exist(config, required_keys)

    if args.input_list:
        if args.input_type == 's':
            sample_names = {name: [name] for name in args.input_list}
            condition_names = {}
        elif args.input_type == 'c':
            condition_names = {name: [name] for name in args.input_list}
            sample_names = {}
        else:
            files_only = True
    else:
        files_only = False
        sample_names = config['riboseq_samples']
        condition_names = ribo_utils.get_riboseq_replicates(config)

    # check output path
    if os.path.exists(args.out):
        args.out = os.path.join(args.out, '')
    else:
        msg = "Invalid output path or wrong permission: {}. Quitting.".format(args.out)
        raise OSError(msg)

    use_config_fields = {}
    # generate an AutoSql format (.as) file describing the fields
    if args.configure_fields:
        fields_to_keep = []
        extra_fields = []
        f = open(args.configure_fields, 'r')
        lines = f.readlines()
        f.close()
        as_file = args.out + 'SelectedFields.as'
        f = open(str(as_file), 'w')
        f.write("{} {}\n".format("table", "bedSourceSelectedFields"))
        f.write('''"{}"\n'''.format("Browser extensible data selected fields."))
        f.write("{}\n".format("("))
        n_fields = 0
        for line_no, line in enumerate(lines):
            l = line.strip()
            if not l:
                n_fields = line_no
                break
            fields = l.split(',')
            fields_to_keep.append(fields[0])
            f.write("{}\t{};\t{}\n".format(fields[1],fields[2],fields[3]))
        bed_type = "bed" + str(len(fields_to_keep))
        if n_fields:
            for line in lines[n_fields+1:]:
                l = line.strip()
                fields = l.split(',')
                extra_fields.append(fields[0])
                f.write("{}\t{};\t{}\n".format(fields[1],fields[2],fields[3]))
            bed_type += "+" + str(len(extra_fields))
            fields_to_keep += extra_fields
        f.write("{}\n".format(")"))
        f.close()
        use_config_fields['as_file'] = as_file
        use_config_fields['bed_type'] = bed_type
    else:
        fields_to_keep = bed_utils.bed12_field_names

    if files_only:
        for bed_file in args.input_list:
            if not os.path.exists(bed_file):
                msg = "Could not find the bed file: {}. Quitting.".format(bed_file)
                raise FileNotFoundError(msg)
            bed_file_name = _get_bed(bed_file, fields_to_keep, args)
            bed = bed_file_name + '.tmp.bed'
            bb = bed_file_name + '.bb'
            _convert(bed, bb, use_config_fields, args)

        return

    note_str = config.get('note', None)
    # keep multimappers?
    is_unique = not ('keep_riboseq_multimappers' in config)
    # and the smoothing parameters if present
    fraction = config.get('smoothing_fraction', None)
    reweighting_iterations = config.get('smoothing_reweighting_iterations', None)

    # first the samples
    for name in sorted(sample_names.keys()):

        msg = "Processing sample: {}".format(name)
        logger.info(msg)

        lengths, offsets = ribo_utils.get_periodic_lengths_and_offsets(config, name, is_unique=is_unique)

        predicted_orfs = filenames.get_riboseq_predicted_orfs(config['riboseq_data'], name,
            length=lengths, offset=offsets, is_unique=is_unique, note=note_str, fraction=fraction,
            reweighting_iterations=reweighting_iterations, is_filtered=True)

        if not os.path.exists(predicted_orfs):
            msg = ("Could not find the predictions bed file for {}. ({}). Quitting.".
               format(name, predicted_orfs))
            raise FileNotFoundError(msg)

        bed_file_name = _get_bed(predicted_orfs, fields_to_keep, args)
        bed = bed_file_name + '.tmp.bed'
        bb = bed_file_name + '.bb'
        _convert(bed, bb, use_config_fields, args)

    # then conditions
    lengths = None
    offsets = None
    for name in sorted(condition_names.keys()):

        msg = "Processing condition: {}".format(name)
        logger.info(msg)

        predicted_orfs = filenames.get_riboseq_predicted_orfs(config['riboseq_data'], name,
            length=lengths, offset=offsets, is_unique=is_unique, note=note_str,
            fraction=fraction, reweighting_iterations=reweighting_iterations, is_filtered=True)

        if not os.path.exists(predicted_orfs):
            msg = ("Could not find the predictions bed file for {}. ({}). Quitting.".
                   format(name, predicted_orfs))
            raise FileNotFoundError(msg)

        bed_file_name = _get_bed(predicted_orfs, fields_to_keep, args)
        bed = bed_file_name + '.tmp.bed'
        bb = bed_file_name + '.bb'
        _convert(bed, bb, use_config_fields, args)

if __name__ == '__main__':
    main()
