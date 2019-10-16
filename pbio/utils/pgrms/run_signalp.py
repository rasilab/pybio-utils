#! /usr/bin/env python3

import os

import argparse
import logging

import multiprocessing

import pbio.misc.logging_utils as logging_utils
import pbio.misc.shell_utils as shell_utils

from Bio import SeqIO

logger = logging.getLogger(__name__)


def _split_fasta(args):
    fasta = list(SeqIO.parse(args.fasta, 'fasta'))

    # limit the max number of sequences to pass to SignalP to ~1000
    fasta_len = len(fasta)
    fasta_chunk_size = int(fasta_len / args.num_cpus) + 1
    if fasta_len >= 1000 * args.num_cpus:
        fasta_chunk_size = 1000
    fasta_chunks = []
    for i in range(0, fasta_len, fasta_chunk_size):
        seqs = [s for s in fasta[i:i + fasta_chunk_size]]
        fasta_chunks.append(seqs)

    # write each chunk to disk
    msg = "Writing fasta files to disk."
    logger.info(msg)

    base_name = os.path.split(args.fasta)[1]
    base_name = os.path.splitext(base_name)[0]
    fasta_chunk_ids = []
    for idx, seqs in enumerate(fasta_chunks):
        file_name = os.path.join(args.tmp, '{}.{}.fasta'.format(base_name, idx))
        fasta_chunk_ids.append(file_name)
        SeqIO.write(seqs, file_name, 'fasta')

    return fasta_chunk_ids


def _run_signalp(fasta_chunk_ids, signalp_opt_str, args):
    # create list of cmd so that we can spawn processes asynchronously
    # keep track of output
    all_cmds = []
    all_in_files = []
    all_out_files = []
    for fasta_in in fasta_chunk_ids:
        out = os.path.splitext(fasta_in)[0] + '.out'
        all_out_files.append(out)
        cmd = "signalp {} {} > {}".format(signalp_opt_str, fasta_in, out)
        all_cmds.append(cmd)
        all_in_files.append(fasta_in)

    msg = "Initialise process pool. Calling SignalP..."
    logger.info(msg)

    pool = multiprocessing.Pool(processes=args.num_cpus)
    pool_results = [pool.apply_async(shell_utils.call_if_not_exists, args=(cmd, out_files),
                    kwds={'in_files': [in_files], 'overwrite': args.overwrite})
                    for cmd, out_files, in_files in zip(all_cmds, all_out_files, all_in_files)]
    pool.close()
    pool.join()

    for fasta_in in fasta_chunk_ids:
        try:
            os.unlink(fasta_in)
            msg = "Removing: {}".format(fasta_in)
            logger.info(msg)
        except OSError:
            msg = "Could not remove: {}".format(fasta_in)
            logger.info(msg)

    return all_out_files


def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
            description="This script wraps calls to SignalP v4.1. It handles arbitrary large"
            "multi-fasta files by splitting them into multiple smaller files. SignalP is run on"
            "each file and the results are concatenated in the report file. The program must"
            "be available on the user's path. The default method is used, i.e. results are"
            "chosen from the best predictions from either SignalP-TM or SignalP-noTM networks.")

    parser.add_argument('fasta', help="The input (fasta) file.")
    parser.add_argument('out', help="The output file name. This uses "
                                    "the 'short format' option from SignalP.")
    parser.add_argument('tmp', help="A temporary output directory for writing files. This is"
            "different form the SignalP temporary directory.")

    parser.add_argument('--num-cpus', type=int)
    parser.add_argument('--overwrite', help="If this flag is present, then existing files"
                                            "will be overwritten.", action='store_true')

    parser.add_argument('--use-tmp', help="Option passed to SignalP. Use tmp instead of the"
            "default 'outputDir' set up at installation.", action='store_true')
    parser.add_argument('-k', '--keep', help="Option passed to SignalP. Keep SignalP outputDir/tmp" 
            "directory. Unless specified, this will be removed.", action='store_true')
    parser.add_argument('-c', '--cut', help="Option passed to SignalP. Truncate the input"
        "sequences to the specified length from the N-terminal. Default is 70 residues, 0 disables"
        "truncation.", type=int)
    parser.add_argument('-m', '--minimal', help="Option passed to SignalP. Set minimal predicted"
        "peptide signal length.", type=int)
    parser.add_argument('--mature-fasta', help="Option passed to SignalP. File name to write "
        "mature sequences based on the predictions.", type=str)
    parser.add_argument('--network-type', help="Option passed to SignalP. Use network and models"
        "trained on specified organism type. We use default D-cutoff for noTM networks for all"
        "organisms.", type=str, choices=['euk', 'gram-', 'gram+'], default='euk')

    logging_utils.add_logging_options(parser)
    args = parser.parse_args()
    logging_utils.update_logging(args)

    # create sub-directory for output, fasta files are also re-written in this directory
    if os.path.exists(args.tmp):
        args.tmp = os.path.join(args.tmp, '') + 'signalp/'
        if not os.path.isdir(args.tmp):
            os.makedirs(args.tmp)
    else:
        msg = "Invalid tmp path or wrong permission: {}".format(args.tmp)
        raise OSError(msg)

    # Sort SignalP options
    # redirect SignalP log
    s_verbose = ""
    if args.log_file:
        s_verbose = "-v"
    s_tmp = ""
    if args.use_tmp:
        s_tmp = "-T {}".format(args.tmp)
    s_keep = ""
    if args.keep:
        s_keep = "-k"
    s_cut = ""
    if args.cut:
        s_cut = "-c {}".format(str(args.cut))
    s_minimal = ""
    if args.minimal:
        s_minimal = "-M {}".format(str(args.minimal))
    s_mature = ""
    if args.mature_fasta:
        s_mature = "-m {}".format(args.mature_fasta.strip())
    s_type = ""
    if args.network_type:
        s_type = "-t {}".format(args.network_type.strip())

    signalp_opts = {'v': s_verbose, 'T': s_tmp, 'k': s_keep, 'c': s_cut, 'M': s_minimal,
                    'm': s_mature, 't': s_type}
    signalp_opt_str = ' '.join(signalp_opts.values()).strip()

    if not args.num_cpus:
        args.num_cpus = multiprocessing.cpu_count()

    # split input fasta and write chunks to disk
    fasta_chunk_ids = _split_fasta(args)
    # call SignalP on each chunk
    all_out_files = _run_signalp(fasta_chunk_ids, signalp_opt_str, args)
    # collate the outputs and clean up
    msg = "Collate results and writing to: {}.".format(args.out)
    logger.info(msg)
    # deal with header first
    cmd = '''sed -i -e '1d' -e 's/#/>/g' {}'''.format(all_out_files[0])
    shell_utils.call_if_not_exists(cmd, all_out_files[0],
                in_files=[all_out_files[0]], overwrite=True, call=True)
    # create a temporary file to write all results
    cat_file = os.path.join(args.tmp, 'all_results.cat')
    cmd = 'cat {} > {}'.format(' '.join(all_out_files), cat_file)
    shell_utils.call_if_not_exists(cmd, cat_file,
                in_files=all_out_files, overwrite=True, call=True)
    # finally create output
    cmd = '''awk '!/^#/' {} > {}'''.format(cat_file, args.out)
    shell_utils.call_if_not_exists(cmd, args.out,
                in_files=[cat_file], overwrite=args.overwrite, call=True)

    msg = "Cleaning tmp."
    logger.info(msg)
    all_out_files.append(cat_file)
    for out_file in all_out_files:
        try:
            os.unlink(out_file)
        except OSError:
            msg = "Could not remove: {}".format(out_file)
            logger.info(msg)

if __name__ == '__main__':
    main()
