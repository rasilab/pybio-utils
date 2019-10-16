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

    # limit the max number of sequences to pass to TMHMM to ~1000
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


def _run_tmhmm(fasta_chunk_ids, args):
    # create list of cmd so that we can spawn processes asynchronously
    # keep track of output
    all_cmds = []
    all_in_files = []
    all_out_files = []
    for fasta_in in fasta_chunk_ids:
        out = os.path.splitext(fasta_in)[0] + '.out'
        all_out_files.append(out)
        cmd = "tmhmm -short -workdir {} {} > {}".format(args.tmp, fasta_in, out)
        all_cmds.append(cmd)
        all_in_files.append(fasta_in)

    msg = "Initialise process pool. Calling TMHMM..."
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
            description="This script wraps calls to TMHMM v2.0. It handles arbitrary large"
            "multi-fasta files by splitting them into multiple smaller files. TMHMM is run on"
            "each file and the results are concatenated in the report file. The program must"
            "be available on the user's path.")

    parser.add_argument('fasta', help="The input (fasta) file.")
    parser.add_argument('out', help="The output file name. This uses "
                                    "the 'short format' option from TMHMM.")
    parser.add_argument('tmp', help="A temporary output directory.")

    parser.add_argument('--num-cpus', type=int)
    parser.add_argument('--overwrite', help="If this flag is present, then existing files"
                                            "will be overwritten.", action='store_true')

    logging_utils.add_logging_options(parser)
    args = parser.parse_args()
    logging_utils.update_logging(args)

    # create sub-directory for output, fasta files are also re-written in this directory
    if os.path.exists(args.tmp):
        args.tmp = os.path.join(args.tmp, '') + 'tmhmm/'
        if not os.path.isdir(args.tmp):
            os.makedirs(args.tmp)
    else:
        msg = "Invalid tmp path or wrong permission: {}".format(args.tmp)
        raise OSError(msg)

    if not args.num_cpus:
        args.num_cpus = multiprocessing.cpu_count()

    # split input fasta and write chunks to disk
    fasta_chunk_ids = _split_fasta(args)
    # call TMHMM on each chunk
    all_out_files = _run_tmhmm(fasta_chunk_ids, args)
    # collate the outputs and clean up
    msg = "Collate results and writing to: {}.".format(args.out)
    logger.info(msg)
    cmd = 'cat {} > {}'.format(' '.join(all_out_files), args.out)
    shell_utils.call_if_not_exists(cmd, args.out, in_files=all_out_files,
                                   overwrite=args.overwrite, call=True)
    msg = "Cleaning tmp."
    logger.info(msg)
    for out_file in all_out_files:
        try:
            os.unlink(out_file)
        except OSError:
            msg = "Could not remove: {}".format(out_file)
            logger.info(msg)


if __name__ == '__main__':
    main()
