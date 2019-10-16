###
#   This file contains helper functions for working with BAM files. It mostly
#   wraps calls to pysam and often returns results of interest as pandas data
#   frames.
###

import os
import sys
import logging

logger = logging.getLogger(__name__)

def check_bam_file(filename, check_index=False, raise_on_error=True, logger=logger):
    """ This function wraps a call to "samtools view" and pipes the output to 
        /dev/null. Additionally, it optionally checks that the index is present.
        
        Optionally, it raises an exception if the return code is 
        not 0. Otherwise, it writes a "critical" warning message.

        Args:
            filename (str): a path to the bam file

            check_index (bool): check whether the index is present

            raise_on_error (bool): whether to raise an OSError (if True) or log
                a "critical" message (if false)

            logger (logging.Logger): a logger for writing the message if an
                error is not raised

        Returns:
            bool: whether the file was valid

        Raises:
            OSError: if quickcheck does not return 0 and raise_on_error is True

        Imports:
            misc.utils
    """
    import pbio.misc.utils as utils
    import pbio.misc.shell_utils as shell_utils
    
    programs = ['samtools']
    shell_utils.check_programs_exist(programs)

    dev_null = utils.abspath("dev", "null")

    cmd = "samtools view -h {} > {}".format(filename, dev_null)

    ret = shell_utils.check_call_step(cmd, raise_on_error=False)

    if ret != 0:
        msg = "The bam/sam file does not appear to be valid: {}".format(filename)

        if raise_on_error:
            raise OSError(msg)

        logger.critical(msg)
        return False

    # now look for the index
    if not check_index:
        return True

    cmd = "samtools idxstats {}".format(filename)

    ret = shell_utils.check_call_step(cmd, raise_on_error=False)

    if ret != 0:
        msg = "The bam/sam file does not appear to have a valid index: {}".format(filename)

        if raise_on_error:
            raise OSError(msg)

        logger.critical(msg)
        return False


    # then the file and the index was okay
    return True


def sort_bam_file(bam, sorted_bam, args):
    """ Sort bam file, wrapping a call to "samtools sort" via
    shell_utils.call_if_not_exists.

    Args:
        bam (str): path to bam file
        sorted_bam (str): path to sorted bam file
        args (namespace): calling arguments
    """

    import pbio.misc.shell_utils as shell_utils

    call = not args.do_not_call
    keep_delete_files = args.keep_intermediate_files or args.do_not_call

    sam_tmp_str = ""
    if args.tmp is not None:
        sam_tmp_dir = os.path.join(args.tmp, "{}_samtools".format(args.name))
        if not os.path.exists(sam_tmp_dir):
            os.makedirs(sam_tmp_dir)
        sam_tmp_str = "-T {}".format(sam_tmp_dir)

    cmd = "samtools sort {} -@{} -o {} {}".format(
        bam,
        args.num_cpus,
        sorted_bam,
        sam_tmp_str
    )

    in_files = [bam]
    out_files = [sorted_bam]
    to_delete = [bam]
    file_checkers = {
        sorted_bam: check_bam_file
    }
    shell_utils.call_if_not_exists(cmd, out_files, in_files=in_files,
        file_checkers=file_checkers, overwrite=args.overwrite, call=call,
        keep_delete_files=keep_delete_files, to_delete=to_delete)


def index_bam_file(bam, args):
    """ Index bam file, wrapping a call to "samtools index" via
    shell_utils.call_if_not_exists.

    Args:
        bam (str): path to bam file
        args (namespace): calling arguments
    """

    import pbio.misc.shell_utils as shell_utils

    call = not args.do_not_call

    bai = bam + ".bai"
    cmd = "samtools index -b {} {}".format(bam, bai)
    in_files = [bam]
    out_files = [bai]
    shell_utils.call_if_not_exists(cmd, out_files, in_files=in_files,
        overwrite=args.overwrite, call=call)


def get_pysam_alignment_file(f, mode=None, **kwargs):
    """ This function checks the type of a given object and returns a pysam
        AlignmentFile object. If the object is already an AlignmentFile, it
        is simply returned. If it is a string, then it is treated as a file
        path and opened as an AlignmentFile. Otherwise, an Exception is raised.

        Args:
            f (obj): either an existing AlignmentFile or a path to a bam/sam file

            mode (str): the mode for opening the file, if f is a file path

            **kwargs : other keyword arguments to pass to the AlignmentFile
                constructor

        Returns:
            pysam.AlignmentFile: the AlignmentFile for f

        Raises:
            TypeError: if f is neither an AlignmentFile nor string
        Import:
            pysam
    """

    import pysam

    if isinstance(f, pysam.AlignmentFile):
        return f

    elif isinstance(f, str):
        return pysam.AlignmentFile(f, mode=mode, **kwargs)

    else:
        msg = "Could not interpret value as pysam.AlignmentFile: {}".format(f)
        raise ValueError(msg)


def count_unique_reads(bam, **kwargs):
    """ Count the number of distinct reads, based on the read names.

    Parameters
    ----------
    bam: string or pysam AlignmentFile
        Either an existing AlignmentFile or a path to a bam/sam file

    \*\*kwargs: key-value pairs
        Other arguments for the pysam.AlignmentFile constructor. These are
        ignored if bam is already an AlignmentFile.

    Returns
    -------
    reads_and_counts: dict of string -> int
        A map from the name of each unique read to the number of times it
        occurs.The read names are derived from 
        pysam.AlignedSegnment.query_name.
    """
    import collections

    alignments = get_pysam_alignment_file(bam, **kwargs)

    reads_and_counts = collections.defaultdict(int)
    for a in alignments:
        reads_and_counts[a.query_name] += 1

    return reads_and_counts

def count_aligned_reads(bam, **kwargs):
    """ Count the number of aligned reads in a bam file.

    Internally, this function uses get_length_distribution, so please see its
    documentation for more details. In particular, though, it skips secondary
    and supplementary alignments to avoid double-counting multimappers.

    Parameters
    ----------
    bam: string or pysam.AlignmentFile
        Either an open AlignmentFile or a path to a bam file

    \*\*kwargs: key-value pairs
        Other arguments for the pysam.AlignmentFile constructor. These are
       
    Returns
    -------
    num_aligned_reads: int
        the number of reads with at least one alignment in the file
    """
    import numpy as np

    length_distribution_df = get_length_distribution(bam, **kwargs)
    num_aligned_reads = np.sum(length_distribution_df['count'])
    return num_aligned_reads

def get_length_distribution(bam, progress_bar=False, **kwargs):
    """ Count the number of aligned reads of each length is the bam file.

    This function excludes "unmapped", "secondary", and "supplementary"
    alignments, so each aligned read is only counted once.

    Parameters
    ----------
    bam: string or pysam.AlignmentFile
        Either an open AlignmentFile or a path to a bam file

    progress_bar: bool
        Whether to show a progress bar

    \*\*kwargs: key-value pairs
        Other arguments for the pysam.AlignmentFile constructor. These are
        ignored if bam is already an AlignmentFile.

    Returns
    -------
    length_distribution_df: pd.DataFrame
        A data frame with the columns: length, count
    """
    import collections
    import tqdm
    import pbio.misc.pandas_utils as pandas_utils

    msg = "Opening alignment file"
    logger.debug(msg)
    
    alignment_file = get_pysam_alignment_file(bam, **kwargs)
    length_distribution = collections.defaultdict(int)

    if progress_bar:
        msg = "Counting number of alignments"
        logger.debug(msg)

        num_alignments = alignment_file.count(until_eof=True)
        alignment_file.reset()

    alignments = alignment_file.fetch(until_eof=True)

    if progress_bar:
        alignments = tqdm.tqdm(
            alignments,
            leave=True,
            file=sys.stdout,
            total=num_alignments
        )

    msg = "Collecting read length distribution"
    logger.debug(msg)

    for alignment in alignments:
        if (alignment.is_unmapped or 
            alignment.is_secondary or 
            alignment.is_supplementary):
            continue
        length_distribution[alignment.qlen] += 1

    msg = "Converting distribution to data frame"
    logger.debug(msg)

    length_distribution_df = pandas_utils.dict_to_dataframe(
        length_distribution,
        key_name='length',
        value_name='count'
    )

    return length_distribution_df

def count_uniquely_mapping_reads(bam, logger=logger):
    """ Count the number of reads in the bam file which map uniquely.

    Specifically, this function looks for the "NH:i:1" flag in the alignments.

    Parameters
    ----------
    bam : string
        The filename of the bam file

    logger : logging.logger
        A logger to use for writing status information

    Returns
    -------
    num_uniquely_mapping_reads : int
        The number of reads which map uniquely (i.e., have the "NH:i:1" flag)
    """
    import pbio.misc.shell_utils as shell_utils

    cmd = "samtools view -h {} | grep '^@\|NH:i:1	'| wc -l".format(bam)
    res = shell_utils.check_output(cmd)
    num_uniquely_aligned_reads = int(res.strip())
    return num_uniquely_aligned_reads


def filter_by_unique_tag_value_pair(alignments, unique_tag, unique_value):
    """Return alignment satisfying a given tag from the SAM optional
    fields. Standard tags are displayed as TAG:TYPE:VALUE. There is no
    check on the validity and/or presence of the pair TAG:VALUE.
        Arguments:
            alignments: the iterator over the collection of reads
            unique_tag (str): the SAM TAG
            unique_value: the SAM VALUE (type must match TYPE)
        Return:
            generator
    """

    for a in alignments:
        tag = a.get_tag(unique_tag)
        if tag == unique_value:
            yield a
        else:
            yield None


def filter_stranded_by_flag(alignments, strand):
    """Return alignment satisfying a given orientation, based
    on the SAM flags. See `filter_stranded_library`.
        Arguments:
            alignments: the iterator over the collection of reads
            strand (str): library strandedness (fr or rf)
        Return:
            generator
        """

    for a in alignments:
        keep = a.is_reverse
        if strand == 'fr':
            keep = not keep
        if keep:
            yield a
        else:
            yield None


def filter_stranded_library(align_in, align_out, strand, logger=logger):
    """Remove alignments based on library strandedness. Alignments
    must be in transcript coordinates, otherwise filtering cannot
    be based only on the SAM flag. This does not filter secondary
    alignments, if present.

        Arguments:
            align_in (string) : the "stranded" BAM filename

            align_out (string) : the output BAM filename

            strand (str) : either "fr" (secondstrand) or "rf" (firststrand)

            logger (logging.Logger): a logger to which status messages will
                be written
        Return:
            None
    """

    import tqdm

    msg = "Filter alignments based on library strandedness."
    logger.debug(msg)

    strand_choices = ['fr', 'rf']
    if strand not in strand_choices:
        msg = "Invalid strand flag, expected one of {}".format(strand_choices)
        raise ValueError(msg)

    bam = get_pysam_alignment_file(align_in)
    alignments = bam.fetch()
    num_alignments = bam.count()

    out_bam = get_pysam_alignment_file(align_out, mode="wb", template=bam)

    for a in tqdm.tqdm(filter_stranded_by_flag(alignments, strand),
                       leave=True, file=sys.stdout, total=num_alignments):
        if a is not None:
            out_bam.write(a)
    out_bam.close()


def remove_multimappers(align_in, align_out, n_map=1, index=True,
                        logger=logger):
    """Remove multimappers from align_in and write uniquely mapped reads
    to align_out, and index that file. This is solely based on the SAM NH
    attribute, which enumerates multiple alignments of a read starting with 1.

    N.B. If the file is sorted, the sort order is preserved. The file will
    be indexed only if it is already sorted.

    N.B. This requires the SAM NH tag. If using STAR, default output include
    NH HI AS nM attributes. See also STAR options --outSAMprimaryFlag and
    --outSAMmultNmax, which may affect the behaviour of this function.

        Arguments:
            align_in (string) : the BAM filename which may
                contain multimappers

            align_out (string) : the BAM filename which contains
                only the uniquely mapped reads.

            n_map (int) : the SAM NH tag value corresponding to the number
            of loci a read maps to.

            index (bool) : index the output.

            logger (logging.Logger): a logger to which status messages will
                be written
        Return:
            None
    """

    import tqdm
    import pbio.misc.shell_utils as shell_utils

    msg = "Removing multimappers."
    logger.debug(msg)

    bam = get_pysam_alignment_file(align_in)
    alignments = bam.fetch()
    num_alignments = bam.count()

    out_bam = get_pysam_alignment_file(align_out, mode="wb", template=bam)

    for a in tqdm.tqdm(filter_by_unique_tag_value_pair(alignments, 'NH', n_map),
                       leave=True, file=sys.stdout, total=num_alignments):
        if a is not None:
            out_bam.write(a)
    out_bam.close()

    if index:
        msg = "Attempt to index uniquely mapped reads."
        logger.debug(msg)
        try:
            cmd = "samtools index -b {}".format(align_out)
            shell_utils.check_call(cmd, call=True)
        except:
            msg = ("Failed to index uniquely mapped reads."
                   "Check if BAM file is sorted!")
            logger.debug(msg)


# temporary keep this function ** need to adjust calls in Rp-Bp, B-tea
def remove_multimapping_reads(align_in, align_out, call=True, tmp=None, 
            logger=logger):
    """ This functions wraps calls to samtools which remove the multimapping
        reads from the align_in file and writes them to align_out. It also
        indexes the output file.

        N.B. This function *does not* attempt to keep the "best" alignment.
        Rather, it discards all reads with more than one alignment.

        Args:
            align_in (string) : the filename of the BAM file which may 
                contain multimappers

            align_out (string) : the filename of the BAM file which contains
                only the uniquely mapping reads.

            call (boolean) : If desired, the function will only print out
                the calls to samtools, but not actually call them.

            tmp (string) : the path where temporary files for samtools sort
                will be stored. If not given, then the samtools default tmp
                choice will be used.

            logger (logging.Logger): a logger to which status messages will
                be written

        Returns:
            None
    """
    import pbio.misc.shell_utils as shell_utils

    programs = ['samtools']
    shell_utils.check_programs_exist(programs)

    sam_tmp_str = ""
    if tmp is not None:
        sam_tmp_str = "-T {}".format(tmp)

    msg = "Removing multimappers and sorting the remaining reads"
    logger.debug(msg)

    cmd = ("samtools view -h {} | grep '^@\|NH:i:1	' | samtools view -bS  - "
            "| samtools sort -  -o {} {}".format(align_in, align_out, sam_tmp_str))
    shell_utils.check_call(cmd, call=call)

    # index
    msg = "Indexing the sorted reads"
    logger.debug(msg)

    cmd = "samtools index -b {}".format(align_out)
    shell_utils.check_call(cmd, call=call)

def get_five_prime_ends(bam, progress_bar=True, count=True, logger=logger):
    """ This function extracts the 5' ends of all reads in the bam file. It
        returns the results as a BED6+1 data frame. The additional field is
        the length of the read of the alignment. A tqdm progress bar is used
        to show the progress.

        Args:
            bam (string or pysam.AlignmentFile): either the path to a bam file,
                or an open pysam.AlignmentFile. See get_pysam_alignment_file
                for more details about how this is interpreted.

            count (bool): if True, then the number of reads in the bam file
                will first be counted. This allows the progress bar to be more
                accurate.

        Imports:
            numpy
            pandas
            pysam
            tqdm
    """
    import numpy as np
    import pandas as pd
    import pysam
    import tqdm

    # first, make sure we have an alignment file
    bam = get_pysam_alignment_file(bam)

    if count:
        msg = "Counting the number of alignments"
        logger.debug(msg)
        num_alignments = bam.count()
    else:
        num_alignments = None

    alignments = bam.fetch()

    lengths = np.zeros(num_alignments)
    five_prime_ends = np.zeros(num_alignments)
    seqnames = np.full(num_alignments, None, dtype=object)
    strands = np.full(num_alignments, None, dtype=object)

    msg = "Extracting 5' ends of reads from alignments"
    logger.debug(msg)

    for i, a in enumerate(tqdm.tqdm(alignments, leave=True, file=sys.stdout, total=num_alignments)):
        five_prime_ends[i] = a.reference_start
        if a.is_reverse:
            five_prime_ends[i] = a.reference_end

        lengths[i] = a.qlen
        seqnames[i] = a.reference_name
        strands[i] = "+"
        
        if a.is_reverse:
            strands[i] = "-"

    msg = "Constructing data frame with 5' ends of reads from alignments"
    logger.debug(msg)

    alignment_df = pd.DataFrame()
    alignment_df['seqname'] = seqnames
    alignment_df['start'] = five_prime_ends
    alignment_df['end'] = five_prime_ends + 1
    alignment_df['id'] = "."
    alignment_df['score'] = 0
    alignment_df['strand'] = strands
    alignment_df['length'] = lengths

    return alignment_df

def get_read_identifiers(bam):
    """ Collect the read identifiers from the bam file in a set.

    Parameters
    ----------
    bam : string or pysam.AlignmentFile
        The bam file. If it is a pysam.AlignmentFile, it will be "consumed"
        during this operation.

    Returns
    -------
    read_ids : set
        The read identifiers from the bam file, in a set
    """
    bam = get_pysam_alignment_file(bam)

    ids = {
        alignment.query_name for alignment in bam
    }

    return ids


