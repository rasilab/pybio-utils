
"""Module provides functions to interact with and handle options
for STAR and Flexbar.
"""

import logging
logger = logging.getLogger(__name__)


# generic functions


def get_final_args(default_args, args):
    """ Extract the final flags and options from args to construct the
    command line call to programs. Arguments args override defaults.

        Parameters
        ---------
        default_args: dict of default arguments
        args: list of string arguments

        Returns
        -------
        final_options_str: string
            a string containing the final options suitable to pass to another command
        """

    from collections import defaultdict

    # create dict from list of strings
    final_options = defaultdict(list)
    final_options_str = ''
    if args is not None:
        for opt in args:
            final_options['{}'.format(opt.rsplit()[0].strip('--'))].append(' '.join(opt.rsplit()[1:]))

        final_options_str = ' '.join(['--{} {}'.format(key, val) for (key, values) in
                                      final_options.items() for val in values])

    # now search if key exist, else use default
    for key, val in default_args.items():
        if key not in final_options:
            if isinstance(val, (str,)) and (len(str(val)) > 0):
                val = val
            elif isinstance(val, (int, float)) and (len(str(val)) > 0):
                val = str(val)
            elif len(val) > 0:
                # assume this is a list
                val = ' '.join(str(v) for v in val)
            else:
                msg = "Default argument type not supported, or no argument given"
                logger.critical(msg)

            key_val_str = '--{} {}'.format(key, val)
            final_options_str = '{}'.format(' '.join([final_options_str, key_val_str]))

    return final_options_str


# STAR

def create_star_tmp(tmp_path:str, tmp_name:str='STAR'):
    """ Ensure the specified directory is ready for use as the temp directory
    for STAR.
    
    In particular, it ensures that tmp_path exists and that tmp_path/tmp_name 
    does not exist (so STAR can create it).

    N.B. This function *does not* use sophisticated file locking mechanisms,
         but it should suffice for "typical" use cases.

    Parameters
    ----------
    tmp_path: string
        the path to the STAR tmp directory NOT including the STAR tmp 
        directory itself

    tmp_name: string
        the name of the STAR tmp directory

    Returns
    -------
    star_temp_path: string
        the complete path to the STAR tmp directory (tmp_path/tmp_name)

    Side effects
    ------------
    tmp_path will be created, if it does not already exist.

    tmp_path/tmp_name will be removed if it already existed.
    """
    import os
    import shutil

    # the path to the star directory must exist
    if not os.path.exists(tmp_path):
        os.makedirs(tmp_path)
    
    star_tmp_dir = os.path.join(tmp_path, tmp_name)
    # remove the actual star directory, if it exists
    if os.path.exists(star_tmp_dir):
        shutil.rmtree(star_tmp_dir)

    return star_tmp_dir

def get_star_index_files(path:str, with_sjdb_files:bool=False):
    """ Find the file names necessary for a STAR index rooted at path.

        Parameters
        ----------
        path: string
            the path to the directory containing the STAR index

        with_sjdb_files: bool
            whether to include the splice junction database files in the 
            returned list

        Returns
        -------
        star_index_files: list of strings
            the complete paths to all of the files expected for a STAR index 
            rooted at path (include sjdb/exon files, if specified)
    """
    import os

    star_files = [
        "chrLength.txt",
        "chrNameLength.txt",
        "chrName.txt",
        "chrStart.txt",
        "Genome",
        "genomeParameters.txt",
        "SA",
        "SAindex"
    ]

    sjdb_files = [
        "exonGeTrInfo.tab",
        "exonInfo.tab",
        "geneInfo.tab",
        "sjdbInfo.txt",
        "sjdbList.fromGTF.out.tab",
        "sjdbList.out.tab",
        "transcriptInfo.tab"
    ]

    if with_sjdb_files:
        sf = star_files + sjdb_files
    else:
        sf = star_files

    star_index_files = [ os.path.join(path, star_file) for star_file in sf]
    return star_index_files

def read_star_tr_file(filename:str):
    """ Read a STAR transcript info file into a data frame. 
    
    The file is assumed to contain the number of transcripts on the first row.
    Each row then contains information about each transcript. The columns are 
    assumed to have the following meanings:

        - [ID]      transcript id
        - [S]       transcript (genomic) start
        - [E]       transcript (genomic) end
        - [Emax]    transcript (genomic) end minus "max" (???)
        - [Str]     strand
        - [ExN]     number of exons
        - [ExI]     index of first exon (in STAR exon database)

    N.B. These semantics come from the STAR 2.5.1b source code. They appear 
    stable but could change for new (or old) versions of STAR. The column names
    (in square brackets) follow the internal STAR naming conventions.

    Parameters
    ----------
    filename: string
        the complete path to the STAR transcriptInfo.tab file

    Returns
    -------
    transcript_info: pd.DataFrame
        a data frame containing the transcripts info. The order of the 
        transcripts will be the same as the order in the file.
    """
    import pandas as pd

    column_names = ['ID', 'S', 'E', 'Emax', 'Str', 'ExN', 'ExI']
    transcript_info = pd.read_csv(
        filename, 
        header=None, 
        names=column_names, 
        skiprows=1, 
        sep='\t'
    )

    return transcript_info


def add_star_options(parser, star_executable:str="STAR"):
    """ Add options to a cmd parser to call STAR.

    N.B. This is primarily intended for use with the rp-bp and b-tea projects.

    Parameters
    ----------
    parser: argparse.ArgumentParser
        The parser to which the options will be added

    star_executable: string
        The name of the star executable. For example, "STARlong" is typically
        used for mapping longer reads, while "STAR" is for most HTS data.
    star_options: list of strings
        Additional options to pass to STAR
    """

    star_options = parser.add_argument_group("STAR options")

    star_options.add_argument('--star-executable', help="The name of the STAR "
        "executable", default=star_executable)

    star_options.add_argument('--star-options', help="""A space-delimited
        list of options to pass to star (for the mapping step only). Each option
        must be quoted separately as in "--starOption value", using soft 
        quotes, where '--starOption' is the long parameter name from star, and 'value'
        is the value given to this parameter. If specified, star options will override 
        default settings.""", nargs='*', type=str)


def get_star_options_string(args):
    """ Extract the flags and options specified for STAR added with 
    add_star_options.

    Parameters
    ---------
    args: argparse.Namespace
        The parsed arguments

    Returns
    -------
    star_options: string
        a string containing STAR options suitable to pass to another command
    """
    import shlex

    args_dict = vars(args)

    star_options = ['star_executable']

    # create a new dictionary mapping from the flag to the value
    star_options = {'--{}'.format(o.replace('_', '-')): args_dict[o]
        for o in star_options if args_dict[o] is not None}

    s = ' '.join("{} {}".format(k, shlex.quote(v))
                    for k, v in star_options.items())

    # if additional options
    if args_dict['star_options']:
        star_additional_options_str = "--star-options {}".format(
            ' '.join(shlex.quote(star_option) for star_option
                     in args_dict['star_options']))
        s = "{}".format(' '.join([s, star_additional_options_str]))

    return s


# Flexbar


def add_flexbar_options(parser):
    """ Add options to a cmd parser to call flexbar.

    N.B. This is primarily intended for use with the rp-bp and b-tea projects.

    Parameters
    ----------
    parser: argparse.ArgumentParser
        The parser to which the options will be added

    flexbaroptions: list of strings
        Additional options to pass to flexbar
    """

    flexbar_options = parser.add_argument_group("Flexbar options")

    flexbar_options.add_argument('--flexbar-options', help="""Optional argument: a space-delimited 
        list of options to pass to flexbar. Each option must be quoted separately as in 
        "--flexbarOption value", using soft quotes, where '--flexbarOption'
        is the long parameter name from flexbar and 'value' is the value given to this parameter. 
        If specified, flexbar options will override default settings.""", nargs='*', type=str)


def get_flexbar_options_string(args):
    """ Extract the flags and options specified for flexbar added with
    add_flexbar_options.

    Parameters
    ---------
    args: argparse.Namespace
        The parsed arguments

    Returns
    -------
    flexbar_options: string
        a string containing flexbar options suitable to pass to another command
    """
    import shlex

    args_dict = vars(args)

    s = ""
    if args_dict['flexbar_options']:
        flexbar_option_str = "--flexbar-options {}".format(
            ' '.join(shlex.quote(flx_op) for flx_op in args_dict['flexbar_options']))
        s = "{}".format(' '.join([s, flexbar_option_str]))

    return s


# Bowtie 2


def get_bowtie2_index_files(base_index_name):
    """ This function returns a list of all of the files necessary for a Bowtie2 index
        that was created with the given base_index_name.

        Args:
            base_index_name (string): the path and base name used to create the bowtie2 index

        Returns:
            list of strings: the paths to all of the files expected for a Bowtie2 index
                based on the provided index_name
    """
    bowtie_extensions = ['.1.bt2', '.2.bt2', '.3.bt2', '.4.bt2', '.rev.1.bt2', '.rev.2.bt2']

    bowtie_files = ['{}{}'.format(base_index_name, bowtie_extension)
                    for bowtie_extension in bowtie_extensions]

    return bowtie_files