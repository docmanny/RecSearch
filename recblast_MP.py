import logging
import re
import subprocess
import sqlite3
from itertools import product, repeat
from functools import partial, reduce
from collections import OrderedDict
from contextlib import redirect_stdout
from datetime import datetime as dt
from inspect import isgenerator
from io import StringIO
from operator import itemgetter
from pathlib import Path
from time import sleep
from copy import deepcopy

import multiprocess as multiprocessing
from Bio import SearchIO, SeqIO, SeqFeature
from Bio import __version__ as bp_version
from Bio.Blast import NCBIXML, NCBIWWW
from Bio.Blast.Record import Blast as BioBlastRecord
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
# import multiprocessing
from BioSQL import BioSeqDatabase
from BioSQL.BioSeq import DBSeqRecord

from Auxilliary import print, ProgressBar, merge_ranges, translate_annotation


def _percent_identity_searchio(hit, is_protein=True):
    """Calculates percent identity based on entire hit. Adapted from UCSC BLAT FAQ and Biopython."""
    from math import log
    size_mul = 3 if is_protein else 1
    qali_size = size_mul * sum([i[-1] - i[0] for i in merge_ranges([(hsp.query_start, hsp.query_end) for hsp in hit])])
    tali_size = sum([i[-1] - i[0] for i in merge_ranges([(hsp.hit_start, hsp.hit_end) for hsp in hit])])
    ali_size = min(qali_size, tali_size)
    if ali_size <= 0:
        return 0
    size_dif = qali_size - tali_size
    size_dif = 0 if size_dif < 0 else size_dif
    sum_match = sum([i.match_num for i in hit])
    sum_rep = sum([i.match_rep_num for i in hit])
    sum_mismatch = sum([i.mismatch_num for i in hit])
    total = size_mul * (sum_match + sum_rep + sum_mismatch)
    if total != 0:
        millibad = (1000 * (sum([i.mismatch_num for i in hit]) * size_mul + sum([i.query_gap_num for i in hit]) +
                            round(3 * log(1 + size_dif)))) / total
    else:
        raise Exception('Somehow your total in the percent_identity function was 0, so you broke the script!')
    perc_ident = 100 - (millibad*0.1)
    return perc_ident


def get_searchdb(search_type, species, db_loc, verbose=1, indent=0):
    """Finds and returns the appropriate search database for the given species and search type.

    This function automates the process of selecting the search database needed by the selected search program,
    like BLAST or BLAT, so that the user does not need to preoccupy themselves with providing said information
    for a large number of species. For BLAST* that depend on protein databases (BLASTP and BLASTX), the function
    searches for files matching the form 'Genus_species_protein.*' in the given directory; for BLAST* that depend
    on DNA databases (BLASTN, TBLASTN, and TBLASTX), it instead looks for files 'Genus_species_genome.*'.
    If '-transcript' is added to the end of any of the DNA-dependent BLAST*, then instead the function will
    search for files in the style of 'Genus_species_transcript.*'. In the case of BLAT searches, the program will
    similarly search for 'Genus_species*.2bit', or for 'Genus_species*transcript.2bit' if '-transcript' is added
    after the search type.
    In all usage cases, if the program does not find files matching the 'Genus_species' format, it will try to
    find the files using a case-insensitive search using the 6-letter abbreviated form of the species name.

    Usage::
    >>> get_searchdb('blastp', 'Homo sapiens', '/path/to/search/files')
    /path/to/search/files/Homo_Sapiens_protein.*
    >>> get_searchdb('tblastn', 'Homo sapiens', '/path/to/search/files')
    /path/to/search/files/HomSap_genome.*
    >>> get_searchdb('blastn-transcript', 'Homo sapiens', '/path/to/search/files')
    /path/to/search/files/HomSap_transcript.*
    >>> get_searchdb('blat', 'Homo sapiens', '/path/to/search/files')
    /path/to/search/files/HomSap.2bit
    >>> get_searchdb('blat-transcript', 'Homo sapiens', '/path/to/search/files')
    /path/to/search/files/HomSap_transcript.2bit

    Arguments::
    :param str search_type: The name of the search method (blast or blat, and sub-type: blastp, blastn, blat, tblat...)
    :param str species: Name of species associated with the database. If there is a space, it will be replaced with an
    underscore.
    :param str db_loc: Path to folder containing collection of search databases.
    :param int verbose: How verbose should the output be. Zero suppresses all output, 2 is max verbosity.
    :param int indent: Indent level for printed output.
    :return str:  Path to the identified search database.
    """
    if verbose:
        print('Blast DB set to auto, choosing blast_db...', indent=indent)
    species = species.replace(' ', '_')
    if verbose > 1:
        print('Blast DB location set to: ', db_loc, indent=indent)
    if search_type.lower() in ['blastp', 'blastx']:
        db_type = 'protein'
    elif search_type.lower() in ['blastn', 'tblastn', 'tblastx']:
        db_type = 'genome'
    elif search_type.lower() in ['blastn-transcript', 'tblastn-transcript', 'tblastx-transcript']:
        db_type = 'transcript'
    elif search_type.lower() in ['blat', 'tblat', 'translated_blat',
                                 'untranslated_blat', 'oneshot blat', 'oneshot tblat']:
        db_type = 'blat'
    elif search_type.lower() in ['blat-transcript']:
        db_type = 'blat-transcript'
    elif search_type.lower() in ['tblat-transcript']:
        db_type = 'tblat-transcript'
    else:
        print('Unable to determing blast db type!', indent=indent)
        raise Exception('Improper search type given: ', search_type)
    if verbose > 1:
        print('DB type: ', db_type, indent=indent)
    db_path = Path(db_loc).absolute()
    if not db_path.exists():
        db_path = Path(db_loc)
    if db_path.exists() and db_path.is_dir():
        if db_type == 'blat':
            glob_path = [i for i in db_path.glob('{0}*.2bit'.format(species.replace(' ', '_')))]
        elif db_type in ['blat-transcript', 'tblat-transcript']:
            glob_path = [i for i in db_path.glob('{0}*transcript.2bit'.format(species.replace(' ', '_')))]
        else:
            glob_path = [i for i in db_path.glob('{0}_{1}*'.format(species.replace(' ', '_'), db_type))]
        if glob_path == []:
            if verbose:
                print('No DB found! Trying again with abbreviated species name')
            species_abbv = ''.join([i[0:3] for i in species.title().split(' ')])
            # making it insensitive to case for Glob
            species_abbv_insensitive = ''.join(['[{0}{1}]'.format(c.lower(), c.upper()) for c in species_abbv if c.isalpha()])
            if verbose:
                print('Abbreviated species name: ', species_abbv, indent=indent)
                print('RegEx species abbreviation: ', species_abbv_insensitive, indent=indent)
            if db_type == 'blat':
                glob_path = [i for i in db_path.glob('{0}*.2bit'.format(species_abbv_insensitive))]
            elif db_type in ['blat-transcript', 'tblat-transcript']:
                glob_path = [i for i in db_path.glob('{0}*transcript.2bit'.format(species_abbv_insensitive))]
            else:
                glob_path = [i for i in db_path.glob('{0}_{1}*'.format(species_abbv_insensitive, db_type))]
        try:
            if verbose:
                print(glob_path, indent=indent)
            if isinstance(glob_path, list):
                search_db = sorted(glob_path, reverse=True)[0].stem
            else:
                search_db = glob_path.stem
        except IndexError:
            print('WARNING: COULD NOT FIND DATABASE! ABORTING!', indent=indent)
            raise Exception('DatabaseError:', 'No databases were found!')
    else:
        raise Exception('DB_Path {} does not exist!'.format(str(db_path)))
    if verbose:
        print('{0} DB chosen: {1}'.format(search_type, search_db), indent=indent)
    return search_db


def blat_server(twobit, order='start', host='localhost', port=20000, type='blat', log='/dev/null', species=None,
                BLASTDB='/usr/db/blat', verbose = 1, indent=0, try_limit=10, **gfserver_kwargs):
    """Convenience function that controls a gfServer. Still in alpha.

    This function serves as a python wrapper for the Bash gfServer command. The user can either provide a .2bit file,
    or else can provide a species and set 'twobit="auto"' to have the function use 'get_searchdb()' to find a .2bit file
    automatically. By default, the function is set to start up a new gfServer instance, but using the 'order' parameter,
    the user can execute any of the standard gfServer commands such as 'stop' and 'status'.
    To start a gfServer, the function first probes the selected port (default is 20000) to ensure its unused; if it is
    currently in use, the program then goes port-by-port in ascending order until it finds an empty port to use for the
    server. Then, it simply calls the gfServer command with all the keyword arguments required, as well as with any
    extra arguments provided by the user.

    Usage::
    >>>blat_server(twobit='hg38.2bit', port=20000, verbose=3)
    gfServer start localhost 20001 -canStop -stepSize=5 hg38.2bit
    # Waits 30 seconds, then starts calling 'gfServer status localhost 20001' every 30 seconds for 5 minutes
    # If at any point 'gfServer status' returns something that is not an error or "Couldn't connect...", it
    # returns the port where the server was opened.
    20001
    >>>blat_server(twobit='auto', port=20000, species='Homo sapiens', verbose=3)
    # Calls get_searchdb('blat', 'Homo sapiens', db_loc=BLATDB)
    # Internally, will return a .2bit file such as 'Homo_sapiens.2bit'
    20001
    >>>blat_server(twobit='hg38.2bit', port=20000, order='status', verbose=3)
    # If the server is active:
    1
    >>>blat_server(twobit='hg38.2bit', port=20000, order='status', verbose=3)
    # If the server either has not been started or is not yet active:
    0
    >>>blat_server(twobit='hg38.2bit', port=20000, order='status', verbose=3)
    # If the server returns an error
    Exception(...)


    :param str twobit: A path to the .2bit file to be used for the server. Can also be set to 'auto'.
    :param str order: A command for gfServer. Can be one of the following: start, stop, status, files, query (requires
    a nucleotide sequence in fasta format), protQuery (requires a protein sequence in fasta format), transQuery
    (requires a nucleotide sequence in fasta format), pcr (requires arguments fPrimer, rPrimer, maxDistance), direct
    (requires probe.fa, file(s).nib), or pcrDirect (requires fPrimer, rPrimer, file(s).nib).
    :param str host: Address at which to host the server.
    :param int port: Port number that will be assigned to server. If in use, will test new port number in increments of
    1 until a free port is found.
    :param str type: Type of server to be hosted. 'blat' will start a DNA server, 'tblat' will start a DNAX server for
    protein queries.
    :param str log: Path and name of log file to be written.
    :param str species: Species name that get_searchdb() will use to find .2bit file when twobit='auto'.
    :param str BLASTDB: Path to the folder containing .2bit file.
    :param int verbose: Level of verbosity of function output. 0 suppresses all output, 3 is max verbosity.
    :param int indent: Indentation level of print output.
    :param int try_limit: Number of tries at 30-second intervals that function should probe the gfServer before timeout.
    :param gfserver_kwargs: keyword arguments to be passed on to gfServer.
    :return: if order='start', returns the port of the new gfServer; if order='status', returns 0 if there was no
    connection, or 1 if the server is active and responding.
    """
    # Regular: gfServer start localhost portX -stepSize=5 -log=untrans.log database.2bit
    # Prot>DNAX:  gfServer start localhost portY -trans -mask -log=trans.log database.2bit
    gfserver_suppl_args = list()
    if twobit == 'auto' and order != 'stop':
        if verbose:
            print('2bit set to auto: searching for 2bit file for species ', species, indent=indent)
        blat_2bit = get_searchdb(search_type='blat', species=species, db_loc=BLASTDB,
                                 verbose=verbose, indent=indent + 1)
        twobit = Path(BLASTDB, blat_2bit + '.2bit').absolute()
        twobit = str(twobit) if twobit.is_file() else None
        if twobit is None:
            raise Exception('Invalid 2bit file!')
    for key, item in gfserver_kwargs.items():
        if key == 'order':
            order = item
        elif key == 'host':
            host = item
        elif key == 'port':
            port = item
        else:
            gfserver_suppl_args.append('-{0}={1}'.format(key, item))
    if order == 'status':
        gfcheck = subprocess.Popen('gfServer status {0} {1}'.format(str(host), str(port)), stdout=subprocess.PIPE,
                                   stderr=subprocess.STDOUT, universal_newlines=True, shell=True,
                                   executable='/bin/bash')
        out, _ = gfcheck.communicate()
        if "couldn't connect to localhost" in out.lower():
            return 0
        elif "error" in out.lower():
            raise Exception(out)
        else:
            return 1
    elif order == 'stop':
        subprocess.check_call('gfServer stop {0} {1}'.format(str(host), str(port)), stdout=subprocess.PIPE,
                              stderr=subprocess.STDOUT, universal_newlines=True, shell=True,
                              executable='/bin/bash')
        return
    else:
        print(order)
        #Todo: make the portsniffer its own function and make sure it works properly.
        portfinder = subprocess.check_output('/home/manny/Scripts/oneshot/checkifportisopen.sh {}'.format(str(port)),
                                             universal_newlines=True, shell=True, executable='/bin/bash')
        port = portfinder.rstrip()

        gfserver_cmd = ['gfServer', str(order), str(host), str(port), '-canStop']
        if type == 'blat':
            gfserver_cmd.append('-stepSize=5')
        elif type == 'tblat':
            gfserver_cmd += ['-trans', '-mask']
        if gfserver_suppl_args != []:
            gfserver_cmd += gfserver_suppl_args
        gfserver_cmd_str = ' '. join(gfserver_cmd+[twobit])
        if verbose > 2:
            print(gfserver_cmd_str, indent=indent)
        subprocess.Popen(gfserver_cmd_str, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                         universal_newlines=True, shell=True, executable='/bin/bash')
        tries = 0
        while tries <= try_limit:
            sleep(30)
            gfcheck = subprocess.Popen('gfServer status {0} {1}'.format(str(host), str(port)), stdout=subprocess.PIPE,
                                        stderr=subprocess.STDOUT, universal_newlines=True, shell=True,
                                        executable='/bin/bash')
            out, _ = gfcheck.communicate()
            if verbose > 2:
                print(out)
            if "couldn't connect to localhost" in out.lower():
                tries += 1
            elif "error" in out.lower():
                raise Exception(out)
            else:
                if verbose:
                    print(out)
                return port
        if tries > try_limit:
            raise Exception('Timed out!')


def id_ranker(record, perc_score, expect, perc_length, perc_ident, perc_span=0.1, min_hsps=1, hsp_cumu_score=True,
              seq_method = 'whole', align_scorelist=list(), indent=0, verbose=1, method='all', samestrand=True):
    """Filters results based on score, expectation value, length, percent identity, and span; returns a sorted list.

    :param record: Either a SearchIO.QueryResult or a Bio.Blast.Record.
    :param float perc_score: Minimum percentage of top score for a hit.
    :param float expect: Maximum e-value for a hit.
    :param float perc_length: Minimum percent of the longest hit by query coverage for a hit.
    :param int perc_ident: Minimum percent identity of a hit.
    :param float perc_span: Minimum percentage of the longest length-span of an HSP within a hit.
    :param int min_hsps: Minimum number of hsps needed per hit.
    :param bool hsp_cumu_score: Use cumulative HSP scores rather than individual scores. (LEAVE AS DEFAULT)
    :param str seq_method: How should id_ranker compile the sequence ranges? [Default: 'whole']
    :param list align_scorelist: (LEAVE AS DEFAULT)
    :param int indent: Indent level for pretty print. [Default: 0]
    :param int verbose: Level of verbose output? [Default: 1]
    :param bool samestrand: Should the function filter hits with HSPs on different strands? [Default:True]
    :return:
    """
    id_list = []
    if verbose:
        print('Beginning ID_Ranker...', indent=indent)
    if isinstance(record, SearchIO.QueryResult):
        print('SearchIO detected.', indent=indent)
        # Figure out what kind of search program was run:
        if record.program == 'blat':
            if verbose > 2:
                print('Results obtained from BLAT run.', indent=indent+1)
        elif 'blast' in record.program:
            if verbose > 2:
                print('Results obtained from BLAST run.', indent=indent+1)
        else:
            raise Exception('Sorry, your program {} is not yet implemented for RecBlast!'.format(record.program))

        # Create filter functions:
        def hit_minhsps(hit):
            return len(hit.hsps) >= min_hsps

        def hit_minscores(hit):
            return sum([hsp.score for hsp in hit.hsps]) >= int(perc_score * top_score)

        def hit_minlength(hit):
            return sum([i[-1] - i[0] for i in merge_ranges([(hsp.query_start, hsp.query_end) for hsp in hit])
                        ]) >= perc_length * top_length

        def hit_perc_id(hit):
            #return _percent_identity_searchio(hit) >= perc_ident
            return True

        def hit_hsp_span(hit):
            top_span = max([hsp.hit_span for hsp in hit])
            hit = hit.filter(lambda hsp: hsp.hit_span >= int(perc_span * top_span))
            return hit

        def flatten(list):
            return [item for sublist in list for item in sublist]

        def hit_same_strand(hit):
            x = [bla.hit_strand_all for bla in hit.hsps]
            y = all(s > 0 for s in flatten(x)) or all(s < 0 for s in flatten(x)) or \
                all(s == 0 for s in flatten(x)) or None
            return y

        # Two more functions for to make great progress
        def sort_scores(hit):
            return sum([hsp.score for hsp in hit.hsps])

        def hit_target_span(hit):
            return list(merge_ranges([(hsp.hit_start, hsp.hit_end) for hsp in hit]))

        # Get top stats:
        top_score = max([sum([hsp.score for hsp in hit.hsps]) for hit in record])
        if verbose > 1:
            print('Top score for {}:\t'.format(record.id), top_score, indent=indent)
        top_length = max([sum([i[-1] - i[0] for i in merge_ranges([(hsp.query_start, hsp.query_end)
                                                                   for hsp in hit])
                               ]) for hit in record])
        if verbose > 1:
            print('Longest hit for {}:\t'.format(record.id), top_length, indent=indent)

        if verbose >2:
            print("ALL HITS STATS:")
            print('Hit Name:\t|\t# HSPs\t|\tScore:\t|\tLength:\t|\tP.Ident\t|\thsp_span_list\t|')
            for hit in record:
                name = hit.id
                n_hsp = len(hit.hsps)
                score = sum([hsp.score for hsp in hit.hsps])
                length = sum([i[-1] - i[0] for i in merge_ranges([(hsp.query_start, hsp.query_end)
                                                                   for hsp in hit])])
                ident = _percent_identity_searchio(hit)
                span = [hsp.hit_span for hsp in hit]
                print('{HitName}\t|\t{HSP}\t|\t{Score}\t|\t{Length}\t|\t{PIdent}\t|\t{span_list}\t|'.format(HitName=name,
                                                                                                      HSP=n_hsp,
                                                                                                      Score=score,
                                                                                                      Length=length,
                                                                                                      PIdent=ident,
                                                                                                      span_list=span))
        # Execute filters:
        # HSP
        if verbose > 1:
            print('Number of hits for {}:\t'.format(record.id), len(record), indent=indent)
            print('Filtering based on min. number of HSPs...', indent=indent)
        record1 = record.hit_filter(hit_minhsps)
        if verbose > 1:
            print('Number of hits for {}:\t'.format(record1.id), len(record1), indent=indent)
        if not record1:
            raise Exception('No hits in Query Results have {} or more HSPs!'.format(min_hsps))
        # HSP.span
        if verbose > 1:
            print('Filtering out all HSPs with a span below {} of the longest HSP...'.format(perc_span), indent=indent)
        record2 = record1.hit_map(hit_hsp_span)
        if verbose > 1:
            print('Number of hits for {}:\t'.format(record2.id), len(record2), indent=indent)
        # Length
        if verbose > 1:
            print('Filtering out all hits shorter than {}...'.format(perc_length*top_length), indent=indent)
        record3 = record2.hit_filter(hit_minlength)
        if verbose > 1:
            print('Number of hits for {}:\t'.format(record3.id), len(record3), indent=indent)
        if not record3:
            raise Exception('No hits in Query Results above min_length {0}!'.format((top_length * perc_length)))
        # Score
        if verbose > 1:
            print('Filtering out all hits with scores less than {}...'.format(top_score * perc_score), indent=indent)
        record4 = record3.hit_filter(hit_minscores)
        if verbose > 1:
            print('Number of hits for {}:\t'.format(record4.id), len(record4), indent=indent)
        if not record4:
            raise Exception('No hits in Query Results above minimum score {0}!'.format((top_score * perc_score)))
        # Percent Identity
        if verbose > 1:
            print('Filtering out all hits with a percent identity below {}...'.format(perc_ident),
                  indent=indent)
        record5 = record4.hit_filter(hit_perc_id)
        if verbose > 1:
            print('Number of hits for {}:\t'.format(record5.id), len(record5), indent=indent)
        if not record5:
            raise Exception('No hits in Query Results above minimum score {0}!'.format((top_score * perc_score)))
        # If strand is set to strict, it will filter out hits that have different strandedness in HSPs
        if samestrand:
            if verbose > 1:
                print('Filtering out all hits whose HSPs are not on the same strand...', indent=indent)
            record6 = record5.hit_filter(hit_same_strand)
            if verbose > 1:
                print('Number of hits for {}:\t'.format(record6.id), len(record6), indent=indent)
            if not record6:
                raise Exception('No hits in Query Results with all HSPs on the same strand!')
        # Sorting them for good measure
        if verbose > 1:
            print('Sorting all hits by descending scores!', indent=indent)
        record6.sort(key=sort_scores, reverse=True, in_place=True)

        # Add items to id_list
        for hit in record6:
            seq_name = hit.id
            if seq_method == 'whole':
                hts = hit_target_span(hit)
                seq_range = '[:{0}-{1}]'.format(hts[0][0], hts[-1][-1])
            elif seq_method == 'strict':
                seq_range = '[:' + ''.join(['{0}-{1};'.format(i[0], i[-1])
                                            for i in hit_target_span(hit)]).rstrip(';') + ']'
            else:
                seq_range = ''
            strands = flatten([bla.hit_strand_all for bla in hit.hsps])
            if all(s > 0 for s in strands):
                seq_range += '(+)'
            elif all(s < 0 for s in strands):
                seq_range += '(-)'
            elif all(s == 0 for s in strands):
                seq_range += '(0)'
            else:
                seq_range += '(N)'
            seq_score = sum([hsp.score for hsp in hit.hsps])
            if verbose > 2:
                print("Adding hit {} to id list".format(seq_name+seq_range), indent=indent)
            id_list.append((seq_name, seq_range, seq_score))
            if method == 'best hit':
                print('Best Hit Reciprocal BLAST was selected, ending Reverse BLASTS after first annotation!',
                      indent=indent)
                break
    else:
        truthple = []
        subject_range = []
        query_start_end = []
        if verbose > 1:
            print('Sorting through alignment\'s HSPs to get top scores of all alignments...', indent=indent)
        for alignment in record.alignments:
            subject_range_hsp = []
            query_start_end_hsp = []
            hsp_scorelist = []
            if verbose > 3:
                print('Number of HSPs: ', len(alignment.hsps), indent=indent+1)
            for hsp in alignment.hsps:
                hsp_scorelist.append(hsp.score)
                subject_range_hsp.append(hsp.sbjct_start)
                subject_range_hsp.append(hsp.sbjct_end)
                query_start_end_hsp.append((hsp.query_start, hsp.query_end))
            hsp_scorelist.sort(reverse=True)
            query_start_end.append([i for i in merge_ranges(query_start_end_hsp)])
            subject_range.append((subject_range_hsp[0], subject_range_hsp[-1]))
            if verbose > 3:
                print("HSP Score List:", indent=indent+1)
                print(hsp_scorelist, indent=indent+2)
            if hsp_cumu_score:
                align_scorelist.append(sum(hsp_scorelist))
            else:
                align_scorelist.append(hsp_scorelist[0])
        if verbose > 1:
            print('Done with first-round sorting!', indent=indent)
        if verbose > 3:
            for align_index, alignment in enumerate(record.alignments):
                print(alignment.title, indent=indent)
                print("\tAlignment Score List:", indent=indent+1)
                print(align_scorelist[align_index], indent=indent+2)
                print("\tQuery_start_end:", indent=indent+1)
                print(query_start_end[align_index], indent=indent+2)
                print("\tSubject Range:", indent=indent+1)
                print(subject_range[align_index], indent=indent+2)
        if verbose > 1:
            print('Sorting through alignments to get all hits above threshold...', indent=indent)

        for align_index, alignment in enumerate(record.alignments):
            score_threshold = (perc_score * align_scorelist[0])
            length_alignment = sum([i[-1] - i[0] for i in query_start_end[align_index]])
            align_len_threshold = record.query_length * perc_length
            if hsp_cumu_score:
                hsp_scoretotal = sum([hsp.score for hsp in alignment.hsps])
                truthple.append((align_index, hsp_scoretotal >= score_threshold, hsp.expect <= expect,
                                 length_alignment >= align_len_threshold, alignment.title))
            else:
                truthple.append((align_index, hsp.score >= score_threshold, hsp.expect <= expect,
                                 length_alignment >= align_len_threshold, alignment.title))
        if verbose > 3:
            print('List of Alignments and Criteria Status:', indent=indent)
            print('i\t', 'Score\t', 'Expect\t', 'Length\t', 'Alignment.Title', indent=indent)
            for i in range(len(record.alignments)):
                print("{0}\t{1}\t{2}\t{3}\t{4}".format(truthple[i][0], truthple[i][1], truthple[i][2],
                      truthple[i][3], truthple[i][4]), indent=indent)
        for i in truthple:
            if i[1] and i[2] and i[3]:
                id_list.append((record.alignments[i[0]].title,
                                '[:{0}-{1}]'.format(subject_range[i[0]][0],
                                                    subject_range[i[0]][1]),
                                align_scorelist[i[0]]))
                if verbose > 2:
                    print("Alignment {} added to id_list!".format(i[4]), indent=indent)
            else:
                if verbose > 2:
                    print("WARNING: ALIGNMENT {} FAILED TO MEET CRITERIA!".format(i[4]), indent=indent)
                    if not i[1]:
                        print('Score was below threshold!', indent=indent)
                    if not i[2]:
                        print('Expect was below threshold!', indent=indent)
                    if not i[3]:
                        print('Length was below threshold!', indent=indent)
        sorted(id_list, reverse=True, key=itemgetter(2))
    return id_list


def biosql_get_sub_db_names(passwd, db="bioseqdb", driver="psycopg2", user="postgres", host="localhost"):
    """A convenience wrapper for getting all the sub-database names in a BioSQL-formatted database.

    :param str passwd: The password for the database.
    :param str db: The name of the database.
    :param str driver: The driver BioSQL will use to access the database.
    :param str user: The username used to access the database.
    :param str host: The host of the database.
    :return: a list of sub-database names.
    """
    from BioSQL import BioSeqDatabase
    server = BioSeqDatabase.open_database(driver=driver, user=user, passwd=passwd, host=host, db=db)
    sub_db_name_list = [i for i in server.keys()]
    return sub_db_name_list


def biosql_DBSeqRecord_to_SeqRecord(DBSeqRecord_, off=False):
    """Converts a DBSeqRecord object into a SeqRecord object.

    Motivation of this function was two-fold: first, it makes type testing simpler; and second, DBSeqRecord does
    not have a functional implementation of the translate method.
    :param DBSeqRecord DBSeqRecord_: The DBSeqRecord object to be converted.
    :param bool off: Don't actually convert the DBSeqRecord. [Default: False]
    :return:
    """
    assert isinstance(DBSeqRecord_, DBSeqRecord), 'Input must be a DBSeqRecord!'
    if off:
        return DBSeqRecord_
    else:
        return SeqRecord(seq=Seq(str(DBSeqRecord_.seq)), id=DBSeqRecord_.id, name=DBSeqRecord_.name,
                         description=DBSeqRecord_.description, dbxrefs=DBSeqRecord_.dbxrefs,
                         features=DBSeqRecord_.features, annotations=DBSeqRecord_.annotations,
                         letter_annotations=DBSeqRecord_.letter_annotations)


def format_range(seqrange, addlength, indent, verbose):
    try:
        lextend = -int(addlength[0])
        rextend = int(addlength[1])
    except Exception as err:
        print(err)
        lextend = 0
        rextend = 0
    try:
        lrange = int(seqrange[0])
        rrange = int(seqrange[1])
    except Exception as err:
        print(err)
        lrange = 0
        rrange = -1
    try:
        strand = seqrange[2]
    except KeyError:
        strand = '(0)'
    if verbose > 1:
        print('Original range: {0}-{1}{2}'.format(lrange, rrange, strand), indent=indent)
        print('Adding {0} steps to the beginning and {1} steps to the end of the sequence!'.format(lextend,
                                                                                                   rextend),
              indent=indent)
    if lrange > rrange:
        if strand == '(-)':
            strand = '(+)'
        else:
            strand = '(-)'
        lrange = seqrange[1]
        rrange = seqrange[0]
    newrange = tuple(map(lambda x, y: int(x) + y, (lrange, rrange), (lextend, rextend)))
    if verbose > 2:
        print('New range: {0}-{1}{2}'.format(lrange, rrange, strand), indent=indent)
    return (newrange[0], newrange[1], strand)


class RecBlastContainer(dict):
    """RecBlastContainer class containing all intermediary and final RecBlast outputs."""
    def __init__(self, target_species, query_record, **kwargs):
        super(dict, self).__init__()
        if isinstance(query_record, SeqIO.SeqRecord):
            self[target_species] = {query_record.id: dict(proc_id=kwargs.pop('proc_id', str()),
                                                            query_record=query_record,
                                                            query_species=kwargs.pop('query_species',str()),
                                                            forward_blast=dict(blast_results=kwargs.pop('blast_results',
                                                                                                        BioBlastRecord),
                                                                               blast_errors=kwargs.pop('blast_errors', '')),
                                                            forward_ids=dict(ids=kwargs.pop('ids', list()),
                                                                             missing_ids=kwargs.pop('missing_ids', list()),
                                                                             pretty_ids=kwargs.pop('pretty_ids', list())),
                                                            recblast_unanno=kwargs.pop('recblast_unanno', list()),
                                                            reverse_blast=dict(blast_results=kwargs.pop('blast_results',
                                                                                                        BioBlastRecord),
                                                                               blast_errors=kwargs.pop('blast_errors', '')),
                                                            reverse_ids=dict(ids=kwargs.pop('ids', list()),
                                                                             missing_ids=kwargs.pop('missing_ids', list()),
                                                                             pretty_ids=kwargs.pop('pretty_ids', list())),
                                                            recblast_results=kwargs.pop('recblast_results', list()),
                                                            output_paths=dict(forward_blast_output=kwargs.pop('forward_\
                                                                                                               blast_output',
                                                                                                              Path()),
                                                                              forward_id_score_output=kwargs.pop('forward_\
                                                                                                                  id_score_\
                                                                                                                  output',
                                                                                                                 Path()),
                                                                              recblast_output_unanno=kwargs.pop('recblast_\
                                                                                                                 output_\
                                                                                                                 unanno',
                                                                                                                Path()),
                                                                              reverse_blast_output=kwargs.pop('reverse_\
                                                                                                               blast_\
                                                                                                               output',
                                                                                                              list()),
                                                                              recblast_output=kwargs.pop('recblast_output',
                                                                                                         Path()),
                                                                              blast_nohits=kwargs.pop('blast_nohits',
                                                                                                      Path())
                                                                              )
                                                            )}
        elif isinstance(query_record, dict):
            self[target_species] = query_record
        else:
            raise AssertionError('query_record must either be a query dict from another RecBlastContainer '
                                 'or a SeqRecord Object!')

    def __getattr__(self, attr):
        return self.get(attr)

    def __setattr__(self, key, value):
        self.__setitem__(key, value)

    def __setitem__(self, key, value):
        super(RecBlastContainer, self).__setitem__(key, value)
        self.__dict__.update({key: value})

    def __delattr__(self, item):
        self.__delitem__(item)

    def __delitem__(self, key):
        super(RecBlastContainer, self).__delitem__(key)
        del self.__dict__[key]

    def __getstate__(self):
        return self

    def __setstate__(self, state):
        self.update(state)
        self.__dict__ = self

    def generate_stats(self):
        self.stats = RecBlastStats(self)

    def result_filter(self, FUN, *args, species=None, query=None, **kwargs):
        replace_internal = kwargs.pop('replace_internal', True)
        n_proc = kwargs.pop('n_proc', 1)

        if '__dict__' in self.keys():
            del self['__dict__']

        if species is None:
            rc=[]
            species = list(self.keys())
            n_species = len(species)
            if n_proc > n_species:
                n_pool = n_species
                n_sub_pool = int(n_proc/n_species)
            else:
                n_pool = n_proc
                n_sub_pool = 0
            """
            Pool = multiprocessing.Pool(n_pool)
            tmp = dict(replace_internal=False, n_proc=n_sub_pool)
            tmp.update(**kwargs)
            rc_handle = Pool.starmap_async(partial(self.result_filter, FUN), product(species, query, repeat(tmp)))
            rc = rc_handle.get()
            Pool.close()
            """

            for spec in species:
                rc.append(self.result_filter(FUN=FUN, *args, species=spec, query=query, replace_internal=False,
                                             n_proc=n_sub_pool, **kwargs))

            return sum(rc)
        elif query is None:
            rc=[]
            query = list(self[species].keys())
            n_pool=n_proc
            n_sub_pool=1
            """
            Pool = multiprocessing.Pool(n_pool)
            rc_handle = Pool.map_async(partial(self.result_filter, FUN=FUN, *args, species=species,
                                               replace_internal=False, n_proc=n_sub_pool, **kwargs), query)
            rc = rc_handle.get()
            Pool.close()
            """
            for q in query:
                rc.append(self.result_filter(FUN=FUN, *args, species=species, query=q, replace_internal=False,
                                             **kwargs))
            return sum(rc)
        else:
            if replace_internal:
                self[species][query]['recblast_results'] = list(filter(partial(FUN, *args, **kwargs),
                                                                       self[species][query]['recblast_results']))
                return self
            else:
                query_record = {}
                query_record[query] = deepcopy(self[species][query])
                summary_statistic = kwargs.pop('summary_statistic', None)
                recblast_object = kwargs.pop('recblast_object', 'recblast_results')
                assert isinstance(recblast_object, str), 'recblast_object must be a str!'
                if summary_statistic:
                    try:
                        stat = summary_statistic(self[species][query][recblast_object], *args, **kwargs)
                    except KeyError:
                        raise KeyError('Record {0}, {1} has no key {2}'.format(species, query, recblast_object))
                    query_record[query]['recblast_results'] = list(filter(partial(FUN, *args, stat=stat, **kwargs),
                                                                          self[species][query]['recblast_results']))
                else:
                    query_record[query]['recblast_results'] = list(filter(partial(FUN, *args, **kwargs),
                                                                          self[species][query]['recblast_results']))
                return RecBlastContainer(target_species=species, query_record=query_record)




    def result_map(self, func, *args, species=None, query=None, **kwargs):
        pass

    def write(self, file_loc=None, filetype='fasta', **kwargs):
        if file_loc is None:
            date_str = dt.now().strftime('%y-%m-%d_%I-%M-%p')
            file_loc = Path('./RecBlast_output/{0}/'.format(date_str)).absolute()
            try:
                file_loc.mkdir(parents=True)
            except FileExistsError:
                pass
        elif isinstance(file_loc, Path):
            pass
        elif isinstance(file_loc, str):
            file_loc = Path(file_loc)
        else:
            raise TypeError('file_loc must be either a string or a Path object!')
        if self == dict():
            print('rc_container is empty!')
            return 0
        else:
            if filetype.lower() in 'sqlite3':
                tbn = kwargs.pop('table_name', 'RecBlastOutput')
                odb = kwargs.pop('outdb', None)
                nwrite = self._write_sqlite(outdb=odb, sqlfile=file_loc, table_name=tbn, **kwargs)
            elif 'bed' in filetype.lower():
                filename = kwargs.pop('filename', 'RecBlastOutput.bed').replace(' ', '_')
                filename += '' if filename.endswith('.bed') else '.bed'
                if filetype.lower() == 'bed-min':
                    nwrite = self._write_bed(file_loc=file_loc, filename=filename, col=4, **kwargs)
                else:
                    nwrite = self._write_bed(file_loc=file_loc, filename=filename, **kwargs)
            elif filetype.lower() in 'gff3':
                nwrite = self._write_gff3(file_loc=file_loc, **kwargs)
            else:
                nwrite = self._write_files(file_loc=file_loc, filetype=filetype, **kwargs)
        return nwrite

    def _write_bed(self, file_loc, filename, **kwargs):
        col = kwargs.pop('col', 12)
        nwrite = 0
        for species in self.keys():
            bed = []
            recblast_output = file_loc.joinpath(species.replace(' ', '_')+'_'+str(filename))
            try:
                recblast_output.parent.mkdir(parents=True)
            except FileExistsError:
                pass
            print('Output Location for bed file of {0}:\t{1}'.format(species, str(recblast_output)))
            for query, record in self[species].items():
                if record['recblast_results'] == []:
                    continue
                for index, hit in enumerate(record['recblast_results']):
                    try:
                        feat = hit.features[0]
                    except IndexError:
                        feat = SeqFeature.SeqFeature()
                    try:
                        loc = feat.location
                        try:
                            start = str(loc.start)
                        except Exception:
                            start = '0'
                        try:
                            end = str(loc.end)
                        except Exception:
                            end = '0'
                        try:
                            strand = str(loc.strand)
                            if strand is None or strand == '0':
                                strand = '.'
                        except Exception:
                            strand = '.'
                        loc = (str(start), str(end), str(strand))
                    except AttributeError:
                        loc = ('0','0','.')
                    try:
                        score = str(feat.qualifiers['score'])
                    except KeyError:
                        score = '.'
                    try:
                        thickStart = str(feat.qualifiers['thickStart'])
                    except KeyError:
                        thickStart = '.'
                    try:
                        thickEnd = str(feat.qualifiers['thickEnd'])
                    except KeyError:
                        thickEnd = '.'
                    try:
                        itemRGB = str(feat.qualifiers['itemRGB'])
                    except KeyError:
                        itemRGB = '.'
                    try:
                        blockCount = str(feat.qualifiers['blockCount'])
                    except KeyError:
                        blockCount = '.'
                    try:
                        blockSizes = str(feat.qualifiers['blockSizes'])
                    except KeyError:
                        blockSizes = '.'
                    try:
                        blockStarts = str(feat.qualifiers['blockStarts'])
                    except KeyError:
                        blockStarts = '.'

                    items = [hit.name,
                             loc[0],
                             loc[1],
                             str(query+'_'+str(index)),
                             score,
                             loc[2],
                             thickStart,
                             thickEnd,
                             itemRGB,
                             blockCount,
                             blockSizes,
                             blockStarts]
                    items = [str(i) for i in items][0:col]
                    line = '\t'.join(items) + '\n'
                    bed.append(line)
                nwrite += 1
            with recblast_output.open('a+') as rc_out:
                rc_out.writelines(bed)
        return nwrite

    def _write_sqlite(self, outdb, sqlfile, table_name, **kwargs):
        nwrite = 0
        max_hit_col = kwargs.pop('max_hit_col', 1)
        col_to_make = kwargs.pop('col_to_make', 0)
        row_number = kwargs.pop('row_number', 1)
        if outdb is None:
            outdb = sqlite3.connect(str(sqlfile))
        cursor = outdb.cursor()
        # Create the table if it doesn't exist
        command = 'CREATE TABLE IF NOT EXISTS {tn} (id INT PRIMARY KEY , target_species TEXT, query_species TEXT, ' \
                  'query_name TEXT, query_record TEXT, hit_1 TEXT)'.format(tn=table_name)
        command = command.replace('-', '_')
        cursor.execute(command)
        outdb.commit()
        species = [i for i in self.keys() if i != '__dict__']
        for spc in species:
            targets = list(self[spc].keys())
            for target in targets:
                rc = self[spc][target]
                ts = spc
                qn = target
                qs = rc['query_species']
                qr = rc['query_record'].format('fasta')
                hits = rc['recblast_results']
                if hits == list():
                    print('No RecBlast hits in species {0} for sequence {1}!'.format(spc, rc['query_record'].id))
                    continue
                else:
                    try:
                        n_hits = len(hits)
                        if n_hits == 0:
                            hits = ['NA']
                        # Check if there's more hits than columns, and if so, add more columns to the table
                        elif n_hits > max_hit_col:
                            col_to_make = n_hits - max_hit_col
                            while col_to_make > 0:
                                col_to_make -= 1
                                max_hit_col += 1
                                hc_name = 'hit_{}'.format(max_hit_col)
                                try:
                                    cursor.execute("ALTER TABLE {tn} ADD COLUMN '{cn}' TEXT".format(tn=table_name,
                                                                                                    cn=hc_name))
                                except sqlite3.OperationalError as err:
                                    if 'duplicate column name' in str(err).lower():
                                        continue
                                    else:
                                        raise
                        # Construct value string to pass to db
                        colvals = "{id}, '{ts}', '{qs}', '{qn}', '{qr}'".format(id=row_number,
                                                                                ts=ts,
                                                                                qs=qs,
                                                                                qn=qn,
                                                                                qr=qr)
                        colvals += ', ' + ', '.join(["'{}'".format(hit.format('fasta')) for hit in hits])
                        colnames = 'id, target_species, query_species, query_name, query_record'
                        if n_hits > 0:
                            colnames += ', ' + ', '.join(['hit_{}'.format(n + 1) for n in range(n_hits)])
                        else:
                            colnames += ', hit_1'
                        # Add row
                        cursor.execute("INSERT OR IGNORE INTO {tn} ({cols}) VALUES ({vals})".format(tn=table_name,
                                                                                                    cols=colnames,
                                                                                                    vals=colvals))
                        nwrite += 1
                        row_number += 1
                        outdb.commit()
                    except AttributeError as err:
                        print('WARNING! Could not write output of RecBlastContainer[{0}, {1}]'.format(ts, qn))
                        print('Error:\n', err, indent=1)
                        continue
        outdb.commit()
        outdb.close()
        return nwrite

    def _write_gff3(self, file_loc, **kwargs):
        # Example row
        # seqid\tRecBlast\tduplicated_pseudogene,supported_by_sequence_similarity\tstart\tend\t0\t+\t0\t
        nwrite = 0

        return nwrite

    def _write_files(self, file_loc, filetype, **kwargs):
        """

        :param file_loc:
        :param filetype:
        :return:
        """
        nwrite = 0
        species = [i for i in self.keys() if i != '__dict__']
        for spc in species:
            targets = list(self[spc].keys())
            for target in targets:
                rc_local = self[spc][target]
                recblast_sequence = rc_local['recblast_results']
                if recblast_sequence == list():
                    print('No RecBlast hits in species {0} for sequence {1}!'.format(spc,
                                                                                     rc_local[
                                                                                         'query_record'].id))
                    continue
                else:
                    recblast_output = file_loc.joinpath(rc_local['output_paths']['recblast_output'])
                    print('Output Location:\t', str(recblast_output))
                    try:
                        recblast_output.parent.mkdir(parents=True)
                    except FileExistsError:
                        pass
                with recblast_output.open('a+') as rc_out:
                    if isinstance(recblast_sequence, list):
                        if isinstance(recblast_sequence[0], SeqRecord):
                            nwrite += SeqIO.write(recblast_sequence, rc_out, filetype)
                        else:
                            nwrite += rc_out.write('\n'.join(['>{}'.format(i) for i in recblast_sequence]))
                    elif isinstance(recblast_sequence, SeqRecord):
                        nwrite += SeqIO.write(recblast_sequence, rc_out, filetype)
                    else:
                        nwrite += rc_out.write('\n'.join(['>{}'.format(i) for i in recblast_sequence]))
        return nwrite

    def __str__(self):
        super(RecBlastContainer, self).__str__()
        strobj = ''
        i = '\t'
        for species, species_dict in self.items():
            strobj += species + ':\n'
            for query, query_dict in species_dict.items():
                strobj += 1*i + query + ":\n"
                for key, value in query_dict.items():
                    strobj += 2*i + key + ":\n"
                    if isinstance(value, dict):
                        for subkey, subvalue in value.items():
                            strobj += 3 * i + subkey + ":\n"
                            if isinstance(subvalue, dict):
                                for subsubkey, subsubvalue in value.items():
                                    strobj += 4 * i + subsubkey + ":\n"
                                    val = str(subsubvalue).replace('\n', '\n' + 5*i)
                                    strobj += 5 * i + val + "\n"
                            else:
                                val = str(subvalue).replace('\n', '\n' + 4*i)
                                strobj += 4 * i + val + "\n"
                    else:
                        val = str(value).replace('\n', '\n' + 3*i)
                        strobj += 3 * i + val + "\n"
        return strobj

    def __add__(self, other):
        assert isinstance(other, RecBlastContainer), "Cannot add a non-RecBlastContainer object to a RecBlastContainer!"
        for species in other.keys():
            if species in self.keys():
                other_sub = other[species]
                self_sub = self[species]
                for query in other_sub.keys():
                    if query in self_sub.keys():
                        pass  # Note: will never overwrite pre-existing query results for a species.
                    else:
                        self[species].update({query: other_sub[query]})
            else:
                self.update({species: other[species]})
        return self

    def __radd__(self, other):
        if isinstance(other, int):
            pass
        else:
            for species in other.keys():
                if species in self.keys():
                    other_sub = other[species]
                    self_sub = self[species]
                    for query in other_sub.keys():
                        if query in self_sub.keys():
                            pass  # Note: will never overwrite pre-existing query results for a species.
                        else:
                            self[species].update({query: other_sub[query]})
                else:
                    self.update({species: other[species]})
        return self


class RecBlastStats(object):
    def __init__(self, RecBlastRecord):
        self.stats = {}
        for species, query_record in RecBlastRecord.items():
            for query, record in query_record.items():
                self.stats[(species, query)] = {}
                self.stats[(species, query)]['query_len'] = len(record['query_record'])
        self.n_searches = sum([len(RecBlastRecord[species]) for species in RecBlastRecord.keys()])



def blast(seq_record, target_species, database, filetype="fasta", blast_type='blastn',
          local_blast=False, expect=0.005, megablast=True, use_index=False, blastoutput_custom=Path(), perc_ident=75,
          verbose=True, indent=0, n_threads=1, write=False, BLASTDB='/usr/db/blastdb/', outtype='pslx', return_raw = False,
          **blast_kwargs):
    """
    A wrapper function for BLAST and BLAT searches.
    :param seq_record: The record containing the query sequence for the search. Can be either a SeqIO.SeqRecord or
                       a string with the file loaction.
    :param str target_species: The species whose sequence database will be queried.
    :param str database: The name of the database to be used in the search.
    :param str filetype: Filetype of seq_record (if seq_record is a SeqRecord object, leave as default.
                         [default: 'fasta']
    :param str blast_type: Type of search to conduct. Can be a BLAST type (blastn, blastp, blastx, tblastn, tblastx) or
                           a BLAT type (blat, tblat). [Default: 'blastn']
    :param bool local_blast: Should the search be conducted locally or on remote servers? (BLAT searches are always
                             local.) [Default: False]
    :param float expect: Highest expect value of BLAST results to be returned. [Default: 0.005]
    :param bool megablast: Should MegaBLAST be used for nucleotide searches? [Default: True]
    :param bool use_index: Should BLAST use indexes associated with the database files? [Default: False]
    :param Path blastoutput_custom: Path() object of file location for output to be written.
    :param int perc_ident: Minimum percent identity required of results to be returned [Default: 75]
    :param bool verbose: Verbose output? [Default: True]
    :param int indent: Indent level for pretty print. [Default: 0]
    :param int n_threads: Number of threads to allocate for BLAST [Default: 1]
    :param bool write: Write output of searches? [Default: False]
    :param str BLASTDB: Path of databases for either BLAST or BLAT. [Default: '/usr/db/blastdb'
    :param str outtype: Output type. (see options for BLAST and BLAT) [Default: pslx]
    :param bool return_raw: Return raw output rather than processed BioBlastRecord? [Default: False]
    :param blast_kwargs: Additional keyword arguments to pass on to BLAST/BLAT.
    :return: blast_record, blast_err
    """
    if isinstance(seq_record, SeqIO.SeqRecord):
        pass
    else:
        seq_record = SeqIO.read(seq_record, filetype)
    args = dict()
    if verbose:
        print("Now starting BLAST...", indent=indent)

    if blast_type.lower() in ['blat', 'tblat', 'blat-transcript', 'tblat-transcript']:
        if verbose > 1:
            print('Search Type: ', blast_type, indent=indent)
        args_expanded = ['gfClient', 'localhost', str(database), '/', '/dev/stdin', '/dev/stdout']
        if blast_type.lower() == 'tblat':
            args_expanded += ['-t=dnax', '-q=prot']
        args_expanded += ['minIdentity={}'.format(perc_ident), '-out={}'.format(outtype)]
        try:
            if verbose:
                print('Running BLAT command:', indent=indent)
                print(args_expanded, indent=indent+1)
            blat = subprocess.Popen(args_expanded, stdout=subprocess.PIPE, universal_newlines=True, cwd=BLASTDB,
                                    stdin=subprocess.PIPE, stderr=subprocess.PIPE)

            blat_raw, blat_raw_err = blat.communicate(input=seq_record.format('fasta'))
            if blat_raw_err:
                raise Exception(blat_raw_err)
            head = subprocess.Popen(["head", "-n", "-1"], universal_newlines=True, stdin=subprocess.PIPE,
                                    stdout=subprocess.PIPE)

            blat_handle = head.communicate(input=blat_raw)
            if verbose > 2:
                print(blat_handle[0], indent=indent)
            if verbose:
                print('Done!', indent=indent)
            if isinstance(blat_handle, str):
                blat_result = blat_handle
                blat_err = None
            else:
                blat_result, blat_err = blat_handle
        except subprocess.CalledProcessError:
            raise
        blast_result, blast_err = blat_result, blat_err
        with StringIO(blast_result) as fin:
            try:
                if outtype == 'pslx':
                    blast_record = SearchIO.read(fin, format='blat-psl', pslx=True)
                elif outtype == 'psl':
                    blast_record = SearchIO.read(fin, format='blat-psl')
                elif outtype == 'blast8':
                    blast_record = SearchIO.read(fin, format='blast-tab')
                elif outtype == 'blast9':
                    blast_record = SearchIO.read(fin, format='blast-tab', comments=True)
                elif outtype == 'blast':
                    blast_record = SearchIO.read(fin, format='blast-xml')
                else:
                    raise Exception('Invalid out type')
            except ValueError:
                if verbose:
                    print('No Query Results were found in handle for seq_record {}!'.format(seq_record.id))
                    raise ValueError
            except Exception:
                if write:
                    try:
                        blastoutput_custom.mkdir(parents=True)
                    except FileExistsError:
                        pass
                    with blastoutput_custom.open('w') as pslx:
                        pslx.write(blast_result)
                print('Error reading BLAT results! Aborting!')
                print('Error details:\n')
                raise

    else:
        # search_rec_type = 'blast-xml'
        if local_blast:
            # build up the BLAST arguments:
            args.update({'-db': database, '-evalue': expect,
                         '-outfmt': outtype,
                         '-num_threads': n_threads})
            if blast_type == 'blastn':
                if megablast:
                    args['-task'] = 'megablast'
                if use_index:
                    args['-use_index'] = use_index
                args['-perc_identity'] = perc_ident
            args_expanded = list()
            [(args_expanded.append(j), args_expanded.append(k)) for j, k in args.items()]
            if verbose:
                print('Running Local Blast...', indent=indent)
                print('Options:', indent=indent)
                print(args_expanded, indent=indent+1)
            # TODO: expand behavior here to use other variants
            if blast_type in ["blastn", "blastp", "blastx", "tblastx", "tblastn"]:
                blast_cline = [blast_type] + args_expanded
                try:
                    blast_handle = (subprocess.check_output([str(i) for i in blast_cline],
                                                            input=seq_record.format('fasta'),
                                                            universal_newlines=True, cwd=BLASTDB))
                    if isinstance(blast_handle, str):
                        blast_result = blast_handle
                        blast_err = None
                    else:
                        blast_result, blast_err = blast_handle

                except subprocess.CalledProcessError:
                    raise
            else:
                raise Exception("Invalid blast choice!")

        else:
            args.update(dict(program=str(blast_type), database=str(database), sequence=seq_record.format('fasta'),
                             entrez_query='"{}"[ORGN]'.format(target_species), expect=expect, perc_ident=perc_ident))
            if megablast & (blast_type == 'blastn'):
                args['megablast'] = 'True'
            if blast_kwargs:
                args.update(blast_kwargs)
            if verbose:
                print('Submitting Remote BLAST! Options passed:', indent=indent)
                for k, v in args.items():
                    print('{0}\t=\t{1}'.format(k, v), indent=indent+1)
            try:
                blast_result = NCBIWWW.qblast(**args)
                blast_err = None
            except Exception as err:
                print(err)
                raise err

        print(blast_result)
        print(blast_err)
        if verbose:
            print('Done with Blast!', indent=indent)
        if return_raw:
            return blast_result, blast_err
        if isinstance(blast_result, StringIO):
            blast_record = NCBIXML.read(blast_result)
        else:
            with StringIO(''.join(blast_result)) as fin:
                try:
                    blast_record = NCBIXML.read(fin)
                except Exception as err:
                    print('Error reading Blast Results! Aborting!', indent=indent)
                    print('Error details:\n', err, indent=indent)
                    raise err

    if blast_type in ['blat', 'tblat']:
        pass
        # TODO: once I'm more familiar with SearchIO, fill in some of the unknowns like targetdb, etc
    return blast_record, blast_err


def biosql_seq_lookup_cascade(dtbase, sub_db_name, id_type, identifier, indent=0, verbose=False):
    seqrec = SeqRecord()
    try_get_id = True
    if id_type == 'scaffold':
        lookup_key = 'name'
    else:
        lookup_key = id_type
    if try_get_id:
        try:
            if verbose:
                print("\t\tNow searching database {0} for {1}: {2}".format(sub_db_name, id_type, identifier),
                      indent=indent)
            seqrec = biosql_DBSeqRecord_to_SeqRecord(dtbase.lookup(**{lookup_key: identifier}))
            if verbose:
                print('\tGot sequence for {}!'.format(identifier), indent=indent)
            try_get_id = False
        except KeyError:
            sleep(0.1)
            seqrec = biosql_DBSeqRecord_to_SeqRecord(dtbase.lookup(**{lookup_key: identifier}))
            if verbose:
                print('\tGot sequence for {}!'.format(identifier), indent=indent)
            try_get_id = False
        except IndexError as err:
            if verbose:
                print("WARNING: couldn't find {0} using given ID type... \n Full error: {1}".format(identifier, err),
                      indent=indent)

    if try_get_id:
        identifier_sans_subnumber = identifier.split('.')[0]
        if verbose:
            print('\t\tSeeing if removing any sub-numbers (acc: xxxxxx.1 for example) helps...', indent=indent)
            print('\t\tIdentifier: ', identifier_sans_subnumber, indent=indent)
        try:
            if verbose:
                print("\t\tNow searching database {0} for {1}: {2}".format(sub_db_name, id_type,
                                                                           identifier_sans_subnumber), indent=indent)
            seqrec = biosql_DBSeqRecord_to_SeqRecord(dtbase.lookup(**{lookup_key: identifier_sans_subnumber}))
            if verbose:
                print('\tGot sequence for {}!'.format(identifier), indent=indent)
            try_get_id = False
        except KeyError:
            sleep(0.1)
            seqrec = biosql_DBSeqRecord_to_SeqRecord(dtbase.lookup(**{lookup_key: identifier_sans_subnumber}))
            if verbose:
                print('\tGot sequence for {}!'.format(identifier), indent=indent)
            try_get_id = False
        except IndexError as err1:
            if verbose:
                print("WARNING: couldn't find {0} using abbreviated ID... \n Full error: {1}"
                      .format(identifier_sans_subnumber, err1), indent=indent)
    if try_get_id:
        try:
            if verbose:
                print('\t\tAttempting to search using Primary ID instead of declared type:', indent=indent)
            seqrec = biosql_DBSeqRecord_to_SeqRecord(dtbase.lookup(primary_id=identifier))
            if verbose:
                print('\tGot sequence for {}!'.format(identifier), indent=indent)
            try_get_id = False
        except KeyError:
            sleep(0.1)
            seqrec = biosql_DBSeqRecord_to_SeqRecord(dtbase.lookup(primary_id=identifier))
            if verbose:
                print('\tGot sequence for {}!'.format(identifier), indent=indent)
            try_get_id = False
        except IndexError as err2:
            if verbose:
                print("WARNING: couldn't find {0} using Primary ID... \n full error: {1}".format(identifier, err2),
                      indent=indent)
    if try_get_id:
        try:
            if verbose:
                print('\t\tAttempting to search using name instead of declared type:', indent=indent)
            seqrec = biosql_DBSeqRecord_to_SeqRecord(dtbase.lookup(name=identifier))
            if verbose:
                print('\tGot sequence for {}!'.format(identifier), indent=indent)
            try_get_id = False
        except KeyError:
            sleep(0.1)
            seqrec = biosql_DBSeqRecord_to_SeqRecord(dtbase.lookup(name=identifier))
            if verbose:
                print('\tGot sequence for {}!'.format(identifier), indent=indent)
            try_get_id = False
        except IndexError as err3:
            if verbose:
                print("WARNING: Still couldn't find {0} using name search: \n full error: {1}".format(identifier, err3),
                      indent=indent)

    if try_get_id:
        try:
            lookup_key = input('Last shot, chose an ID type: '
                               '[accession, primary_id, gi, version, display_id, name]')
            if lookup_key == 'exit':
                exit(exit())
            seqrec = biosql_DBSeqRecord_to_SeqRecord(dtbase.lookup(**{lookup_key: identifier}))
            if verbose:
                print('\tGot sequence for {}!'.format(identifier), indent=indent)
        except IndexError as err5:
            if verbose:
                print("WARNING: COULD NOT FIND SEQUENCES FOR ID:{0}: \n full error: {1}".format(identifier, err5),
                      indent=indent)
    return seqrec


class GetSeqMP(multiprocessing.Process):
    def __init__(self, task_queue, result_queue, db, host, driver, user, passwd, sub_db_name, verbose, server=None):
        multiprocessing.Process.__init__(self)
        self.task_queue = task_queue
        self.result_queue = result_queue
        self.db = db
        self.host = host
        self.driver = driver
        self.user = user
        self.password = passwd
        self.sub_db_name = sub_db_name
        self.verbose = verbose
        if server is None:
            self.server = BioSeqDatabase.open_database(driver=self.driver, user=self.user, passwd=self.password,
                                                       host=self.host, db=self.db)
        else:
            self.server = server

    def run(self):
        while True:
            next_task = self.task_queue.get()
            if next_task is None:
                if self.verbose:
                    print('\tAll GetSeq tasks in this process are complete')
                self.task_queue.task_done()
                break
            answer = next_task(server=self.server, sub_db_name=self.sub_db_name)
            self.task_queue.task_done()
            self.result_queue.put(answer)


class BioSeqLookupCascade(object):
    def __init__(self, id_type, identifier, verbose, indent):
        self.id_type = id_type
        self.identifier = identifier
        self.verbose = verbose
        self.indent = indent

    def __call__(self, sub_db_name, server):
        if self.verbose:
            print('\tFetching sequence: ', self.identifier, indent=self.indent)
        try:
            dtbase = server[sub_db_name]
        except KeyError as err:
            print('Woah! KeyError!', err, indent=self.indent)
            print('Waiting for 0.1 second and rerunning in case it was a mistake!', indent=self.indent)
            sleep(0.1)
            try:
                dtbase = server[sub_db_name]
            except KeyError:
                raise

        seqrec = biosql_seq_lookup_cascade(dtbase=dtbase, sub_db_name=sub_db_name, id_type=self.id_type,
                                           indent=self.indent, identifier=self.identifier, verbose=self.verbose)
        return self.identifier, seqrec


def biosql_get_record_mp(sub_db_name, passwd='', id_list=list(), id_type='accession', driver="psycopg2", indent=0,
                         user="postgres", host="localhost", db="bioseqdb", num_proc=2, verbose=True, server=None):
    """
    
    :param sub_db_name: 
    :param passwd: 
    :param id_list: 
    :param id_type: 
    :param driver: 
    :param indent: 
    :param user: 
    :param host: 
    :param db: 
    :param num_proc: 
    :param verbose: 
    :param server: 
    :return: 
    if __name__ == '__main__':
        biosql_get_record_mp(sub_db_name='MyoLuc2.0', passwd='',
                             id_list=['NW_005871148', 'NW_005871300', 'NW_005871148'], id_type='accession',
                             driver="psycopg2", user="postgres",
                             host="localhost", db="bioseqdb", verbose=True)
    """
    idents = multiprocessing.JoinableQueue()
    results = multiprocessing.Queue()

    # num = multiprocessing.cpu_count() * 2
    if verbose > 2:
        print('\tStarting BioSQL_get_record_mp', indent=indent)
    num_jobs = len(id_list)
    seqdict = dict()
    getseqs = [GetSeqMP(idents, results, db=db, host=host, driver=driver, user=user, passwd=passwd,
                        sub_db_name=sub_db_name, verbose=verbose, server=server) for i in range(num_proc)]
    for gs in getseqs:
        gs.start()

    for item in id_list:
        idents.put(BioSeqLookupCascade(id_type=id_type, identifier=item, verbose=verbose, indent=indent))

    for i in range(num_proc):
        idents.put(None)

    while num_jobs:
        temp = results.get()
        print(temp, indent=indent)
        temp[1].name = temp[0]
        seqdict[temp[0]] = temp[1]
        num_jobs -= 1
    if verbose:
        print('Done with Biosql_get_record_mp!', indent=indent)
        print('Closing processes!', indent=indent)
    for gs in getseqs:
        if gs.is_alive():
            gs.join()

    return seqdict


def id_search(id_rec, id_type='brute', verbose=True, indent=0):
    # Define the regex functions
    p = [re.compile('(gi)([| :_]+)(\d\d+\.?\d*)(.*)'),  # regex for gi
         re.compile('([AXNYZ][MWRCPGTZ]|ref)([| _:]+)(\d\d+\.?\d*)(.*)'),  # regex for refseq accession
         re.compile('(scaffold)([| _:]+)(\d+\.?\d*)(.*)'),    # regex for scaffolds
         re.compile('(id)([| :_]+)(\d\d+\.?\d*)(.*)'),  # regex for generic ID
         re.compile('(chr)([| :_]?)(\w*\d*[.|_\-]?\d*)\s*(.*)'),    # regex for chr
         re.compile(':(\d+)-(\d+)'),  # regex for sequence range
         re.compile('(\w+)([| :_]?)(\[?:?\d+\]?)(.*)'),     # regex for assembly
         re.compile('(\S+)(.*)'), # regex for gene symbol
         re.compile('\([-+0N]\)'),  # regex for strand
         ]

    id_list_ids = []  # Initialized list of IDs
    seq_range = {}  # Initialized dict of sequence ranges

    # Begin search:
    if verbose > 1:
        print('ID File Loaded, performing regex search for identifiers...', indent=indent)
        print('ID Specified as: ', id_type, indent=indent)
    if id_type == 'brute':
        if bool(p[1].findall(id_rec)):
            id_type = 'accession'
            if verbose > 1:
                print(p[1].findall(id_rec), indent=indent)
        elif bool(p[0].findall(id_rec)):
            id_type = 'gi'
            if verbose > 1:
                print(p[0].findall(id_rec), indent=indent)
        elif bool(p[2].findall(id_rec)):
            id_type = 'scaffold'
            if verbose > 1:
                print(p[2].findall(id_rec), indent=indent)
        elif bool(p[3].findall(id_rec)):
            id_type = 'id'
            if verbose > 1:
                print(p[3].findall(id_rec))
        elif bool(p[4].findall(id_rec)):
            id_type = 'chr'
            if verbose > 1:
                print(p[4].findall(id_rec))
        elif bool(p[6].findall(id_rec)):
            id_type = 'assembly'
            if verbose > 1:
                print(p[6].findall(id_rec))
        elif bool(p[7].findall(id_rec)):
            id_type = 'symbol'
            if verbose > 1:
                print(p[7].findall(id_rec))
        else:
            raise Exception('Couldn\'t identify the id!')
        if verbose > 1:
            print('Brute Force was set, tested strings for all pre-registered IDs. ID was selected as type ',
                  id_type, indent=indent)
    if id_type == 'gi':
        if bool(p[0].findall(id_rec)):
            found_id = True
            if verbose > 1:
                print('Successfully found GI numbers, compiling list!', indent=indent)
            item_parts = p[0].findall(id_rec)
            if verbose > 1:
                print('Item:\t', item_parts, indent=indent)
            id_list_ids.append(item_parts[0][0:3])
            if bool(p[5].findall(id_rec)):
                # Seq_range will be a list of tuples where the second element is the range, and the first
                # is the ID. This way, the function accommodates sequences with a subrange and sequences without a
                # subrange.
                sr_tuple = p[5].findall(id_rec)[0]
                if bool(p[8].findall(id_rec)[0]):
                    strand = p[8].findall(id_rec)[0]
                else:
                    strand = '(N)'
                seq_range[''.join(p[0].findall(id_rec)[0][0:3])] = (sr_tuple[0], sr_tuple[1], strand)
                if verbose > 1:
                    print('Found sequence delimiters in IDs!', indent=indent)
        else:
            found_id = False
    elif id_type == 'accession':
        if bool(p[1].findall(id_rec)):
            found_id = True
            if verbose > 1:
                print('Successfully found accession numbers, compiling list!', indent=indent)
            item_parts = p[1].findall(id_rec)
            if verbose > 1:
                print('Item:\t', item_parts, indent=indent)
            id_list_ids.append(item_parts[0][0:3])
            if bool(p[5].findall(id_rec)):
                sr_tuple = p[5].findall(id_rec)[0]
                if bool(p[8].findall(id_rec)[0]):
                    strand = p[8].findall(id_rec)[0]
                else:
                    strand = '(N)'
                seq_range[''.join(p[1].findall(id_rec)[0][0:3])] = (sr_tuple[0], sr_tuple[1], strand)
                if verbose > 1:
                    print('Found sequence delimiters in IDs!', indent=indent)
        else:
            found_id = False
    elif id_type == 'id':
        if bool(p[3].findall(id_rec)):
            found_id = True
            if verbose > 1:
                print('Successfully found ID numbers, compiling list!', indent=indent)
            item_parts = p[3].findall(id_rec)
            if verbose > 1:
                print('Item:\t', item_parts, indent=indent)
            id_list_ids.append(item_parts[0][0:3])
            if bool(p[5].findall(id_rec)):
                sr_tuple = p[5].findall(id_rec)[0]
                if bool(p[8].findall(id_rec)[0]):
                    strand = p[8].findall(id_rec)[0]
                else:
                    strand = '(N)'
                seq_range[''.join(p[3].findall(id_rec)[0][0:3])] = (sr_tuple[0], sr_tuple[1], strand)
                if verbose > 1:
                    print('Found sequence delimiters in IDs!', indent=indent)
        else:
            found_id = False
    elif id_type == 'scaffold':
        if bool(p[2].findall(id_rec)):
            found_id = True
            if verbose > 1:
                print('Successfully found ID numbers, compiling list!', indent=indent)
            item_parts = p[2].findall(id_rec)
            if verbose > 1:
                print('Item:\t', item_parts, indent=indent)
            if verbose > 2:
                print('Appending {} to id list!'.format(item_parts[0][0:3]))
            id_list_ids.append(item_parts[0][0:3])
            if bool(p[5].findall(id_rec)):
                if verbose > 1:
                    print('Found sequence delimiters in IDs!', indent=indent)
                item_id = ''.join(p[2].findall(id_rec)[0][0:3])
                sr_tuple = p[5].findall(id_rec)[0]
                if bool(p[8].findall(id_rec)[0]):
                    strand = p[8].findall(id_rec)[0]
                else:
                    strand = '(N)'
                seq_range[item_id] = (sr_tuple[0], sr_tuple[1], strand)
                if verbose > 1:
                    print('Seq range: ', seq_range[item_id], indent=indent)
        else:
            found_id = False
    elif id_type == 'chr':
        if bool(p[4].findall(id_rec)):
            found_id = True
            if verbose > 1:
                print('Successfully found ID numbers, compiling list!', indent=indent)
            item_parts = p[4].findall(id_rec)
            if verbose > 1:
                print('Item:\t', item_parts, indent=indent)
                print(item_parts[0])
                for i, item in enumerate(item_parts[0]):
                    print(i, item)
            id_list_ids.append(item_parts[0][0:3])
            if bool(p[5].findall(id_rec)):
                sr_tuple = p[5].findall(id_rec)[0]
                if bool(p[8].findall(id_rec)[0]):
                    strand = p[8].findall(id_rec)[0]
                else:
                    strand = '(N)'
                seq_range[''.join(p[4].findall(id_rec)[0][0:3])] = (sr_tuple[0], sr_tuple[1], strand)
                if verbose > 1:
                    print('Found sequence delimiters in IDs!', indent=indent)
        else:
            found_id = False
    elif id_type == 'assembly':
        if bool(p[6].findall(id_rec)):
            found_id = True
            if verbose > 1:
                print('Successfully found ID numbers, compiling list!', indent=indent)
            item_parts = p[6].findall(id_rec)
            if verbose > 1:
                print('Item:\t', item_parts, indent=indent)
            if verbose > 2:
                print('Appending {} to id list!'.format(item_parts[0][0:3]))
            id_list_ids.append(item_parts[0][0:3])
            if bool(p[5].findall(id_rec)):
                sr_tuple = p[5].findall(id_rec)[0]
                if bool(p[8].findall(id_rec)[0]):
                    strand = p[8].findall(id_rec)[0]
                else:
                    strand = '(N)'
                seq_range[''.join(p[6].findall(id_rec)[0][0:3])] = (sr_tuple[0], sr_tuple[1], strand)
                if verbose > 1:
                    print('Found sequence delimiters in IDs!', indent=indent)
        else:
            found_id = False
    elif id_type == 'symbol':
        if bool(p[7].findall(id_rec)):
            found_id = True
            if verbose > 1:
                print('Successfully found ID numbers, compiling list!', indent=indent)
            item_parts = p[7].findall(id_rec)
            if verbose > 1:
                print('Item:\t', item_parts, indent=indent)
                print(item_parts[0])
                for i, item in enumerate(item_parts[0]):
                    print(i, item)
            id_list_ids.append(item_parts[0][0:3])
            if bool(p[5].findall(id_rec)):
                sr_tuple = p[5].findall(id_rec)[0]
                if bool(p[8].findall(id_rec)[0]):
                    strand = p[8].findall(id_rec)[0]
                else:
                    strand = '(N)'
                seq_range[''.join(p[7].findall(id_rec)[0][0:3])] = (sr_tuple[0], sr_tuple[1], strand)
                if verbose > 1:
                    print('Found sequence delimiters in IDs!', indent=indent)
    else:
        found_id = False
    if found_id:
        return p, id_list_ids, seq_range, id_type
    else:
        raise Exception('ID Error', 'Could not get ID!')


class FetchSeqMP(multiprocessing.Process):
    def __init__(self, id_queue, seq_out_queue, delim, id_type, server, species, source, db, add_length, indent,
                 host, driver, version, user, passwd, email, output_type, batch_size, verbose, n_subthreads):
        multiprocessing.Process.__init__(self)
        # Queues
        self.id_queue = id_queue
        self.seq_out_queue = seq_out_queue
        # ID parsing
        self.delim = delim
        self.id_type = id_type
        # Search Items
        self.server = server  # Declaring the server here lets me pass it on to everything else by inheritance
        self.species = species
        self.source = source
        self.db = db
        self.add_length = add_length
        # SQL items
        self.host = host
        self.driver = driver
        self.version = version
        self.user = user
        self.passwd = passwd
        # Entrez items
        self.email = email
        self.batch_size = batch_size
        # Function attributes
        self.output_type = output_type
        self.verbose = verbose
        self.indent = indent
        self.n_subthreads = n_subthreads

    def run(self):
        while True:
            fs_instance = self.id_queue.get()
            if fs_instance is None:
                self.id_queue.task_done()
                print('All FetchSeqs in Queue completed!', indent=self.indent)
                break
            try:
                seq_dict, miss_items = fs_instance(passwd=self.passwd, id_type=self.id_type, driver=self.driver,
                                                   user=self.user, host=self.host, db=self.db, delim=self.delim,
                                                   server=self.server, version=self.version, add_length=self.add_length,
                                                   species=self.species, source=self.source, verbose=self.verbose,
                                                   n_threads=self.n_subthreads, indent=self.indent)
            except Exception as err:
                print('ERROR!')
                print(err)
                seq_dict = dict()
                miss_items = list()
            self.id_queue.task_done()
            self.seq_out_queue.put((seq_dict, miss_items))
        return


class FetchSeq(object):  # The meat of the script
    def __init__(self, id_rec):
        self.id_rec = id_rec

    def __call__(self, delim, species, version, source, passwd, id_type, driver, user, host, db, n_threads, server,
                 verbose, add_length, indent):

        # out_file = Path(output_name + '.' + output_type)
        if verbose > 1:
            print('Full header for Entry:', indent=indent)
            print(self.id_rec, indent=indent)

        p, id_list_ids, seq_range, id_type = id_search(self.id_rec, id_type=id_type, verbose=verbose)
        id_attr_dict = {}
        for i in id_list_ids:
            try:
                id_attr_dict[''.join(i[0:3])] =  (''.join(i[0:3]), i[2])
            except IndexError:
                id_attr_dict[''.join(i[0:3])] = (''.join(i[0:3]), 0)
                
        #print(seq_range)

        if verbose > 1:
            print('ID list: ', indent=indent)
            for index, ID_item in enumerate(id_list_ids):
                try:
                    print(index + 1, ': {0}    {1}'.format(''.join(ID_item),
                                                           '-'.join(seq_range[''.join(ID_item)])), indent=indent)
                except KeyError:
                    print(index + 1, ': {0}    {1}'.format(''.join(ID_item), '(No Range)'), indent=indent)
        # Armed with the ID list, we fetch the sequences from the appropriate source
        if source.lower() == "entrez":
            raise Exception('Not yet implemented, sorry!!!')
        elif source.lower() == "sql":
            if verbose > 1:
                print('Searching for sequences in local SQL db...', indent=indent)
            if verbose > 2:
                print('Please note the sub_databases of server:\n\t', [str(i) for i in server.keys()], indent=indent)
            if version.lower() == 'auto':
                sub_db_list = []
                sub_db_name = ''.join([i[0:3] for i in species.title().split(' ')])
                for sub_db in server.keys():
                    if sub_db_name in sub_db:
                        sub_db_list.append(sub_db)
                if len(sub_db_list) < 1:
                    raise NameError('sub_db does not exist!')
                elif len(sub_db_list) == 1:
                    sub_db_name = sub_db_list[0]
                else:
                    if verbose:
                        print('Multiple database versions found!', indent=indent)
                        print(sub_db_list, indent=indent)
                        print('Selecting highest DB', indent=indent)
                    sub_db_name = sorted(sub_db_list, reverse=True)[0]
                if verbose:
                    print('Sub-DB chosen was ', sub_db_name, indent=indent)
            else:
                sub_db_name = ''.join([i[0:3] for i in species.title().split(' ')]) + version
            id_list_search = [''.join(i[0:3]) for i in id_list_ids]
            try:
                seqdict = biosql_get_record_mp(sub_db_name=sub_db_name, passwd=passwd, id_list=id_list_search,
                                               id_type=id_type, driver=driver, user=user,
                                               host=host, db=db, num_proc=n_threads, server=server, verbose=True)
            except Exception as err:
                print('Please note the sub_databases of server:\n\t', [str(i) for i in server.keys()], indent=indent)
                raise Exception('Exception!', err)
            itemsnotfound = [''.join(x) for x in id_list_ids if ''.join(x) not in seqdict.keys()]
            if itemsnotfound is not None:
                if verbose > 1:
                    print('Some items were not found. List of items will be saved to the file '
                          'items_not_found.output', indent=indent)
                    for item in itemsnotfound:
                        print(item, indent=indent)
                        # with open(str(out_file.cwd()) + 'items_not_found.output', 'w') as missingitems:
                        #     missingitems.writelines(itemsnotfound)
            else:
                itemsnotfound = None
            keys = [k for k in seqdict.keys()]
            if verbose > 1:
                print("Sequence Dictionary keys:", indent=indent)
                print(keys, indent=indent)
            if bool(seq_range):
                seqrange_ids = [ids for ids in seq_range.keys()]
                if verbose > 1:
                    print('Sequence Range IDs:', indent=indent)
                    print(seqrange_ids, indent=indent)
                for k in keys:
                    seqdict[k].features.append(SeqFeature.SeqFeature(type='duplicate'))
                    if seqdict[k].id in seqrange_ids:
                        if verbose > 1:
                            print('For sequence {}, found a sequence range!'.format(str(seqdict[k].id)), indent=indent)
                            print('Full length of sequence: {}'.format(len(seqdict[k])), indent=indent)
                        if id_type == 'gi':
                            seq_description_full = p[0].findall(seqdict[k].description)[0]
                        elif id_type == 'accession':
                            seq_description_full = p[1].findall(seqdict[k].description)[0]
                        elif id_type == 'scaffold':
                            seq_description_full = p[2].findall(seqdict[k].description)[0]
                        elif id_type == 'id':
                            seq_description_full = p[3].findall(seqdict[k].description)[0]
                        else:
                            seq_description_full = p[5].findall(seqdict[k].description)[0]
                    else:
                        if verbose > 1:
                            print('No sequence range found, continuing...', indent=indent)
                        continue
                    seq_range[k] = format_range(seqrange=seq_range[k], addlength=add_length, indent=indent+1, verbose=verbose)
                    id_range = ':' + '-'.join([str(i) for i in seq_range[k]])
                    if verbose > 1:
                        print('Sequence range: ', seq_range, indent=indent)
                    if int(seq_range[k][0]) > int(seq_range[k][1]):
                        if verbose > 1:
                            print('Starting seq_range is larger than ending seq_range - sequence is '
                                  'in the (-) direction! Reversing...', indent=indent)
                        seqdict[k].seq = seqdict[k][
                                         int(seq_range[k][1]):int(seq_range[k][0])].seq.reverse_complement()
                        seqdict[k].features[0].location = SeqFeature.FeatureLocation(int(seq_range[k][1]),
                                                                                     int(seq_range[k][0]),
                                                                                     strand=-1)
                    else:
                        try:
                            strand = seq_range[k][2]
                            if strand == '(-)':
                                if verbose > 1:
                                    print('Sequence was labeled as being in the (-) direction! Reversing...',
                                          indent=indent)
                                seqdict[k].seq = seqdict[k][
                                                 int(seq_range[k][1]):int(seq_range[k][0])].seq.reverse_complement()
                            else:
                                seqdict[k] = seqdict[k][int(seq_range[k][0]):int(seq_range[k][1])]
                        except KeyError:
                            hasstrand = False
                            seqdict[k] = seqdict[k][int(seq_range[k][0]):int(seq_range[k][1])]
                            seqdict[k].features[0].location = SeqFeature.FeatureLocation(int(seq_range[k][0]),
                                                                                         int(seq_range[k][1]),
                                                                                         strand=1)
                        else:
                            hasstrand = True
                    if verbose > 1:
                        print('Seq_description_full: ', seq_description_full, indent=indent)
                        print('id_range: ', id_range[1:], indent=indent)
                    if hasstrand:
                        seqdict[k].description= ''.join(seq_description_full[0:3]) + id_range + strand + \
                                                     str(seq_description_full[3])
                        if strand == '(-)':
                            if int(seq_range[k][0]) > int(seq_range[k][0]):
                                seqdict[k].features[0].location = SeqFeature.FeatureLocation(int(seq_range[k][1]),
                                                                                             int(seq_range[k][0]),
                                                                                             strand=-1)
                            else:
                                seqdict[k].features[0].location = SeqFeature.FeatureLocation(int(seq_range[k][0]),
                                                                                             int(seq_range[k][1]),
                                                                                             strand=-1)
                        else:
                            seqdict[k].features[0].location = SeqFeature.FeatureLocation(int(seq_range[k][0]),
                                                                                         int(seq_range[k][1]),
                                                                                         strand=1)
                    else:
                        if int(seq_range[k][0]) > int(seq_range[k][1]):
                            seqdict[k].description = ''.join(seq_description_full[0:3]) + id_range + '(-)' + \
                                                     str(seq_description_full[3])
                        else:
                            seqdict[k].description = ''.join(seq_description_full[0:3]) + id_range + '(+)' + \
                                                     str(seq_description_full[3])
                    if verbose > 1:
                        print('Sequence Description: \n\t', seqdict[k].description, indent=indent)
                    seqdict[k].id += id_range
                    try:
                        seqdict[k].features[0].qualifiers['score'] = id_attr_dict[k][1]
                        seqdict[k].name = id_attr_dict[k][0]
                    except KeyError:
                        seqdict[k].features[0].qualifiers['score'] = 0
                        seqdict[k].name = id_attr_dict[k][0]
                    if verbose > 1:
                        print('Sequence ID: \n\t', seqdict[k].id, indent=indent)
                        if id_range:
                            print('Length of subsequence with range {0}: {1}'.format(id_range, len(seqdict[k])),
                                  indent=indent)

            if verbose > 1:
                print('Sequence Record post-processing, to be saved:', indent=indent)
                print(seqdict, indent=indent)
            if verbose > 1:
                print('Finished getting sequence!', indent=indent)
            return seqdict, itemsnotfound
        elif source == "fasta":  # Note: anecdotally, this doesn't run terribly fast - try to avoid.
            # TODO: have this work like SQL does.

            seqdict = SeqIO.index(db, source, key_function=lambda identifier: \
                                                                   p[0].search(p[2].search(identifier).group()).group())
            itemsnotfound = [''.join(x) for x in id_list_ids if ''.join(x) not in seqdict.keys()]
            if itemsnotfound is not None:
                if verbose > 1:
                    print('Some items were not found. List of items will be saved to the file '
                          'items_not_found.output', indent=indent)
                    for item in itemsnotfound:
                        print(item, indent=indent)
                        # with open(str(out_file.cwd()) + 'items_not_found.output', 'w') as missingitems:
                        #     missingitems.writelines(itemsnotfound)
            else:
                itemsnotfound = None
            keys = [k for k in seqdict.keys()]
            if verbose > 1:
                print("Sequence Dictionary keys:", indent=indent)
                print(keys, indent=indent)
            if bool(seq_range):
                seqrange_ids = [ids for ids in seq_range.keys()]
                if verbose > 1:
                    print('Sequence Range IDs:', indent=indent)
                    print(seqrange_ids, indent=indent)
                for k in keys:
                    seqdict[k].features.append(SeqFeature.SeqFeature(type='duplicate'))
                    if seqdict[k].id in seqrange_ids:
                        if verbose > 1:
                            print('For sequence {}, found a sequence range!'.format(str(seqdict[k].id)), indent=indent)
                            print('Full length of sequence: {}'.format(len(seqdict[k])), indent=indent)
                        if id_type == 'gi':
                            seq_description_full = p[0].findall(seqdict[k].description)[0]
                        elif id_type == 'accession':
                            seq_description_full = p[1].findall(seqdict[k].description)[0]
                        elif id_type == 'scaffold':
                            seq_description_full = p[2].findall(seqdict[k].description)[0]
                        elif id_type == 'id':
                            seq_description_full = p[3].findall(seqdict[k].description)[0]
                        else:
                            seq_description_full = p[5].findall(seqdict[k].description)[0]
                    else:
                        if verbose > 1:
                            print('No sequence range found, continuing...', indent=indent)
                        continue
                    seq_range[k] = format_range(seqrange=seq_range[k], addlength=add_length, indent=indent+1, verbose=verbose)
                    id_range = ':' + '-'.join([str(i) for i in seq_range[k]])
                    if verbose > 1:
                        print('Sequence range: ', seq_range, indent=indent)
                    if int(seq_range[k][0]) > int(seq_range[k][1]):
                        if verbose > 1:
                            print('Starting seq_range is larger than ending seq_range - sequence is '
                                  'in the (-) direction! Reversing...', indent=indent)
                        seqdict[k].seq = seqdict[k][
                                         int(seq_range[k][1]):int(seq_range[k][0])].seq.reverse_complement()
                        seqdict[k].features[0].location = SeqFeature.FeatureLocation(int(seq_range[k][1]),
                                                                                     int(seq_range[k][0]),
                                                                                     strand=-1)
                    else:
                        try:
                            strand = seq_range[k][2]
                            if strand == '(-)':
                                if verbose > 1:
                                    print('Sequence was labeled as being in the (-) direction! Reversing...',
                                          indent=indent)
                                seqdict[k].seq = seqdict[k][
                                                 int(seq_range[k][1]):int(seq_range[k][0])].seq.reverse_complement()
                            else:
                                seqdict[k] = seqdict[k][int(seq_range[k][0]):int(seq_range[k][1])]
                        except KeyError:
                            hasstrand = False
                            seqdict[k] = seqdict[k][int(seq_range[k][0]):int(seq_range[k][1])]
                            seqdict[k].features[0].location = SeqFeature.FeatureLocation(int(seq_range[k][0]),
                                                                                         int(seq_range[k][1]),
                                                                                         strand=1)
                        else:
                            hasstrand = True
                    if verbose > 1:
                        print('Seq_description_full: ', seq_description_full, indent=indent)
                        print('id_range: ', id_range[1:], indent=indent)
                    if hasstrand:
                        seqdict[k].description = ''.join(seq_description_full[0:3]) + id_range + strand + \
                                                 str(seq_description_full[3])
                        if strand == '(-)':
                            if int(seq_range[k][0]) > int(seq_range[k][0]):
                                seqdict[k].features[0].location = SeqFeature.FeatureLocation(int(seq_range[k][1]),
                                                                                             int(seq_range[k][0]),
                                                                                             strand=-1)
                            else:
                                seqdict[k].features[0].location = SeqFeature.FeatureLocation(int(seq_range[k][0]),
                                                                                             int(seq_range[k][1]),
                                                                                             strand=-1)
                        else:
                            seqdict[k].features[0].location = SeqFeature.FeatureLocation(int(seq_range[k][0]),
                                                                                         int(seq_range[k][1]),
                                                                                         strand=1)
                    else:
                        if int(seq_range[k][0]) > int(seq_range[k][1]):
                            seqdict[k].description = ''.join(seq_description_full[0:3]) + id_range + '(-)' + \
                                                     str(seq_description_full[3])
                        else:
                            seqdict[k].description = ''.join(seq_description_full[0:3]) + id_range + '(+)' + \
                                                     str(seq_description_full[3])
                    if verbose > 1:
                        print('Sequence Description: \n\t', seqdict[k].description, indent=indent)
                    seqdict[k].id += id_range
                    if verbose > 1:
                        print('Sequence ID: \n\t', seqdict[k].id, indent=indent)
                        if id_range:
                            print('Length of subsequence with range {0}: {1}'.format(id_range, len(seqdict[k])),
                                  indent=indent)
            if verbose > 1:
                print('Sequence Record post-processing, to be saved:', indent=indent)
                print(seqdict, indent=indent)
            if verbose > 1:
                print('Finished getting sequence!', indent=indent)
            return seqdict, itemsnotfound
            # SeqIO.write([seqdict[key] for key in seqdict.keys()], str(out_file), output_type)
        elif source == "2bit":
            seqdict = {}
            itemsnotfound = []
            id_list = [''.join(i[0:3]) for i in id_list_ids]
            for id in id_list:
                was_reversed = False
                print(seq_range[id])
                try:
                    seq_range[id] = format_range(seqrange=seq_range[id], addlength=add_length, indent=indent+1,
                                                 verbose=verbose)
                    if seq_range[id][0] < seq_range[id][1]:
                        id_full = '{0}:{1}-{2}'.format(id, seq_range[id][0], seq_range[id][1])
                    elif seq_range[id][0] > seq_range[id][1]:
                        id_full = '{0}:{2}-{1}'.format(id, seq_range[id][0], seq_range[id][1])
                        was_reversed = True
                    else:
                        id_full = '{0}'.format(id)
                except KeyError:
                    print(id, "does not have a SeqRange, continuing!", indent=indent)
                    id_full = '{0}'.format(id)
                else:
                    if verbose:
                        print('Found SeqRange, full id:\t', id_full, indent=indent)
                    command = ["twoBitToFa", '{0}:{1}'.format(db, id_full), '/dev/stdout']
                    if verbose > 1:
                        print('Command:', indent=indent)
                        print(' '.join(command), indent=indent+1)
                twoBitToFa_handle = subprocess.check_output(command, universal_newlines=True, stdin=subprocess.PIPE,
                                                            stderr = subprocess.PIPE)
                if type(twoBitToFa_handle) is str:
                    seq_out = twoBitToFa_handle
                else:
                    seq_out, seq_err = twoBitToFa_handle
                    raise Exception(seq_err)

                if seq_out is not None:
                    if verbose:
                        print('Got sequence for ', id_full, indent=indent)
                    if verbose > 3:
                        print(seq_out, indent=indent+1)
                    with StringIO(seq_out) as output:
                        seqdict[id] = SeqIO.read(output, 'fasta')
                    try:
                        strand = seq_range[id][2]
                        if strand == '(-)':
                            if verbose > 1:
                                print('Sequence was labeled as being in the (-) direction! Reversing...',
                                      indent=indent)
                            seqdict[id].seq = seqdict[id].seq.reverse_complement()
                    except KeyError:
                        hasstrand = False
                        seqdict[id] = seqdict[id][int(seq_range[id][0]):int(seq_range[id][1])]
                        seqdict[id].features.append(SeqFeature.SeqFeature(type='duplicate'))
                        seqdict[id].features[0].location = SeqFeature.FeatureLocation(int(seq_range[id][0]),
                                                                                      int(seq_range[id][1]),
                                                                                      strand=1)
                    else:
                        hasstrand = True
                    seqdict[id].features.append(SeqFeature.SeqFeature(type='duplicate'))
                    if hasstrand:
                        seqdict[id].description += strand
                        if strand == '(-)':

                            if int(seq_range[id][0]) > int(seq_range[id][0]):
                                seqdict[id].features[0].location = SeqFeature.FeatureLocation(int(seq_range[id][1]),
                                                                                             int(seq_range[id][0]),
                                                                                             strand=-1)
                            else:
                                seqdict[id].features[0].location = SeqFeature.FeatureLocation(int(seq_range[id][0]),
                                                                                             int(seq_range[id][1]),
                                                                                             strand=-1)
                        else:
                            seqdict[id].features[0].location = SeqFeature.FeatureLocation(int(seq_range[id][0]),
                                                                                         int(seq_range[id][1]),
                                                                                         strand=1)
                    else:
                        if was_reversed:
                            seqdict[id].description += '(-)'
                        else:
                            seqdict[id].description += '(+)'
                else:
                    itemsnotfound.append(id)
                try:
                    seqdict[id].features[0].qualifiers['score'] = id_attr_dict[id][1]
                    seqdict[id].name = id_attr_dict[id][0]
                except KeyError:
                    seqdict[id].name = id_attr_dict[id][0]
                    seqdict[id].features[0].qualifiers['score'] = 0
            return seqdict, itemsnotfound
        else:
            raise Exception('Not a valid database source: {}'.format(source))


def fetchseqMP(ids, species, write=False, output_name='', delim='\t', id_type='brute', server=None, source="SQL",
               db="bioseqdb", host='localhost', driver='psycopg2', version='1.0', user='postgres', passwd='', email='',
               batch_size=50, output_type="fasta", verbose=1, n_threads=1, n_subthreads=1, add_length=(0, 0), indent=0):

    if isgenerator(ids):
        if verbose > 1:
            print('Received generator!', indent=indent)
    elif isinstance(ids, list):
        if verbose > 1:
            print('Received list!', indent=indent)
    else:
        if verbose > 1:
            print('Reading ID File... ', indent=indent)
        with ids.open('w') as in_handle:
            id_prelist = [line.strip() for line in in_handle]  # list of each line in the file
            print('Done!', indent=indent)
        ids = [id_item for id_item in filter(None, id_prelist) if id_item]
        if not id_prelist or id_prelist is None:
            if verbose:
                print('id_prelist is empty!', indent=indent)
            return 'None'

    if verbose > 1:
        print('Readied ids!', indent=indent)

    id_list = multiprocessing.JoinableQueue()
    results = multiprocessing.Queue()

    if server is None and 'sql' in source.lower():
        try:
            if verbose > 1:
                print('No server received, opening server...', indent=indent)
            server = BioSeqDatabase.open_database(driver=driver, user=user, passwd=passwd, host=host, db=db)
            if verbose > 1:
                print('Done!', indent=indent)
        except:
            if verbose > 1:
                print('FAILED!', indent=indent)
            raise
    elif 'sql' in source.lower():
        if verbose > 1:
            print('Received server handle:', indent=indent)
            print(server, indent=indent)
        if verbose > 2:
            print('Please note the sub_databases of server:\n\t', [str(i) for i in server.keys()], indent=indent)
    elif source.lower() in ['fasta','2bit','twobit']:
        print('Search type: ', source, indent=indent)
    if verbose > 1:
        print('Creating FecSeq Processes...', indent=indent)
    fs_instances = [FetchSeqMP(id_queue=id_list, seq_out_queue=results,
                               delim=delim, id_type=id_type, server=server, species=species, source=source, db=db,
                               host=host, driver=driver, version=version, user=user, passwd=passwd, email=email,
                               output_type=output_type, batch_size=batch_size, verbose=verbose,
                               n_subthreads=n_subthreads, add_length=add_length, indent=indent+1)
                    for i in range(n_threads)]
    if verbose > 1:
        print('Done! Starting processes...', indent=indent)
    for fs in fs_instances:
        fs.start()
    if verbose > 1:
        print('Done!', indent=indent)
        print('Assigning FetchSeq records to queue... ', indent=indent)
    for id_rec in ids:
        id_list.put(FetchSeq(id_rec=id_rec))
    for fs in fs_instances:
        id_list.put(None)
    if verbose > 1:
        print('Done!', indent=indent)
    output_dict = dict()
    missing_items_list = list()
    if verbose > 1:
        print('Getting sequences from processes... ', indent=indent)
    print('--------------', n_threads, indent=indent)
    n_jobs = len(ids)
    # for i in range(len(ids)):
    while n_jobs:
        seq, missing = results.get()
        output_dict.update(seq)
        missing_items_list.append(missing)
        n_jobs -= 1
    if verbose > 1:
        print('Done! Finished fetching sequences!', indent=indent)
        print('Closing processes!', indent=indent)
    for fs in fs_instances:
        if fs.is_alive():
            fs.join()
    if write:
        SeqIO.write([output_dict[i] for i in output_dict.keys()], output_name, output_type)
        return
    else:
        return output_dict, missing_items_list


class RecBlastMP_Thread(multiprocessing.Process):
    """
    RecBlast_MP_Thread_Handle is the first branch to be made. It will perform the actual RecBlast.
    """

    def __init__(self, proc_id, rb_queue, rb_results_queue, fw_blast_db, infile_type, output_type, BLASTDB,
                 query_species, blast_type_1, blast_type_2, local_blast_1, local_blast_2, rv_blast_db, expect,
                 perc_score, perc_span, outfolder, indent, reciprocal_method, hit_name_only,
                 perc_ident, perc_length, megablast, email, id_type, fw_source, fw_id_db, fetch_batch_size, passwd,
                 host, user, driver, fw_id_db_version, verbose, n_threads, fw_blast_kwargs, rv_blast_kwargs,
                 write_intermediates):
        multiprocessing.Process.__init__(self)
        self.name = proc_id
        self.outfolder = outfolder
        self.rb_queue = rb_queue
        self.rb_results_queue = rb_results_queue
        self.fw_blast_db = fw_blast_db
        self.infile_type = infile_type
        self.output_type = output_type
        self.query_species = query_species
        self.blast_type_1 = blast_type_1
        self.blast_type_2 = blast_type_2
        self.local_blast_1 = local_blast_1
        self.local_blast_2 = local_blast_2
        self.rv_blast_db = rv_blast_db
        self.expect = expect
        self.perc_score = perc_score
        self.perc_span = perc_span
        self.perc_ident = perc_ident
        self.perc_length = perc_length
        self.megablast = megablast
        self.email = email
        self.id_type = id_type
        self.fw_source = fw_source
        self.host = host
        self.user = user
        self.driver = driver
        self.fw_id_db = fw_id_db
        self.batch_size = fetch_batch_size
        self.passwd = passwd
        self.fw_id_db_version = fw_id_db_version
        self.verbose = verbose
        self.n_threads = n_threads
        self.fw_blast_kwargs = fw_blast_kwargs
        self.rv_blast_kwargs = rv_blast_kwargs
        self.BLASTDB = BLASTDB
        self.write_intermediates = write_intermediates
        self.indent = indent
        self.reciprocal_method = reciprocal_method
        self.hit_name_only = hit_name_only


    def run(self):  # The meat of the script
        master_out = self.outfolder.joinpath('Proc-{0}.log'.format(self.name)).absolute()
        if self.verbose > 1:
            print(master_out, indent=self.indent)
        master_out_handle = master_out.open('w')
        with redirect_stdout(master_out_handle):
            while True:
                rb_instance = self.rb_queue.get()
                if rb_instance is None:
                    self.rb_queue.task_done()
                    break
                try:
                    output = rb_instance(fw_blast_db=self.fw_blast_db, BLASTDB=self.BLASTDB,
                                         infile_type=self.infile_type, output_type=self.output_type,
                                         query_species=self.query_species,
                                         blast_type_1=self.blast_type_1, blast_type_2=self.blast_type_2,
                                         local_blast_1=self.local_blast_1,
                                         local_blast_2=self.local_blast_2,
                                         rv_blast_db=self.rv_blast_db, expect=self.expect, perc_score=self.perc_score,
                                         perc_ident=self.perc_ident, perc_span=self.perc_span,
                                         perc_length=self.perc_length, megablast=self.megablast, email=self.email,
                                         id_type=self.id_type,
                                         fw_source=self.fw_source, fw_id_db=self.fw_id_db, fetch_batch_size=self.batch_size,
                                         passwd=self.passwd,
                                         fw_id_db_version=self.fw_id_db_version, verbose=self.verbose, indent=self.indent,
                                         n_threads=self.n_threads,
                                         host=self.host, reciprocal_method = self.reciprocal_method,
                                         user=self.user, driver=self.driver,
                                         fw_blast_kwargs=self.fw_blast_kwargs, rv_blast_kwargs=self.rv_blast_kwargs,
                                         proc_id=self.name, write_intermediates=self.write_intermediates,
                                         hit_name_only=self.hit_name_only
                                         )
                    self.rb_queue.task_done()
                    self.rb_results_queue.put(output)
                except Exception as err:
                    print('Woah! Something went wrong! Aborting!')
                    print('Here\'s the error:\n', err)
                    self.rb_results_queue.put(dict(error=err, proc_id=self.name))
        master_out_handle.close()
        return


class RecBlast(object):
    def __init__(self, seq_record, target_species):
        self.seq_record = seq_record
        transtab = str.maketrans('!@#$%^&*();:.,\'\"/\\?<>|[]{}-=+', '_____________________________')
        self.seq_record.id = self.seq_record.id.translate(transtab)
        self.target_species = target_species

    def __call__(self, fw_blast_db, infile_type, output_type, BLASTDB, reciprocal_method,
                 query_species, blast_type_1, blast_type_2, local_blast_1, local_blast_2, rv_blast_db, expect,
                 perc_score, indent, hit_name_only, perc_span,
                 perc_ident, perc_length, megablast, email, id_type, fw_source, fw_id_db, fetch_batch_size, passwd,
                 host, user, driver, fw_id_db_version, verbose, n_threads, fw_blast_kwargs, rv_blast_kwargs,
                 write_intermediates, proc_id):
        # Simple shunt to minimize having to rewrite code.
        target_species = self.target_species

        # Creating the RecBlast Container
        rc_container_full = RecBlastContainer(target_species=target_species, query_record=self.seq_record,
                                              query_species=query_species)
        # Shorthand reference for full container
        rc_container = rc_container_full[target_species][self.seq_record.id]
        rc_container.update(proc_id=proc_id)
        if verbose:
            print('[Proc: {0}] [Seq.Name: {1}] [Target: {2}]'.format(proc_id, self.seq_record.id, target_species),
                  indent=indent)
            print('Parameters:')
            for kwarg, value in locals().items():
                print(kwarg, ':\t', value, indent=indent+1)
        indent += 1
        if verbose > 1:
            print('Creating handles for intermediary outputs...', indent=indent)
        # Handle for output paths:
        output_paths = rc_container['output_paths']
        output_paths.update(
            forward_blast_output=Path("{0}_recblast_out".format(target_species).replace(' ', '_') + '/' +
                                      "{0}_{1}_tmp".format(blast_type_1,
                                                           self.seq_record.id
                                                           ).replace(' ', '_') + '/' +
                                      "{0}_{1}_{2}_to_{3}.xml".format(blast_type_1,
                                                                      self.seq_record.id,
                                                                      query_species,
                                                                      target_species
                                                                      ).replace(' ', '_')),
            forward_id_score_output=Path("{0}_recblast_out".format(target_species).replace(' ', '_') + '/' +
                                         "{0}_{1}_tmp".format(blast_type_1,
                                                              self.seq_record.id
                                                              ).replace(' ', '_') + '/' +
                                         "{0}_{1}_{2}_to_{3}.ID_Scores".format(blast_type_1,
                                                                               self.seq_record.id,
                                                                               query_species,
                                                                               target_species
                                                                               ).replace(' ', '_')),
            recblast_output_unanno=Path("{0}_recblast_out".format(target_species).replace(' ', '_') + '/' +
                                        "{0}_{1}_tmp".format(blast_type_1,
                                                             self.seq_record.id
                                                             ).replace(' ', '_') + '/' +
                                        "unannotated_{0}_{1}.fasta".format(blast_type_1,
                                                                           self.seq_record.id
                                                                           ).replace(' ', '_')),
            recblast_output=Path("{0}_recblast_out".format(target_species).replace(' ', '_') + '/' +
                                 "{0}_{1}.{2}".format(blast_type_1,
                                                      self.seq_record.id,
                                                      output_type
                                                      ).replace(' ', '_')),
            blast_nohits=Path("{0}_recblast_out".format(target_species).replace(' ', '_') + '/' +
                              "{0}_{1}.no-hits".format(blast_type_1,
                                                       self.seq_record.id
                                                       ).replace(' ', '_'))
        )
        if blast_type_1 in ['blat', 'tblat']:
            if verbose > 1:
                print(blast_type_1, 'was selected. Setting fw_blast_db as port:', indent=indent, end='')
            try:
                fw_blast_db = fw_blast_db[target_species]
                if verbose > 1:
                    print(str(fw_blast_db), indent=0)
            except KeyError:
                print('No port found for species ', target_species)
                return rc_container_full
            except TypeError:
                assert isinstance(fw_blast_db, Path), "Expected Path object in lieu of a subscriptable dictionary!"
        elif 'oneshot' in blast_type_1.lower():
            if verbose > 1:
                print('Since blat was selecting, searching for appropriate .2bit file for blat...', indent=indent)
            blat_2bit = get_searchdb(search_type=blast_type_1, species=target_species, db_loc=BLASTDB,
                                    verbose=verbose, indent=indent+1)
            fw_blast_db = Path(BLASTDB, blat_2bit+'.2bit').absolute()
            fw_blast_db = str(fw_blast_db) if fw_blast_db.is_file() else None
            if fw_blast_db is None:
                raise Exception('Invalid 2bit file!')
            if verbose > 1:
                print('fw_blast_db: ', fw_blast_db, indent=indent+1)

        if fw_blast_db == 'auto':
            if verbose:
                print('Blast type set to auto!', indent=indent)
            try:
                fw_blast_db = get_searchdb(search_type=blast_type_1, species=target_species, db_loc=BLASTDB,
                                           verbose=verbose, indent=indent+1)
            except Exception as Err:
                print(Err)
                return rc_container_full
        fwblast = rc_container['forward_blast']

        if fw_blast_db == 'skip':
            if verbose:
                print("Skipping Forward Blast! Using local file instead: ", indent=indent+1)
                print('Opening Forward blast output located at ',
                      str(output_paths['forward_blast_output'].absolute()), indent=indent+1)
            with output_paths['forward_blast_output'].open("r") as forward_blasthits:
                if 'blat' in blast_type_1:
                    fwblastrecord = SearchIO.read(forward_blasthits, 'blat-psl')
                else:
                    fwblastrecord = NCBIXML.read(forward_blasthits)
        elif isinstance(fw_blast_db, Path):
            with fw_blast_db.open('r') as forward_blasthits:
                if fw_blast_db.suffix == '.psl':
                    fwblastrecord = SearchIO.read(forward_blasthits, 'blat-psl')
                elif fw_blast_db.suffix == '.pslx':
                    fwblastrecord = SearchIO.read(forward_blasthits, 'blat-psl', pslx=True)
                elif fw_blast_db.suffix == '.xml':
                    fwblastrecord = SearchIO.read(forward_blasthits, 'blast-xml')
                else:
                    raise Exception('FileTypeError: Invalid File Type!')
        else:
            if verbose:
                print("Performing forward BLAST for {}... ".format(self.seq_record.id), indent=indent+1)
            try:
                if fw_blast_kwargs:
                    if verbose:
                        for key, item in fw_blast_kwargs.items():
                            print('{0}\t=\t{1}'.format(key, item), indent=indent+2)
                    try:
                        fwblastrecord, blast_err = blast(seq_record=self.seq_record, target_species=target_species,
                                                         database=fw_blast_db,
                                                         filetype=infile_type, BLASTDB=BLASTDB,
                                                         blast_type=blast_type_1, local_blast=local_blast_1, expect=expect,
                                                         megablast=megablast, n_threads=n_threads,
                                                         blastoutput_custom=output_paths['forward_blast_output'],
                                                         perc_ident=perc_ident, verbose=verbose, write=write_intermediates,
                                                         **fw_blast_kwargs)
                    except ValueError:
                        rc_container['recblast_unanno'] = [SeqRecord('')]
                        return rc_container_full
                else:
                    try:
                        fwblastrecord, blast_err = blast(seq_record=self.seq_record, target_species=target_species,
                                                         database=fw_blast_db,
                                                         filetype=infile_type, BLASTDB=BLASTDB,
                                                         blast_type=blast_type_1, local_blast=local_blast_1, expect=expect,
                                                         megablast=megablast, n_threads=n_threads,
                                                         blastoutput_custom=output_paths['forward_blast_output'],
                                                         perc_ident=perc_ident, verbose=verbose, write=write_intermediates)
                    except ValueError:
                        rc_container['recblast_unanno'] = [SeqRecord('')]
                        return rc_container_full
            except Exception as err:
                print('WARNING! UNCATCHED EXCEPTION OCCURED!')
                print(err)
                return rc_container_full
            if blast_err:
                fwblast['blast_errors'] = blast_err
                print('Received Error after FW-BLAST:', blast_err)
                return rc_container_full
        fwblast['blast_results'] = fwblastrecord
        if not fwblastrecord:
            print('Forward BLAST record was empty!')
            return rc_container_full
        if verbose:
            print('Forward blast done!', indent=indent)
            print('Culling Results based on given criteria...', indent=indent)

        try:
            f_id_ranked = id_ranker(fwblastrecord, perc_ident=perc_ident, perc_score=perc_score, perc_span=perc_span,
                                    expect=expect, perc_length=perc_length, verbose=verbose, indent=indent+1)
        except Exception as err:
            print('WARNING! UNCATCHED ERROR IN ID_RANKER!')
            print(err)
            return rc_container_full
        if verbose:
            print('Done!', indent=indent)
        f_id_out_list = ['{0}\t{1}\t{2}\n'.format(id_i[0], id_i[1], id_i[2]) for id_i in f_id_ranked]
        fw_ids = rc_container['forward_ids']
        fw_ids['ids'] = f_id_ranked
        fw_ids['pretty_ids'] = f_id_out_list

        if not f_id_out_list:
            print('Forward Blast yielded no hits, continuing to next sequence!')
            return rc_container_full
        if blast_type_1.lower() in ['blat', 'tblat'] and fw_id_db == 'auto':
            if verbose > 1:
                print('Since blat was selecting, setting fw_id_db equal to fw_blast_db', indent=indent)
            blat_2bit = get_searchdb(search_type=blast_type_1, species=target_species, db_loc=BLASTDB,
                                    verbose=verbose, indent=indent+1)
            fw_id_db = Path(BLASTDB, blat_2bit+'.2bit').absolute()
            fw_id_db = str(fw_id_db) if fw_id_db.is_file() else None
            if fw_id_db is None:
                raise Exception('Invalid 2bit file!')
            if verbose > 1:
                print(fw_id_db)
        try:
            if 'sql' in fw_source.lower():
                server = BioSeqDatabase.open_database(driver=driver, user=user, passwd=passwd,
                                                      host=host, db=fw_id_db)
            if fw_source == 'strict':
                if verbose:
                    print('Fetching Sequence from hit results!', indent=indent)
                seq_dict = {}
                missing_items = []
                assert isinstance(fwblastrecord, SearchIO.QueryResult), 'Strict fetching only implemented for SearchIO!'
                for hit in fwblastrecord:
                    seq_frags = [rec.hit for rec in hit.fragments]
                    first = True
                    for i in seq_frags:
                        if first:
                            seq = i
                            first = False
                        else:
                            seq += i
                        seq.description = ''.join([str(i) for i in seq.id])
                    seq_dict[seq.id] = seq
            else:
                server = None
                if verbose:
                    print('Beginning Fetchseq!', indent=indent)
            # Note: BioSQL is NOT thread-safe! Throws tons of errors if executed with more than one thread!
                seq_dict, missing_items = fetchseqMP(ids=f_id_out_list, species=target_species, delim='\t',
                                                     id_type='brute', server=server, source=fw_source,
                                                     db=fw_id_db, host=host, driver=driver,
                                                     version=fw_id_db_version, user=user,
                                                     passwd=passwd, email=email, batch_size=fetch_batch_size,
                                                     output_type=output_type, output_name='', write=write_intermediates,
                                                     verbose=verbose, n_threads=1, indent=indent+1,
                                                     n_subthreads=1)
            if verbose:
                print('Done with fetching!', indent=indent)
        except IndexError:
            print('WARNING! FETCHSEQ FAILED! SKIPPING THIS SEQUENCE!')
            return rc_container_full

        if missing_items == []:
            fw_ids['missing_items'] = missing_items
            print('Items were missing!', indent=indent)
            for i in missing_items:
                print(i, indent=indent+1)

        recblast_sequence = []

        if seq_dict:
            for item in f_id_ranked:
                seq_dict[item[0]].features[0].qualifiers['score'] = item[2]
                allitems = ''.join([str(i) for i in item])
                if verbose >2:
                    print(allitems)
                id_list_ids = id_search(allitems, verbose=verbose, indent=indent+1)[1]
                if verbose > 3:
                    print('f-id-ranked: ', indent=indent)
                    print(id_list_ids, indent=indent+1)
                # Todo: find a more efficient way to do this:
                if verbose:
                    print('Sorting Sequence Dict into ordered list...', indent=indent)
                for otherkey, otheritem in seq_dict.items():
                    print(id_list_ids[0], indent=indent+1)
                    print(otheritem.description, indent=indent+1)
                    if ''.join([str(i) for i in id_list_ids[0]]) in otheritem.description:
                        recblast_sequence.append(otheritem)

        else:
            err = 'No SeqDict was returned for record {0} in process {1}!'.format(''.join((self.target_species,
                                                                                           self.seq_record.id)),
                                                                                   proc_id)
            print(err)
            raise Exception(err)
        if not isinstance(recblast_sequence, list):
            recblast_sequence = [recblast_sequence]
        rc_container['recblast_unanno'] = recblast_sequence if recblast_sequence != list() else [SeqRecord('')]
        if recblast_sequence == list():
            return rc_container_full
        if verbose:
            print('Preparing for Reverse BLAST...', indent=indent)
        for index, entry_record in enumerate(recblast_sequence):
            if verbose:
                print("Entry {} in unannotated RecBlast Hits:\n".format(entry_record.id), indent=indent+1)
                for item in [entry_record.id, entry_record.description, entry_record.seq[0:10] + '...' +
                             entry_record.seq[-1]]:
                    print(item, indent=indent+2)
            output_paths['reverse_blast_output'].append(
                Path("{0}_recblast_out".format(target_species).replace(' ', '_') +
                     '/' +
                     "{0}_{1}_tmp".format(blast_type_2,
                                          self.seq_record.id
                                          ).replace(' ', '_') +
                     '/' +
                     "{0}_{1}_{3}_to_{2}_{4}.xml".format(blast_type_2,
                                                         self.seq_record.id,
                                                         query_species,
                                                         target_species,
                                                         entry_record.name
                                                         ).replace(' ', '_')
                     ))
            if verbose:
                print('Performing Reverse Blast:', indent=indent)
            if blast_type_2 in ['blat', 'tblat', 'blat-transcript', 'tblat-transcript']:
                try:
                    print(rv_blast_db)
                    rv_blast_db_i = rv_blast_db[query_species]
                except KeyError:
                    print('No port found for species ', query_species)
                    return rc_container_full
            if (rv_blast_db == 'auto') | (rv_blast_db == 'auto-transcript'):
                try:
                    rv_blast_db_i = get_searchdb(search_type=blast_type_2, species=query_species, db=BLASTDB,
                                                 verbose=verbose, indent=indent+1)
                except Exception as Err:
                    print(Err)
                    return rc_container_full
            rv_blast = rc_container['reverse_blast']
            if rv_blast_db == 'skip':
                pass
            elif rv_blast_db == 'stop':
                if verbose:
                    print('Not performing reverse blast!', indent=0)
                return rc_container_full
            elif 'oneshot' in blast_type_2.lower():
                if verbose > 1:
                    print('Since a oneshot blat was selected, searching for appropriate .2bit file for blat...',
                          indent=indent)
                rv_blat_2bit = get_searchdb(search_type=blast_type_2, species=query_species, db_loc=BLASTDB,
                                            verbose=verbose, indent=indent + 1)
                rv_blast_db_i = Path(BLASTDB, rv_blat_2bit + '.2bit').absolute()
                rv_blast_db_i = str(rv_blast_db_i) if rv_blast_db_i.is_file() else None
                if rv_blast_db_i is None:
                    raise Exception('Invalid 2bit file!')
                if verbose > 1:
                    print('rv_blast_db: ', rv_blast_db_i, indent=indent + 1)
            else:
                try:
                    if rv_blast_kwargs:
                        if verbose:
                            for key, item in rv_blast_kwargs.items():
                                print("Reverse BLAST keywords:", indent=indent+1)
                                print(key, '\t', item, indent=indent+2)
                        rvblastrecord, blast_err = blast(seq_record=entry_record, target_species=query_species,
                                                         database=rv_blast_db_i,
                                                         filetype=infile_type, BLASTDB=BLASTDB,
                                                         blast_type=blast_type_2, local_blast=local_blast_2,
                                                         expect=expect, n_threads=n_threads,
                                                         query_length=len(self.seq_record),
                                                         megablast=megablast, indent=indent+1,
                                                         blastoutput_custom=output_paths['reverse_blast_output'][index],
                                                         perc_ident=perc_ident, verbose=verbose,
                                                         write=write_intermediates,
                                                         **rv_blast_kwargs)
                    else:
                        rvblastrecord, blast_err = blast(seq_record=entry_record, target_species=query_species,
                                                         database=rv_blast_db_i,
                                                         filetype=infile_type, BLASTDB=BLASTDB,
                                                         blast_type=blast_type_2, local_blast=local_blast_2,
                                                         expect=expect, n_threads=n_threads, indent=indent+1,
                                                         query_length=len(self.seq_record),
                                                         megablast=megablast,
                                                         blastoutput_custom=output_paths['reverse_blast_output'][index],
                                                         perc_ident=perc_ident, verbose=verbose,
                                                         write=write_intermediates)
                except ValueError as err:
                    print('No reverse hits were found for seq_record {}, continuing!'.format(entry_record.name))
                    continue
                except Exception as err:
                    print('Warning! Uncaught exception!')
                    print(err)
                    return rc_container_full
                if blast_err:
                    print('Received blast error file back from Reverse BLAST!')
                    rv_blast['blast_errors'] = blast_err
                    return rc_container_full
            rv_blast['blast_results'] = rvblastrecord
            if verbose:
                print('Done with Reverse Blast!', indent=indent)
                print('Culling results using given criteria...', indent=indent)
            try:
                reverse_hits = id_ranker(rvblastrecord, perc_ident=perc_ident, perc_score=perc_score,
                                         perc_span=perc_span, perc_length=perc_length, expect=expect, verbose=verbose,
                                         indent=indent+1, method=reciprocal_method)
            except Exception as err:
                print('No Reverse Blast Hits were found for this hit!', indent=indent + 1)
                print('Continuing to next Sequence!', indent=indent + 1)
                continue
            print('Reverse BLAST hits:', indent=indent+1)
            print(reverse_hits, indent=indent+2)
            reverse_blast_annotations = []
            for anno in reverse_hits:
                try:
                    new_anno = translate_annotation(anno[0])
                except Exception:
                    new_anno = anno[0]
                finally:
                    reverse_blast_annotations.append('\t |[ {0} {1} ({2}) ]|'.format(new_anno, anno[1], anno[2]))
            if not reverse_blast_annotations:
                print('No Reverse Blast Hits were found for this hit!', indent=indent+1)
                print('Continuing to next Sequence!', indent=indent+1)
                continue
            else:
                if verbose > 1:
                    print('Done. Annotating RecBlast Hits:', indent=indent+1)
            rv_ids = rc_container['reverse_ids']
            rv_ids['ids'] = reverse_blast_annotations
            if reciprocal_method == 'best hit':
                print('Best Hit Reciprocal BLAST was selected, ending Reverse BLASTS after first annotation!',
                      indent=indent)
                entry_record.description += '|-|' + reverse_blast_annotations[0] if \
                                                    isinstance(reverse_blast_annotations, list) \
                                                    else reverse_blast_annotations
            if verbose > 3:
                print(entry_record, indent=indent+2)
        if hit_name_only:
            recblast_sequence = [entry_record.id + '\t' +  entry_record.description for entry_record in recblast_sequence]
        if not isinstance(recblast_sequence, list):
            recblast_sequence = [recblast_sequence]
        if recblast_sequence == list():
            recblast_sequence.append(SeqRecord(''))
        rc_container['recblast_results'] = recblast_sequence
        print('PID', 'SeqName', sep='\t')
        print(proc_id, self.seq_record.id, sep='\t')
        return rc_container_full
        # def run(self):


def recblastMP(seqfile, target_species, fw_blast_db='auto', rv_blast_db='auto-transcript', infile_type='fasta',
               output_type='fasta',
               host='localhost', user='postgres', driver='psycopg2',
               query_species='Homo sapiens', blast_type_1='blastn', blast_type_2='blastn', local_blast_1=False,
               local_blast_2=False,
               expect=10, perc_score=0.5, perc_span=0.1, perc_ident=0.50, perc_length=0.5, megablast=True,
               email='', run_name='default',
               id_type='brute', fw_source='sql', fw_id_db='bioseqdb', fetch_batch_size=50,
               passwd='', hit_name_only=False, min_mem=False,
               fw_id_db_version='auto', BLASTDB='/usr/db/blastdb', indent=0,
               verbose='v', max_n_processes='auto', n_threads=2, write_intermediates=False, write_final=True,
               reciprocal_method = 'best hit', fw_blast_kwargs=None, rv_blast_kwargs=None):
    """

    >>>rc_out = recblastMP('nr_Homo_sapiens_protein_GRCh38p9.fa', ['Loxodonta africana','Procavia capensis','Trichechus manatus'],
                           fw_blast_db={'Loxodonta africana':20007,
                                         'Procavia capensis':20008,
                                         'Trichechus manatus':20009}, rv_blast_db={'Homo sapiens':20005},
                           infile_type='fasta',
                           output_type='database', host='localhost', user='postgres', driver='psycopg2',
                           query_species='Homo sapiens', blast_type_1='tblat', blast_type_2='blat-transcript',
                           local_blast_1=True, local_blast_2=True,
                           expect=10, perc_score=0.009, perc_span=0.1, perc_ident=0.69, perc_length=0.001,
                           megablast=True, email='', run_name='EHM_AvA',
                           id_type='brute', fw_source='2bit', fw_id_db='auto', fetch_batch_size=50,
                           passwd='', hit_name_only=True, min_mem=True,
                           fw_id_db_version='auto', BLASTDB='/usr/db/BLAT', indent=0,
                           verbose='vvv', max_n_processes='auto', n_threads=2, write_intermediates=False,
                           write_final=True, reciprocal_method = 'best hit', fw_blast_kwargs=None, rv_blast_kwargs=None)
    :param seqfile:
    :param target_species:
    :param fw_blast_db:
    :param rv_blast_db:
    :param infile_type:
    :param output_type:
    :param host:
    :param user:
    :param driver:
    :param query_species:
    :param blast_type_1:
    :param blast_type_2:
    :param local_blast_1:
    :param local_blast_2:
    :param expect:
    :param perc_score:
    :param perc_ident:
    :param perc_length:
    :param megablast:
    :param email:
    :param id_type:
    :param fw_source:
    :param fw_id_db:
    :param fetch_batch_size:
    :param passwd:
    :param hit_name_only:
    :param min_mem:
    :param fw_id_db_version:
    :param BLASTDB:
    :param indent:
    :param verbose:
    :param max_n_processes:
    :param n_threads:
    :param write_intermediates:
    :param write_final:
    :param reciprocal_method:
    :param fw_blast_kwargs:
    :param rv_blast_kwargs:
    :return:
    """
    # Verbose-ometer
    if isinstance(verbose, str):
        verbose = verbose.lower().count('v')
    elif isinstance(verbose, int) and verbose > 0:
        pass
    else:
        raise TypeError('Verbose must be either be an integer greater than or equal to zero, or a number of v\'s equal '
                        'to the desired level of verbosity')
    if verbose:
        print('BioPython version: ', bp_version)
        print('Beginning RecBlastMP!')
    if verbose == 1:
        print('Basic verbose mode active. Will print only essential commentary.')
    elif verbose == 2:
        print('Verbose mode was set to 2. Will elaborate considerably about the script.')
    elif verbose == 3:
        print('Debugging-level verbose mode set. You will be innunadated by text. Brace yourself, and hold on to your '
              'console.')
    elif (verbose > 10) and (verbose < 20):
        print('Listen, I only have so much to say. Clearly you want me to say a lot, and I get it. This world can be '
              'chaotic, and in times of coding need, you might be seeking a little more comfort. But I\'m just a '
              'script. There\'s only so much I can say about running Reciprocal-Best-Hit-Blasts. I can run quietly, '
              'I can tell you basics, I can even tell you debug-level stuff that will drive you crazy with code - but '
              'after that, I don\'t know what else to tell you, buddy. \n\n Anyways, now that that\'s out of the way,'
              'here\'s some verbose output:\n')
    elif verbose == 50:
        print("V FOR VERBOSE: \n"
              "\"Voilà! In view, a humble vaudevillian veteran cast vicariously as both victim and villain by the \n"
              "vicissitudes of Fate. This visage, no mere veneer of vanity, is a vestige of the vox populi, now \n"
              "vacant, vanished. However, this valourous visitation of a bygone vexation stands vivified and has \n"
              "vowed to vanquish these venal and virulent vermin vanguarding vice and vouchsafing the violently \n"
              "vicious and voracious violation of volition! The only verdict is vengeance; a vendetta held as a \n"
              "votive, not in vain, for the value and veracity of such shall one day vindicate the vigilant and \n"
              "the virtuous. \n"
              "Verily, this vichyssoise of verbiage veers most verbose, so let me simply add that it's my very good \n"
              "honour to meet you and you may call me [Reciprocal-Best-Hit-BLAST Script].\" \n"
              "\t - V \n"
              "Moore, Alan, David Lloyd, Steve Whitaker, and Siobhan Dodds. V for Vendetta. New York: DC Comics, 2005.")
    if verbose > 3:
        multiprocessing.log_to_stderr(logging.DEBUG)
    #########################################################################

    # Converting perc_ident to integer because that's how these programs roll
    perc_ident = perc_ident*100
    #########################################################################

    # Multiprocessing set-up
    if verbose > 1:
        print('Creating queues... ', end='')
    rb_queue = multiprocessing.JoinableQueue()
    rb_results = multiprocessing.Queue()
    if verbose > 1:
        print('Done!')
    #########################################################################

    # Check Seqfile to make sure its real
    if verbose >= 1:
        print('Loading SeqFile records... ')
    seqfile_path = ''
    if isinstance(seqfile, str):
        seqfile_path = Path(seqfile)
    elif isinstance(seqfile, Path):
        seqfile_path = seqfile
    #########################################################################

    # Loads rec_handle when seqfile_path is defined:
    try:
        if seqfile_path == '':
            if isinstance(seqfile, dict):
                rec_handle = list()
                for key, item in seqfile.items():
                    rec_handle.append(item)
            else:
                rec_handle = list(seqfile)
        else:
            rec_handle = [i for i in SeqIO.parse(str(seqfile_path.absolute()), infile_type)]
    except FileNotFoundError:
        raise
    else:
        if verbose >= 1:
            print('Done!')
    #########################################################################

    # Calculation of processes to run
    if verbose > 1:
        print('Automatically calculating n_processes:')
    if isinstance(target_species, list):
        n_species = len(target_species)
    else:
        n_species = 1
    n_rec = len(rec_handle)
    if n_rec < 1:
        raise Exception('SeqFile was empty!')
    n_jobs = n_rec * n_species
    if isinstance(max_n_processes, str):
        max_n_processes = multiprocessing.cpu_count()
        print('CPU count: ', max_n_processes)
    if n_jobs * n_threads > max_n_processes:
        print('Optimal number of processes would be above max limit, using limit instead!')
        n_processes = int(max_n_processes / n_threads)
    else:
        n_processes = n_jobs
        if verbose:
            print('Number of processes to be made: ', n_processes)
    #########################################################################

    # Make output folder
    if run_name == 'default':
        date_str = dt.now().strftime('%y-%m-%d_%I-%M-%p')
        outfolder = Path('./RecBlast_output/{0}/'.format(date_str))
    else:
        outfolder = Path('./RecBlast_output/{0}/'.format(run_name))
    try:
        outfolder.mkdir(parents=True)
    except FileExistsError:
        pass
    #########################################################################


    # Blat stuff
    # Check if BLAT, then if so, make sure that fw_blast_db is a dictionary with a port for each species:
    if blast_type_1.lower() in ['blat', 'tblat']:
        assert isinstance(fw_blast_db, dict) or isinstance(fw_blast_db, Path), "For BLAT searches, fw_blast_db must be " \
                                                                               "a dictionary with valid species-port " \
                                                                               "key pairs; OR a Path object"
    if blast_type_2.lower() in ['blat', 'tblat', 'blat-transcript', 'tblat-transcript']:
        assert isinstance(rv_blast_db, dict), "For BLAT searches, rv_blast_db must be a dictionary with " \
                                              "valid species-port key pairs"

    if isinstance(target_species, str):
        if target_species.lower() == 'all':
            if verbose:
                print('Target species set to all. Searching server for all available species:')
            target_species = biosql_get_sub_db_names(passwd=passwd, db=fw_id_db, driver=driver, user=user,
                                                     host=host)
            if verbose:
                print(target_species, indent=1)
    server_activated = {}
    if isinstance(rv_blast_db, dict):
        if rv_blast_db[query_species] == 'auto':
            rv_server_online = False
        else:
            print('Checking status of rv_server')
            rv_server_online = blat_server('auto', order='status', host=host, port=rv_blast_db[query_species],
                                           species=query_species, BLASTDB=BLASTDB, verbose=verbose, indent=1,
                                           type=blast_type_1)
        if not rv_server_online and 'blat' in blast_type_2.lower():
            rv_blast_db[query_species] = blat_server('auto', 'start', host=host, port=30000, species=query_species,
                                                     BLASTDB=BLASTDB, verbose=verbose, indent=1, type=blast_type_2)
            server_activated[query_species] = rv_blast_db[query_species]
        elif not rv_server_online and 'tblat' in blast_type_2.lower():
            rv_blast_db[query_species] = blat_server('auto', 'start', host=host, port=30000, species=query_species,
                                                     BLASTDB=BLASTDB, verbose=verbose, indent=1, type=blast_type_2)
            server_activated[query_species] = rv_blast_db[query_species]
    if isinstance(fw_blast_db, dict):
        if isinstance(target_species, list):
            for species in target_species:
                if fw_blast_db[species] == 'auto':
                    print('fw_blast_db was set to "auto!"')
                    fw_server_online = False
                else:
                    fw_server_online = blat_server('auto','status', host=host, port=fw_blast_db[species],
                                                   species=species, BLASTDB=BLASTDB, verbose=verbose, indent=1,
                                                   type=blast_type_1)
                    print('Status of the forward server: ', fw_server_online)
                if not fw_server_online and 'blat' in blast_type_1.lower():
                    print('Forward server was not online, starting!')
                    fw_blast_db[species] = blat_server('auto','start', host=host, port=20000, species=species,
                                                       BLASTDB=BLASTDB, verbose=verbose, indent=1, type=blast_type_1)
                    server_activated[species] = fw_blast_db[species]
                elif not fw_server_online and 'tblat' in blast_type_1.lower():
                    print('Forward server was not online, starting!')
                    fw_blast_db[species] = blat_server('auto', 'start', host=host, port=20000, species=species,
                                                       BLASTDB=BLASTDB, verbose=verbose, indent=1, type=blast_type_1)
                    server_activated[species] = fw_blast_db[species]
            print(fw_blast_db)
        else:
            if fw_blast_db[target_species] == 'auto':
                print('fw_blast_db was set to "auto!"')
                fw_server_online = False
            else:
                fw_server_online = blat_server('auto', 'status', host=host, port=fw_blast_db[target_species],
                                               species=target_species, BLASTDB=BLASTDB, verbose=verbose, indent=1,
                                               type=blast_type_1)
                print('Status of the forward server: ', fw_server_online)
            if not fw_server_online and 'blat' in blast_type_1.lower():
                print('Forward server was not online, starting!')
                fw_blast_db[target_species] = blat_server('auto', 'start', host=host, port=20000,
                                                          species=target_species,
                                                          BLASTDB=BLASTDB, verbose=verbose, indent=1, type=blast_type_1)
                server_activated[target_species] = fw_blast_db[target_species]
            elif not fw_server_online and 'tblat' in blast_type_1.lower():
                print('Forward server was not online, starting!')
                fw_blast_db[target_species] = blat_server('auto', 'start', host=host, port=20000,
                                                          species=target_species,
                                                          BLASTDB=BLASTDB, verbose=verbose, indent=1, type=blast_type_1)
                server_activated[target_species] = fw_blast_db[target_species]
            print(fw_blast_db)
    #########################################################################

    # RecBlast Thread init
    if verbose >= 1:
        print('Creating RecBlast Threads... ')
    rec_blast_instances = [RecBlastMP_Thread(proc_id=str(i + 1), rb_queue=rb_queue, rb_results_queue=rb_results,
                                             fw_blast_db=fw_blast_db, BLASTDB=BLASTDB,
                                             infile_type=infile_type, output_type=output_type,
                                             query_species=query_species,
                                             blast_type_1=blast_type_1, blast_type_2=blast_type_2,
                                             local_blast_1=local_blast_1, local_blast_2=local_blast_2,
                                             rv_blast_db=rv_blast_db, expect=expect, perc_score=perc_score,
                                             perc_ident=perc_ident, perc_span=perc_span,
                                             perc_length=perc_length, megablast=megablast, email=email, id_type=id_type,
                                             fw_source=fw_source, fw_id_db=fw_id_db, fetch_batch_size=fetch_batch_size,
                                             passwd=passwd, reciprocal_method=reciprocal_method,
                                             fw_id_db_version=fw_id_db_version, verbose=verbose, n_threads=n_threads,
                                             host=host, outfolder=outfolder, indent=indent+1,
                                             hit_name_only = hit_name_only,
                                             user=user, driver=driver, write_intermediates=write_intermediates,
                                             fw_blast_kwargs=fw_blast_kwargs, rv_blast_kwargs=rv_blast_kwargs)
                           for i in range(n_processes)]
    for rcb in rec_blast_instances:
        rcb.start()
    #########################################################################

    # Progress bar
    progbar = ProgressBar(n_jobs, fmt=ProgressBar.FULL)
    if verbose:
        progbar()
    #########################################################################

    # Load species list and add RecBlasts to queue
    if isinstance(target_species, list):
        for species in target_species:
            # Initiate RecBlast Tasks
            for rec in rec_handle:
                if verbose > 2:
                    print('Sequence: ', rec.name)
                rb_queue.put(RecBlast(seq_record=rec, target_species=species))
    else:
        for rec in rec_handle:
            if verbose > 2:
                print('Sequence: ', rec.name)
            rb_queue.put(RecBlast(seq_record=rec, target_species=target_species))
    #########################################################################

    # Drop poison pills in queue
    for rcb in rec_blast_instances:
        rb_queue.put(None)
    #########################################################################

    # Collect results
    recblast_out = list()
    write_args = {}
    if output_type == 'database':
        output_type = 'sqlite3'
    if output_type in ['sqlite3']:
        sql_file = outfolder.joinpath(''.join(run_name if run_name != 'default'
                                                            else 'RecBlastOutput') + '.sqlite')
        write_args['table_name'] = ''.join(run_name if run_name != 'default'
                                           else 'RecBlastOutput').replace('-', '_')
        write_args['max_hit_col'] = 1
        write_args['col_to_make'] = 0
    write_args['row_number'] = 1
    while n_jobs:
        try:
            recblast_out.append(rb_results.get())
        except Exception as err:
            print(err)
        else:
            progbar.current += 1
            progbar()
            n_jobs -= 1
            # Get rid of any blank dicts that like to pop up inside the RBC every once in a while
            _ = recblast_out[-1].pop('__dict__', None)
            if write_final:
                try:
                    write_args['file_loc'] = sql_file if output_type in ['sqlite3'] else outfolder
                    if recblast_out[-1] != {}:
                        if 'error' in recblast_out[-1].keys():
                            print('RBC with an error was returned!')
                            print('Proc_ID: ', recblast_out[-1]['proc_id'])
                            print(recblast_out[-1]['error'])
                        else:
                            write_args['row_number'] += recblast_out[-1].write(filetype=output_type, **write_args)
                    else:
                        del recblast_out[-1]
                except Exception as err:
                    print('WARNING! Could not write output of RecBlast #{}'.format(len(recblast_out)))
                    print(err)
            if min_mem:
                recblast_out[-1] = 1
    """if output_type in ['database', 'sqlite3']:
        outdb.commit()
        outdb.close()"""
    #########################################################################

    # Kill any living threads at this point, for good measure
    for rcb in rec_blast_instances:
        if rcb.is_alive():
            rcb.join()
    #########################################################################

    """if server_activated != {}:
        if verbose:
            print('List of servers opened for search:', indent=indent)
        for species, port in server_activated.items():
            print(species, '\t', port, indent=indent+1, end='')
            try:
                blat_server('auto', 'stop', host=host, port=port, species=species,
                            BLASTDB=BLASTDB, verbose=verbose, indent=1, type=blast_type_1)
            except Exception as err:
                print('Couldn\'t stop server! Exception:')
                print(err, indent=indent+1)
            else:
                print('\tStopped!')
    """
    # Combine all the rc_out records into one big one and return it:
    recblast_out = sum([rc for rc in recblast_out if rc != {}])
    return recblast_out
    #########################################################################



def recblast_write(rc_container, verbose=1, outfolder=None):
    if not isinstance(outfolder, Path):
        date_str = dt.now().strftime('%y-%m-%d_%I-%M-%p')
        outfolder = Path('./RecBlast_output/{0}/'.format(date_str)).absolute()
        try:
            outfolder.mkdir(parents=True)
        except FileExistsError:
            pass
    if isinstance(rc_container, list):
        for rc in rc_container:
            recblast_write(rc)
    else:
        if rc_container == dict():
            print('rc_container is empty!')
            return
        else:
            species = [i for i in rc_container.keys() if i != '__dict__']
            for spc in species:
                targets = list(rc_container[spc].keys())
                for target in targets:
                    rc_local = rc_container[spc][target]
                    recblast_sequence = rc_local['recblast_results']
                    if recblast_sequence == list():
                        print('No RecBlast hits in species {0} for sequence {1}!'.format(spc,
                                                                                         rc_local['query_record'].name))
                        continue
                    else:
                        recblast_output = outfolder.joinpath(rc_local['output_paths']['recblast_output'])
                        print('Output Location:\t', str(recblast_output))
                        try:
                            recblast_output.parent.mkdir(parents=True)
                            if verbose > 1:
                                print('\t\tCreated directory \"{}\"!'.format(str(recblast_output.parent)))
                        except FileExistsError:
                            if verbose > 1:
                                print('\t\tDirectory \"{}\" already exists! Continuing!'.format(
                                    str(recblast_output.parent)))
                        with recblast_output.open('w') as rc_out:
                            SeqIO.write(recblast_sequence, rc_out, 'fasta')

"""
def blat(seq_record, target_species, database, blat_type, indent, BLASTDB, verbose, out='pslx'):
    if verbose > 1:
        print('Search Type: ', blat_type, indent=indent)
    settings = ['-t=dnax -q=prot -repMatch=2253',
                '-stepSize=5 -repMatch=1024 -minScore=0 -minIdentity=0'][blat_type == 'blat']
    args_expanded = ['blat', settings, str(database), '/dev/stdin', '/dev/stdout', '-out={}'.format(out)]
    try:
        if verbose:
            print('Running BLAT command:', indent=indent)
            print(''.join(args_expanded), indent=indent+1)
        blat = subprocess.Popen(args_expanded, stdout=subprocess.PIPE, universal_newlines=True, cwd=BLASTDB,
                                stdin=subprocess.PIPE)
        blat_raw, blat_raw_err = blat.communicate(input=seq_record.format('fasta'))
        if blat_raw_err:
            raise Exception(blat_raw_err)
        head = subprocess.Popen(["head", "-n", "-1"], universal_newlines=True, stdin=subprocess.PIPE,
                                stdout=subprocess.PIPE)
        blat_handle = head.communicate(input=blat_raw)
        if verbose > 2:
            print(blat_handle[0], indent=indent)
        if verbose:
            print('Done!', indent=indent)
        if isinstance(blat_handle, str):
            blat_result = blat_handle
            blat_err = None
        else:
            blat_result, blat_err = blat_handle
    except subprocess.CalledProcessError:
        raise
    with StringIO(blat_result) as fin:
        try:
            if out == 'pslx':
                blat_record = SearchIO.read(fin, format='blat-psl', pslx=True)
            elif out == 'psl':
                blat_record = SearchIO.read(fin, format='blat-psl')
            elif out == 'blast8':
                blat_record = SearchIO.read(fin, format='blast-tab')
            elif out == 'blast9':
                blat_record = SearchIO.read(fin, format='blast-tab', comments=True)
            elif out == 'blast':
                blat_record = SearchIO.read(fin, format='blast-xml')
            else:
                raise Exception('Invalid out type')
        except Exception:
            with Path('./{0}_{1}.pslx'.format(target_species, seq_record.name)).open('w') as pslx:
                pslx.write(blat_result)
            print('Error reading forward BLAT results! Aborting!')
            print('Error details:\n')
            raise
"""

class RecBlastControl(object):
    def __init__(self, file, header=True):
        self.options_list = []
        self.file = Path(file)
        self.header = header
        try:
            self.file.exists()
            self.file.is_file()
        except Exception:
            raise
        if self.header == True:
            skip_first_line = True
        else:
            skip_first_line = False
        with self.file.open() as ctl_file:
            tmplines = [line.rstrip() for line in ctl_file.readlines()]
        if tmplines == []:
            raise Exception('No lines in temlines!')
        for line in tmplines:
            if skip_first_line:
                skip_first_line = False
                continue
            elif line[0] == "#":
                continue
            options = [opt.split(',') if len(opt.split(',')) > 1 else opt for opt in line.split('\t')]
            self.options_list.append(OrderedDict(seqfile=options[0],
                                                 target_species=options[1],
                                                 fw_blast_db=options[2] if options[2] != '*' else 'auto',
                                                 rv_blast_db=options[3] if options[3] != '*' else 'auto',
                                                 infile_type=options[4] if options[4] != '*' else 'fasta',
                                                 output_type=options[5] if options[5] != '*' else 'fasta',
                                                 host=options[6] if options[6] != '*' else 'localhost',
                                                 user=options[7] if options[7] != '*' else 'postgres',
                                                 driver=options[8] if options[8] != '*' else 'psycopg2',
                                                 query_species=options[9] if options[9] != '*' else 'Homo sapiens',
                                                 blast_type_1=options[10] if options[10] != '*' else 'blastn',
                                                 blast_type_2=options[11] if options[11] != '*' else 'blastn',
                                                 local_blast_1=options[12] if options[12] != '*' else False,
                                                 local_blast_2=options[13] if options[13] != '*' else False,
                                                 expect=options[14] if options[14] != '*' else 10,
                                                 perc_score=options[15] if options[15] != '*' else 0.5,
                                                 perc_ident=options[16] if options[16] != '*' else 50,
                                                 perc_length=options[17] if options[17] != '*' else 0.5,
                                                 megablast=options[18] if options[18] != '*' else True,
                                                 email=options[19] if options[19] != '*' else '',
                                                 id_type=options[20] if options[20] != '*' else 'brute',
                                                 fw_source=options[21] if options[21] != '*' else 'sql',
                                                 fw_id_db=options[22] if options[22] != '*' else 'bioseqdb',
                                                 fetch_batch_size=options[23] if options[23] != '*' else 50,
                                                 passwd=options[24] if options[24] != '*' else '',
                                                 fw_id_db_version=options[25] if options[25] != '*' else 'auto',
                                                 BLASTDB=options[26] if options[26] != '*' else '/usr/db/blastdb',
                                                 indent=options[27] if options[27] != '*' else 0,
                                                 verbose=options[28] if options[28] != '*' else 'v',
                                                 max_n_processes=options[29] if options[29] != '*' else 'auto',
                                                 n_threads=options[30] if options[30] != '*' else 2,
                                                 write_intermediates=options[31] if options[31] != '*' else False,
                                                 write_final=options[32] if options[32] != '*' else True,
                                                 fw_blast_kwargs=options[33] if options[33] != '*' else None,
                                                 rv_blast_kwargs=options[34] if options[34] != '*' else None,
                                                 ))

    def __str__(self):
        tmp = ''
        for index, options in enumerate(self.options_list):
            tmp += 'Run {}:\n'.format(index+1)
            tmp += '\tKeyword\tValue\n\t'
            tmp += '_'*7 + '\t' + '_'*5 + '\n'
            for key, item in options.items():
                tmp += '\t{0}\t{1}'.format(key, item)
                tmp += '\n'
        return tmp

    def __call__(self):
        rc_all = []
        for options in self.options_list:
            try:
                rc_all += recblastMP(**options)
            except Exception as err:
                print(err)
                continue
        return rc_all