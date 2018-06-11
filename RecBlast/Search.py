import re
from math import log
from time import sleep
import subprocess
from io import StringIO
from pathlib import Path
from operator import itemgetter
from Bio import SearchIO, SeqIO
from Bio.Blast import NCBIXML, NCBIWWW
from RecBlast import print, merge_ranges
from RecBlast.WarningsExceptions import *


class Search(object):
    def __init__(self, search_type):
        self.search_type = search_type

    def __call__(self, seq_record, species, database, database_path, local,
                 indent, perc_ident, verbose, database_port=None,
                 expect=None, megablast=True, n_threads=1, write=False, filetype=None,
                 **kwargs):
        # query_length = len(seq_record)
        if isinstance(database, Path):
            return self.load(database)
        elif isinstance(database, str) and database != 'stop':
            return self.load(Path(database))
        elif database == 'stop':
            raise StopRecBlast()
        elif self.search_type in ["blastn", "blastp", "blastx", "tblastx", "tblastn"]:
            if verbose > 1:
                print(self.search_type, 'was selected.', indent=indent)
            dt = self.blast_prep(search_type=self.search_type, db_loc=database_path, database=database,
                                 species=species, verbose=verbose, indent=indent)
            return self.blast_run(seq_record=seq_record, species=species, database=dt.name, filetype=filetype,
                                  blast_type=self.search_type, local_blast=local, expect=expect, megablast=megablast,
                                  use_index=False, perc_ident=perc_ident, verbose=verbose, indent=indent,
                                  n_threads=n_threads, blastdb=database_path, outtype=5, return_raw=False,
                                  **kwargs)
        elif self.search_type in ['blat', 'tblat', 'blat-transcript', 'tblat-transcript']:
            if verbose > 1:
                print(self.search_type, 'was selected.', indent=indent)
            port = self.blat_prep(database_port=database_port, species=species, verbose=verbose, indent=indent)
            return self.blat_run(seq_record=seq_record, local=local, port=port,
                                 filetype=filetype, blat_type=self.search_type, perc_ident=perc_ident,
                                 verbose=verbose, indent=indent, blatdb=database_path, outtype='pslx')
        else:
            raise SearchEngineNotImplementedError('Invalid selection for search type!')

    @staticmethod
    def blast_run(seq_record, species, database, blast_type, filetype="fasta",
                  local_blast=False, expect=0.005, megablast=True, use_index=False, perc_ident=75,
                  verbose=True, indent=0, n_threads=1, blastdb='/usr/db/blastdb/', outtype=5,
                  return_raw=False, **kwargs):
        """A wrapper function for BLAST searches.
        :param seq_record: The record containing the query sequence for the search. Can be either a SeqIO.SeqRecord or
                           a string with the file loaction.
        :param str species: The species whose sequence database will be queried.
        :param Union[dict, str, Path] database: The name of the database to be used in the search.
        :param str blast_type: Type of BLAST search being performed
        :param str filetype: Filetype of seq_record (if seq_record is a SeqRecord object, leave as default.
                             [default: 'fasta']
        :param bool local_blast: Should the search be conducted locally or on remote servers? (BLAT searches are always
                                 local.) [Default: False]
        :param float expect: Highest expect value of BLAST results to be returned. [Default: 0.005]
        :param bool megablast: Should MegaBLAST be used for nucleotide searches? [Default: True]
        :param bool use_index: Should BLAST use indexes associated with the database files? [Default: False]
        :param int perc_ident: Minimum percent identity required of results to be returned [Default: 75]
        :param bool verbose: Verbose output? [Default: True]
        :param int indent: Indent level for pretty print. [Default: 0]
        :param int n_threads: Number of threads to allocate for BLAST [Default: 1]
        :param str blastdb: Path of databases for either BLAST or BLAT. [Default: '/usr/db/blastdb'
        :param int outtype: Output type. (see options for BLAST and BLAT) [Default: pslx]
        :param bool return_raw: Return raw output rather than processed BioBlastRecord? [Default: False]
        :param kwargs: Additional keyword arguments to pass on to BLAST/BLAT.
        :return: blast_record, blast_err
        """

        if isinstance(seq_record, SeqIO.SeqRecord):
            pass
        else:
            seq_record = SeqIO.read(seq_record, filetype)
        args = dict()
        if verbose:
            print("Now starting BLAST...", indent=indent)
        if local_blast:
            # build up the BLAST arguments:
            args.update({'-db': str(database), '-evalue': expect,
                         '-outfmt': str(outtype),
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
                print('Running BLAST locally...', indent=indent)
                print('Options:', indent=indent)
                print(args_expanded, indent=indent + 1)
            if blast_type in ["blastn", "blastp", "blastx", "tblastx", "tblastn"]:
                blast_cline = [blast_type] + args_expanded
                try:
                    blast_handle = subprocess.check_output([str(i) for i in blast_cline],
                                                           input=seq_record.format('fasta'),
                                                           universal_newlines=True, cwd=blastdb)
                    if isinstance(blast_handle, str):
                        blast_result = blast_handle
                        blast_err = None
                    else:
                        blast_result, blast_err = blast_handle

                except subprocess.CalledProcessError:
                    raise
            else:
                raise SearchError("Invalid blast choice!")
        else:
            args.update(dict(program=str(blast_type), database=str(database), sequence=seq_record.format('fasta'),
                             entrez_query='"{}"[ORGN]'.format(species), expect=expect, perc_ident=perc_ident))
            if megablast & (blast_type == 'blastn'):
                args['megablast'] = 'True'
            if kwargs:
                args.update(**kwargs)
            if verbose:
                print('Submitting Remote BLAST! Options passed:', indent=indent)
                for k, v in args.items():
                    print('{0}\t=\t{1}'.format(k, v), indent=indent + 1)
            try:
                blast_result = NCBIWWW.qblast(**args)
                blast_err = None
            except Exception as err:
                print(type(err), err)
                raise err

        if verbose:
            print('Done with Blast!', indent=indent)
        if return_raw:
            return blast_result, blast_err
        else:
            if isinstance(blast_result, StringIO):
                blast_record = NCBIXML.read(blast_result)
            else:
                try:
                    with StringIO(''.join(blast_result)) as fin:
                        blast_record = NCBIXML.read(fin)
                except Exception as err:
                    print('Error reading Blast Results! Aborting!', indent=indent)
                    print('Error details:\n', err, indent=indent)
                    raise err
            return blast_record, blast_err

    @staticmethod
    def blat_run(seq_record, port, local="localhost", filetype="fasta", blat_type='blat', perc_ident=None,
                 verbose=True, indent=0, blatdb='/usr/db/blastdb/', outtype='pslx'):
        """A wrapper function for BLAT searches.
        :param seq_record: The record containing the query sequence for the search. Can be either a SeqIO.SeqRecord or
                           a string with the file loaction.
        :param int port: Port of the gfServer to be queried
        :param str local: Host address.
        :param str filetype: Filetype of seq_record (if seq_record is a SeqRecord object, leave as default.
                             [default: 'fasta']
        :param str blat_type: Type of search to conduct. Can be a BLAST type (blastn, blastp, blastx, tblastn, tblastx) 
                              or a BLAT type (blat, tblat). [Default: 'blastn']
        :param int perc_ident: Minimum percent identity required of results to be returned [Default: 75]
        :param bool verbose: Verbose output? [Default: True]
        :param int indent: Indent level for pretty print. [Default: 0]
        :param str blatdb: Path of databases for either BLAST or BLAT. [Default: '/usr/db/blastdb'
        :param str outtype: Output type. (see options for BLAST and BLAT) [Default: pslx]
        :return: blat_record, blat_err
        """
        if isinstance(seq_record, SeqIO.SeqRecord):
            pass
        elif isinstance(seq_record, str):
            seq_record = SeqIO.read(seq_record, filetype)
        else:
            raise TypeError('seq_record was of type {}, must be either '
                            'a str with filepath or a SeqRecord object!'.format(type(seq_record)))

        if verbose:
            print("Now starting BLAT...", indent=indent)
        if verbose > 1:
            print('Search Type: ', blat_type, indent=indent)

        args_expanded = ['gfClient', local, str(port), '/', '/dev/stdin', '/dev/stdout']
        args_expanded += ['-t=dnax', '-q=prot'] if blat_type.lower() == 'tblat' else []
        args_expanded += ['minIdentity={}'.format(perc_ident if perc_ident else 0),
                          '-out={}'.format(outtype)]

        try:
            if verbose > 1:
                print('BLAT command:', indent=indent)
                print(' '.join(args_expanded), indent=indent + 1)

            blat = subprocess.Popen(args_expanded, stdout=subprocess.PIPE,
                                    universal_newlines=True, cwd=blatdb,
                                    stdin=subprocess.PIPE, stderr=subprocess.PIPE)

            blat_raw, blat_raw_err = blat.communicate(input=seq_record.format('fasta'))
            if blat_raw_err:
                raise SearchError(blat_raw_err)

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
        blat_result, blast_err = blat_result, blat_err
        blat_record = None
        with StringIO(blat_result) as fin:
            try:
                if outtype == 'pslx':
                    blat_record = SearchIO.read(fin, format='blat-psl', pslx=True)
                elif outtype == 'psl':
                    blat_record = SearchIO.read(fin, format='blat-psl')
                elif outtype == 'blast8':
                    blat_record = SearchIO.read(fin, format='blast-tab')
                elif outtype == 'blast9':
                    blat_record = SearchIO.read(fin, format='blast-tab', comments=True)
                elif outtype == 'blast':
                    blat_record = SearchIO.read(fin, format='blast-xml')
                else:
                    raise SearchError('Invalid out type')
            except ValueError:
                if verbose:
                    print('No Query Results were found in handle for seq_record {}!'.format(seq_record.id))
                raise NoHitsError('No Query Results were found in handle for seq_record {}!'.format(seq_record.id))
            except Exception as err:
                print('Error reading BLAT results! Aborting!')
                print('Error details:\n')
                raise err
        return blat_record, blat_err

    @staticmethod
    def blat_prep(database_port, species, verbose, indent):
        if isinstance(database_port, dict):
            try:
                blat_port = database_port[species]
                if verbose > 1:
                    print('Using port {0} for gfServer of species {1}.'.format(blat_port, species), indent=indent)
            except KeyError:
                raise SearchError('No 2bit found for species {}!'.format(species))
        elif isinstance(database_port, int):
            blat_port = database_port
        elif isinstance(database_port, str):
            try:
                blat_port = int(database_port)
            except ValueError:
                raise SearchError('Invalid option "{}" was passed to database_port! database_port must be '
                                  'either a dictionary of species-port pairs or an integer!'.format(database_port))
        else:
            raise SearchError('Invalid option of type "{}" was passed to database_port! database_port must be '
                              'either a dictionary of species-port pairs or an '
                              'integer!'.format(str(type(database_port))))
        return blat_port

    @staticmethod
    def blast_prep(search_type, database, species, verbose, indent, db_loc):
        if database == 'auto' or database == 'auto-transcript':
            if verbose > 1:
                print('Blast type set to auto!', indent=indent)
            try:
                blast_db = get_searchdb(search_type=search_type, species=species, db_loc=db_loc,
                                        verbose=verbose, indent=indent + 1)
            except Exception:
                raise SearchError('No BLAST database was found for species {}!'.format(species))
        elif isinstance(database, dict):
            try:
                blast_db = database[species]
                if verbose > 1:
                    print('Using {} as BLAST database!'.format(blast_db), indent=indent)
            except KeyError:
                raise SearchError('No BLAST database was found for species {}!'.format(species))
        elif isinstance(database, str) or isinstance(database, Path):
            blast_db = database
        else:
            raise SearchError('Invalid type given for database!')
        return blast_db

    @staticmethod
    def load(database):
        try:
            if database.exists() and database.is_file():
                rec = None
                with database.open('r') as forward_blasthits:
                    if database.suffix == '.psl':
                        rec = SearchIO.read(forward_blasthits, 'blat-psl')
                    elif database.suffix == '.pslx':
                        rec = SearchIO.read(forward_blasthits, 'blat-psl', pslx=True)
                    elif database.suffix == '.xml':
                        rec = SearchIO.read(forward_blasthits, 'blast-xml')
                    else:
                        raise SearchError('Database file "{}" could not be loaded - '
                                          'Must be either a PSL, PSLX, or BLAST-XML file!'.format(str(database)))
            else:
                raise FileNotFoundError()
        except FileNotFoundError:
            raise SearchError('Database file "{}" was not found!'.format(str(database)))
        return rec


def id_search(id_rec, id_type='brute', verbose=2, indent=0, custom_regex=None, regex_only=False):
    """

    EX:
    gi =
    refseq_accession = 'XP_010883249.1'
    scaffold = 'scaffold_145\t[:1033526-1034566](-)\t190
    id =
    chr = 'chrX[:3047971-3259961](-)119'
    seq_range =
    assembly1 = 'KN678312.1	[:9787-29116](+)	478'
    assembly2 = 'KN678312.1	[:9787-29116](+)	478'
    symbol = 'TP53'
    symbol = 'INS [:259-568](+) (161)'
    strand = '+'

    :param id_rec:
    :param id_type:
    :param custom_regex:
    :param regex_only:
    :param verbose:
    :param indent:
    :return:
    """
    # Define the regex functions
    p = dict(gi=re.compile('(\Agi[| _:]+[0-9.]+)'
                           '([| \t:_])?\[?(:?\d+-?\d+)?\]?([| \t:_])?(.*)'),
             accession=re.compile('(\A[AXNYZ][MWRCPGTZ][| _:]+[0-9.]+|\Aref[| _:]+[0-9.]+)'
                                  '([| \t:_])?\[?(:?\d+-?\d+)?\]?([| \t:_])?(.*)'),
             scaffold=re.compile('(\Ascaffold[| _:]+[0-9.]+)'
                                 '([| \t:_])?\[?(:?\d+-?\d+)?\]?([| \t:_])?(.*)'),
             id=re.compile('(\Aid[| _:]*[0-9.]+)'
                           '([| \t:_])?\[?(:?\d+-?\d+)?\]?([| \t:_])?(.*)'),
             chr=re.compile('(\Achr[| _:]*[A-Za-z0-9.]+)'
                            '([| \t:_])??\[?(:?\d+-?\d+)?\]?([| \t:_])?(.*)'),
             assembly=re.compile('(\A[A-Za-z]+[0-9.]+)'
                                 '([| \t:_])?\[?(:?\d+-?\d+)?\]?([| \t:_])?(.*)'),
             assembly_broad=re.compile('(\b[ALYB]+[0-9.]+)'
                                       '([| \t:_])?\[?(:?\d+-?\d+)?\]?([| \t:_])?(.*)'),
             symbol=re.compile('(\A\S+)([| \t:_])?\[?(:?\d+-?\d+)?\]?([| \t:_])?(.*)'),
             seq_range=re.compile(':?(\d+)-(\d+)'),
             strand=re.compile('(\([-+0N]\))'),
             score=re.compile('\d\d*')
             )
    if custom_regex is not None:
        p = {'custom': custom_regex}
        id_type = 'custom'

    # Begin search:
    if verbose > 1:
        print('ID Loaded, performing regex search for identifiers...', indent=indent)
        print('ID type: ', id_type, indent=indent)
    if id_type == 'brute':
        for tmp_type in ['accession', 'gi', 'scaffold', 'id', 'chr', 'assembly', 'assembly_broad', 'symbol']:
            if bool(p[tmp_type].findall(id_rec)):
                if verbose > 1:
                    print('Brute Force was set, tested strings for all pre-registered IDs.', indent=indent)
                    print('ID was selected as type {0}!'.format(tmp_type), indent=indent + 1)
                if regex_only:
                    return p[tmp_type]
                else:
                    return id_search(id_rec=id_rec, id_type=tmp_type, verbose=verbose, indent=indent)
        raise IDError('Couldn\'t identify the id type of line: {}!'.format(id_rec))
    else:
        try:
            item_parts = p[id_type].findall(id_rec)[0]
            if verbose > 1:
                print('Successfully found {0}, compiling list!'.format(id_type), indent=indent)
                print('Item:\t', '\t'.join(item_parts), indent=indent + 1)
        except IndexError:
            raise IDError('Could not identify patterns in {0} with id_type={1}, '
                          'is the id_search sequence correct?'.format(id_rec, id_type))
        try:

            item_parts = list(item_parts)
            item_parts[0] = item_parts[0] if not isinstance(item_parts[0], str) else ''.join(item_parts[0])

            if item_parts[2]:
                try:
                    sr_tuple = p['seq_range'].findall(item_parts[2])[0]
                    if verbose > 1:
                        print('Found sequence delimiters in IDs!', indent=indent)
                        print(sr_tuple, indent=indent + 1)
                except IndexError:
                    raise IDError('A positive match for a sequence range was found '
                                  '({0}), yet no hits were identified! Confirm that '
                                  'the regex is correct and try again!'.format(item_parts[2]))
            else:
                sr_tuple = (0, -1)
            if item_parts[4]:
                try:
                    strand = p['strand'].findall(item_parts[4])[0]
                except IndexError:
                    strand = '(N)'
                try:
                    score = p['score'].findall(item_parts[4])[0]
                except IndexError:
                    score = 0
            else:
                strand = '(N)'
                score = '0'
            if verbose > 1:
                if strand != '(N)':
                    print('Strand info found: {0}'.format(strand), indent=indent)
                if score != '0':
                    print('Score info found: {0}'.format(score), indent=indent)

            seq_range = (int(sr_tuple[0]), int(sr_tuple[1]), strand, int(score))
            return p, item_parts[0], seq_range, id_type

        except IndexError:
            raise IDError('Could not identify patterns in {0} with id_type={1}, '
                          'is the id_search sequence correct?'.format(id_rec, id_type))


def percent_identity_searchio(hit, is_protein=True):
    """Calculates percent identity based on entire hit. Adapted from UCSC BLAT FAQ and Biopython."""

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
    perc_ident = 100 - (millibad * 0.1)
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
        print('Search DB set to auto, choosing search_db...', indent=indent)
    species = species.replace(' ', '_')
    if verbose > 1:
        print('Search DB location set to: ', db_loc, indent=indent)
    db_type_dict = {
        'blastx': "protein",
        'blastp': "protein",
        'blastn': "genome",
        'tblastn': "genome",
        'tblastx': "genome",
        'blastn-transcript': "transcript",
        'tblastn-transcript': "transcript",
        'tblastx-transcript': "transcript",
        'blat': "blat",
        'tblat': "blat",
        'blat-transcript': 'blat-transcript',
        'tblat-transcript': 'tblat-transcript'
    }
    try:
        db_type = db_type_dict[search_type]
    except KeyError:
        print('Unable to determine search db type!', indent=indent)
        raise SearchError('Improper search type given ({})!'.format(search_type))
    if verbose > 1:
        print('DB type: ', db_type, indent=indent)
    db_path = Path(db_loc).absolute()
    if not db_path.exists():
        db_path = Path(db_loc)
    if db_path.exists() and db_path.is_dir():
        if db_type == 'blat':
            glob_path = [i for i in
                         db_path.glob('{0}*.2bit'.format(species.replace(' ', '_')))]  # Todo: generalize extension
        elif db_type in ['blat-transcript', 'tblat-transcript']:
            glob_path = [i for i in db_path.glob('{0}*transcript.2bit'.format(species.replace(' ', '_')))]
        else:
            glob_path = [i for i in db_path.glob('{0}_{1}*'.format(species.replace(' ', '_'), db_type))]
        if not glob_path:
            if verbose:
                print('No DB found! Trying again with abbreviated species name', indent=indent)
            species_abbv = ''.join([i[0:3] for i in species.title().split('_')])
            # making it insensitive to case for Glob
            species_abbv_insensitive = ''.join(['[{0}{1}]'.format(c.lower(),
                                                                  c.upper()) for c in species_abbv if c.isalpha()])
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
                search_db = sorted(glob_path, reverse=True)[0]
            else:
                search_db = glob_path
        except IndexError:
            print('WARNING: COULD NOT FIND DATABASE! ABORTING!', indent=indent)
            raise DatabaseNotFoundError('', 'No databases were found!')
    else:
        raise DatabaseNotFoundError('DB_Path {} does not exist!'.format(str(db_path)))
    if verbose:
        print('{0} DB chosen: {1}'.format(search_type, str(search_db)), indent=indent)
    return search_db


def blat_server(twobit, order='start', host='localhost', port=20000, type='blat', log='/dev/null', species=None,
                search_db_loc='/usr/db/blat', verbose=1, indent=0, try_limit=10, **kwargs):
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
    :param str search_db_loc: Path to the folder containing .2bit file.
    :param int verbose: Level of verbosity of function output. 0 suppresses all output, 3 is max verbosity.
    :param int indent: Indentation level of print output.
    :param int try_limit: Number of tries at 30-second intervals that function should probe the gfServer before timeout.
    :param kwargs: keyword arguments to be passed on to gfServer.
    :return: if order='start', returns the port of the new gfServer; if order='status', returns 0 if there was no
    connection, or 1 if the server is active and responding.
    """
    # Regular: gfServer start localhost portX -stepSize=5 -log=untrans.log database.2bit
    # Prot>DNAX:  gfServer start localhost portY -trans -mask -log=trans.log database.2bit
    gfserver_suppl_args = list()
    if twobit == 'auto' and order != 'stop':
        if verbose:
            print('2bit set to auto: searching for 2bit file for species ', species, indent=indent)
        twobit = get_searchdb(search_type='blat', species=species, db_loc=search_db_loc,
                              verbose=verbose, indent=indent + 1)
        if twobit.exists() and twobit.is_file():
            twobit = twobit.name
        else:
            raise BLATServerError('Invalid 2bit file!')
    for key, item in kwargs.items():
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
            raise BLATServerError(out)
        else:
            return 1
    elif order == 'stop':
        subprocess.check_call('gfServer stop {0} {1}'.format(str(host), str(port)), stdout=subprocess.PIPE,
                              stderr=subprocess.STDOUT, universal_newlines=True, shell=True,
                              executable='/bin/bash')
        return
    else:
        print(order)
        # Todo: make the portsniffer its own function and make sure it works properly.
        portfinder = subprocess.check_output('/home/manny/Scripts/oneshot/checkifportisopen.sh {}'.format(str(port)),
                                             universal_newlines=True, shell=True, executable='/bin/bash')
        port = portfinder.rstrip()

        gfserver_cmd = ['gfServer', str(order), str(host), str(port), '-canStop']
        if type == 'blat':
            gfserver_cmd.append('-stepSize=5')
        elif type == 'tblat':
            gfserver_cmd += ['-trans', '-mask']
        if gfserver_suppl_args:
            gfserver_cmd += gfserver_suppl_args
        gfserver_cmd_str = ' '.join(gfserver_cmd + [twobit])
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
                raise BLATServerError(out)
            else:
                if verbose:
                    print(out)
                return port
        if tries > try_limit:
            raise TimeoutError('Timed out!')


def id_ranker(record, perc_score, perc_query_span, perc_ident, expect=None,
              indent=0, verbose=1, same_strand=True, return_only=None):
    """Filters results based on score, expectation value, length, percent identity, and span; returns a sorted list.

    :param query_record record: Either a SearchIO.QueryResult or a Bio.Blast.Record.
    :param float perc_score: Minimum percentage of top score for a hit.
    :param float expect: Maximum e-value for a hit (BLAST-only).
    :param float perc_query_span: Minimum percent of the longest hit by query coverage for a hit.
    :param int perc_ident: Minimum percent identity of a hit.
    :param int indent: Indent level for pretty print. [Default: 0]
    :param int verbose: Level of verbose output? [Default: 1]
    :param bool same_strand: Should the function filter hits with HSPs on different strands? [Default:True]
    :param return_only: Should all or only one id be returned?
    :return list: Returns a list of tuples containing the final hit data in BED6 format.  
    """
    id_list = []
    if verbose:
        print('Beginning ID_Ranker...', indent=indent)
    if record.program == 'blat':
        if verbose > 2:
            print('Results obtained from BLAT run.', indent=indent + 1)
    elif 'blast' in record.program:
        if verbose > 2:
            print('Results obtained from BLAST run.', indent=indent + 1)
    else:
        raise NotImplementedError('Sorry, your program {} is not yet '
                                  'implemented for RecBlast!'.format(record.program))

    # Create filter functions:
    def hsp_minscores(hsp):
        return hsp.score >= int(perc_score * top_score)

    def hsp_min_query_span(hsp):
        return hsp.query_span >= perc_query_span * top_length

    def hsp_perc_ident(hsp):
        return hsp.ident_pct >= perc_ident

    def hsp_same_strand(hsp):
        if same_strand:
            return all([i == hsp.hit_strand_all[0] for i in hsp.hit_strand_all])
        else:
            return True

    def hit_sort_scores(hit):
        return sum([hsp.score for hsp in hit.hsps])

    def hsp_sort_scores(hsp):
        return hsp.score

    # Get top stats:
    top_score = max([max([hsp.score for hsp in hit.hsps]) for hit in record])
    if verbose > 1:
        print('Top score for {}:\t'.format(record.id), top_score, indent=indent)
    top_length = max([max([hsp.query_span for hsp in hit]) for hit in record])
    if verbose > 1:
        print('Longest hit for {}:\t'.format(record.id), top_length, indent=indent)

    if verbose > 2:
        print("ALL HITS STATS:")
        print('|\tHit Name:\t|\t# HSPs\t|\tScore:\t|\tLength:\t|\tP.Ident\t|')
        print("==========================================================")
        for hit in record:
            name = hit.id
            n_hsp = len(hit.hsps)
            print('|\t{HitName}\t|\t{HSP}\t|'.format(HitName=name, HSP=n_hsp))
            print("------------------------------------------------------")
            for hsp in hit:
                print('|\t{id}\t|\t{hf}\t|\t{score}\t|\t{length}\t|\t{ident}\t|'.format(id=hsp.hit_id,
                                                                                        hf=len(hsp),
                                                                                        score=hsp.score,
                                                                                        length=hsp.hit_span,
                                                                                        ident=hsp.ident_pct))

    # Execute filters:
    # query_span
    if verbose > 1:
        print('Number of HSPs for {}:\t'.format(record.id), sum([len(i.hsps) for i in record]), indent=indent)
        print('Filtering out all HSPs shorter than {}...'.format(perc_query_span * top_length), indent=indent)
    record = record.hsp_filter(hsp_min_query_span) if perc_query_span else record
    if not record:
        text = ('No hits in Query Results match a stretch of the query sequence longer than '
                '{0}!').format((top_length * perc_query_span))
        raise NoHitsError(text)
    # Score
    if verbose > 1:
        print('Number of HSPs for {}:\t'.format(record.id), sum([len(i.hsps) for i in record]), indent=indent)
        print('Filtering out all HSPs with scores less than {}...'.format(top_score * perc_score), indent=indent)
    record = record.hsp_filter(hsp_minscores) if perc_score else record
    if not record:
        text = 'No hits in Query Results have a score above the minimum of {0}!'.format((top_score * perc_score))
        raise NoHitsError(text)
    if verbose > 1:
        print('Number of HSPs for {}:\t'.format(record.id), sum([len(i.hsps) for i in record]), indent=indent)
        print('Filtering out all HSPs with percent identity below {}...'.format(perc_ident), indent=indent)
    record = record.hsp_filter(hsp_perc_ident) if perc_ident else record
    if not record:
        text = 'No hits in Query Results have a percent identity above {}%!'.format(round(perc_ident * 100, 2))
        raise NoHitsError(text)
    if verbose > 1:
        print('Number of HSPs for {}:\t'.format(record.id), sum([len(i.hsps) for i in record]), indent=indent)
        if same_strand:
            print('Filtering out all HSPs that have fragments on opposite strands...')
    record = record.hsp_filter(hsp_same_strand) if same_strand else record
    if not record:
        text = 'No hits in Query Results are on the same strand!'
        raise NoHitsError(text)
    # Sorting them for good measure
    if verbose > 1:
        print('Sorting all hits by descending scores!', indent=indent)
    record.sort(key=hit_sort_scores, reverse=True, in_place=True)
    for hit in record:
        hit.sort(key=hsp_sort_scores, reverse=True, in_place=True)
    if verbose > 1:
        print('Done!', indent=indent)
    # Add items to id_list
    # Big note: think in HSPs, not Hits
    n = 1
    for hit in record:
        for hsp in hit:
            # some quick strand math:
            if hsp._has_hit_strand:
                strands = set(hsp.hit_strand_all)
                if len(strands) == 1:
                    strand = "+" if strands == {1} else "-"
                else:
                    strand = "."
            else:
                strand = "."
            if verbose > 2:
                print("Adding hit {chr}:{s}-{e}({st}) to id list".format(chr=hsp.hit_id,
                                                                         s=str(hsp.hit_range[0]),
                                                                         e=str(hsp.hit_range[1]),
                                                                         st=strand),
                      indent=indent)
            # A little witchcraft before we do though
            # turns out hsp.hit_start_all won't necessarily start with the starting point of the hit...
            # That means we need to zip hit_start_all and hit_span_all, sort by the first one, then de-zip.
            block_starts, block_spans = zip(*sorted(zip(hsp.hit_start_all, hsp.hit_span_all), key=itemgetter(0)))

            # chr (start,end) id score strand thickStart thickEnd rgb blockcount blockspans blockstarts query_span
            id_list.append([hsp.hit_id, hsp.hit_range, hsp.query_id, hsp.score, strand, hsp.hit_range[0],
                            hsp.hit_range[1], "255,0,0", len(hsp.hit_start_all),
                            ",".join([str(i) for i in block_spans]),
                            ",".join([str(i - hsp.hit_range[0]) for i in block_starts]), hsp.query_range])

            if return_only and n == return_only:
                print('Returning only the top {} hits, ending here!'.format(return_only),
                      indent=indent)
                return id_list
            n += 1
    return id_list
