import re
from time import sleep
import subprocess
from io import StringIO
from pathlib import Path
from operator import itemgetter
from Bio import SearchIO, SeqIO
from Bio.Blast import NCBIXML, NCBIWWW
from RecBlast import print, merge_ranges, flatten
from RecBlast.WarningsExceptions import *


class Search(object):
    def __init__(self, search_type):
        self.search_type = search_type

    def __call__(self, seq_record, species, database, filetype, search_db_path, search_type, local, expect,
                 n_threads, megablast, indent, perc_ident, verbose, write,
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
                print(search_type, 'was selected.', indent=indent)
            dt = self.blast_prep(search_type=self.search_type, db_loc=search_db_path, database=database,
                                 species=species, verbose=verbose, indent=indent)
            return self.blast_run(seq_record=seq_record, species=species, database=dt, filetype=filetype,
                                  blast_type=search_type, local_blast=local, expect=expect, megablast=megablast,
                                  use_index=False, perc_ident=perc_ident, verbose=verbose, indent=indent,
                                  n_threads=n_threads, blastdb=search_db_path, outtype=5, return_raw=False,
                                  **kwargs)
        elif self.search_type in ['blat', 'tblat', 'blat-transcript', 'tblat-transcript']:
            if verbose > 1:
                print(search_type, 'was selected.', indent=indent)
            dt = self.blat_prep(database=database, species=species, verbose=verbose, indent=indent)
            return self.blat_run(seq_record=seq_record, species=species, database=dt, filetype=filetype,
                                 blat_type=search_type, perc_ident=perc_ident, verbose=verbose, indent=indent,
                                 blatdb=search_db_path, outtype='pslx', **kwargs)
        else:
            raise SearchEngineNotImplementedError('Invalid selection for search type!')

    @staticmethod
    def blast_run(seq_record, species, database, filetype="fasta", blast_type='blastn',
                  local_blast=False, expect=0.005, megablast=True, use_index=False, perc_ident=75,
                  verbose=True, indent=0, n_threads=1, blastdb='/usr/db/blastdb/', outtype=5,
                  return_raw=False, **kwargs):
        """A wrapper function for BLAST searches.
        :param seq_record: The record containing the query sequence for the search. Can be either a SeqIO.SeqRecord or
                           a string with the file loaction.
        :param str species: The species whose sequence database will be queried.
        :param Union[dict, str, Path] database: The name of the database to be used in the search.
        :param str filetype: Filetype of seq_record (if seq_record is a SeqRecord object, leave as default.
                             [default: 'fasta']
        :param str blast_type: Type of search to conduct. Can be a BLAST type (blastn, blastp, blastx, tblastn, tblastx)
                               or a BLAT type (blat, tblat). [Default: 'blastn']
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
                print(args_expanded, indent=indent+1)
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
                    print('{0}\t=\t{1}'.format(k, v), indent=indent+1)
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
    def blat_run(seq_record, database, filetype="fasta", blat_type='blat', perc_ident=75, verbose=True,
                 indent=0, blatdb='/usr/db/blastdb/', outtype='pslx', **kwargs):
        """A wrapper function for BLAT searches.
        :param seq_record: The record containing the query sequence for the search. Can be either a SeqIO.SeqRecord or
                           a string with the file loaction.
        :param Union[dict, int] database: The name of the database to be used in the search.
        :param str filetype: Filetype of seq_record (if seq_record is a SeqRecord object, leave as default.
                             [default: 'fasta']
        :param str blat_type: Type of search to conduct. Can be a BLAST type (blastn, blastp, blastx, tblastn, tblastx) 
                              or a BLAT type (blat, tblat). [Default: 'blastn']
        :param int perc_ident: Minimum percent identity required of results to be returned [Default: 75]
        :param bool verbose: Verbose output? [Default: True]
        :param int indent: Indent level for pretty print. [Default: 0]
        :param str blatdb: Path of databases for either BLAST or BLAT. [Default: '/usr/db/blastdb'
        :param str outtype: Output type. (see options for BLAST and BLAT) [Default: pslx]
        :param kwargs: Additional keyword arguments to pass on to BLAST/BLAT.
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

        args_expanded = ['gfClient', 'localhost', str(database), '/', '/dev/stdin', '/dev/stdout']
        args_expanded += ['-t=dnax', '-q=prot'] if blat_type.lower() == 'tblat' else []
        args_expanded += ['minIdentity={}'.format(perc_ident), '-out={}'.format(outtype)]

        try:
            if verbose > 1:
                print('BLAT command:', indent=indent)
                print(' '.join(args_expanded), indent=indent+1)

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
                    raise ValueError
            except Exception as err:
                print('Error reading BLAT results! Aborting!')
                print('Error details:\n')
                raise err
        return blat_record, blat_err

    @staticmethod
    def blat_prep(database, species, verbose, indent):
        if isinstance(database, dict):
            try:
                blat_port = database[species]
                if verbose > 1:
                    print('Using port {} as blat port.'.format(blat_port), indent=indent)
            except KeyError:
                raise SearchError('No port found for species {}!'.format(species))
        elif isinstance(database, int):
            blat_port = database
        elif isinstance(database, str):
            try:
                blat_port = int(database)
            except ValueError:
                raise SearchError('Invalid option "{}" was passed to database! Database must be '
                                  'either a dictionary of species-port pairs or an integer!'.format(database))
        else:
            raise SearchError('Invalid option of type "{}" was passed to database! Database must be '
                              'either a dictionary of species-port pairs or an integer!'.format(str(type(database))))
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
    strand = '(+)'

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
                    print('ID was selected as type {0}!'.format(tmp_type), indent=indent+1)
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
                print('Item:\t', '\t'.join(item_parts), indent=indent+1)
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
                        print(sr_tuple, indent=indent+1)
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
        raise SearchError('Improper search type given ({})!'.format(search_type))
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
        if not glob_path:
            if verbose:
                print('No DB found! Trying again with abbreviated species name')
            species_abbv = ''.join([i[0:3] for i in species.title().split(' ')])
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
                search_db = sorted(glob_path, reverse=True)[0].stem
            else:
                search_db = glob_path.stem
        except IndexError:
            print('WARNING: COULD NOT FIND DATABASE! ABORTING!', indent=indent)
            raise DatabaseNotFoundError('', 'No databases were found!')
    else:
        raise DatabaseNotFoundError('DB_Path {} does not exist!'.format(str(db_path)))
    if verbose:
        print('{0} DB chosen: {1}'.format(search_type, search_db), indent=indent)
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
        blat_2bit = get_searchdb(search_type='blat', species=species, db_loc=search_db_loc,
                                 verbose=verbose, indent=indent + 1)
        twobit = Path(search_db_loc, blat_2bit + '.2bit').absolute()
        twobit = str(twobit) if twobit.is_file() else None
        if twobit is None:
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
                raise BLATServerError(out)
            else:
                if verbose:
                    print(out)
                return port
        if tries > try_limit:
            raise TimeoutError('Timed out!')


def id_ranker(record, perc_score, expect, perc_length, perc_ident, perc_span=0.1, min_hsps=1, hsp_cumu_score=True,
              seq_method='whole', indent=0, verbose=1, method='all', samestrand=True):
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
    :param int indent: Indent level for pretty print. [Default: 0]
    :param int verbose: Level of verbose output? [Default: 1]
    :param str method: Return all ranked hits ('all'), or only the top hit ('best-hit')? [Default: 'all']
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
            raise NotImplementedError('Sorry, your program {} is not yet '
                                      'implemented for RecBlast!'.format(record.program))

        # Create filter functions:
        def hit_minhsps(hit):
            return len(hit.hsps) >= min_hsps

        def hit_minscores(hit):
            return sum([hsp.score for hsp in hit.hsps]) >= int(perc_score * top_score)

        def hit_minlength(hit):
            return sum([i[-1] - i[0] for i in merge_ranges([(hsp.query_start, hsp.query_end) for hsp in hit])
                        ]) >= perc_length * top_length

        def hit_perc_id(hit):
            # return _percent_identity_searchio(hit) >= perc_ident
            return True

        def hit_hsp_span(hit):
            top_span = max([hsp.hit_span for hsp in hit])
            hit = hit.filter(lambda hsp: hsp.hit_span >= int(perc_span * top_span))
            return hit

        def hit_same_strand(hit):
            x = [bla.hit_strand_all for bla in hit.hsps]
            y = all(s > 0 for s in flatten(x)) or all(s < 0 for s in flatten(x)) or \
                all(s == 0 for s in flatten(x)) or None
            return y

        # Three more functions for to make great progress
        def sort_scores(hit):
            return sum([hsp.score for hsp in hit.hsps])

        def hit_target_span(hit):
            return list(merge_ranges([(hsp.hit_start, hsp.hit_end) for hsp in hit]))

        def hit_query_coverage(hit):
            return list(merge_ranges(flatten([list(merge_ranges(hsp.query_range_all)) for hsp in hit])))

        # Get top stats:
        top_score = max([sum([hsp.score for hsp in hit.hsps]) for hit in record])
        if verbose > 1:
            print('Top score for {}:\t'.format(record.id), top_score, indent=indent)
        top_length = max([sum([i[-1] - i[0] for i in merge_ranges([(hsp.query_start, hsp.query_end)
                                                                   for hsp in hit])
                               ]) for hit in record])
        if verbose > 1:
            print('Longest hit for {}:\t'.format(record.id), top_length, indent=indent)

        if verbose > 2:
            print("ALL HITS STATS:")
            print('Hit Name:\t|\t# HSPs\t|\tScore:\t|\tLength:\t|\tP.Ident\t|\thsp_span_list\t|')
            for hit in record:
                name = hit.id
                n_hsp = len(hit.hsps)
                score = sum([hsp.score for hsp in hit.hsps])
                length = sum([i[-1] - i[0] for i in merge_ranges([(hsp.query_start, hsp.query_end) for hsp in hit])])
                ident = percent_identity_searchio(hit)
                span = [hsp.hit_span for hsp in hit]
                print('{HitName}\t|\t{HSP}\t|\t{Score}\t|'
                      '\t{Length}\t|\t{PIdent}\t|\t{span_list}\t|'.format(HitName=name, HSP=n_hsp, Score=score,
                                                                          Length=length, PIdent=ident, span_list=span))
        # Execute filters:
        # HSP
        if verbose > 1:
            print('Number of hits for {}:\t'.format(record.id), len(record), indent=indent)
            print('Filtering based on min. number of HSPs...', indent=indent)
        record1 = record.hit_filter(hit_minhsps)
        if verbose > 1:
            print('Number of hits for {}:\t'.format(record1.id), len(record1), indent=indent)
        if not record1:
            raise NoHitsError('No hits in Query Results have {} or more HSPs!'.format(min_hsps))
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
            raise NoHitsError('No hits in Query Results above min_length {0}!'.format((top_length * perc_length)))
        # Score
        if verbose > 1:
            print('Filtering out all hits with scores less than {}...'.format(top_score * perc_score), indent=indent)
        record4 = record3.hit_filter(hit_minscores)
        if verbose > 1:
            print('Number of hits for {}:\t'.format(record4.id), len(record4), indent=indent)
        if not record4:
            raise NoHitsError('No hits in Query Results above minimum score {0}!'.format((top_score * perc_score)))
        # Percent Identity
        if verbose > 1:
            print('Filtering out all hits with a percent identity below {}...'.format(perc_ident),
                  indent=indent)
        record5 = record4.hit_filter(hit_perc_id)
        if verbose > 1:
            print('Number of hits for {}:\t'.format(record5.id), len(record5), indent=indent)
        if not record5:
            raise NoHitsError('No hits in Query Results above minimum score {0}!'.format((top_score * perc_score)))
        # If strand is set to strict, it will filter out hits that have different strandedness in HSPs
        if samestrand:
            if verbose > 1:
                print('Filtering out all hits whose HSPs are not on the same strand...', indent=indent)
            record6 = record5.hit_filter(hit_same_strand)
            if verbose > 1:
                print('Number of hits for {}:\t'.format(record6.id), len(record6), indent=indent)
            if not record6:
                raise NoHitsError('No hits in Query Results with all HSPs on the same strand!')
        else:
            record6 = record5
        # Sorting them for good measure
        if verbose > 1:
            print('Sorting all hits by descending scores!', indent=indent)
        record6.sort(key=sort_scores, reverse=True, in_place=True)
        if verbose > 1:
            print('Done!', indent=indent)
        # Add items to id_list
        for hit in record6:
            seq_name = hit.id
            hts = hit_target_span(hit)
            seq_cov = hit_query_coverage(hit)
            if seq_method == 'whole':
                seq_range = [hts[0][0], hts[-1][-1]]
                seq_coverage = [seq_cov[0][0], seq_cov[-1][-1]]
            elif seq_method == 'strict':
                seq_range = [(i[0], i[-1]) for i in hts]
                seq_coverage = [(s[0], s[-1]) for s in seq_cov]
            else:
                seq_range = ''
                seq_coverage = ''
            strands = flatten([bla.hit_strand_all for bla in hit.hsps])
            if all(s > 0 for s in strands):
                strand = '(+)'
            elif all(s < 0 for s in strands):
                strand = '(-)'
            elif all(s == 0 for s in strands):
                strand = '(0)'
            else:
                strand = '(N)'
            seq_score = sum([hsp.score for hsp in hit.hsps])
            if verbose > 2:
                print("Adding hit {} to id list".format(seq_name + ':' + '-'.join([str(i) for i in seq_range[0:2]])),
                      indent=indent)
            id_list.append((seq_name, seq_range, strand, seq_score, seq_coverage))
            if method == 'best hit':
                print('Best Hit Reciprocal BLAST was selected, ending Reverse BLASTS after first annotation!',
                      indent=indent)
                break
    else:
        RecBlastWarning('No guarantees that this is going to work as of commit 81d3d36')
        align_scorelist = []
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
