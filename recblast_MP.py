import logging
import multiprocessing
import re
import subprocess
import sys
from datetime import datetime as dt
from inspect import isgenerator
from io import StringIO
from operator import itemgetter
from pathlib import Path
from time import sleep

from Bio import SearchIO
from Bio import SeqIO
from Bio import __version__ as bp_version
from Bio.Blast import NCBIXML
from Bio.Blast.Record import Blast as BioBlastRecord
from BioSQL import BioSeqDatabase


class ProgressBar(object):
    """Adapted from Romuald Brunet at StackExchange"""
    DEFAULT = 'Progress: %(bar)s %(percent)3d%%'
    FULL = '%(bar)s %(current)d/%(total)d (%(percent)3d%%) %(remaining)d to go'

    def __init__(self, total, width=100, fmt=DEFAULT, symbol='=',
                 output=sys.stderr):
        assert len(symbol) == 1

        self.total = total
        self.width = width
        self.symbol = symbol
        self.output = output
        self.fmt = re.sub(r'(?P<name>%\(.+?\))d', r'\g<name>%dd' % len(str(total)), fmt)

        self.current = 0

    def __call__(self):
        percent = self.current / float(self.total)
        size = int(self.width * percent)
        remaining = self.total - self.current
        bar = '[' + self.symbol * size + ' ' * (self.width - size) + ']'

        args = {
            'total': self.total,
            'bar': bar,
            'current': self.current,
            'percent': percent * 100,
            'remaining': remaining
        }
        print('\r' + self.fmt % args, file=self.output, end='')

    def done(self):
        self.current = self.total
        self()
        print('', file=self.output)


def get_searchdb(search_type, species, db, verbose=True, level=1):
    if verbose:
        print('\t' * level, 'Forward Blast DB set to auto, chosing fw_blast_db...')
    if verbose > 1:
        print('\t' * level, 'Blast DB location set to: ', db)
    if search_type.lower() in ['blastp', 'blastx', 'tblastx']:
        db_type = 'protein'
    elif search_type.lower() in ['blastn', 'tblastn']:
        db_type = 'genome'
    elif search_type.lower() in ['blat', 'tblat', 'translated_blat', 'untranslated_blat']:
        db_type = 'blat'
    else:
        print('\t' * level, 'Unable to determing blast db type!')
        raise Exception('Improper search type given: ', search_type)
    if verbose > 1:
        print('\t' * level, 'DB type: ', db_type)
    db_path = Path(db).absolute()
    if db_path.exists() and db_path.is_dir():
        if db_type == 'blat':
            glob_path = db_path.glob('{0}*.2bit'.format(species.replace(' ', '_')))
        else:
            glob_path = db_path.glob('{0}_{1}*'.format(species.replace(' ', '_'), db_type))
        try:
            if isgenerator(glob_path):
                glob_path = [i for i in glob_path]
            if verbose:
                print('\t' * level, glob_path)
            if isinstance(glob_path, list):
                search_db = sorted(glob_path, reverse=True)[0].stem
            else:
                search_db = glob_path.stem
        except IndexError:
            print('\t' * level, 'WARNING: COULD NOT FIND DATABASE! ABORTING!')
            raise Exception('DatabaseError:', 'No databases were found!')
    if verbose:
        print('\t' * level, '{0} DB chosen: {1}'.format(search_type, search_db))
    return search_db


def nuc_to_prot_compare(prot_sub, nuc_query, perc_ident, perc_length, trans_table=1, verbose=True):
    from Bio.SeqRecord import SeqRecord
    from Bio import AlignIO
    from Bio.Align.Applications import MuscleCommandline
    tnuc_list = []
    for strand, nuc in [(+1, nuc_query), (-1, nuc_query.reverse_complement())]:
        for frame in range(3):
            tnuc_list += [SeqRecord(seq=p, id='{0}-({1})-{2}.{3}'.format(nuc_query.id, strand, frame + 1, n),
                                    description="Translation of {0} ORF#{1} in frame {2} on {3} strand. Original Description:".format(
                                        nuc_query.id, n, frame + 1, strand) + nuc_query.description
                                    ) for n, p in enumerate(nuc[frame:].seq.translate().split("*")) if
                          len(p) >= (len(prot_sub) * perc_length)]
    temp_fasta = StringIO()
    SeqIO.write([prot_sub] + tnuc_list, temp_fasta, 'fasta')
    muscle_cline = MuscleCommandline()
    stdout, stderr = muscle_cline(stdin=temp_fasta.getvalue())
    align = AlignIO.read(StringIO(stdout), "fasta")
    for record in align:
        print("%s - %s" % (str(record.seq), record.id))
    return align


def blat_server(twobit, status, host='blatMachine', port='2000', stepSize=5, log='blatMachine.log', **gfserver_kwargs):
    # Regular: gfServer start blatMachine portX -stepSize=5 -log=untrans.log database.2bit
    # Prot>DNA:  gfServer start blatMachine portY -trans -mask -log=trans.log database.2bit
    gfserver_suppl_args = list()
    for key, item in gfserver_kwargs.items():
        if key == 'status':
            status = item
        elif key == 'host':
            host = item
        elif key == 'port':
            port = item
        else:
            gfserver_suppl_args.append('-{0}={1}'.format(key, item))
    gfserver_cmd = ['gfServer', status, host, port, '-canStop', '-stepSize={}'.format(stepSize), '-log={}'.format(log),
                    twobit]
    if gfserver_suppl_args != list():
        gfserver_cmd += gfserver_suppl_args
    with subprocess.Popen(gfserver_cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT) as proc:
        for line in proc.stdout:
            if 'Server ready for queries!' in line:
                print(line)
                return
    return




def merge_ranges(ranges):
    """
    Merge overlapping and adjacent ranges and yield the merged ranges in order.
    The argument must be an iterable of pairs (start, stop).
    (Source: Gareth Rees, StackExchange)

    >>> list(merge_ranges([(5,7), (3,5), (-1,3)]))
    [(-1, 7)]
    >>> list(merge_ranges([(5,6), (3,4), (1,2)]))
    [(1, 2), (3, 4), (5, 6)]
    >>> list(merge_ranges([]))
    []
    """
    ranges = iter(sorted(ranges))
    try:
        current_start, current_stop = next(ranges)
    except StopIteration:  # ranges is empty
        return
    for start, stop in ranges:
        if start > current_stop:
            # Gap between segments: output current segment and start a new one.
            yield current_start, current_stop
            current_start, current_stop = start, stop
        else:
            # Segments adjacent or overlapping: merge.
            current_stop = max(current_stop, stop)
    yield current_start, current_stop


def id_ranker_searchio(qresults, min_hsps=1):
    assert isinstance(qresults, SearchIO.QueryResult), 'input record must be a SearchIO.QueryResult object!!!'

    # lambda functions for filtering:
    hit_minhsps = lambda hit: len(hit.hsps) >= min_hsps
    qresults = qresults.hit_filter(hit_minhsps)
    if not qresults:
        raise Exception('No hits in Query Results above min_hsps!')
    for hit in qresults:
        pass



def id_ranker(blastrecord, perc_score, expect, perc_length, hsp_cumu_score=True, align_scorelist=list(),
              indent='', verbose=True):
    indent_internal = indent + '\t'
    truthple = []
    id_list = []
    subject_range = []
    query_start_end = []
    if verbose > 1:
        print(indent, 'Sorting through alignment\'s HSPs to get top scores of all alignments...')
    for alignment in blastrecord.alignments:
        subject_range_hsp = []
        query_start_end_hsp = []
        hsp_scorelist = []
        if verbose > 3:
            print('Number of HSPs: ', len(alignment.hsps))
        for hsp in alignment.hsps:
            hsp_scorelist.append(hsp.score)
            subject_range_hsp.append(hsp.sbjct_start)
            subject_range_hsp.append(hsp.sbjct_end)
            query_start_end_hsp.append((hsp.query_start, hsp.query_end))
        hsp_scorelist.sort(reverse=True)
        query_start_end.append([i for i in merge_ranges(query_start_end_hsp)])
        subject_range.append((subject_range_hsp[0], subject_range_hsp[-1]))
        if verbose > 3:
            print(indent_internal, "HSP Score List: \n", indent_internal+'\t', hsp_scorelist)
        if hsp_cumu_score:
            align_scorelist.append(sum(hsp_scorelist))
        else:
            align_scorelist.append(hsp_scorelist[0])
    if verbose > 1:
        print(indent, 'Done with first-round sorting!')
    if verbose > 3:
        for align_index, alignment in enumerate(blastrecord.alignments):
            print(indent_internal, alignment.title)
            print(indent_internal, "\tAlignment Score List: \n\t", indent_internal, align_scorelist[align_index])
            print(indent_internal, "\tQuery_start_end: \n\t", indent_internal, query_start_end[align_index])
            print(indent_internal, "\tSubject Range: \n\t", indent_internal, subject_range[align_index])
    if verbose > 1:
        print(indent, 'Sorting through alignments to get all hits above threshold...')

    for align_index, alignment in enumerate(blastrecord.alignments):
        '''if verbose >3:
            print(indent_internal, 'Alignment title: ', alignment.title)'''
        score_threshold = (perc_score * align_scorelist[0])
        length_alignment = sum([i[-1] - i[0] for i in query_start_end[align_index]])
        align_len_threshold = blastrecord.query_length * perc_length
        if hsp_cumu_score:
            hsp_scoretotal = sum([hsp.score for hsp in alignment.hsps])
            truthple.append((align_index, hsp_scoretotal >= score_threshold, hsp.expect <= expect,
                             length_alignment >= align_len_threshold, alignment.title))
        else:
            truthple.append((align_index, hsp.score >= score_threshold, hsp.expect <= expect,
                             length_alignment >= align_len_threshold, alignment.title))
    if verbose > 3:
        print(indent, 'List of Alignments and Criteria Status:')
        print(indent_internal, 'i\t', 'Score\t', 'Expect\t', 'Length\t', 'Alignment.Title')
        for i in range(len(blastrecord.alignments)):
            print(indent_internal, "{0}\t{1}\t{2}\t{3}\t{4}".format(truthple[i][0], truthple[i][1], truthple[i][2],
                                                                    truthple[i][3], truthple[i][4]))
    for i in truthple:
        if i[1] and i[2] and i[3]:
            id_list.append((blastrecord.alignments[i[0]].title,
                            '[:{0}-{1}]'.format(subject_range[i[0]][0],
                                                subject_range[i[0]][1]),
                            align_scorelist[i[0]]))
            if verbose > 2:
                print(indent, "Alignment {} added to id_list!".format(i[4]))
        else:
            if verbose > 2:
                print(indent, "WARNING: ALIGNMENT {} FAILED TO MEET CRITERIA!".format(i[4]))
                if not i[1]:
                    print(indent_internal, 'Score was below threshold!')
                if not i[2]:
                    print(indent_internal, 'Expect was below threshold!')
                if not i[3]:
                    print(indent_internal, 'Length was below threshold!')
    sorted(id_list, reverse=True, key=itemgetter(2))
    return id_list


def biosql_get_sub_db_names(passwd, db="bioseqdb", driver="psycopg2", user="postgres", host="localhost"):
    from BioSQL import BioSeqDatabase
    server = BioSeqDatabase.open_database(driver=driver, user=user, passwd=passwd, host=host, db=db)
    sub_db_name_list = [i for i in server.keys()]
    return sub_db_name_list


def biosql_DBSeqRecord_to_SeqRecord(DBSeqRecord_, off=False):
    """
    As I wrote this script I realized very quickly that there are a great many issues with DBSeqRecords throwing bugs
    left and right, so I got tired of trying to debug a blackbox and instead decided to convert it to something that
    works.
    :param DBSeqRecord_:
    :return:
    """
    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import Seq
    if off:
        return DBSeqRecord_
    else:
        return SeqRecord(seq=Seq(str(DBSeqRecord_.seq)), id=DBSeqRecord_.id, name=DBSeqRecord_.name,
                         description=DBSeqRecord_.description, dbxrefs=DBSeqRecord_.dbxrefs,
                         features=DBSeqRecord_.features, annotations=DBSeqRecord_.annotations,
                         letter_annotations=DBSeqRecord_.letter_annotations)


class RecBlastContainer(dict):
    """
    dict(proc_id = '', query_record = '', forward_blast = dict(forward_blast_results = '', forward_blast_errors = ''),
         forward_ids = dict(forward_ids = list, missing_ids = list), recblast_unanno = list,
         reverse_blast = dict(reverse_blast_results = '',
         )
    """

    def __init__(self, target_species, query_record):
        super(dict, self).__init__()
        assert isinstance(query_record, SeqIO.SeqRecord), "Query Record MUST be a SeqRecord!"
        self[target_species] = {query_record.name: dict(proc_id=str(), query_record=query_record,
                                                        forward_blast=dict(blast_results=BioBlastRecord,
                                                                           blast_errors=''),
                                                        forward_ids=dict(ids=list(), missing_ids=list(),
                                                                         pretty_ids=list()),
                                                        recblast_unanno=list(),
                                                        reverse_blast=dict(blast_results=BioBlastRecord,
                                                                           blast_errors=''),
                                                        reverse_ids=dict(ids=list(), missing_ids=list(),
                                                                         pretty_ids=list()),
                                                        recblast_results=list(),
                                                        output_paths=dict(forward_blast_output=Path(),
                                                                          forward_id_score_output=Path(),
                                                                          recblast_output_unanno=Path(),
                                                                          reverse_blast_output=list(),
                                                                          recblast_output=Path(),
                                                                          blast_nohits=Path()
                                                                          )
                                                        )}

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

    """
    def __str__(self, indent=''):
        super(RecBlastContainer, self).__str__()
        strobj = ''
        if isinstance(self, dict):
            for k, v in self.items():
                if indent:
                    strobj += ''.join([indent, '|---- Key "', str(k), '":', '\n'])
                else:
                    strobj += ''.join([indent, '| Key "', str(k), '":', '\n'])
                if isinstance(v, dict):
                    indent += '\t'
                    strobj += v.__str__(indent)
                else:
                    strobj += ''.join(['\t', indent, '|---', str(v).replace('\n', '\n\t' + indent + '|--- '), '\n'])
        return strobj
    """

def blast(seq_record, target_species, database, query_species="Homo sapiens", filetype="fasta", blast_type='blastn',
          local_blast=False, expect=0.005, megablast=True, use_index=False, blastoutput_custom="", perc_ident=75,
          verbose=True, n_threads=1, write=False, BLASTDB='/usr/db/blastdb/', **blast_kwargs):
    from pathlib import Path

    from Bio.Blast import NCBIWWW
    if isinstance(seq_record, SeqIO.SeqRecord):
        pass
    else:
        seq_record = SeqIO.read(seq_record, filetype)
    args = dict()
    if verbose:
        print("Now starting BLAST...")

    if blast_type.lower() in ['blat', 'tblat']:
        blat_db = get_searchdb(search_type=blast_type, species=target_species, db=BLASTDB, verbose=bool(verbose),
                               level=1)
        args_expanded = ['gfClient', 'localhost', str(database), BLASTDB, seq_record.format('fasta'), '/dev/stdout']
        if blast_type.lower() == 'tblat':
            args_expanded += ['-t=dnax', '-q=prot']
        args_expanded += ['minIdentity={}'.format(perc_ident), '-out=pslx']
        try:
            blat_handle = (subprocess.check_output(args_expanded, universal_newlines=True, cwd=BLASTDB))
            if isinstance(blat_handle, str):
                blat_result = blat_handle
                blat_err = None
            else:
                blat_result, blat_err = blat_handle
        except subprocess.CalledProcessError:
            raise
        blast_result, blast_err = blat_result, blat_err
        search_rec_type = 'blat-psl'

    else:
        search_rec_type = 'blast-xml'
        if local_blast:
            # build up the BLAST arguments:
            args.update({'-db': database, '-evalue': expect,
                         '-outfmt': '5',
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
                print('Running Local Blast...')
                print('Options:')
                print('\t', args_expanded)
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
                args['megablast'] = True
            if blast_kwargs:
                args.update(blast_kwargs)
            if verbose:
                print('Submitting Remote BLAST! Options passed:')
                for k, v in args.items():
                    print('\t {0}\t=\t{1}'.format(k, v))
            try:
                blast_result = NCBIWWW.qblast(**args)
                blast_err = None
            except Exception as err:
                print(err)
                raise err

        if verbose:
            print('Done with Blast!')
    """
    with StringIO(blast_result) as fin:
        try:
            blast_record = SearchIO.read(fin, search_rec_type)
        except Exception as err:
            print('Error reading Forward Blast Results! Aborting!')
            print('Error details:\n', err)
            raise err
    """
    with StringIO(blast_result) as fin:
        try:
            blast_record = NCBIXML.read(fin)
        except Exception as err:
            print('Error reading Forward Blast Results! Aborting!')
            print('Error details:\n', err)
            raise err
    if blast_type in ['blat', 'tblat']:
        pass
        # TODO: once I'm more familiar with SearchIO, fill in some of the unknowns like targetdb, etc
    if write:
        if blastoutput_custom == '':
            blastoutput_custom = Path("{0}_blast".format(target_species),
                                      "{0}_{1}_{2}_to_{3}.xml".format(blast_type, seq_record.name,
                                                                      query_species, target_species)).absolute()
        else:
            blastoutput_custom = Path(blastoutput_custom).absolute()
        try:
            blastoutput_custom.parent.mkdir(parents=True)
        except FileExistsError:
            pass
        with blastoutput_custom.open("w") as fxml:
            fxml.write(blast_record)
    else:
        return blast_record, blast_err


def biosql_seq_lookup_cascade(dtbase, sub_db_name, id_type, identifier, verbose=False):
    seqrec = ''
    try_get_id = True
    if id_type == 'scaffold':
        lookup_key = 'name'
    else:
        lookup_key = id_type
    if try_get_id:
        try:
            if verbose:
                print("\t\tNow searching database {0} for {1}: {2}".format(sub_db_name, id_type, identifier))
            seqrec = biosql_DBSeqRecord_to_SeqRecord(dtbase.lookup(**{lookup_key: identifier}))
            if verbose:
                print('\tGot sequence for {}!'.format(identifier))
            try_get_id = False
        except KeyError:
            sleep(0.1)
            seqrec = biosql_DBSeqRecord_to_SeqRecord(dtbase.lookup(**{lookup_key: identifier}))
            if verbose:
                print('\tGot sequence for {}!'.format(identifier))
            try_get_id = False
        except IndexError as err:
            if verbose:
                print("WARNING: couldn't find {0} using given ID type... \n Full error: {1}".format(identifier, err))

    if try_get_id:
        identifier_sans_subnumber = identifier.split('.')[0]
        if verbose:
            print('\t\tSeeing if removing any sub-numbers (acc: xxxxxx.1 for example) helps...')
            print('\t\tIdentifier: ', identifier_sans_subnumber)
        try:
            if verbose:
                print("\t\tNow searching database {0} for {1}: {2}".format(sub_db_name, id_type,
                                                                           identifier_sans_subnumber))
            seqrec = biosql_DBSeqRecord_to_SeqRecord(dtbase.lookup(**{lookup_key: identifier_sans_subnumber}))
            if verbose:
                print('\tGot sequence for {}!'.format(identifier))
            try_get_id = False
        except KeyError:
            sleep(0.1)
            seqrec = biosql_DBSeqRecord_to_SeqRecord(dtbase.lookup(**{lookup_key: identifier_sans_subnumber}))
            if verbose:
                print('\tGot sequence for {}!'.format(identifier))
            try_get_id = False
        except IndexError as err1:
            if verbose:
                print("WARNING: couldn't find {0} using abbreviated ID... \n Full error: {1}"
                      .format(identifier_sans_subnumber,
                              err1))
    if try_get_id:
        try:
            if verbose:
                print('\t\tAttempting to search using Primary ID instead of declared type:')
            seqrec = biosql_DBSeqRecord_to_SeqRecord(dtbase.lookup(primary_id=identifier))
            if verbose:
                print('\tGot sequence for {}!'.format(identifier))
            try_get_id = False
        except KeyError:
            sleep(0.1)
            seqrec = biosql_DBSeqRecord_to_SeqRecord(dtbase.lookup(primary_id=identifier))
            if verbose:
                print('\tGot sequence for {}!'.format(identifier))
            try_get_id = False
        except IndexError as err2:
            if verbose:
                print("WARNING: couldn't find {0} using Primary ID... \n full error: {1}".format(identifier, err2))
    if try_get_id:
        try:
            if verbose:
                print('\t\tAttempting to search using name instead of declared type:')
            seqrec = biosql_DBSeqRecord_to_SeqRecord(dtbase.lookup(name=identifier))
            if verbose:
                print('\tGot sequence for {}!'.format(identifier))
            try_get_id = False
        except KeyError:
            sleep(0.1)
            seqrec = biosql_DBSeqRecord_to_SeqRecord(dtbase.lookup(name=identifier))
            if verbose:
                print('\tGot sequence for {}!'.format(identifier))
            try_get_id = False
        except IndexError as err3:
            if verbose:
                print("WARNING: Still couldn't find {0} using name search: \n full error: {1}".format(identifier, err3))

    if try_get_id:
        try:
            lookup_key = input('Last shot, chose an ID type: '
                            '[accession, primary_id, gi, version, display_id, name]')
            if lookup_key == 'exit':
                exit(exit(), 'Script ended!')
            seqrec = biosql_DBSeqRecord_to_SeqRecord(dtbase.lookup(**{lookup_key: identifier}))
            if verbose:
                print('\tGot sequence for {}!'.format(identifier))
        except IndexError as err5:
            if verbose:
                print("WARNING: COULD NOT FIND SEQUENCES FOR ID:{0}: \n full error: {1}".format(identifier, err5))
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
    def __init__(self, id_type, identifier, verbose):
        self.id_type = id_type
        self.identifier = identifier
        self.verbose = verbose

    def __call__(self, sub_db_name, server):
        if self.verbose:
            print('\tFetching sequence: ', self.identifier)
        try:
            dtbase = server[sub_db_name]
        except KeyError as err:
            print('Woah! KeyError!', err)
            print('Waiting for 0.1 second and rerunning in case it was a mistake!')
            sleep(0.1)
            try:
                dtbase = server[sub_db_name]
            except KeyError:
                raise

        seqrec = biosql_seq_lookup_cascade(dtbase=dtbase, sub_db_name=sub_db_name, id_type=self.id_type,
                                           identifier=self.identifier, verbose=self.verbose)
        return self.identifier, seqrec


def biosql_get_record_mp(sub_db_name, passwd='', id_list=list(), id_type='accession', driver="psycopg2",
                         user="postgres", host="localhost", db="bioseqdb", num_proc=2, verbose=True, server=None):
    """

    :param sub_db_name:
    :param passwd:
    :param id_list:
    :param id_type:
    :param driver:
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
        print('\tStarting BioSQL_get_record_mp')
    num_jobs = len(id_list)
    seqdict = dict()
    getseqs = [GetSeqMP(idents, results, db=db, host=host, driver=driver, user=user, passwd=passwd,
                        sub_db_name=sub_db_name, verbose=verbose, server=server) for i in range(num_proc)]
    for gs in getseqs:
        gs.start()

    for item in id_list:
        idents.put(BioSeqLookupCascade(id_type=id_type, identifier=item, verbose=verbose))

    for i in range(num_proc):
        idents.put(None)

    while num_jobs:
        temp = results.get()
        print(temp)
        seqdict[temp[0]] = temp[1]
        num_jobs -= 1
    if verbose:
        print('Done with Biosql_get_record_mp!')
        print('Closing processes!')
    for gs in getseqs:
        if gs.is_alive():
            gs.join()

    return seqdict


def id_search(id_rec, id_type='brute', verbose=True):
    # Define the regex functions
    p = [re.compile('(gi)([| :_]+)(\d\d+\.?\d*)(.*)'),  # regex for gi
         re.compile('([AXNYZ][MWRCPGTZ]|ref)([| _:]+)(\d\d+\.?\d*)(.*)'),  # regex for accession
         re.compile('(scaffold)([| _:]+)(\d+\.?\d*)(.*)'),    # regex for scaffolds
         re.compile('(id)([| :_]+)(\d\d+\.?\d*)(.*)'),  # regex for generic ID
         re.compile(':(\d+)-(\d+)'),  # regex for sequence range
         ]

    id_list_ids = []  # Initialized list of IDs
    seq_range = {}  # Initialized dict of sequence ranges

    # Begin search:
    if verbose > 1:
        print('ID File Loaded, performing regex search for identifiers...')
        print('ID Specified as: ', id_type)
    if id_type == 'brute':
        if bool(p[1].findall(id_rec)):
            id_type = 'accession'
            if verbose > 1:
                print(p[1].findall(id_rec))
        elif bool(p[0].findall(id_rec)):
            id_type = 'gi'
            if verbose > 1:
                print(p[0].findall(id_rec))
        elif bool(p[2].findall(id_rec)):
            id_type = 'scaffold'
            if verbose > 1:
                print(p[2].findall(id_rec))
        elif bool(p[3].findall(id_rec)):
            id_type = 'id'
            if verbose > 1:
                print(p[3].findall(id_rec))
        else:
            raise Exception('Couldn\'t identify the id!')
        if verbose > 1:
            print('Brute Force was set, tested strings for all pre-registered IDs. ID was selected as type ',
                  id_type)
    if id_type == 'gi':
        if bool(p[0].findall(id_rec)):
            found_id = True
            if verbose > 1:
                print('Successfully found GI numbers, compiling list!')
            item_parts = p[0].findall(id_rec)
            if verbose > 1:
                print('Item:\t', item_parts)
            id_list_ids.append(item_parts[0][0:3])
            if bool(p[4].findall(id_rec)):
                # Seq_range will be a list of tuples where the second element is the range, and the first
                # is the ID. This way, the function accommodates sequences with a subrange and sequences without a
                # subrange.
                seq_range[''.join(p[0].findall(id_rec)[0][0:3])] = p[4].findall(id_rec)[0]
                if verbose > 1:
                    print('Found sequence delimiters in IDs!')
        else:
            found_id = False
    elif id_type == 'accession':
        if bool(p[1].findall(id_rec)):
            found_id = True
            if verbose > 1:
                print('Successfully found accession numbers, compiling list!')
            item_parts = p[1].findall(id_rec)
            if verbose > 1:
                print('Item:\t', item_parts)
            id_list_ids.append(item_parts[0][0:3])
            if bool(p[4].findall(id_rec)):
                seq_range[''.join(p[1].findall(id_rec)[0][0:3])] = p[4].findall(id_rec)[0]
                if verbose > 1:
                    print('Found sequence delimiters in IDs!')
        else:
            found_id = False
    elif id_type == 'id':
        if bool(p[3].findall(id_rec)):
            found_id = True
            if verbose > 1:
                print('Successfully found ID numbers, compiling list!')
            item_parts = p[3].findall(id_rec)
            if verbose > 1:
                print('Item:\t', item_parts)
            id_list_ids.append(item_parts[0][0:3])
            if bool(p[4].findall(id_rec)):
                seq_range[''.join(p[2].findall(id_rec)[0][0:3])] = p[4].findall(id_rec)[0]
                if verbose > 1:
                    print('Found sequence delimiters in IDs!')
        else:
            found_id = False
    elif id_type == 'scaffold':
        if bool(p[2].findall(id_rec)):
            found_id = True
            if verbose > 1:
                print('Successfully found ID numbers, compiling list!')
            item_parts = p[2].findall(id_rec)
            if verbose > 1:
                print('Item:\t', item_parts)
            id_list_ids.append(item_parts[0][0:3])
            if bool(p[4].findall(id_rec)):
                seq_range[''.join(p[2].findall(id_rec)[0][0:3])] = p[4].findall(id_rec)[0]
                if verbose > 1:
                    print('Found sequence delimiters in IDs!')
        else:
            found_id = False
    else:
        found_id = False
    if found_id:
        return p, id_list_ids, seq_range, id_type
    else:
        raise Exception('ID Error', 'Could not get ID!')


class FetchSeqMP(multiprocessing.Process):
    def __init__(self, id_queue, seq_out_queue, delim, id_type, server, species, source, db, add_length,
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
        self.n_subthreads = n_subthreads

    def run(self):
        while True:
            fs_instance = self.id_queue.get()
            if fs_instance is None:
                self.id_queue.task_done()
                print('All FetchSeqs in Queue completed!')
                break
            try:
                seq_dict, miss_items = fs_instance(passwd=self.passwd, id_type=self.id_type, driver=self.driver,
                                                   user=self.user, host=self.host, db=self.db, delim=self.delim,
                                                   server=self.server, version=self.version, add_length=self.add_length,
                                                   species=self.species, source=self.source, verbose=self.verbose,
                                                   n_threads=self.n_subthreads)
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
                 verbose, add_length):
        # out_file = Path(output_name + '.' + output_type)
        if verbose > 1:
            print('Full header for Entry:')
            print(self.id_rec)

        p, id_list_ids, seq_range, id_type = id_search(self.id_rec, id_type=id_type, verbose=verbose)

        if verbose > 1:
            print('ID list: ')
            for index, ID_item in enumerate(id_list_ids):
                print(index + 1, ': ', ''.join(ID_item))

        # Armed with the ID list, we fetch the sequences from the appropriate source
        if source.lower() == "entrez":
            raise Exception('Not yet implemented, sorry!!!')
        elif source.lower() == "sql":
            if verbose > 1:
                print('Searching for sequences in local SQL db...')
            if verbose > 2:
                print('Please note the sub_databases of server:\n\t', [str(i) for i in server.keys()])
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
                        print('Multiple database versions found!')
                        print(sub_db_list)
                        print('Selecting highest DB')
                    sub_db_name = sorted(sub_db_list, reverse=True)[0]
                if verbose:
                    print('Sub-DB chosen was ', sub_db_name)
            else:
                sub_db_name = ''.join([i[0:3] for i in species.title().split(' ')]) + version
            id_list_search = [''.join(i[0:3]) for i in id_list_ids]
            try:
                seqdict = biosql_get_record_mp(sub_db_name=sub_db_name, passwd=passwd, id_list=id_list_search,
                                               id_type=id_type, driver=driver, user=user,
                                               host=host, db=db, num_proc=n_threads, server=server, verbose=True)
            except Exception as err:
                print('Please note the sub_databases of server:\n\t', [str(i) for i in server.keys()])
                raise Exception('Exception!', err)
            itemsnotfound = [''.join(x) for x in id_list_ids if ''.join(x) not in seqdict.keys()]
            if itemsnotfound:
                if verbose > 1:
                    print('Some items were not found. List of items will be saved to the file '
                          'items_not_found.output')
                    for item in itemsnotfound:
                        print(item)
                        # with open(str(out_file.cwd()) + 'items_not_found.output', 'w') as missingitems:
                        #     missingitems.writelines(itemsnotfound)
            else:
                itemsnotfound = None
            keys = [k for k in seqdict.keys()]
            if verbose > 1:
                print("Sequence Dictionary keys:")
                print(keys)
            if bool(seq_range):
                seqrange_ids = [ids for ids in seq_range.keys()]
                if verbose > 1:
                    print('Sequence Range IDs:')
                    print(seqrange_ids)
                for k in keys:
                    if seqdict[k].id in seqrange_ids:
                        if verbose > 1:
                            print('For sequence {}, found a sequence range!'.format(str(seqdict[k].id)))
                            print('Full length of sequence: {}'.format(len(seqdict[k])))
                        if id_type == 'gi':
                            seq_description_full = p[0].findall(seqdict[k].description)[0]
                        elif id_type == 'accession':
                            seq_description_full = p[1].findall(seqdict[k].description)[0]
                        elif id_type == 'scaffold':
                            seq_description_full = p[2].findall(seqdict[k].description)[0]
                        elif id_type == 'id':
                            seq_description_full = p[3].findall(seqdict[k].description)[0]
                        else:
                            seq_description_full = p[4].findall(seqdict[k].description)[0]
                    else:
                        if verbose > 1:
                            print('No sequence range found, continuing...')
                        continue
                    if (add_length[0] != 0) | (add_length[1] != 0):
                        if verbose > 1:
                            print('Adding {0} steps to the beginning and {1} steps '
                                  'to the end of the sequence!'.format(add_length[0], add_length[1]))
                        if verbose > 2:
                            print(seq_range[k][0], seq_range[k][1], sep='\t')
                            print(-add_length[0], add_length[1], sep='\t')
                            print('_'*len(seq_range[k][0]+seq_range[k][1]))
                        add_length = (-int(add_length[0]), add_length[1])
                        seq_range[k] = tuple(map(lambda x, y: int(x) + y, seq_range[k], add_length))
                        if verbose > 2:
                            print(seq_range[k][0], seq_range[k][1], sep='\t')
                    id_range = ':' + '-'.join([str(i) for i in seq_range[k]])
                    if verbose > 1:
                        print('Sequence range: ',seq_range)
                    if int(seq_range[k][0]) > int(seq_range[k][1]):
                        """tmp_id = seqdict[k].id
                        tmp_name = seqdict[k].name
                        tmp_desc = seqdict[k].description
                        tmp_dbxrefs = seqdict[k].dbxrefs
                        tmp_feat = seqdict[k].features
                        tmp_annotations = seqdict[k].annotations
                        tmp_let_anno = seqdict[k].letter_annotations"""
                        if verbose > 1:
                            print('Starting seq_range is larger than ending seq_range - sequence is '
                                  'in the (-) direction! Reversing...')
                        seqdict[k].seq = seqdict[k][
                                         int(seq_range[k][1]):int(seq_range[k][0])].seq.reverse_complement()

                    else:
                        seqdict[k] = seqdict[k][int(seq_range[k][0]):int(seq_range[k][1])]
                    if verbose > 1:
                        print('Seq_description_full: ', seq_description_full)
                        print('id_range: ', id_range[1:])
                    if int(seq_range[k][0]) > int(seq_range[k][1]):
                        seqdict[k].description = ''.join(seq_description_full[0:3]) + id_range + '(-)' + \
                                                 str(seq_description_full[3])
                    else:
                        seqdict[k].description = ''.join(seq_description_full[0:3]) + id_range + '(+)' + \
                                                 str(seq_description_full[3])
                    if verbose > 1:
                        print('Sequence Description: \n\t', seqdict[k].description)
                    seqdict[k].id += id_range
                    if verbose > 1:
                        print('Sequence ID: \n\t', seqdict[k].id)
                        if id_range:
                            print('Length of subsequence with range {0}: {1}'.format(id_range, len(seqdict[k])))
            if verbose > 1:
                print('Sequence Record post-processing, to be saved:')
                print(seqdict)
            if verbose > 1:
                print('Finished getting sequence!')
            return seqdict, itemsnotfound
        elif source == "fasta":  # Note: anecdotally, this doesn't run terribly fast - try to avoid.
            # TODO: have this work like SQL does.
            
            seqdict = SeqIO.index(db, source,
                                  key_function=lambda identifier: p[0].search(
                                      p[2].search(identifier).group()).group())
            itemsnotfound = [x for x in id_list_ids if x not in seqdict.keys()]
            if itemsnotfound:
                if verbose > 1:
                    print('Some items were not found. List of items will be saved to the file items_not_found.output')
                    for item in itemsnotfound:
                        print(item)
                        # with open(str(out_file.cwd()) + 'items_not_found.output', 'w') as missingitems:
                        #    missingitems.writelines(itemsnotfound)
            else:
                itemsnotfound = None
            keys = [k for k in seqdict.keys()]
            if verbose > 1:
                print("Sequence Dictionary keys:")
                print(keys)
            if bool(seq_range):
                seqrange_ids = [ids for ids in seq_range.keys()]
                if verbose > 1:
                    print('Sequence Range IDs:')
                    print(seqrange_ids)
                for k in keys:
                    if seqdict[k].id in seqrange_ids:
                        if verbose > 1:
                            print('For sequence {}, found a sequence range!'.format(str(seqdict[k].id)))
                            print('\tFull length of sequence: {}'.format(len(seqdict[k])))
                        if id_type == 'gi':
                            seq_description_full = p[0].findall(seqdict[k].description)[0]
                        elif id_type == 'accession':
                            seq_description_full = p[1].findall(seqdict[k].description)[0]
                        elif id_type == 'id':
                            seq_description_full = p[2].findall(seqdict[k].description)[0]
                        else:
                            seq_description_full = p[4].findall(seqdict[k].description)[0]
                    if verbose > 1:
                        print(int(seq_range[k][0]))
                        print(int(seq_range[k][1]))
                    id_range = ':' + '-'.join(seq_range[k])
                    seqdict[k] = seqdict[k][int(seq_range[k][0]):int(seq_range[k][1])]
                    seqdict[k].description = ''.join(seq_description_full[0:3]) + id_range + \
                                             str(seq_description_full[3])
                    seqdict[k].id += id_range
                    if verbose > 1:
                        print('\tLength of subsequence with range{0}: {1}'.format(id_range, len(seqdict[k])))
                    else:
                        if verbose > 1:
                            print('No sequence range found, continuing...')
                        continue
            if verbose > 1:
                print('Done!')
            return seqdict, itemsnotfound
            # SeqIO.write([seqdict[key] for key in seqdict.keys()], str(out_file), output_type)
        else:
            raise Exception('Not a valid database source!')


def fetchseqMP(ids, species, write=False, output_name='', delim='\t', id_type='brute', server=None, source="SQL",
               db="bioseqdb", host='localhost', driver='psycopg2', version='1.0', user='postgres', passwd='', email='',
               batch_size=50, output_type="fasta", verbose=1, n_threads=1, n_subthreads=1, add_length=(0,0)):

    if isgenerator(ids):
        if verbose > 1:
            print('Received generator!')
    elif isinstance(ids, list):
        if verbose > 1:
            print('Received list!')
    else:
        if verbose > 1:
            print('Reading ID File... ')
        with ids.open('w') as in_handle:
            id_prelist = [line.strip() for line in in_handle]  # list of each line in the file
            print('Done!')
        ids = [id_item for id_item in filter(None, id_prelist) if id_item]
        if not id_prelist or id_prelist is None:
            if verbose:
                print('id_prelist is empty!')
            return 'None'

    if verbose > 1:
        print('Readied ids!')

    id_list = multiprocessing.JoinableQueue()
    results = multiprocessing.Queue()

    if server is None:
        try:
            if verbose > 1:
                print('No server received, opening server...')
            server = BioSeqDatabase.open_database(driver=driver, user=user, passwd=passwd, host=host, db=db)
            if verbose > 1:
                print('Done!')
        except:
            if verbose > 1:
                print('FAILED!')
            raise
    else:
        if verbose > 1:
            print('Received server handle:')
            print(server)
        if verbose > 2:
            print('Please note the sub_databases of server:\n\t', [str(i) for i in server.keys()])
    if verbose > 1:
        print('Creating FecSeq Processes...')
    fs_instances = [FetchSeqMP(id_queue=id_list, seq_out_queue=results,
                               delim=delim, id_type=id_type, server=server, species=species, source=source, db=db,
                               host=host, driver=driver, version=version, user=user, passwd=passwd, email=email,
                               output_type=output_type, batch_size=batch_size, verbose=verbose,
                               n_subthreads=n_subthreads, add_length=add_length)
                    for i in range(n_threads)]
    if verbose > 1:
        print('Done! Starting processes...')
    for fs in fs_instances:
        fs.start()
    if verbose > 1:
        print('Done!')
        print('Assigning FetchSeq records to queue... ')
    for id_rec in ids:
        id_list.put(FetchSeq(id_rec=id_rec))
    for fs in fs_instances:
        id_list.put(None)
    if verbose > 1:
        print('Done!')
    output_dict = dict()
    missing_items_list = list()
    if verbose > 1:
        print('Getting sequences from processes... ')
    print('--------------', n_threads)
    n_jobs = len(ids)
    # for i in range(len(ids)):
    while n_jobs:
        seq, missing = results.get()
        output_dict.update(seq)
        missing_items_list.append(missing)
        n_jobs -= 1
    if verbose > 1:
        print('Done! Finished fetching sequences!')
        print('Closing processes!')
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
                 perc_score, outfolder,
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

    def run(self):  # The meat of the script
        master_out = self.outfolder.joinpath('Proc-{0}.log'.format(self.name)).absolute()
        if self.verbose >1:
            print(master_out)
        master_out_handle = master_out.open('w')
        old_sysout = sys.__stdout__
        sys.stdout = master_out_handle
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
                                     perc_ident=self.perc_ident,
                                     perc_length=self.perc_length, megablast=self.megablast, email=self.email,
                                     id_type=self.id_type,
                                     fw_source=self.fw_source, fw_id_db=self.fw_id_db, fetch_batch_size=self.batch_size,
                                     passwd=self.passwd,
                                     fw_id_db_version=self.fw_id_db_version, verbose=self.verbose,
                                     n_threads=self.n_threads,
                                     host=self.host,
                                     user=self.user, driver=self.driver,
                                     fw_blast_kwargs=self.fw_blast_kwargs, rv_blast_kwargs=self.rv_blast_kwargs,
                                     proc_id=self.name, write_intermediates=self.write_intermediates,
                                     )
                self.rb_queue.task_done()
                self.rb_results_queue.put(output)
            except Exception as err:
                print('Woah! Something went wrong! Aborting!')
                print('Here\'s the error:\n', err)
                self.rb_results_queue.put(dict(bla='bla'))
        sys.stdout = old_sysout
        master_out_handle.close()
        return


class RecBlast(object):
    def __init__(self, seq_record, target_species):
        self.seq_record = seq_record
        transtab = str.maketrans('!@#$%^&*();:.,\'\"/\\?<>|[]{}-=+', '_____________________________')
        self.seq_record.name = self.seq_record.name.translate(transtab)
        self.target_species = target_species

    def __call__(self, fw_blast_db, infile_type, output_type, BLASTDB,
                 query_species, blast_type_1, blast_type_2, local_blast_1, local_blast_2, rv_blast_db, expect,
                 perc_score,
                 perc_ident, perc_length, megablast, email, id_type, fw_source, fw_id_db, fetch_batch_size, passwd,
                 host, user, driver, fw_id_db_version, verbose, n_threads, fw_blast_kwargs, rv_blast_kwargs,
                 write_intermediates, proc_id):
        # Simple shunt to minimize having to rewrite code.
        target_species = self.target_species
        # Creating the RecBlast Container
        rc_container_full = RecBlastContainer(target_species=target_species, query_record=self.seq_record)
        # Shorthand reference for full container
        rc_container = rc_container_full[target_species][self.seq_record.name]
        rc_container.update(proc_id=proc_id)
        # Indentation level for print()
        level = 0
        if verbose:
            print('\t' * level, '[Proc: {0}] [Seq.Name: {1}] [Target: {2}]'.format(proc_id,
                                                                                   self.seq_record.name,
                                                                                   target_species))
        level += 1
        if verbose > 1:
            print('\t' * level, '\t Creating handles for intermediary outputs...')
        # Handle for output paths:
        output_paths = rc_container['output_paths']
        output_paths.update(
            forward_blast_output=Path("{0}_recblast_out".format(target_species).replace(' ', '_') + '/' +
                                      "{0}_{1}_tmp".format(blast_type_1,
                                                           self.seq_record.name
                                                           ).replace(' ', '_') + '/' +
                                      "{0}_{1}_{2}_to_{3}.xml".format(blast_type_1,
                                                                      self.seq_record.name,
                                                                      query_species,
                                                                      target_species
                                                                      ).replace(' ', '_')),
            forward_id_score_output=Path("{0}_recblast_out".format(target_species).replace(' ', '_') + '/' +
                                         "{0}_{1}_tmp".format(blast_type_1,
                                                              self.seq_record.name
                                                              ).replace(' ', '_') + '/' +
                                         "{0}_{1}_{2}_to_{3}.ID_Scores".format(blast_type_1,
                                                                               self.seq_record.name,
                                                                               query_species,
                                                                               target_species
                                                                               ).replace(' ', '_')),
            recblast_output_unanno=Path("{0}_recblast_out".format(target_species).replace(' ', '_') + '/' +
                                        "{0}_{1}_tmp".format(blast_type_1,
                                                             self.seq_record.name
                                                             ).replace(' ', '_') + '/' +
                                        "unannotated_{0}_{1}.fasta".format(blast_type_1,
                                                                           self.seq_record.name
                                                                           ).replace(' ', '_')),
            recblast_output=Path("{0}_recblast_out".format(target_species).replace(' ', '_') + '/' +
                                 "{0}_{1}.{2}".format(blast_type_1,
                                                      self.seq_record.name,
                                                      output_type
                                                      ).replace(' ', '_')),
            blast_nohits=Path("{0}_recblast_out".format(target_species).replace(' ', '_') + '/' +
                              "{0}_{1}.no-hits".format(blast_type_1,
                                                       self.seq_record.name
                                                       ).replace(' ', '_'))
        )
        """if write_intermediates:
            for path_type, path in output_paths.items():
                try:
                    path.absolute().parent.mkdir(parents=True)
                    if verbose > 1:
                        print('\t' * level, '\t\tCreated directory \"{}\"!'.format(str(path.absolute().parent)))
                except FileExistsError:
                    if verbose > 1:
                        print('\t' * level, '\t\tDirectory \"{}\" already exists! Continuing!'.format(
                            str(path.absolute().parent)))
        """
        if blast_type_1 in ['blat', 'tblat']:
            try:
                fw_blast_db = fw_blast_db[target_species]
            except KeyError:
                print('No port found for species ', target_species)
                return rc_container_full
        if fw_blast_db == 'auto':
            try:
                fw_blast_db = get_searchdb(search_type=blast_type_1, species=target_species, db=BLASTDB,
                                           verbose=verbose, level=level + 1)
            except Exception as Err:
                print(Err)
                return rc_container_full
        fwblast = rc_container['forward_blast']

        if fw_blast_db == 'skip':
            if verbose:
                print('\t' * level, "\tSkipping Forward Blast! Using local file instead: ")
                print('\t' * (level + 1), '\tOpening Forward blast output located at ',
                      str(output_paths['forward_blast_output'].absolute()))
            with output_paths['forward_blast_output'].open("r") as forward_blasthits:
                fwblastrecord = NCBIXML.read(forward_blasthits)
        else:
            if verbose:
                print('\t' * level, "\tPerforming forward BLAST for {}... ".format(self.seq_record.name))
            level += 1
            try:
                if fw_blast_kwargs:
                    if verbose:
                        for key, item in fw_blast_kwargs.items():
                            print('\t' * level, '{0}\t=\t{1}'.format(key, item))
                    fwblastrecord, blast_err = blast(seq_record=self.seq_record, target_species=target_species,
                                                     database=fw_blast_db, query_species=query_species,
                                                     filetype=infile_type, BLASTDB=BLASTDB,
                                                     blast_type=blast_type_1, local_blast=local_blast_1, expect=expect,
                                                     megablast=megablast, n_threads=n_threads,
                                                     blastoutput_custom=str(output_paths['forward_blast_output']),
                                                     perc_ident=perc_ident, verbose=verbose, write=write_intermediates,
                                                     **fw_blast_kwargs)
                else:
                    fwblastrecord, blast_err = blast(seq_record=self.seq_record, target_species=target_species,
                                                     database=fw_blast_db, query_species=query_species,
                                                     filetype=infile_type, BLASTDB=BLASTDB,
                                                     blast_type=blast_type_1, local_blast=local_blast_1, expect=expect,
                                                     megablast=megablast, n_threads=n_threads,
                                                     blastoutput_custom=str(output_paths['forward_blast_output']),
                                                     perc_ident=perc_ident, verbose=verbose, write=write_intermediates)
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
        level -= 1
        if verbose:
            print('\t' * level, 'Forward blast done!')
            print('\t' * level, 'Culling Results based on given criteria...')
        try:
            f_id_ranked = id_ranker(fwblastrecord, perc_score=perc_score, expect=expect, perc_length=perc_length,
                                    verbose=verbose)
        except Exception as err:
            print('WARNING! UNCATCHED ERROR IN ID_RANKER!')
            print(err)
            return rc_container_full
        if verbose:
            print('\t' * (level + 1), 'Done!')
        f_id_out_list = ['{0}\t{1}\t{2}\n'.format(id_i[0], id_i[1], id_i[2]) for id_i in f_id_ranked]
        fw_ids = rc_container['forward_ids']
        fw_ids['ids'] = f_id_ranked
        fw_ids['pretty_ids'] = f_id_out_list
        """
        if write_intermediates:
            if verbose:
                print('\t' * (level + 1), 'Writing Forward ID Hits to output!')
            with output_paths['forward_id_score_output'].open('w') as fidout:
                fidout.writelines(f_id_out_list)
        """
        if not f_id_out_list:
            print('Forward Blast yielded no hits, continuing to next sequence!')
            return rc_container_full
        try:
            server = BioSeqDatabase.open_database(driver=driver, user=user, passwd=passwd,
                                                  host=host, db=fw_id_db)
            # Note: BioSQL is NOT thread-safe! Throws tons of errors if executed with more than one thread!
            seq_dict, missing_items = fetchseqMP(ids=f_id_out_list, species=target_species, delim='\t',
                                                 id_type='brute', server=server, source=fw_source,
                                                 db=fw_id_db, host=host, driver=driver,
                                                 version=fw_id_db_version, user=user,
                                                 passwd=passwd, email=email, batch_size=fetch_batch_size,
                                                 output_type=output_type, output_name='', write=write_intermediates,
                                                 verbose=verbose, n_threads=1,
                                                 n_subthreads=1)
            if verbose:
                print('\t' * level, 'Done with fetching!')
        except IndexError:
            print('WARNING! FETCHSEQ FAILED! SKIPPING THIS SEQUENCE!')
            return rc_container_full

        if missing_items:
            fw_ids['missing_items'] = missing_items
            print('\t' * level, 'Items were missing!')
            for i in missing_items:
                print('\t' * level, i)

        recblast_sequence = []
        if seq_dict:
            for item in f_id_ranked:
                junk1, id_list_ids, seq_range, junk2 = id_search(''.join([str(i) for i in item]), verbose=verbose)
                if verbose > 3:
                    print('\t' * level, 'f-id-ranked: ', id_list_ids)
                # Todo: find a more efficient way to do this:
                for otherkey, otheritem in seq_dict.items():
                    print('\t' * level, id_list_ids[0])
                    print('\t' * level, otheritem.description)
                    if ''.join([str(i) for i in id_list_ids[0]]) in otheritem.description:
                        recblast_sequence.append(otheritem)
        rc_container['recblast_unanno'] = recblast_sequence
        if verbose:
            print('\t' * level, 'Preparing for Reverse BLAST...')
        for index, entry_record in enumerate(recblast_sequence):
            if verbose:
                print('\t' * level, "Entry {} in unannotated RecBlast Hits:\n".format(entry_record.name))
                for item in [entry_record.id, entry_record.description, entry_record.seq]:
                    print('\t' * level, '\t', item)
            level += 1
            output_paths['reverse_blast_output'].append(
                Path("{0}_recblast_out".format(target_species).replace(' ', '_') +
                     '/' +
                     "{0}_{1}_tmp".format(blast_type_2,
                                          self.seq_record.name
                                          ).replace(' ', '_') +
                     '/' +
                     "{0}_{1}_{3}_to_{2}_{4}.xml".format(blast_type_2,
                                                         self.seq_record.name,
                                                         query_species,
                                                         target_species,
                                                         entry_record.name
                                                         ).replace(' ', '_')
                     ))
            if verbose:
                print('\t' * level, 'Performing Reverse Blast:')
            if blast_type_1 in ['blat', 'tblat']:
                try:
                    rv_blast_db = rv_blast_db[query_species]
                except KeyError:
                    print('No port found for species ', query_species)
                    return rc_container_full
            if (rv_blast_db == 'auto') | (rv_blast_db == 'auto-transcript'):
                try:
                    rv_blast_db = get_searchdb(search_type=blast_type_2, species=target_species, db=BLASTDB,
                                               verbose=verbose, level=level + 1)
                except Exception as Err:
                    print(Err)
                    return rc_container_full
            rv_blast = rc_container['reverse_blast']
            if rv_blast_db == 'skip':
                pass
            elif rv_blast_db == 'stop':
                if verbose:
                    print('\t' * level, 'Not performing reverse blast!')
                return rc_container_full
            else:
                try:
                    if rv_blast_kwargs:
                        if verbose:
                            for key, item in rv_blast_kwargs.items():
                                print('\t' * level, key, '\t', item)
                        rvblastrecord, blast_err = blast(seq_record=entry_record, target_species=target_species,
                                                         database=rv_blast_db, query_species=query_species,
                                                         filetype=infile_type, BLASTDB=BLASTDB,
                                                         blast_type=blast_type_2, local_blast=local_blast_2,
                                                         expect=expect, n_threads=n_threads,
                                                         query_length=len(self.seq_record),
                                                         megablast=megablast,
                                                         blastoutput_custom=str(
                                                                           output_paths['reverse_blast_output'][index]),
                                                         perc_ident=perc_ident, verbose=verbose,
                                                         write=write_intermediates,
                                                         **rv_blast_kwargs)
                    else:
                        rvblastrecord, blast_err = blast(seq_record=entry_record, target_species=target_species,
                                                         database=rv_blast_db, query_species=query_species,
                                                         filetype=infile_type, BLASTDB=BLASTDB,
                                                         blast_type=blast_type_2, local_blast=local_blast_2,
                                                         expect=expect, n_threads=n_threads,
                                                         query_length=len(self.seq_record),
                                                         megablast=megablast, blastoutput_custom=str(
                                output_paths['reverse_blast_output'][index]),
                                                         perc_ident=perc_ident, verbose=verbose,
                                                         write=write_intermediates)
                except Exception as err:
                    print('Warning! Uncaught exception!')
                    print(err)
                    return rc_container_full
                if blast_err:
                    print('Received blast error file back from Reverse BLAST!')
                    rv_blast['blast_errors'] = blast_err
                    return rc_container_full
            rv_blast['blast_results'] = rvblastrecord
            level -= 1
            if verbose:
                print('\t' * level, 'Done with Reverse Blast!')
                print('\t' * level, 'Culling results using given criteria...')
            reverse_hits = id_ranker(rvblastrecord, perc_score=perc_score, perc_length=perc_length,
                                     expect=expect, verbose=verbose)
            print(reverse_hits)
            reverse_blast_annotations = ['\t |[ {0} {1} ({2}) ]|'.format(anno[0], anno[1], anno[2]) for anno in
                                         reverse_hits]
            if not reverse_blast_annotations:
                print('\t' * level, 'No Reverse Blast Hits were found for this hit!')
                print('\t' * level, 'Continuing to next Sequence!')
            else:
                if verbose > 1:
                    print('\t' * level, 'Done. Annotating RecBlast Hits:')
            rv_ids = rc_container['reverse_ids']
            rv_ids['ids'] = reverse_blast_annotations
            entry_record.description += '|-|' + ''.join(reverse_blast_annotations)
            if verbose > 3:
                print('\t' * level, entry_record)
        rc_container['recblast_results'] = recblast_sequence
        print('\t' * level, proc_id, self.seq_record.name, sep='\t')
        return rc_container_full
        # def run(self):


class SearchMaster(multiprocessing.Process):
    # TODO: Configure this to search using either BLAST or BLAT
    # TODO: Use SearchIO to have this read results
    def __init__(self, search_type='blast'):
        multiprocessing.Process.__init__(self)
        self.search_type = search_type

    def run(self):  # The meat of the script
        while True:
            (search_proc_id, search_instance) = self.search_queue.get()
            if search_instance is None:
                self.search_queue.task_done()
                break
            output = search_instance()
            self.search_queue.task_done()
            self.result_queue.put(output)
        return


class SearchInstance(object):
    def __init__(self):
        pass

    def __call__(self):
        return


def recblastMP(seqfile, target_species, fw_blast_db='auto', rv_blast_db='auto-transcript', infile_type='fasta',
               output_type='fasta',
               host='localhost', user='postgres', driver='psycopg2',
               query_species='Homo sapiens', blast_type_1='blastn', blast_type_2='blastn', local_blast_1=False,
               local_blast_2=False,
               expect=10, perc_score=0.5, perc_ident=50, perc_length=0.5, megablast=True,
               email='',
               id_type='brute', fw_source='sql', fw_id_db='bioseqdb', fetch_batch_size=50, passwd='',
               fw_id_db_version='auto', BLASTDB='/usr/db/blastdb',
               verbose='v', max_n_processes='auto', n_threads=2, write_intermediates=False, write_final=True,
               fw_blast_kwargs=None, rv_blast_kwargs=None):

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

    if verbose > 1:
        print('Creating queues... ')
    rb_queue = multiprocessing.JoinableQueue()
    rb_results = multiprocessing.Queue()
    if verbose > 1:
        print('Done!')
    # Check Seqfile to make sure its real
    if verbose >= 1:
        print('Loading SeqFile records... ')
    seqfile_path = ''
    if isinstance(seqfile, str):
        seqfile_path = Path(seqfile)
    elif isinstance(seqfile, Path):
        seqfile_path = seqfile

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
            rec_handle = [i for i in SeqIO.parse(str(seqfile_path.absolute()), output_type)]
    except FileNotFoundError:
        raise
    else:
        if verbose >= 1:
            print('Done!')

    # Calculation of processes to run
    if verbose > 1:
        print('Automatically calculating n_processes:')
    if isinstance(target_species, list):
        n_species = len(target_species)
    else:
        n_species = 1
    n_jobs = len(rec_handle) * n_species
    if isinstance(max_n_processes, str):
        max_n_processes = multiprocessing.cpu_count()
        print('CPU count: ', max_n_processes)
    if n_jobs * n_threads > max_n_processes:
        print('Optimal number of processes would be above max limit, using limit instead!')
        n_processes = max_n_processes / n_threads
    else:
        n_processes = n_jobs
        if verbose:
            print('Number of processes to be made: ', n_processes)

    # Get date-time and make output folder
    date_str = dt.now().strftime('%y-%m-%d_%I-%M-%p')
    outfolder = Path('./RecBlast_output/{0}/'.format(date_str))
    try:
        outfolder.mkdir(parents=True)
    except FileExistsError:
        pass


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
                                             perc_ident=perc_ident,
                                             perc_length=perc_length, megablast=megablast, email=email, id_type=id_type,
                                             fw_source=fw_source, fw_id_db=fw_id_db, fetch_batch_size=fetch_batch_size,
                                             passwd=passwd,
                                             fw_id_db_version=fw_id_db_version, verbose=verbose, n_threads=n_threads,
                                             host=host, outfolder=outfolder,
                                             user=user, driver=driver, write_intermediates=write_intermediates,
                                             fw_blast_kwargs=fw_blast_kwargs, rv_blast_kwargs=rv_blast_kwargs)
                           for i in range(n_processes)]
    for rcb in rec_blast_instances:
        rcb.start()
    progbar = ProgressBar(n_jobs, fmt=ProgressBar.FULL)
    if verbose:
        progbar()
    # Check if BLAT, then if so, make sure that fw_blast_db is a dictionary with a port for each species:
    if blast_type_1.lower() in ['blat', 'tblat']:
        assert isinstance(fw_blast_db, dict), "For BLAT searches, fw_blast_db must be a dictionary with " \
                                              "valid species-port key pairs"
    if blast_type_2.lower() in ['blat', 'tblat']:
        assert isinstance(rv_blast_db, dict), "For BLAT searches, rv_blast_db must be a dictionary with " \
                                              "valid species-port key pairs"

    if isinstance(target_species, str):
        if target_species.lower() == 'all':
            if verbose:
                print('Target species set to all. Searching server for all available species:')
            target_species = biosql_get_sub_db_names(passwd=passwd, db=fw_id_db, driver=driver, user=user, host=host)
            if verbose:
                print(target_species)

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

    for rcb in rec_blast_instances:
        rb_queue.put(None)
    recblast_out = list()
    while n_jobs:
        try:
            recblast_out.append(rb_results.get())
        except TypeError:
            print('Behold the stupid TypeError Bug!')
            return recblast_out
        else:
            progbar.current += 1
            progbar()
            n_jobs -= 1

    for rcb in rec_blast_instances:
        if rcb.is_alive():
            rcb.join()
    '''
    if 'blat' in blast_type_1.lower():
        if isinstance(target_species, list):
            for species in target_species:
                blat_db = get_searchdb(search_type=blast_type_1, species=species, db=BLASTDB, verbose=bool(verbose),
                                       level=1)
                blat_server(blat_db, 'stop', **fw_blast_kwargs)

        else:
            blat_db = get_searchdb(search_type=blast_type_1, species=target_species, db=BLASTDB,
                                   verbose=bool(verbose), level=1)
            blat_server(blat_db, 'stop', **fw_blast_kwargs)
    '''

    if write_final:
        recblast_write(recblast_out, verbose=verbose, outfolder=outfolder)
    return recblast_out


def recblast_write(rc_container, verbose=True, outfolder=None):
    if outfolder is None:
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


def cleanup_fasta_input(handle, filetype='fasta', write=True):
    from Bio import SeqIO
    oldlist = [i for i in SeqIO.parse(handle, filetype)]
    names = set([i.name for i in oldlist])
    newlist = list()
    for name in names:
        x = [i for i in oldlist if i.name == str(name) and 'Sequenceunavailable' not in i.seq]
        for j in x:
            j.name += '_' + str(j.description).split('|')[2]
            newlist += x
    if write:
        with open(handle + '.clean', 'w') as outf:
            SeqIO.write(newlist, outf, filetype)
    return newlist


def count_dups(recblast_out):
    master_list = list()
    pat = [re.compile('(gi)([| :_]+)(\d\d+\.?\d*)(.*)'),  # regex for gi
           re.compile('([AXNYZ][MWRCPGTZ]|ref)([| _:]+)(\d\d+\.?\d*)(.*)'),  # regex for accession
           re.compile('(id)([| :_]+)(\d\d+\.?\d*)(.*)'),  # regex for generic ID
           re.compile(':(\d+)-(\d+)'),  # regex for sequence range
           re.compile('\[\|(\d*)\|\]')]  # regex for items in annotation
    for rc in recblast_out:
        for species, rc_spec_rec in rc.items():
            for gene, rc_rec in rc_spec_rec.items():
                rc_out = rc_rec['recblast_output']
                for record in rc_out:
                    target_id, annotations = record.description.split('[+]')
                    for p in pat:
                        id_lst = p.findall(target_id)
                        if id_lst:
                            try:
                                master_list += ''.join(id_lst[2:3])
                            except IndexError:
                                try:
                                    master_list += ''.join(id_lst[2])
                                except IndexError:
                                    raise
                            break
