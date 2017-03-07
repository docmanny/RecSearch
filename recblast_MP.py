import multiprocessing
from operator import itemgetter

from Bio import SeqIO
from Bio.Blast import Record as BioBlastRecord
from BioSQL import BioSeqDatabase


def id_ranker(blastrecord, perc_score, expect, perc_length, hsp_cumu_score=True, align_scorelist=list(),
              indent='', verbose=True):
    assert isinstance(blastrecord, BioBlastRecord), "Parameter 'blastrecord' must be a BioBlastRecord!!!"
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


"""
def print(*objects, sep=' ', end='\n', file=sys.stdout, flush=False, lock=None):
    import builtins
    if lock is None:
        builtins.print(*objects, sep=' ', end='\n', file=sys.stdout, flush=False)
    else:
        with lock:
            builtins.print(*objects, sep=' ', end='\n', file=sys.stdout, flush=False)
"""


class RecBlastContainer(dict):
    def __init__(self, *args, **kwargs):
        super(RecBlastContainer, self).__init__(*args, **kwargs)
        for arg in args:
            if isinstance(arg, dict):
                for k, v in arg.items():
                    k.replace('.', '_')
                    if type(v) is dict:
                        v = RecBlastContainer(v)
                    else:
                        pass
                    self[k] = v

        if kwargs:
            for k, v in kwargs.items():
                k.replace('.', '_')
                if k == 'proc_id':
                    if isinstance(v, int):
                        v = 'x' + str(v)
                    self.proc_id = v
                elif k == 'seq_record':
                    self.seq_record = v
                if isinstance(v, dict):
                    v = RecBlastContainer(v)
                    self[k] = v
                else:
                    self[k] = v

    def __getattr__(self, attr):
        try:
            self[attr]
        except KeyError:
            raise
        except AssertionError:
            raise
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

    def __add__(self, other):
        assert isinstance(other, RecBlastContainer), "Other is not a RecBlastContainer!"
        try:
            self.proc_id
        except KeyError:
            raise Exception('Attribute proc_id must be defined!!!')
        try:
            self.seq_record
        except KeyError:
            raise Exception('Attribute seq_record must be defined!!!')
        if self.proc_id == other.proc_id:
            return RecBlastContainer({'proc_id': self.proc_id, self.seq_record.name.replace('.', '_'): self,
                                      other.seq_record.name.replace('.', '_'): other})
        else:
            return RecBlastContainer({str(self.proc_id): self, str(other.proc_id): other},
                                     proc_id=[self.proc_id, other.proc_id])

    def __str__(self, indent=''):
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


def blast(seq_record, target_species, database, query_species="Homo sapiens", filetype="fasta", blast_type='blastn',
          local_blast=False, expect=0.005, megablast=True, blastoutput_custom="", perc_ident=75,
          verbose=True, n_threads=1, use_index=True, write=False, BLASTDB='/usr/db/blastdb/', **kwargs):
    from pathlib import Path
    
    from Bio.Blast import NCBIWWW
    if isinstance(seq_record, SeqIO.SeqRecord):
        pass
    else:
        seq_record = SeqIO.read(seq_record, filetype)
    args = dict()
    if verbose:
        print("Now starting BLAST...")
    if kwargs:
        args.update(**kwargs)
    # Begin by opening recblast_out, and then start with the primary BLAST
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

    if local_blast:
        import subprocess
        args.update({'-db': database, '-evalue': expect,
                     '-outfmt': '5',
                     '-num_threads': n_threads})
        if blast_type == 'blastn':
            if megablast:
                args['-task'] = 'megablast'
            # args['-use_index'] = use_index
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
                blast_handle = (subprocess.check_output([str(i) for i in blast_cline], input=seq_record.format('fasta'),
                                                        universal_newlines=True, cwd=BLASTDB))
                if isinstance(blast_handle, str):
                    blast_result = blast_handle
                    blast_err = None
                else:
                    blast_result, blast_err = blast_handle

                # = blast_call.communicate(input=)
                if verbose > 3:
                    print('Blast Result: ', blast_result)
            except subprocess.CalledProcessError:
                raise
        else:
            raise Exception("Invalid blast choice!")
    else:
        args.update(dict(program=str(blast_type), database=str(database), sequence=seq_record.format('fasta'),
                         entrez_query='"{}"[ORGN]'.format(target_species), expect=expect, perc_ident=perc_ident))
        if megablast & (blast_type == 'blastn'):
            args['megablast'] = True
        if verbose:
            print('Submitting Remote BLAST! Options passed:')
            for k, v in args.items():
                print('\t {0}\t=\t{1}'.format(k, v))
        blast_handle = NCBIWWW.qblast(**args)
        blast_result = blast_handle.read()
        blast_err = None
    if verbose:
        print('Done with Blast!')
    if write:
        with blastoutput_custom.open("w") as fxml:
            fxml.write(blast_result)
    else:
        return blast_result, blast_err


def biosql_seq_lookup_cascade(dtbase, sub_db_name, id_type, identifier, verbose=False):
    seqrec = ''
    try_get_id = True
    if try_get_id:
        try:
            if verbose:
                print("\t\tNow searching database {0} for {1}: {2}".format(sub_db_name, id_type, identifier))
            seqrec = biosql_DBSeqRecord_to_SeqRecord(dtbase.lookup(**{id_type: identifier}))
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
            seqrec = biosql_DBSeqRecord_to_SeqRecord(dtbase.lookup(**{id_type: identifier_sans_subnumber}))
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
        except IndexError as err3:
            if verbose:
                print("WARNING: Still couldn't find {0} using name search: \n full error: {1}".format(identifier, err3))

    if try_get_id:
        try:
            id_type = input('Last shot, chose an ID type: '
                            '[accession, primary_id, gi, version, display_id, name]')
            if id_type == 'exit':
                exit(exit(), 'Script ended!')
            seqrec = biosql_DBSeqRecord_to_SeqRecord(dtbase.lookup(**{id_type: identifier}))
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
        dtbase = server[sub_db_name]
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


class FetchSeqMP(multiprocessing.Process):
    def __init__(self, id_queue, seq_out_queue, missing_items_queue, delim, id_type, server, species, source, db,
                 host, driver, version, user, passwd, email, output_type, batch_size, verbose, lock, n_subthreads):
        multiprocessing.Process.__init__(self)
        # Queues
        self.id_queue = id_queue
        self.seq_out_queue = seq_out_queue
        self.missing_items_queue = missing_items_queue
        # ID parsing
        self.delim = delim
        self.id_type = id_type
        # Search Items
        self.server = server  # Declaring the server here lets me pass it on to everything else by inheritance
        self.species = species
        self.source = source
        self.db = db
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
        self.lock = lock
        self.n_subthreads = n_subthreads

    def run(self):  # runs FetSeq with everything EXCEPT the
        while True:
            fs_instance = self.id_queue.get()
            if fs_instance is None:
                self.id_queue.task_done()
                print('All FetchSeqs in Queue completed!')
                break
            seq_dict, miss_items = fs_instance(passwd=self.passwd, id_type=self.id_type, driver=self.driver,
                                               user=self.user, host=self.host, db=self.db, delim=self.delim,
                                               server=self.server, version=self.version,
                                               species=self.species, source=self.source,
                                               lock=self.lock, verbose=self.verbose, n_threads=self.n_subthreads)
            if self.verbose:
                print('*****************************************1')
            self.id_queue.task_done()
            if self.verbose:
                print('*****************************************2')
            self.seq_out_queue.put(seq_dict)
            if self.verbose:
                print('*****************************************3')
            self.missing_items_queue.put(miss_items)
            if self.verbose:
                print('*****************************************4')
        return


class FetchSeq(object):  # The meat of the script
    def __init__(self, id_rec):
        self.id_rec = id_rec

    def __call__(self, delim, species, version, source, passwd, id_type, driver, user, host, db, n_threads, server,
                 lock, verbose):
        # out_file = Path(output_name + '.' + output_type)
        import re
        if verbose > 1:
            with lock:
                print('Full header for Entry:')
                print(self.id_rec)
        # id_list = [str(item.split(delim)) for item in self.id_rec]  # Breaks the tab sep in the lines into strings
        # Define the regex functions
        p = [re.compile('(gi)([| :_]+)(\d\d+\.?\d*)(.*)'),  # regex for gi
             re.compile('([AXNYZ][MWRCPGTZ]|ref)([| _:]+)(\d\d+\.?\d*)(.*)'),  # regex for accession
             re.compile('(id)([| :_]+)(\d\d+\.?\d*)(.*)'),  # regex for generic ID
             re.compile(':(\d+)-(\d+)'),  # regex for sequence range
             ]
        id_list_ids = []  # Initialized list of IDs
        seq_range = {}  # Initialized dict of sequence ranges

        # Begin search:
        if verbose > 1:
            with lock:
                print('ID File Loaded, performing regex search for identifiers...')
                print('ID Specified as: ', id_type)
        if id_type == 'brute':
            if bool(p[1].findall(self.id_rec)):
                id_type = 'accession'
                if verbose > 1:
                    with lock:
                        print(p[1].findall(self.id_rec))
            elif bool(p[0].findall(self.id_rec)):
                id_type = 'gi'
                if verbose > 1:
                    with lock:
                        print(p[1].findall(self.id_rec))
            elif bool(p[2].findall(self.id_rec)):
                id_type = 'id'
                if verbose > 1:
                    with lock:
                        print(p[1].findall(self.id_rec))
            else:
                raise Exception('Couldn\'t identify the id!')
            if verbose > 1:
                with lock:
                    print(
                        'Brute Force was set, tested strings for all pre-registered IDs. ID was selected as type ',
                        id_type)
        if id_type == 'gi':
            if bool(p[0].findall(self.id_rec)):
                found_id = True
                if verbose > 1:
                    print('Successfully found GI numbers, compiling list!')
                item_parts = p[0].findall(self.id_rec)
                if verbose > 1:
                    print('Item:\t', item_parts)
                id_list_ids.append(item_parts[0][0:3])
                if bool(p[3].findall(self.id_rec)):
                    # Seq_range will be a list of tuples where the second element is the range, and the first
                    # is the ID. This way, the function accommodates sequences with a subrange and sequences without a
                    # subrange.
                    seq_range[''.join(p[0].findall(self.id_rec)[0][0:3])] = p[3].findall(self.id_rec)[0]
                    if verbose > 1:
                        print('Found sequence delimiters in IDs!')
            else:
                found_id = False
        elif id_type == 'accession':
            if bool(p[1].findall(self.id_rec)):
                found_id = True
                if verbose > 1:
                    print('Successfully found accession numbers, compiling list!')
                item_parts = p[1].findall(self.id_rec)
                if verbose > 1:
                    print('Item:\t', item_parts)
                id_list_ids.append(item_parts[0][0:3])
                if bool(p[3].findall(self.id_rec)):
                    seq_range[''.join(p[1].findall(self.id_rec)[0][0:3])] = p[3].findall(self.id_rec)[0]
                    if verbose > 1:
                        print('Found sequence delimiters in IDs!')
            else:
                found_id = False
        elif id_type == 'id':
            if bool(p[2].findall(self.id_rec)):
                found_id = True
                if verbose > 1:
                    print('Successfully found ID numbers, compiling list!')
                item_parts = p[2].findall(self.id_rec)
                if verbose > 1:
                    print('Item:\t', item_parts)
                id_list_ids.append(item_parts[0][0:3])
                if bool(p[3].findall(self.id_rec)):
                    seq_range[''.join(p[2].findall(self.id_rec)[0][0:3])] = p[3].findall(self.id_rec)[0]
                    if verbose > 1:
                        print('Found sequence delimiters in IDs!')
            else:
                found_id = False
        else:
            found_id = False
        """
        while not found_id:
            with lock:
                print('Header identified for first sequence ID:', self.id_rec, sep='\n')
                custom_regex = input(
                    'Couldn\'t find ID using preset patterns... Please enter ID pattern for regex search:')
                if custom_regex[0].lower() == 'q':
                    exit()
                print('Will try again...')
                p.append(re.compile(custom_regex))
                if bool(p[4].findall(self.id_rec)):
                    id_type = input('ID name:')
                    found_id = True
                    print('Successfully found custom ID numbers, compiling list!')
                    item_parts = p[4].findall(self.id_rec)
                    if verbose > 1:
                        print('Item:\t', item_parts)
                    id_list_ids.append(item_parts[0][0:3])
                    if bool(p[3].findall(str(self.id_rec))):
                        seq_range[''.join(p[4].findall(self.id_rec)[0][0:3])] = p[3].findall(self.id_rec)[0]
                        if verbose > 1:
                            print('Found sequence delimiters in IDs!')
                else:
                    print('Sorry, still can\'t find it...')
        """
        if verbose > 1:
            with lock:
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
                    with lock:
                        print('Some items were not found. List of items will be saved to the file '
                              'items_not_found.output')
                        for item in itemsnotfound:
                            print(item)
                            # with open(str(out_file.cwd()) + 'items_not_found.output', 'w') as missingitems:
                            #     missingitems.writelines(itemsnotfound)
            keys = [k for k in seqdict.keys()]
            if verbose > 1:
                with lock:
                    print("Sequence Dictionary keys:")
                    print(keys)
            if bool(seq_range):
                seqrange_ids = [ids for ids in seq_range.keys()]
                if verbose > 1:
                    with lock:
                        print('Sequence Range IDs:')
                        print(seqrange_ids)
                for k in keys:
                    if seqdict[k].id in seqrange_ids:
                        if verbose > 1:
                            with lock:
                                print('For sequence {}, found a sequence range!'.format(str(seqdict[k].id)))
                                print('Full length of sequence: {}'.format(len(seqdict[k])))
                        if id_type == 'gi':
                            seq_description_full = p[0].findall(seqdict[k].description)[0]
                        elif id_type == 'accession':
                            seq_description_full = p[1].findall(seqdict[k].description)[0]
                        elif id_type == 'id':
                            seq_description_full = p[2].findall(seqdict[k].description)[0]
                        else:
                            seq_description_full = p[4].findall(seqdict[k].description)[0]
                    else:
                        if verbose > 1:
                            print('No sequence range found, continuing...')
                        continue
                    id_range = ':' + '-'.join(seq_range[k])
                    if int(seq_range[k][0]) > int(seq_range[k][1]):
                        """tmp_id = seqdict[k].id
                        tmp_name = seqdict[k].name
                        tmp_desc = seqdict[k].description
                        tmp_dbxrefs = seqdict[k].dbxrefs
                        tmp_feat = seqdict[k].features
                        tmp_annotations = seqdict[k].annotations
                        tmp_let_anno = seqdict[k].letter_annotations"""
                        seqdict[k].seq = seqdict[k][
                                         int(seq_range[k][1]):int(seq_range[k][0])].seq.reverse_complement()

                    else:
                        seqdict[k] = seqdict[k][int(seq_range[k][0]):int(seq_range[k][1])]
                    if verbose > 1:
                        with lock:
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
                        with lock:
                            print('Sequence ID: \n\t', seqdict[k].id)
                            if id_range:
                                print('Length of subsequence with range {0}: {1}'.format(id_range, len(seqdict[k])))
            if verbose > 1:
                with lock:
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
                    with lock:
                        print(
                            'Some items were not found. List of items will be saved to the file items_not_found.output')
                        for item in itemsnotfound:
                            print(item)
                            # with open(str(out_file.cwd()) + 'items_not_found.output', 'w') as missingitems:
                            #    missingitems.writelines(itemsnotfound)
            keys = [k for k in seqdict.keys()]
            if verbose > 1:
                with lock:
                    print("Sequence Dictionary keys:")
                    print(keys)
            if bool(seq_range):
                seqrange_ids = [ids for ids in seq_range.keys()]
                if verbose > 1:
                    with lock:
                        print('Sequence Range IDs:')
                        print(seqrange_ids)
                for k in keys:
                    if seqdict[k].id in seqrange_ids:
                        if verbose > 1:
                            with lock:
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
                        with lock:
                            print(int(seq_range[k][0]))
                            print(int(seq_range[k][1]))
                    id_range = ':' + '-'.join(seq_range[k])
                    seqdict[k] = seqdict[k][int(seq_range[k][0]):int(seq_range[k][1])]
                    seqdict[k].description = ''.join(seq_description_full[0:3]) + id_range + str(
                        seq_description_full[3])
                    seqdict[k].id += id_range
                    if verbose > 1:
                        with lock:
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
            with lock:
                raise Exception('Not a valid database source!')


def fetchseqMP(ids, species, write=False, output_name='', delim='\t', id_type='brute', server=None, source="SQL",
               db="bioseqdb", host='localhost', driver='psycopg2', version='1.0', user='postgres', passwd='', email='',
               batch_size=50, output_type="fasta", verbose=1, lock=None, n_threads=2, n_subthreads=1):

    from inspect import isgenerator
    if isgenerator(ids):
        if verbose > 1:
            if lock is None:
                print('Received generator!')
            else:
                with lock:
                    print('Received generator!')
    elif isinstance(ids, list):
        if verbose > 1:
            if lock is None:
                print('Received list!')
            else:
                with lock:
                    print('Received list!')
    else:
        if verbose > 1:
            with lock:
                print('Reading ID File... ')
        with ids.open('w') as in_handle:
            id_prelist = [line.strip() for line in in_handle]  # list of each line in the file
            print('Done!')
        ids = list(filter(None, id_prelist))
        if not id_prelist or id_prelist is None:
            if verbose:
                with lock:
                    print('id_prelist is empty!')
            return None

    if verbose > 1:
        if lock is None:
            print('Readied ids!')
        else:
            with lock:
                print('Readied ids!')

    id_list = multiprocessing.JoinableQueue()
    results = multiprocessing.Queue()
    miss_items_queue = multiprocessing.Queue()

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
    fs_instances = [FetchSeqMP(id_queue=id_list, seq_out_queue=results, missing_items_queue=miss_items_queue,
                               delim=delim, id_type=id_type, server=server, species=species, source=source, db=db,
                               host=host, driver=driver, version=version, user=user, passwd=passwd, email=email,
                               output_type=output_type, batch_size=batch_size, verbose=verbose, lock=lock,
                               n_subthreads=n_subthreads)
                    for i in range(n_threads)]
    if verbose > 1:
        print('Done! Starting processes...')
    for fs in fs_instances:
        fs.start()
    if verbose > 1:
        print('Done!')
        print('Assigning FetchSeq records to queue... ')
    for id_rec in ids:
        if id_rec:
            id_list.put(FetchSeq(id_rec=id_rec))
    for i in range(n_threads):
        id_list.put(None)
    if verbose > 1:
        print('Done!')
    output_dict = dict()
    missing_items_list = list()
    if verbose > 1:
        print('Getting sequences from processes... ')
    print('--------------', n_threads)
    while n_threads:  # Todo: Figure out why there's only one result coming out the queues...
        print('--------------', n_threads)
        print('AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA')
        output_dict.update(results.get())
        print('BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB')
        missing_items_list.append(miss_items_queue.get())
        print('CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC')
        n_threads -= 1
        print('DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD')
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

    def __init__(self, proc_id, rb_queue, rb_results_queue, target_species, fw_blast_db, infile_type, output_type,
                 query_species, blast_type, local_blast_1, local_blast_2, rv_blast_db, expect, perc_score,
                 perc_ident, perc_length, megablast, email, id_type, fw_source, fw_id_db, fetch_batch_size, passwd,
                 host, user, driver, fw_id_db_version, verbose, n_threads, fw_blast_kwargs, rv_blast_kwargs,
                 write_intermediates):
        multiprocessing.Process.__init__(self)
        self.name = proc_id
        self.rb_queue = rb_queue
        self.rb_results_queue = rb_results_queue
        self.target_species = target_species
        self.fw_blast_db = fw_blast_db
        self.infile_type = infile_type
        self.output_type = output_type
        self.query_species = query_species
        self.blast_type = blast_type
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
        self.rv_blast_kwargs = fw_blast_kwargs
        self.write_intermediates = write_intermediates
        self.rc_container = RecBlastContainer(proc_id=self.name)

    def run(self):  # The meat of the script
        while True:
            rb_instance = self.rb_queue.get()
            if rb_instance is None:
                self.rb_queue.task_done()
                break
            output = rb_instance(target_species=self.target_species, fw_blast_db=self.fw_blast_db,
                                 infile_type=self.infile_type, output_type=self.output_type,
                                 query_species=self.query_species,
                                 blast_type=self.blast_type, local_blast_1=self.local_blast_1,
                                 local_blast_2=self.local_blast_2,
                                 rv_blast_db=self.rv_blast_db, expect=self.expect, perc_score=self.perc_score,
                                 perc_ident=self.perc_ident,
                                 perc_length=self.perc_length, megablast=self.megablast, email=self.email,
                                 id_type=self.id_type,
                                 fw_source=self.fw_source, fw_id_db=self.fw_id_db, fetch_batch_size=self.batch_size,
                                 passwd=self.passwd,
                                 fw_id_db_version=self.fw_id_db_version, verbose=self.verbose, n_threads=self.n_threads,
                                 host=self.host,
                                 user=self.user, driver=self.driver,
                                 fw_blast_kwargs=self.fw_blast_kwargs, rv_blast_kwargs=self.rv_blast_kwargs,
                                 proc_id=self.name, write_intermediates=self.write_intermediates,
                                 rc_container=self.rc_container)
            self.rb_queue.task_done()
            self.rb_results_queue.put(output)
        return


class RecBlast(object):
    def __init__(self, seq_record):
        self.seq_record = seq_record

    def __call__(self, target_species, fw_blast_db, infile_type, output_type,
                 query_species, blast_type, local_blast_1, local_blast_2, rv_blast_db, expect, perc_score,
                 perc_ident, perc_length, megablast, email, id_type, fw_source, fw_id_db, fetch_batch_size, passwd,
                 host, user, driver, fw_id_db_version, verbose, n_threads, fw_blast_kwargs, rv_blast_kwargs,
                 write_intermediates, rc_container, proc_id):
        # TODO: Replace lock with calls to a central print processes (ANOTHER TODO)!

        from pathlib import Path
        from Bio.Blast import NCBIXML
        transtab = str.maketrans('!@#$%^&*();:.,\'\"/\\?<>|[]{}-=+', '_____________________________')
        self.seq_record.name = self.seq_record.name.translate(transtab)
        rc_container['seq_record'] = self.seq_record
        if verbose:
            print('[Proc: {0}] [Seq.Name: {1}]'.format(proc_id, self.seq_record.name))

        if verbose > 1:
            print('\t Creating handles for intermediary outputs...')
        forward_blast_output = Path("{0}_recblast_out".format(target_species).replace(' ', '_') + '/' +
                                    "{0}_{1}_tmp".format(blast_type, self.seq_record.name).replace(' ', '_') + '/' +
                                    "{0}_{1}_{2}_to_{3}.xml".format(blast_type, self.seq_record.name, query_species,
                                                                    target_species
                                                                    ).replace(' ', '_'))

        forward_id_score_output = Path("{0}_recblast_out".format(target_species).replace(' ', '_') + '/' +
                                       "{0}_{1}_tmp".format(blast_type, self.seq_record.name).replace(' ', '_') + '/' +
                                       "{0}_{1}_{2}_to_{3}.ID_Scores.tmp".format(blast_type, self.seq_record.name,
                                                                                 query_species,
                                                                                 target_species
                                                                                 ).replace(' ', '_'))

        recblast_output_unanno = Path("{0}_recblast_out".format(target_species).replace(' ', '_') + '/' +
                                      "{0}_{1}_tmp".format(blast_type, self.seq_record.name).replace(' ', '_') + '/' +
                                      "unannotated_{0}_{1}.fasta".format(blast_type,
                                                                         self.seq_record.name
                                                                         ).replace(' ', '_'))
        recblast_output = Path("{0}_recblast_out".format(target_species).replace(' ', '_') + '/' +
                               "{0}_{1}.{2}".format(blast_type, self.seq_record.name, output_type).replace(' ', '_'))
        if write_intermediates:
            try:
                forward_blast_output.absolute().parent.mkdir(parents=True)
                if verbose > 1:
                    print('\t\tCreated directory \"{}\"!'.format(str(forward_blast_output.absolute().parent)))
            except FileExistsError:
                if verbose > 1:
                    print('\t\tDirectory \"{}\" already exists! Continuing!'.format(
                        str(forward_blast_output.absolute().parent)))
            try:
                forward_id_score_output.absolute().parent.mkdir(parents=True)
                if verbose > 1:
                    print('\t\tCreated directory \"{}\"!'.format(str(forward_id_score_output.absolute().parent)))
            except FileExistsError:
                if verbose > 1:
                    print('\t\tDirectory \"{}\" already exists! Continuing!'.format(
                        str(forward_blast_output.absolute().parent)))
            try:
                recblast_output_unanno.absolute().parent.mkdir(parents=True)
                if verbose > 1:
                    print('\t\tCreated directory \"{}\"!'.format(str(recblast_output_unanno.absolute().parent)))
            except FileExistsError:
                if verbose > 1:
                    print('\t\tDirectory \"{}\" already exists! Continuing!'.format(
                        str(forward_blast_output.absolute().parent)))

        if fw_blast_db == 'skip':
            if verbose:
                print("\tSkipping Forward Blast! Using local file instead: ")
            with forward_blast_output.open("r") as forward_blasthits:
                if verbose:
                    print('\tOpening Forward blast output located at ', str(forward_blast_output.absolute()))
                blastrecord = NCBIXML.read(forward_blasthits)
        else:
            if verbose:
                print("\tPerforming forward BLAST for {}... ".format(self.seq_record.name))
            fwblasthandle, blast_err = blast(seq_record=self.seq_record, target_species=target_species,
                                             database=fw_blast_db, query_species=query_species, filetype=infile_type,
                                             blast_type=blast_type, local_blast=local_blast_1, expect=expect,
                                             megablast=megablast, blastoutput_custom=str(forward_blast_output),
                                             perc_ident=perc_ident, verbose=verbose, write=write_intermediates,
                                             **fw_blast_kwargs)
            if local_blast_1:
                if write_intermediates:
                    with forward_blast_output.open() as fin:
                        blastrecord = NCBIXML.read(fin)
                else:
                    from io import StringIO
                    with StringIO(fwblasthandle) as fin:
                        blastrecord = NCBIXML.read(fin)
            else:
                blastrecord = fwblasthandle
        if verbose:
            print('Forward blast done!')
            print('Culling Results based on given criteria...')
        f_id_ranked = id_ranker(blastrecord, perc_score=perc_score, expect=expect, perc_length=perc_length,
                                   verbose=verbose)
        f_id_list = [(id_i[0], id_i[2]) for id_i in f_id_ranked]
        f_id_out_list = ['{0}\t{1}\t{2}\n'.format(id_i[0], id_i[1], id_i[2]) for id_i in f_id_ranked]
        if write_intermediates:
            if verbose:
                print('Writing Forward ID Hits to output!')
            with forward_id_score_output.open('w') as fidout:
                fidout.writelines(f_id_out_list)
        if not f_id_out_list:
            print(Warning('Forward Blast yielded no hits, continuing to next sequence!'))
            return
        try:
            server = BioSeqDatabase.open_database(driver=driver, user=user, passwd=passwd,
                                                  host=host, db=fw_id_db)
            seq_dict, missing_items = fetchseqMP(ids=f_id_out_list, species=target_species, delim='\t',
                                                 id_type='brute', server=server, source=fw_source,
                                                 db=fw_id_db, host=host, driver=driver,
                                                 version=fw_id_db_version, user=user,
                                                 passwd=passwd, email=email, batch_size=fetch_batch_size,
                                                 output_type=output_type, output_name='', write=write_intermediates,
                                                 verbose=verbose, n_threads=1,
                                                 n_subthreads=1)
            if verbose:
                print('Done with fetching!')
        except IndexError:
            print(Warning('WARNING! FETCHSEQ FAILED! SKIPPING THIS SEQUENCE!'))
            return
        # return seq_dict, missing_items
        if verbose:
            print('Preparing for Reverse BLAST...')
        # Todo: write code to allow me to see unannotated blast hits

        for entry_id, entry_record in seq_dict.items():
            if not entry_record.seq:
                print(Warning('Entry {0} in reverbse blast file came back empty'.format(entry_record.name)))
                continue
            if verbose:
                print("Entry {} in unannotated RecBlast Hits:\n".format(entry_id))
                for item in [entry_record.id, entry_record.description, entry_record.seq]:
                    print('\t', item)
            reverse_blast_output = Path("{0}_recblast_out".format(target_species).replace(' ', '_') + '/' +
                                        "{0}_{1}_tmp".format(blast_type, self.seq_record.name).replace(' ', '_') + '/' +
                                        "{0}_{1}_{3}_to_{2}_{4}.xml".format(blast_type, self.seq_record.name,
                                                                            query_species, target_species,
                                                                            entry_id).replace(' ', '_'))
            if write_intermediates:
                try:
                    reverse_blast_output.absolute().parent.mkdir(parents=True)
                    if verbose > 1:
                        print('\tCreated Directory ', str(reverse_blast_output.absolute().parent))
                except FileExistsError:
                    pass
            if verbose:
                print('Performing Reverse Blast:')
            if rv_blast_db == 'skip':
                pass
            elif rv_blast_db == 'stop':
                if verbose:
                    print('Not performing reverse blast!')
                continue
            else:
                rvblasthandle, blast_err = blast(seq_record=entry_record, target_species=target_species,
                                                 database=rv_blast_db, query_species=query_species,
                                                 filetype=infile_type,
                                                 blast_type=blast_type, local_blast=local_blast_2, expect=expect,
                                                 megablast=megablast, blastoutput_custom=str(reverse_blast_output),
                                                 perc_ident=perc_ident, verbose=verbose, write=write_intermediates,
                                                 **rv_blast_kwargs)
            if local_blast_1:
                if write_intermediates:
                    with reverse_blast_output.open() as fin:
                        rvblastrecord = NCBIXML.read(fin)
                else:
                    from io import StringIO
                    with StringIO(rvblasthandle) as fin:
                        rvblastrecord = NCBIXML.read(fin)
            else:
                rvblastrecord = rvblasthandle
            if verbose:
                print('Done with Reverse Blast!')
                print('Culling results using given criteria...')
            reverse_blast_annotations = ['\t |[ {0} {1} ({2}) ]|'.format(anno[0], anno[1], anno[2]) for anno in 
                                         id_ranker(rvblastrecord, perc_score=perc_score, perc_length=perc_length,
                                                   expect=expect, verbose=verbose)
                                         ]
            if not reverse_blast_annotations:
                print(Warning('No Reverse Blast Hits were found! Continuing to next Sequence!'))
                return
            if verbose:
                print('Done. Annotating RecBlast Hits:')
            entry_record.description += '|-|' + ''.join(reverse_blast_annotations)
            if verbose > 3:
                print(entry_record)
        recblast_sequence = []
        for item in f_id_ranked:
            # Todo: find a more efficient way to do this:
            for otherkey, otheritem in seq_dict.items:
                if item[0] in otheritem.description:
                    recblast_sequence.append(otheritem)
        return recblast_output, recblast_sequence
        # def run(self):


class SearchMaster(multiprocessing.Process):
    # TODO: Figure out a way for this process to communicate independently with various process
    # TODO: Configure this to search using either BLAST or BLAT
    # TODO: Use SearchIO to have this read results
    def __init__(self, search_queue, result_queue, search_type='blast'):
        multiprocessing.Process.__init__(self)
        self.search_queue = search_queue
        self.result_queue = result_queue
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


def recblastMP(seqfile, target_species, fw_blast_db='chromosome', infile_type='fasta', output_type='fasta',
               host='localhost', user='postgres', driver='psycopg2',
               query_species='Homo sapiens', blast_type='blastn', local_blast_1=False, local_blast_2=False,
               rv_blast_db='nt', expect=10, perc_score=0.5, perc_ident=50, perc_length=0.5, megablast=True, email='',
               id_type='brute', fw_source='sql', fw_id_db='bioseqdb', fetch_batch_size=50, passwd='',
               fw_id_db_version='1.0',
               verbose='v', n_processes=10, n_threads=1, write_intermediates=False, fw_blast_kwargs=dict,
               rv_blast_kwargs=dict):
    import multiprocessing
    import logging
    from pathlib import Path
    
    from Bio import __version__ as bp_version
    # TODO: Make a process that acts as a central print server, have all print() functions redirect to that!
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
    if verbose >= 4:
        multiprocessing.log_to_stderr(
            logging.DEBUG)

    if verbose >= 2:
        print('Creating queues... ')
    rb_queue = multiprocessing.JoinableQueue()
    rb_results = multiprocessing.Queue()
    if verbose >= 2:
        print('Done!')

    # Check Seqfile to make sure its real
    if verbose >= 1:
        print('Loading SeqFile records... ')
    seqfile_path = Path(seqfile)

    try:
        rec_handle = SeqIO.parse(str(seqfile_path.absolute()), output_type)
        if verbose >= 1:
            print('Done!')
    except FileNotFoundError:
        raise
    if verbose >= 1:
        print('Creating RecBlast Threads... ')
    rec_blast_instances = [RecBlastMP_Thread(proc_id=str(i + 1), rb_queue=rb_queue, rb_results_queue=rb_results,
                                             target_species=target_species, fw_blast_db=fw_blast_db,
                                             infile_type=infile_type, output_type=output_type,
                                             query_species=query_species,
                                             blast_type=blast_type, local_blast_1=local_blast_1,
                                             local_blast_2=local_blast_2,
                                             rv_blast_db=rv_blast_db, expect=expect, perc_score=perc_score,
                                             perc_ident=perc_ident,
                                             perc_length=perc_length, megablast=megablast, email=email, id_type=id_type,
                                             fw_source=fw_source, fw_id_db=fw_id_db, fetch_batch_size=fetch_batch_size,
                                             passwd=passwd,
                                             fw_id_db_version=fw_id_db_version, verbose=verbose, n_threads=n_threads,
                                             host=host,
                                             user=user, driver=driver, write_intermediates=write_intermediates,
                                             fw_blast_kwargs=fw_blast_kwargs, rv_blast_kwargs=rv_blast_kwargs,
                                             )
                           for i in range(n_processes)]
    for rcb in rec_blast_instances:
        rcb.start()
    for rec in [i for i in rec_handle]:
        rb_queue.put(RecBlast(seq_record=rec))
    for i in range(n_processes):
        rb_queue.put(None)

    recblast_out = list()
    while n_processes:
        recblast_out.append(rb_results.get())
        n_processes -= 1
    for rcb in rec_blast_instances:
        if rcb.is_alive():
            rcb.join()
    for recblast_output, recblast_sequence in recblast_out:
        with recblast_output.open('w') as rc_out:
            SeqIO.write(recblast_sequence, rc_out, 'fasta')
    return recblast_out
