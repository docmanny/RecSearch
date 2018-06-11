import subprocess
from io import StringIO
from time import sleep
import multiprocess as multiprocessing
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from BioSQL import BioSeqDatabase
from BioSQL.BioSeq import DBSeqRecord
from Bio import SeqIO, SeqFeature
from RecBlast.Search import id_search
from inspect import isgenerator
from RecBlast import print
from RecBlast.Search import get_searchdb
from RecBlast.WarningsExceptions import *


def format_range(seqrange, strand, addlength, indent, verbose):
    assert isinstance(addlength, tuple), "addlength was of type {}, must be a tuple!".format(type(addlength))
    assert len(addlength) == 2, "addlength must be a tuple of length 2! Received: {}".format(addlength)
    try:
        lextend = -int(addlength[0])
        rextend = int(addlength[1])
    except Exception as err:
        print(type(err), err)
        lextend = 0
        rextend = 0
    try:
        lrange = int(seqrange[0])
        rrange = int(seqrange[1])
    except Exception as err:
        print(type(err), err)
        lrange = 0
        rrange = -1
    if verbose > 1:
        print('Original range: {0}-{1}{2}'.format(lrange, rrange, strand), indent=indent)
        print('Adding {0} steps to the beginning and {1} steps to the end of the sequence!'.format(lextend, rextend),
              indent=indent)
    if lrange > rrange:
        strand = '+' if strand == '-' else '-'
        lrange = seqrange[1]
        rrange = seqrange[0]
    newrange = tuple(map(lambda x, y: int(x) + y, (lrange, rrange), (lextend, rextend)))
    if verbose > 2:
        print('New range: {0}-{1}{2}'.format(lrange, rrange, strand), indent=indent)
    return newrange, strand


def biosql_get_sub_db_names(passwd, database="bioseqdb", driver="psycopg2", user="postgres", host="localhost"):
    """A convenience wrapper for getting all the sub-database names in a BioSQL-formatted database.

    :param str passwd: The password for the database.
    :param str database: The name of the database.
    :param str driver: The driver BioSQL will use to access the database.
    :param str user: The username used to access the database.
    :param str host: The host of the database.
    :return: a list of sub-database names.
    """
    from BioSQL import BioSeqDatabase
    server = BioSeqDatabase.open_database(driver=driver, user=user, passwd=passwd, host=host, database=database)
    sub_db_name_list = [i for i in server.keys()]
    return sub_db_name_list


def biosql_dbseqrecord_to_seqrecord(dbseqrecord_, off=False):
    """Converts a DBSeqRecord object into a SeqRecord object.

    Motivation of this function was two-fold: first, it makes type testing simpler; and second, DBSeqRecord does
    not have a functional implementation of the translate method.
    :param DBSeqRecord dbseqrecord_: The DBSeqRecord object to be converted.
    :param bool off: Don't actually convert the DBSeqRecord. [Default: False]
    :return:
    """
    assert isinstance(dbseqrecord_, DBSeqRecord), ('Input must be a DBSeqRecord, '
                                                   'was of type {}!').format(type(dbseqrecord_))
    if off:
        return dbseqrecord_
    else:
        return SeqRecord(seq=Seq(data=str(dbseqrecord_.seq)), id=dbseqrecord_.id, name=dbseqrecord_.name,
                         description=dbseqrecord_.description, dbxrefs=dbseqrecord_.dbxrefs,
                         features=dbseqrecord_.features, annotations=dbseqrecord_.annotations,
                         letter_annotations=dbseqrecord_.letter_annotations)


def biosql_seq_lookup_cascade(dtbase, sub_db_name, id_type, identifier, indent=0, verbose=False):
    seqrec = SeqRecord(seq='')
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
            seqrec = biosql_dbseqrecord_to_seqrecord(dtbase.lookup(**{lookup_key: identifier}))
            if verbose:
                print('\tGot sequence for {}!'.format(identifier), indent=indent)
            try_get_id = False
        except KeyError:
            sleep(0.1)
            seqrec = biosql_dbseqrecord_to_seqrecord(dtbase.lookup(**{lookup_key: identifier}))
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
            seqrec = biosql_dbseqrecord_to_seqrecord(dtbase.lookup(**{lookup_key: identifier_sans_subnumber}))
            if verbose:
                print('\tGot sequence for {}!'.format(identifier), indent=indent)
            try_get_id = False
        except KeyError:
            sleep(0.1)
            seqrec = biosql_dbseqrecord_to_seqrecord(dtbase.lookup(**{lookup_key: identifier_sans_subnumber}))
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
            seqrec = biosql_dbseqrecord_to_seqrecord(dtbase.lookup(primary_id=identifier))
            if verbose:
                print('\tGot sequence for {}!'.format(identifier), indent=indent)
            try_get_id = False
        except KeyError:
            sleep(0.1)
            seqrec = biosql_dbseqrecord_to_seqrecord(dtbase.lookup(primary_id=identifier))
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
            seqrec = biosql_dbseqrecord_to_seqrecord(dtbase.lookup(name=identifier))
            if verbose:
                print('\tGot sequence for {}!'.format(identifier), indent=indent)
            try_get_id = False
        except KeyError:
            sleep(0.1)
            seqrec = biosql_dbseqrecord_to_seqrecord(dtbase.lookup(name=identifier))
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
            seqrec = biosql_dbseqrecord_to_seqrecord(dtbase.lookup(**{lookup_key: identifier}))
            if verbose:
                print('\tGot sequence for {}!'.format(identifier), indent=indent)
        except IndexError as err5:
            if verbose:
                print("WARNING: COULD NOT FIND SEQUENCES FOR ID:{0}: \n full error: {1}".format(identifier, err5),
                      indent=indent)
    return seqrec


class GetSeqMP(multiprocessing.Process):
    def __init__(self, task_queue, result_queue, database, host, driver, user, passwd, sub_db_name, verbose, server=None):
        multiprocessing.Process.__init__(self)
        self.task_queue = task_queue
        self.result_queue = result_queue
        self.database = database
        self.host = host
        self.driver = driver
        self.user = user
        self.password = passwd
        self.sub_db_name = sub_db_name
        self.verbose = verbose
        if server is None:
            self.server = BioSeqDatabase.open_database(driver=self.driver, user=self.user, passwd=self.password,
                                                       host=self.host, database=self.database)
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
            print('Waiting for 0.1 second and rerunning in case it was a collision!', indent=self.indent)
            sleep(0.1)
            try:
                dtbase = server[sub_db_name]
            except KeyError:
                raise

        seqrec = biosql_seq_lookup_cascade(dtbase=dtbase, sub_db_name=sub_db_name, id_type=self.id_type,
                                           indent=self.indent, identifier=self.identifier, verbose=self.verbose)
        return self.identifier, seqrec


def biosql_get_record(id_list, sub_db_name, passwd='', id_type='accession', driver="psycopg2", indent=0,
                      user="postgres", host="localhost", database="bioseqdb", num_proc=2, verbose=True, server=None):
    """

    :param sub_db_name:
    :param passwd:
    :param id_list:
    :param id_type:
    :param driver:
    :param indent:
    :param user:
    :param host:
    :param database:
    :param num_proc:
    :param verbose:
    :param server:
    :return:
    if __name__ == '__main__':
        biosql_get_record(sub_db_name='MyoLuc2.0', passwd='',
                             id_list=['NW_005871148', 'NW_005871300', 'NW_005871148'], id_type='accession',
                             driver="psycopg2", user="postgres",
                             host="localhost", database="bioseqdb", verbose=True)
    """
    idents = multiprocessing.JoinableQueue()
    results = multiprocessing.Queue()

    # num = multiprocessing.cpu_count() * 2
    if verbose > 2:
        print('\tStarting biosql_get_record_mp', indent=indent)
    id_list = id_list if isinstance(id_list, list) else [id_list]
    num_jobs = len(id_list)
    seqdict = dict()
    getseqs = [GetSeqMP(idents, results, database=database, host=host, driver=driver, user=user, passwd=passwd,
                        sub_db_name=sub_db_name, verbose=verbose, server=server) for _ in range(num_proc)]
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
        print('Done with biosql_get_record_mp!', indent=indent)
        print('Closing processes!', indent=indent)
    for gs in getseqs:
        if gs.is_alive():
            gs.join()

    itemsnotfound = [i for i in id_list if i not in seqdict.keys()]

    return seqdict, itemsnotfound


class FetchSeqMP(multiprocessing.Process):
    def __init__(self, id_queue, seq_out_queue, delim, id_type, server, species, source, database, database_path,
                 add_length, indent, host, driver, version, user, passwd, email, output_type, batch_size, verbose,
                 n_subthreads):
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
        self.database = database
        self.database_path = database_path
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
                id_item, seq, miss_items = fs_instance(passwd=self.passwd, id_type=self.id_type, driver=self.driver,
                                                       user=self.user, host=self.host, database=self.database,
                                                       database_path=self.database_path, delim=self.delim,
                                                       server=self.server, version=self.version,
                                                       add_length=self.add_length,
                                                       species=self.species, source=self.source, verbose=self.verbose,
                                                       n_threads=self.n_subthreads, indent=self.indent)
            except Exception as err:
                print('FetchSeq Error!')
                print(type(err), err)
                id_item, seq = ('', '')
                miss_items = []
            self.id_queue.task_done()
            self.seq_out_queue.put(((id_item, seq), miss_items))
        return


class FetchSeq(object):  # The meat of the script
    def __init__(self, id_rec):
        assert isinstance(id_rec, tuple) or isinstance(id_rec, list), 'id_rec{0} has type {1}!'.format(id_rec,
                                                                                                       type(id_rec))
        self.id_rec = id_rec

    def __call__(self, delim, species, version, source, passwd, id_type, driver, user, host, database, n_threads,
                 server, verbose, add_length, indent, database_path=None):
        if isinstance(database, dict):
            if species in database:
                database = database[species]
            else:
                raise DatabaseNotFoundError('No sequence source database for species {} '
                                            'was found in the provided dict!'.format(species))
        elif database == "auto" and source in ["2bit", "twobit", "blastdb"]:
            database = get_searchdb(search_type=source, species=species, db_loc=database_path,
                                    verbose=verbose, indent=indent + 1).name
        if database_path:
            database = database_path.rstrip("/") + '/' + database
        if verbose > 1:
            print('Full header for Entry:', indent=indent)
            print(self.id_rec, indent=indent)
        (item_chr, seq_range, item_name, score, strand, thickStart,
         thickEnd, rgb, blockcount, blockspans, blockstarts, query_coverage) = self.id_rec

        try:
            if verbose > 1:
                print('Seq range: ', seq_range, indent=indent)
            assert len(seq_range) == 2, 'Seq_range returned a tuple of length != 2!!!'
            old_strand = strand
            if add_length != (0, 0):
                seq_range, strand = format_range(seqrange=seq_range, strand=strand, addlength=add_length,
                                                 indent=indent + 1, verbose=verbose)
                self.id_rec[1] = seq_range
                self.id_rec[4] = strand
            if -1 in seq_range[0:2]:
                id_full = '{0}'.format(item_chr)
            else:
                id_full = '{0}:{1}-{2}'.format(item_chr, seq_range[0], seq_range[1])
        except KeyError:
            raise KeyError('Sequence {0} lacks a seq_range entry!!!'.format(item_chr))
        if verbose:
            print('ID for query:\t', id_full, indent=indent)

        # Armed with the ID list, we fetch the sequences from the appropriate source

        if source.lower() == "entrez":
            seq, itemnotfound = self.entrez(item_chr, seq_range, indent, add_length, verbose)
        elif source.lower() in ["postgresql", "mysql"]:
            seq, itemnotfound = self.sql(id_item=item_chr, seq_range=seq_range, source=source, species=species,
                                         id_type=id_type, user=user, host=host, passwd=passwd, database=database,
                                         n_threads=n_threads, version=version, server=server,
                                         indent=indent, verbose=verbose)
        elif source == "fasta":  # Note: anecdotally, this doesn't run terribly fast - try to avoid.
            seq, itemnotfound = self.fasta(id_item=item_chr, seq_range=seq_range, database=database, source=source,
                                           indent=indent, verbose=verbose)
        elif source in ["2bit", "twobit"]:
            seq, itemnotfound = self.twobit(id_full=id_full, id_item=item_chr, database=database,
                                            indent=indent, verbose=verbose)
        else:
            raise DatabaseNotFoundError('Not a valid database source: {}'.format(source))
        if itemnotfound is not None:
            if verbose > 1:
                print('Some items were not found:', indent=indent)
                print(itemnotfound, indent=indent)
        if old_strand != strand:
            if verbose > 1:
                print('Sequence was inverted! Reverse complementing now...', indent=indent)
            seq.seq = seq.seq.reverse_complement()
            if verbose > 1:
                print('Done!', indent=indent)
        seq.features.append(SeqFeature.SeqFeature(type='duplicate'))
        if strand == '+':
            s = 1
        elif strand == '-':
            s = -1
        else:
            s = "."
        seq.features[0].location = SeqFeature.FeatureLocation(int(seq_range[0]), int(seq_range[1]), strand=s)
        seq.features[0].qualifiers['score'] = score
        seq.features[0].qualifiers['query_coverage'] = query_coverage
        seq.features[0].qualifiers['thickStart'] = thickStart
        seq.features[0].qualifiers['thickEnd'] = thickEnd
        seq.features[0].qualifiers['itemRGB'] = rgb
        seq.features[0].qualifiers['blockCount'] = blockcount
        seq.features[0].qualifiers['blockSizes'] = blockspans
        seq.features[0].qualifiers['blockStarts'] = blockstarts

        seq.name = item_chr
        return id_full, seq, itemnotfound

    @staticmethod
    def entrez(id_item, seq_range, indent, add_length, verbose):
        raise SearchEngineNotImplementedError('Search type "entrez" is not yet implemented, sorry!!!')

    @staticmethod
    def fasta(id_item, seq_range, database, source, indent, verbose):
        regex = id_search(id_item, indent=indent, verbose=verbose, regex_only=True)
        seqdict = SeqIO.index(database, source,
                              key_function=lambda identifier: regex.search(identifier).groups()[0])
        itemnotfound = id_item if id_item not in seqdict.keys() else None
        seq = seqdict[id_item]
        seq = seq[slice(seq_range[0], seq_range[1])]
        return seq, itemnotfound

    @staticmethod
    def sql(id_item, seq_range, source, species, id_type, user, host, passwd, database,
            n_threads, version, server, indent, verbose):
        driver = "mysql" if source.lower() == 'mysql' else "psycopg2"
        if verbose > 1:
            print('Searching for sequences in local SQL database...', indent=indent)
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
        id_list_search = id_item
        try:
            seq_dict, itemnotfound = biosql_get_record(sub_db_name=sub_db_name, passwd=passwd,
                                                       id_list=id_list_search,
                                                       id_type=id_type, driver=driver, user=user,
                                                       host=host, database=database, num_proc=n_threads, server=server,
                                                       verbose=True)
        except Exception as err:
            print('Please note the sub_databases of server:\n\t', [str(i) for i in server.keys()], indent=indent)
            raise err
        seq_ids = list(seq_dict.keys())
        assert len(seq_ids) == 1, 'Multiple sequences were returned for a single query!'
        seq = seq_dict[seq_ids[0]]
        seq = seq[slice(seq_range[0], seq_range[1])]
        return seq, itemnotfound

    @staticmethod
    def twobit(id_full, id_item, database, indent, verbose):
        seq = None
        itemsnotfound = None

        command = ["twoBitToFa", '{0}:{1}'.format(database, id_full), '/dev/stdout']
        if verbose > 1:
            print('Command:', indent=indent)
            print(' '.join(command), indent=indent + 1)
        twobittofa_handle = subprocess.check_output(command, universal_newlines=True, stdin=subprocess.PIPE,
                                                    stderr=subprocess.PIPE)
        if type(twobittofa_handle) is str:
            seq_out = twobittofa_handle
        else:
            seq_out, seq_err = twobittofa_handle
            raise FetchSeqError(seq_err)

        if seq_out is not None:
            if verbose:
                print('Got sequence for ', id_full, indent=indent)
            if verbose > 3:
                print(str(seq_out).replace('\n', '\n' + '\t' * (indent + 1)), indent=indent + 1)
            with StringIO(seq_out) as output:
                seq = SeqIO.read(output, 'fasta')
        else:
            itemsnotfound = id_item
        return seq, itemsnotfound


def fetchseq(ids, species, write=False, output_name='', delim='\t', id_type='brute', server=None, source="SQL",
             database="bioseqdb", database_path=None, host='localhost', driver='psycopg2', version='1.0',
             user='postgres', passwd='', email='', batch_size=50, output_type="fasta", verbose=1, n_threads=1,
             n_subthreads=1, add_length=(0, 0), indent=0):
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
    for id_item in ids:
        assert len(id_item) == 12, ("Item {0} in id_list has {1} items, not 5!\n"
                                    "Format should be: "
                                    "chr, (start,end), id, score, strand, thickStart, thickEnd, rgb, blockcount,"
                                    " blockspans, blockstarts, query_span"
                                    "!").format(" ".join((" ".join(item) if not isinstance(item, str) else item
                                                        for item in id_item)),
                                                len(id_item))
    if verbose > 1:
        print('Readied ids!', indent=indent)

    id_list = multiprocessing.JoinableQueue()
    results = multiprocessing.Queue()

    if 'sql' in source.lower():
        if server is None:
            try:
                if verbose > 1:
                    print('No server received, opening server...', indent=indent)
                server = BioSeqDatabase.open_database(driver=driver, user=user, passwd=passwd, host=host,
                                                      database=database)
                if verbose > 1:
                    print('Done!', indent=indent)
            except Exception as err:
                if verbose > 1:
                    print('Failed to open server!', indent=indent)
                    print(str(type(err)), err, sep=' ', indent=indent)
                raise
        else:
            if verbose > 1:
                print('Received server handle:', indent=indent)
                print(server, indent=indent)
            if verbose > 2:
                print('Please note the sub_databases of server:\n\t', [str(i) for i in server.keys()], indent=indent)
    elif source.lower() in ['fasta', '2bit', 'twobit']:
        print('Search type: ', source, indent=indent)
    else:
        raise SearchEngineNotImplementedError('Search using source {} has not yet been implemented!'.format(source))
    if verbose > 1:
        print('Creating FecSeq Processes...', indent=indent)
    fs_instances = [FetchSeqMP(id_queue=id_list, seq_out_queue=results,
                               delim=delim, id_type=id_type, server=server, species=species, source=source,
                               database=database, database_path= database_path,
                               host=host, driver=driver, version=version, user=user, passwd=passwd, email=email,
                               output_type=output_type, batch_size=batch_size, verbose=verbose,
                               n_subthreads=n_subthreads, add_length=add_length, indent=indent + 1)
                    for _ in range(n_threads)]
    if verbose > 1:
        print('Done! Starting processes...', indent=indent)
    for fs in fs_instances:
        fs.start()
    if verbose > 1:
        print('Done!', indent=indent)
        print('Assigning FetchSeq records to queue... ', indent=indent)
    id_order = []
    for i, id_rec in enumerate(ids):
        try:
            id_order.append("{0}:{1}-{2}".format(id_rec[0],id_rec[1][0],id_rec[1][1]))
        except IndexError:
            id_order.append("{0}".format(id_rec[0]))
        try:
            id_list.put(FetchSeq(id_rec=id_rec))
        except AssertionError as err:
            print(i, type(err), err, sep=' ')
            break
    for _ in fs_instances:
        id_list.put(None)
    if verbose > 1:
        print('Done!', indent=indent)
    output_dict = dict()
    missing_items_list = list()
    if verbose > 1:
        print('Getting sequences from processes... ', indent=indent)
    n_jobs = len(ids)
    while n_jobs:
        seq, missing = results.get()
        output_dict[seq[0]] = seq[1]
        missing_items_list.append(missing)
        n_jobs -= 1
    if verbose > 1:
        print('Done! Finished fetching sequences!', indent=indent)
        print('Closing processes!', indent=indent)
    for fs in fs_instances:
        if fs.is_alive():
            fs.join()
    output_list = [output_dict[i] for i in id_order if i in output_dict]
    if write:
        SeqIO.write(output_list, output_name, output_type)
        return
    else:
        if missing_items_list == [None]:
            missing_items_list = None
        return output_list, missing_items_list


def seq_from_exon(query_record, forward_search_criteria):
    # TODO: filter sequences as in Search.id_ranker()
    from functools import reduce
    query_seq = []
    for hsp in query_record:
        hsp_frags = [hsp_frag.hit for hsp_frag in hsp]
        for hit, hsp_frag in zip(hsp_frags, hsp):
            hit.features[0] = SeqFeature.SeqFeature(SeqFeature.FeatureLocation(*hsp_frag.hit_range),
                                                    strand=hsp_frag.hit_strand)
        query_seq.append(hsp_frags)
    query_seq = [sorted(hsp,
                        key=lambda x: x.features[0].location.start,
                        reverse=all((hsp_fragment.features[0].location.strand == -1 for hsp_fragment in hsp)))
                 for hsp in query_seq]
    query_seq_joint = [reduce(lambda x, y: x + y, hsp) for hsp in query_seq]
