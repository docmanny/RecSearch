import multiprocessing

from BioSQL import BioSeqDatabase

from modules import biosql_DBSeqRecord_to_SeqRecord


def biosql_seq_lookup_cascade(dtbase, sub_db_name, id_type, identifier, verbose=False):
    seqrec = None
    try_get_id = True
    if try_get_id:
        try:
            if verbose:
                print("Now searching database {0} for {1}: {2}".format(sub_db_name, id_type, identifier))
            seqrec = biosql_DBSeqRecord_to_SeqRecord(dtbase.lookup(**{id_type: identifier}))
            if verbose:
                print('Got sequence for {}!'.format(identifier))
            try_get_id = False
        except IndexError as err:
            if verbose:
                print("WARNING: couldn't find {0} using given ID type... \n Full error: {1}".format(identifier, err))
    if try_get_id:
        identifier_sans_subnumber = identifier.split('.')[0]
        if verbose:
            print('Seeing if removing any sub-numbers (acc: xxxxxx.1 for example) helps...')
            print('Identifier: ', identifier_sans_subnumber)
        try:
            if verbose:
                print("Now searching database {0} for {1}: {2}".format(sub_db_name, id_type,
                                                                       identifier_sans_subnumber))
            seqrec = biosql_DBSeqRecord_to_SeqRecord(dtbase.lookup(**{id_type: identifier_sans_subnumber}))
            if verbose:
                print('Got sequence for {}!'.format(identifier))
            try_get_id = False
        except IndexError as err1:
            if verbose:
                print("WARNING: couldn't find {0} using abbreviated ID... \n Full error: {1}"
                      .format(identifier_sans_subnumber,
                              err1))
    if try_get_id:
        try:
            if verbose:
                print('Attempting to search using Primary ID instead of declared type:')
            seqrec = biosql_DBSeqRecord_to_SeqRecord(dtbase.lookup(primary_id=identifier))
            if verbose:
                print('Got sequence for {}!'.format(identifier))
            try_get_id = False
        except IndexError as err2:
            if verbose:
                print("WARNING: couldn't find {0} using Primary ID... \n full error: {1}".format(identifier, err2))
    if try_get_id:
        try:
            if verbose:
                print('Attempting to search using name instead of declared type:')
            seqrec = biosql_DBSeqRecord_to_SeqRecord(dtbase.lookup(name=identifier))
            if verbose:
                print('Got sequence for {}!'.format(identifier))
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
                print('Got sequence for {}!'.format(identifier))
        except IndexError as err5:
            if verbose:
                print("WARNING: COULD NOT FIND SEQUENCES FOR ID:{0}: \n full error: {1}".format(identifier, err5))
    return seqrec


def biosql_get_record_mp(sub_db_name, passwd, id_list=list(), id_type='accession', driver="psycopg2", user="postgres",
                         host="localhost", db="bioseqdb", num_proc = 2, verbose=True):
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
    :return:
    if __name__ == '__main__':
        biosql_get_record_mp(sub_db_name='MyoLuc2.0', passwd='',
                             id_list=['NW_005871148', 'NW_005871300', 'NW_005871148'], id_type='accession',
                             driver="psycopg2", user="postgres",
                             host="localhost", db="bioseqdb", verbose=True)
    """

    class GetSeqMP(multiprocessing.Process):
        def __init__(self, task_queue, result_queue, db, host, driver, user, passwd, sub_db_name, verbose):
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
            self.server = BioSeqDatabase.open_database(driver=self.driver, user=self.user, passwd=self.password,
                                                       host=self.host, db=self.db)

        # Understand this better
        def run(self):
            proc_name = self.name
            while True:
                next_task = self.task_queue.get()
                if next_task is None:
                    if self.verbose:
                        print('Tasks Complete')
                    self.task_queue.task_done()
                    break
                answer = next_task(connection=self.server, sub_db_name=self.sub_db_name)
                self.task_queue.task_done()
                self.result_queue.put(answer)

    class BioSeqLookupCascade(object):
        def __init__(self, id_type, identifier, verbose):
            self.id_type = id_type
            self.identifier = identifier
            self.verbose = verbose

        def __call__(self, sub_db_name, connection):
            server = connection
            dtbase = server[sub_db_name]
            seqrec = biosql_seq_lookup_cascade(dtbase=dtbase, sub_db_name=sub_db_name, id_type=self.id_type,
                                               identifier=self.identifier, verbose=self.verbose)
            # print('What is self?')
            # print(self.a)
            return self.identifier, seqrec

        # def __str__(self):
        #    return 'ARC'

        def run(self):
            if self.verbose:
                print('Fetching sequence: ', self.identifier)

    idents = multiprocessing.JoinableQueue()
    results = multiprocessing.Queue()

    # num = multiprocessing.cpu_count() * 2
    
    num_jobs = len(id_list)
    seqdict = dict()
    getseqs = [GetSeqMP(idents, results, db=db, host=host, driver=driver, user=user, passwd=passwd,
                            sub_db_name=sub_db_name, verbose=verbose) for i in range(num_proc)]
    for gs in getseqs:
        gs.start()

    for i in range(num_jobs):
        idents.put(BioSeqLookupCascade(id_type=id_type, identifier=id_list[i], verbose=verbose))

    for i in range(num_proc):
        idents.put(None)

    while num_jobs:
        temp = results.get()
        print(temp)
        seqdict[temp[0]] = temp[1]
        num_jobs -= 1
    return seqdict





class RecBlastMP(multiprocessing.Process):
    def __init__(self, seq_queue, rec_queue, target_species, fw_blast_db, infile_type, output_type,
                 query_species, blast_type, local_blast_1, local_blast_2, rv_blast_db, expect, perc_score, perc_ident,
                 perc_length, megablast, email, id_type, fw_source, fw_id_db, fetch_batch_size, passwd, fw_id_db_version,
                 verbose, n_threads):
        multiprocessing.Process.__init__(self)
        self.seq_queue = seq_queue
        self.rec_queue = rec_queue
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
        self. perc_ident = perc_ident
        self.perc_length = perc_length
        self.megablast = megablast
        self.email = email
        self.id_type = id_type
        self.fw_source = fw_source
        self.fw_id_db = fw_id_db
        self.batch_size = fetch_batch_size
        self.passwd = passwd
        self.fw_id_db_version = fw_id_db_version
        self.verbose = verbose
        self.n_threads= n_threads
        self.global_lock = multiprocessing.Lock()  # Stops threads from mixing output

    def run(self):  # The meat of the script
        proc_name = self.name
        while True:
            next_task = self.seq_queue.get()
            if next_task is None:
                self.seq_queue.task_done()
                break
            # TODO: put the step-by-step here
            self.seq_queue.task_done()
            self.rec_queue.put(answer)
        return


class Task(object):
    def __init__(self, a):
        self.a = a

    def __call__(self, connection=None):
        pyConn = connection
        pyCursor1 = pyConn.cursor()

        procQuery = 'UPDATE city SET gid_fkey = gid FROM country  WHERE ST_within((SELECT the_geom FROM city WHERE city_id = %s), country.the_geom) AND city_id = %s' % (
            self.a, self.a)

        pyCursor1.execute(procQuery)
        print
        'What is self?'
        print
        self.a

        return self.a

    def __str__(self):
        return 'ARC'

    def run(self):
        print
        'IN'



def recblastMP(seqfile, target_species, fw_blast_db='chromosome', infile_type='fasta', output_type='fasta',
               query_species='Homo sapiens', blast_type='blastn', local_blast_1=False, local_blast_2=False,
               rv_blast_db='nt', expect=10, perc_score=0.5, perc_ident=50, perc_length=0.5, megablast=True, email='',
               id_type='brute', fw_source='psql', fw_id_db='', fetch_batch_size=50, passwd='', fw_id_db_version='1.0',
               verbose='v', n_threads=10):
    import multiprocessing
    import logging
    from pathlib import Path
    from Bio import SeqIO
    from Bio import __version__ as bp_version

    verbose = verbose.count('v')
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
    if verbose:
        print('BioPython version: ', bp_version)
        print('Beginning RecBlastMP!')
    if verbose >= 3:
        multiprocessing.log_to_stderr(logging.DEBUG)

    if verbose >= 2:
        print('Creating queues... ', end='')
    rb_queue = multiprocessing.JoinableQueue()
    rb_results = multiprocessing.Queue()
    if verbose >= 2:
        print('Done!')

    # Check Seqfile to make sure its real
    if verbose >=2:

    seqfile_path = Path(seqfile)

    try:
        rec_handle = SeqIO.parse(seqfile, output_type)
    except FileNotFoundError:
        raise

    rec_blast_instances = [RecBlastMP(rb_queue, rb_results) for i in range(n_threads)]
    for rcb in rec_blast_instances:
        rcb.start()
    for x in [i for i in rec_handle]:
        rb_queue.put(x)
    for i in range(n_threads):
        rb_queue.put(None)

    # Todo: figure out output.
    while num_jobs:
        result = results.get()
        print result
        num_jobs -= 1


    """
    pool = mp.Pool(n_threads)

    #for seq_batch_gen in n_threads:
    try:
        pool.map_async(partial(recblast, target_species=target_species, fw_blast_db=fw_blast_db,
                                   infile_type=infile_type, output_type=output_type, query_species=query_species,
                                   blasttype=blasttype, localblast1=localblast1, localblast2=localblast2,
                                   rv_blast_db=rv_blast_db, expect=expect, scoreperc=scoreperc,
                                   perc_ident=perc_ident, perc_length=perc_length, megablast=megablast, email=email,
                                   id_type=id_type, fw_source=fw_source, fw_id_db=fw_id_db, batch_size=batch_size,
                                   passwd=passwd, fw_id_db_version=fw_id_db, verbose=verbose, parallel = True),
                       seq_batch_gen)
    except KeyboardInterrupt:
        sys.stdout.write('\033[0m')
        sys.stdout.write('User Interupt\n')
    pool.close()
    pool.join()
    """
