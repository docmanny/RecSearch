import logging
from contextlib import redirect_stdout
from datetime import datetime
from pathlib import Path

import multiprocess as multiprocessing

from Bio import SearchIO, SeqIO
from Bio import __version__ as bp_version
from io import StringIO
from Bio.SeqRecord import SeqRecord
from BioSQL import BioSeqDatabase

from RecBlast import print, ProgressBar, __version__
from RecBlast.FetchSeq import fetchseq, biosql_get_sub_db_names
from RecBlast.Search import Search, get_searchdb, blat_server, id_ranker
from RecBlast.RBC import RecBlastContainer
from RecBlast.WarningsExceptions import *
from RecBlast.Auxilliary import translate_annotation


class RecBlastMPThread(multiprocessing.Process):
    """
    RecBlast_MP_Thread_Handle is the first branch to be made. It will perform the actual RecBlastRun.
    """

    def __init__(self, proc_id, rb_queue, rb_results_queue, fw_search_db, infile_type, output_type, id_db_path,
                 query_species, fw_search_type, rv_search_type, fw_search_local, rv_search_local, rv_search_db, expect,
                 perc_score, perc_span, outfolder, indent, reciprocal_method, hit_name_only, translate_hit_name,
                 perc_ident, perc_length, megablast, email, id_type, id_source, id_db, fetch_batch_size, passwd,
                 host, user, driver, id_db_version, verbose, n_threads, fw_search_kwargs, rv_search_kwargs,
                 write_intermediates):
        multiprocessing.Process.__init__(self)
        self.name = proc_id
        self.outfolder = outfolder
        self.rb_queue = rb_queue
        self.rb_results_queue = rb_results_queue
        self.fw_search_db = fw_search_db
        self.infile_type = infile_type
        self.output_type = output_type
        self.query_species = query_species
        self.fw_search_type = fw_search_type
        self.rv_search_type = rv_search_type
        self.fw_search_local = fw_search_local
        self.rv_search_local = rv_search_local
        self.rv_search_db = rv_search_db
        self.expect = expect
        self.perc_score = perc_score
        self.perc_span = perc_span
        self.perc_ident = perc_ident
        self.perc_length = perc_length
        self.megablast = megablast
        self.email = email
        self.id_type = id_type
        self.id_source = id_source
        self.host = host
        self.user = user
        self.driver = driver
        self.id_db = id_db
        self.batch_size = fetch_batch_size
        self.passwd = passwd
        self.translate_hit_name = translate_hit_name
        self.id_db_version = id_db_version
        self.verbose = verbose
        self.n_threads = n_threads
        self.fw_search_kwargs = fw_search_kwargs
        self.rv_search_kwargs = rv_search_kwargs
        self.id_db_path = id_db_path
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
                    output = rb_instance(fw_search_db=self.fw_search_db, id_db_path=self.id_db_path,
                                         infile_type=self.infile_type, output_type=self.output_type,
                                         query_species=self.query_species,
                                         fw_search_type=self.fw_search_type, rv_search_type=self.rv_search_type,
                                         fw_search_local=self.fw_search_local,
                                         rv_search_local=self.rv_search_local,
                                         rv_search_db=self.rv_search_db, expect=self.expect, perc_score=self.perc_score,
                                         perc_ident=self.perc_ident, perc_span=self.perc_span,
                                         perc_length=self.perc_length, megablast=self.megablast, email=self.email,
                                         id_type=self.id_type,
                                         id_source=self.id_source, id_db=self.id_db, 
                                         fetch_batch_size=self.batch_size,
                                         passwd=self.passwd, translate_hit_name=self.translate_hit_name,
                                         id_db_version=self.id_db_version, verbose=self.verbose, 
                                         indent=self.indent,
                                         n_threads=self.n_threads,
                                         host=self.host, reciprocal_method=self.reciprocal_method,
                                         user=self.user, driver=self.driver,
                                         fw_search_kwargs=self.fw_search_kwargs, rv_search_kwargs=self.rv_search_kwargs,
                                         proc_id=self.name, write_intermediates=self.write_intermediates,
                                         hit_name_only=self.hit_name_only
                                         )
                    self.rb_queue.task_done()
                    self.rb_results_queue.put(output)
                except Exception as err:
                    print('Woah! Something went wrong! Aborting!')
                    print('Here\'s the error:\n', type(err), err)
                    self.rb_results_queue.put(dict(error=err, proc_id=self.name))
        master_out_handle.close()
        return


class RecBlastRun(object):
    def __init__(self, seq_record, target_species):
        self.starttime = datetime.now()
        self.seq_record = seq_record
        transtab = str.maketrans('!@#$%^&*();:.,\'\"/\\?<>|[]{}-=+', '_____________________________')
        self.seq_record.id = self.seq_record.id.translate(transtab)
        self.target_species = target_species

    def __call__(self, fw_search_db, infile_type, output_type, id_db_path, reciprocal_method,
                 query_species, fw_search_type, rv_search_type, fw_search_local, rv_search_local, rv_search_db, expect,
                 perc_score, indent, hit_name_only, perc_span, translate_hit_name,
                 perc_ident, perc_length, megablast, email, id_type, id_source, id_db, fetch_batch_size, passwd,
                 host, user, driver, id_db_version, verbose, n_threads, fw_search_kwargs, rv_search_kwargs,
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
            forward_search_output=Path("{0}_recblast_out".format(target_species).replace(' ', '_') + '/' +
                                       "{0}_{1}_tmp".format(fw_search_type,
                                                            self.seq_record.id
                                                            ).replace(' ', '_') + '/' +
                                       "{0}_{1}_{2}_to_{3}.xml".format(fw_search_type,
                                                                       self.seq_record.id,
                                                                       query_species,
                                                                       target_species
                                                                       ).replace(' ', '_')),
            forward_id_score_output=Path("{0}_recblast_out".format(target_species).replace(' ', '_') + '/' +
                                         "{0}_{1}_tmp".format(fw_search_type,
                                                              self.seq_record.id
                                                              ).replace(' ', '_') + '/' +
                                         "{0}_{1}_{2}_to_{3}.ID_Scores".format(fw_search_type,
                                                                               self.seq_record.id,
                                                                               query_species,
                                                                               target_species
                                                                               ).replace(' ', '_')),
            recblast_output_unanno=Path("{0}_recblast_out".format(target_species).replace(' ', '_') + '/' +
                                        "{0}_{1}_tmp".format(fw_search_type,
                                                             self.seq_record.id
                                                             ).replace(' ', '_') + '/' +
                                        "unannotated_{0}_{1}.fasta".format(fw_search_type,
                                                                           self.seq_record.id
                                                                           ).replace(' ', '_')),
            recblast_output=Path("{0}_recblast_out".format(target_species).replace(' ', '_') + '/' +
                                 "{0}_{1}.{2}".format(fw_search_type,
                                                      self.seq_record.id,
                                                      output_type
                                                      ).replace(' ', '_')),
            search_nohits=Path("{0}_recblast_out".format(target_species).replace(' ', '_') + '/' +
                               "{0}_{1}.no-hits".format(fw_search_type,
                                                        self.seq_record.id
                                                        ).replace(' ', '_'))
        )
        
        fw_search = rc_container['forward_search']

        if verbose:
            print("Performing Forward search for {}... ".format(self.seq_record.id), indent=indent)
        try:
            fw_search_kwargs = fw_search_kwargs if isinstance(fw_search_kwargs, dict) else dict()
            try:
                FWSearch = Search(search_type=fw_search_type)
                fwsearchrecord, search_err = FWSearch(seq_record=self.seq_record, species=target_species,
                                                      database=fw_search_db, filetype=infile_type,
                                                      search_db_path=id_db_path,
                                                      search_type=fw_search_type, local=fw_search_local,
                                                      expect=expect,
                                                      megablast=megablast, n_threads=n_threads,
                                                      blastoutput_custom=output_paths['forward_search_output'],
                                                      perc_ident=perc_ident, verbose=verbose, indent=indent+1,
                                                      write=write_intermediates, **fw_search_kwargs)
            except ValueError:
                rc_container['recblast_unanno'] = [SeqRecord('')]
                return rc_container_full
        except Exception as err:
            print('WARNING! UNCATCHED EXCEPTION OCCURED!')
            print(type(err), err)
            return rc_container_full
        if search_err:
            fw_search['search_errors'] = search_err
            print('Forward search returned with an error!', search_err)
            return rc_container_full
        fw_search['search_results'] = fwsearchrecord
        if not fwsearchrecord:
            print('Forward search record was empty!')
            return rc_container_full
        if verbose:
            print('Forward search done!', indent=indent)
            print('Culling results based on given criteria...', indent=indent)

        try:
            f_id_ranked = id_ranker(fwsearchrecord, perc_ident=perc_ident, perc_score=perc_score, perc_span=perc_span,
                                    expect=expect, perc_length=perc_length, verbose=verbose, indent=indent+1)
        except Exception as err:
            print('WARNING! UNCATCHED ERROR IN ID_RANKER!')
            print(type(err), err)
            return rc_container_full
        if verbose:
            print('Done!', indent=indent)
        # f_id_ranked: [(seq_name, seq_range, strand, seq_score, seq_coverage) per hit]
        f_id_out_list = ['{0}:{1}-{2}{3}\t{4}\n'.format(id_i[0], id_i[1][0],
                                                        id_i[1][1], id_i[2], id_i[3]) for id_i in f_id_ranked]
        fw_ids = rc_container['forward_ids']
        fw_ids['ids'] = f_id_ranked
        fw_ids['pretty_ids'] = f_id_out_list

        if not f_id_ranked:
            print('Forward search yielded no hits, continuing to next sequence!')
            return rc_container_full
        if fw_search_type.lower() in ['blat', 'tblat']:
            if verbose > 1:
                print('Since blat was selecting, setting id_db equal to fw_search_db', indent=indent)
            if isinstance(id_db, dict):
                blat_2bit = id_db[target_species]
            elif id_db == 'auto':
                blat_2bit = get_searchdb(search_type=fw_search_type, species=target_species, db_loc=id_db_path,
                                         verbose=verbose, indent=indent+1)
            elif isinstance(id_db, str):
                blat_2bit = id_db
            else:
                raise FileNotFoundError('Invalid 2bit file designation!')
            id_db = Path(id_db_path, blat_2bit+'.2bit').absolute()
            id_db = str(id_db) if id_db.is_file() else None
            if id_db is None:
                raise FileNotFoundError('Invalid 2bit file!')
            if verbose > 1:
                print(id_db)
        if 'sql' in id_source.lower():
            server = BioSeqDatabase.open_database(driver=driver, user=user, passwd=passwd,
                                                  host=host, db=id_db)
        try:
            if id_source == 'strict':
                if verbose:
                    print('Fetching Sequence from hit results!', indent=indent)
                seq_dict = {}
                missing_items = []
                assert isinstance(fwsearchrecord, SearchIO.QueryResult), 'Strict fetching only ' \
                                                                         'implemented for SearchIO!'
                for hit in fwsearchrecord:
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
                seq_dict, missing_items = fetchseq(ids=f_id_ranked, species=target_species, delim='\t',
                                                   id_type='brute', server=server, source=id_source,
                                                   db=id_db, host=host, driver=driver,
                                                   version=id_db_version, user=user,
                                                   passwd=passwd, email=email, batch_size=fetch_batch_size,
                                                   output_type=output_type, output_name='', write=write_intermediates,
                                                   verbose=verbose, n_threads=1, indent=indent+1,
                                                   n_subthreads=1)
            if verbose:
                print('Done with fetching!', indent=indent)
        except IndexError:
            print('WARNING! FETCHSEQ FAILED! SKIPPING THIS SEQUENCE!')
            return rc_container_full

        if not missing_items:
            fw_ids['missing_items'] = missing_items
            print('Items were missing!', indent=indent)
            for i in missing_items:
                print(i, indent=indent+1)

        recblast_sequence = []

        if seq_dict:
            if verbose > 2:
                print('Ranking SeqDicts:', indent=indent)
            for item in f_id_ranked:
                if verbose > 3:
                    print(''.join((str(i) for i in item)), indent=indent + 1)
                try:
                    recblast_sequence.append(seq_dict[item[0]])
                except KeyError as err:
                    raise RecBlastException(type(err), item, ':', err)
        else:
            err = 'No SeqDict was returned for record {0} in process {1}!'.format(''.join((self.target_species,
                                                                                           self.seq_record.id)),
                                                                                  proc_id)
            print(type(err), err)
            raise RecBlastException(err)
        if not isinstance(recblast_sequence, list):
            recblast_sequence = [recblast_sequence]
        rc_container['recblast_unanno'] = recblast_sequence if recblast_sequence != list() else [SeqRecord('')]
        if not recblast_sequence:
            print('Error! "recblast_sequence" came back empty!')
            return rc_container_full
        if verbose:
            print('Preparing for Reverse search...', indent=indent)
        for index, entry_record in enumerate(recblast_sequence):
            assert isinstance(entry_record, SeqRecord), 'Warning! entry_record is of type {} ' \
                                                        'and not a SeqRecord!'.format(str(type(entry_record)))
            if verbose:
                print("Entry {} in unannotated RecBlast Hits:\n".format(entry_record.id), indent=indent+1)
                for item in [entry_record.id, entry_record.description,
                             entry_record.seq[0:10] + '...' + entry_record.seq[-1]]:
                    print(item, indent=indent+2)
            output_paths['reverse_search_output'].append(
                Path("{0}_recblast_out".format(target_species).replace(' ', '_') +
                     '/' +
                     "{0}_{1}_tmp".format(rv_search_type,
                                          self.seq_record.id
                                          ).replace(' ', '_') +
                     '/' +
                     "{0}_{1}_{3}_to_{2}_{4}.xml".format(rv_search_type,
                                                         self.seq_record.id,
                                                         query_species,
                                                         target_species,
                                                         entry_record.name
                                                         ).replace(' ', '_')
                     ))
            if verbose:
                print("Performing Forward search for {}... ".format(entry_record.name), indent=indent)
            rv_search = rc_container['reverse_search']
            try:
                rv_search_kwargs = rv_search_kwargs if isinstance(rv_search_kwargs, dict) else dict()
                try:
                    RVSearch = Search(search_type=rv_search_type)
                    rvsearchrecord, search_err = RVSearch(seq_record=entry_record, species=query_species,
                                                          database=rv_search_db, filetype=infile_type,
                                                          search_db_path=id_db_path,
                                                          search_type=rv_search_type, local=rv_search_local,
                                                          expect=expect,
                                                          megablast=megablast, n_threads=n_threads,
                                                          blastoutput_custom=output_paths['reverse_search_output'],
                                                          perc_ident=perc_ident, verbose=verbose, indent=indent+1, 
                                                          write=write_intermediates, **rv_search_kwargs)
                except ValueError:
                    rc_container['recblast_unanno'] = [SeqRecord('')]
                    return rc_container_full
            except Exception as err:
                print('WARNING! UNCATCHED EXCEPTION OCCURED!')
                print(type(err), err)
                return rc_container_full
            if search_err:
                print('Reverse search returned with errors!')
                rv_search['search_errors'] = search_err
                return rc_container_full
            rv_search['search_results'] = rvsearchrecord
            if verbose:
                print('Done with Reverse search!', indent=indent)
                print('Culling results using given criteria...', indent=indent)
            try:
                reverse_hits = id_ranker(rvsearchrecord, perc_ident=perc_ident, perc_score=perc_score,
                                         perc_span=perc_span, perc_length=perc_length, expect=expect, verbose=verbose,
                                         indent=indent+1, method=reciprocal_method)
            except Exception as err:
                print('No Reverse search hits were found for this hit!', indent=indent + 1)
                print('Continuing to next Sequence!', indent=indent + 1)
                continue
            print('Reverse search hits:', indent=indent+1)
            print(reverse_hits, indent=indent+2)
            reverse_search_annotations = []
            for anno in reverse_hits:
                try:
                    if translate_hit_name:
                        new_anno = translate_annotation(anno[0])
                except Exception as err:
                    new_anno = anno[0]
                finally:
                    reverse_search_annotations.append('\t |[ {0}:{1}-{2}{3} ]|'.format(new_anno, anno[1][0],
                                                                                       anno[1][1], anno[2]))
            if not reverse_search_annotations:
                print('No Reverse search hits were found for this hit!', indent=indent+1)
                print('Continuing to next Sequence!', indent=indent+1)
                continue
            else:
                if verbose > 1:
                    print('Done. Annotating RecBlast Hits:', indent=indent+1)
            rv_ids = rc_container['reverse_ids']
            rv_ids['ids'] = reverse_search_annotations
            if reciprocal_method == 'best hit':
                print('Best Hit Reciprocal BLAST was selected, ending Reverse searches after first annotation!',
                      indent=indent)
                entry_record.description += '|-|' + reverse_search_annotations[0] if \
                    isinstance(reverse_search_annotations, list) \
                    else reverse_search_annotations
            if verbose > 3:
                print(entry_record, indent=indent+2)

        if hit_name_only:
            recblast_sequence = [rec.id + '\t' + rec.description
                                 for rec in recblast_sequence]
        if not isinstance(recblast_sequence, list):
            recblast_sequence = [recblast_sequence]
        if recblast_sequence == list():
            recblast_sequence.append(SeqRecord(''))
        rc_container['recblast_results'] = recblast_sequence
        run_time = str(datetime.now() - self.starttime)
        rc_container['run_time'] = run_time
        print('PID', 'SeqName', sep='\t')
        print(proc_id, self.seq_record.id, sep='\t')
        print('Run time: ', run_time)
        return rc_container_full
        # def run(self):


def recblast_run(seqfile, target_species, fw_search_db='auto', rv_search_db='auto-transcript', infile_type='fasta',
                 output_type='fasta',
                 host='localhost', user='postgres', driver='psycopg2',
                 query_species='Homo sapiens', fw_search_type='blastn', rv_search_type='blastn', fw_search_local=False,
                 rv_search_local=False,
                 expect=10, perc_score=0.5, perc_span=0.1, perc_ident=0.50, perc_length=0.5, megablast=True,
                 email='', run_name='default', output_loc='./RecBlast_output',
                 id_type='brute', id_source='sql', id_db='bioseqdb', fetch_batch_size=50,
                 passwd='', hit_name_only=False, min_mem=False,
                 id_db_version='auto', id_db_path='/usr/db/id_db_path', indent=0, translate_hit_name=True,
                 verbose='v', max_n_processes='auto', n_threads=2, write_intermediates=False, write_final=True,
                 reciprocal_method='best hit', fw_search_kwargs=None, rv_search_kwargs=None):
    """

    >>>rc_out = recblast_run('nr_Homo_sapiens_protein_GRCh38p9.fa',
                             ['Loxodonta africana','Procavia capensis','Trichechus manatus'],
                             fw_search_db={'Loxodonta africana':20007,
                                           'Procavia capensis':20008,
                                           'Trichechus manatus':20009}, rv_search_db={'Homo sapiens':20005},
                             infile_type='fasta',
                             output_type='database', host='localhost', user='postgres', driver='psycopg2',
                             query_species='Homo sapiens', fw_search_type='tblat', rv_search_type='blat-transcript',
                             fw_search_local=True, rv_search_local=True,
                             expect=10, perc_score=0.009, perc_span=0.1, perc_ident=0.69, perc_length=0.001,
                             megablast=True, email='', run_name='EHM_AvA',
                             id_type='brute', id_source='2bit', id_db='auto', fetch_batch_size=50,
                             passwd='', hit_name_only=True, min_mem=True,
                             id_db_version='auto', id_db_path='/usr/db/BLAT', indent=0,
                             verbose='vvv', max_n_processes='auto', n_threads=2, write_intermediates=False,
                             write_final=True, reciprocal_method = 'best hit', fw_search_kwargs=None,
                             rv_search_kwargs=None)
    :param seqfile:
    :param target_species:
    :param fw_search_db:
    :param rv_search_db:
    :param infile_type:
    :param output_type:
    :param host:
    :param user:
    :param driver:
    :param query_species:
    :param fw_search_type:
    :param rv_search_type:
    :param fw_search_local:
    :param rv_search_local:
    :param expect:
    :param perc_score:
    :param perc_ident:
    :param perc_length:
    :param perc_span:
    :param megablast:
    :param email:
    :param id_type:
    :param id_source:
    :param id_db:
    :param fetch_batch_size:
    :param passwd:
    :param translate_hit_name:
    :param hit_name_only:
    :param min_mem:
    :param id_db_version:
    :param id_db_path:
    :param indent:
    :param verbose:
    :param max_n_processes:
    :param n_threads:
    :param write_intermediates:
    :param write_final:
    :param run_name:
    :param output_loc:
    :param reciprocal_method:
    :param fw_search_kwargs:
    :param rv_search_kwargs:
    :return:
    """
    # Verbose-ometer
    if isinstance(verbose, str):
        verbose = verbose.lower().count('v')
    elif isinstance(verbose, int) and verbose >= 0:
        pass
    else:
        raise TypeError('Verbose was type {}; must be either be an integer greater '
                        'than or equal to zero, or a number of v\'s equal '
                        'to the desired level of verbosity'.format(type(verbose)))
    if verbose:
        print('RecBlast version: ', __version__)
        print('Using BioPython version: ', bp_version)
        print('Beginning RecBlastMP!')
    if verbose == 1:
        print('Basic verbose mode active. Will print only essential commentary.')
    elif verbose == 2:
        print('Verbose mode was set to 2. Will elaborate considerably about the script.')
    elif verbose == 3:
        print('Debugging-level verbose mode set. You will be innunadated by text. Brace yourself, and hold on to your '
              'console.')
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
        raise RecBlastException('SeqFile was empty!')
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
        date_str = datetime.now().strftime('%y-%m-%d_%I-%M-%p')
        outfolder = Path(output_loc.rstrip('/') + '/{0}/'.format(date_str))
    else:
        outfolder = Path(output_loc.rstrip('/') + '/{0}/'.format(run_name))
    try:
        outfolder.mkdir(parents=True)
    except FileExistsError:
        pass
    #########################################################################

    # Blat stuff
    # Check if BLAT, then if so, make sure that fw_search_db is a dictionary with a port for each species:
    if fw_search_type.lower() in ['blat', 'tblat']:
        assert isinstance(fw_search_db, dict) or isinstance(fw_search_db, Path), \
            "For BLAT searches, fw_search_db must be a dictionary with valid species-port key pairs; OR a Path object"
    if rv_search_type.lower() in ['blat', 'tblat', 'blat-transcript', 'tblat-transcript']:
        assert isinstance(rv_search_db, dict), "For BLAT searches, rv_search_db must be a dictionary with " \
                                              "valid species-port key pairs"

    if isinstance(target_species, str):

        if target_species.lower() == 'all':
            if verbose:
                print('Target species set to all. Searching server for all available species:')
            target_species = biosql_get_sub_db_names(passwd=passwd, db=id_db, driver=driver, user=user,
                                                     host=host)
            if verbose:
                print(target_species, indent=1)
    server_activated = {}
    if isinstance(rv_search_db, dict):
        if rv_search_db[query_species] == 'auto':
            rv_server_online = False
        else:
            print('Checking status of rv_server')
            rv_server_online = blat_server('auto', order='status', host=host, port=rv_search_db[query_species],
                                           species=query_species, id_db_path=id_db_path, verbose=verbose,
                                           indent=1, type=fw_search_type)
        if not rv_server_online and 'blat' in rv_search_type.lower():
            rv_search_db[query_species] = blat_server('auto', 'start', host=host, port=30000, species=query_species,
                                                      id_db_path=id_db_path, verbose=verbose,
                                                      indent=1, type=rv_search_type)
            server_activated[query_species] = rv_search_db[query_species]
        elif not rv_server_online and 'tblat' in rv_search_type.lower():
            rv_search_db[query_species] = blat_server('auto', 'start', host=host, port=30000, species=query_species,
                                                      id_db_path=id_db_path, verbose=verbose,
                                                      indent=1, type=rv_search_type)
            server_activated[query_species] = rv_search_db[query_species]
    if isinstance(fw_search_db, dict):
        if isinstance(target_species, list):
            for species in target_species:
                if fw_search_db[species] == 'auto':
                    print('fw_search_db was set to "auto!"')
                    fw_server_online = False
                else:
                    fw_server_online = blat_server('auto', 'status', host=host, port=fw_search_db[species],
                                                   species=species, id_db_path=id_db_path, verbose=verbose,
                                                   indent=1, type=fw_search_type)
                    print('Status of the forward server: ', fw_server_online)
                if not fw_server_online and 'blat' in fw_search_type.lower():
                    print('Forward server was not online, starting!')
                    fw_search_db[species] = blat_server('auto', 'start', host=host, port=20000, species=species,
                                                        id_db_path=id_db_path, verbose=verbose,
                                                        indent=1, type=fw_search_type)
                    server_activated[species] = fw_search_db[species]
                elif not fw_server_online and 'tblat' in fw_search_type.lower():
                    print('Forward server was not online, starting!')
                    fw_search_db[species] = blat_server('auto', 'start', host=host, port=20000, species=species,
                                                        id_db_path=id_db_path, verbose=verbose,
                                                        indent=1, type=fw_search_type)
                    server_activated[species] = fw_search_db[species]
            print(fw_search_db)
        else:
            if fw_search_db[target_species] == 'auto':
                print('fw_search_db was set to "auto!"')
                fw_server_online = False
            else:
                fw_server_online = blat_server('auto', 'status', host=host, port=fw_search_db[target_species],
                                               species=target_species, id_db_path=id_db_path,
                                               verbose=verbose, indent=1, type=fw_search_type)
                print('Status of the forward server: ', fw_server_online)
            if not fw_server_online and 'blat' in fw_search_type.lower():
                print('Forward server was not online, starting!')
                fw_search_db[target_species] = blat_server('auto', 'start', host=host, port=20000,
                                                           species=target_species, id_db_path=id_db_path,
                                                           verbose=verbose, indent=1, type=fw_search_type)
                server_activated[target_species] = fw_search_db[target_species]
            elif not fw_server_online and 'tblat' in fw_search_type.lower():
                print('Forward server was not online, starting!')
                fw_search_db[target_species] = blat_server('auto', 'start', host=host, port=20000,
                                                           species=target_species, id_db_path=id_db_path,
                                                           verbose=verbose, indent=1, type=fw_search_type)
                server_activated[target_species] = fw_search_db[target_species]
            print(fw_search_db)
    #########################################################################

    # RecBlast Thread init
    if verbose >= 1:
        print('Creating RecBlast Threads... ')
    rec_blast_instances = [RecBlastMPThread(proc_id=str(i + 1), rb_queue=rb_queue, rb_results_queue=rb_results,
                                            fw_search_db=fw_search_db, id_db_path=id_db_path,
                                            infile_type=infile_type, output_type=output_type,
                                            query_species=query_species,
                                            fw_search_type=fw_search_type, rv_search_type=rv_search_type,
                                            fw_search_local=fw_search_local, rv_search_local=rv_search_local,
                                            rv_search_db=rv_search_db, expect=expect, perc_score=perc_score,
                                            perc_ident=perc_ident, perc_span=perc_span,
                                            translate_hit_name=translate_hit_name,
                                            perc_length=perc_length, megablast=megablast, email=email, id_type=id_type,
                                            id_source=id_source, id_db=id_db, fetch_batch_size=fetch_batch_size,
                                            passwd=passwd, reciprocal_method=reciprocal_method,
                                            id_db_version=id_db_version, verbose=verbose, n_threads=n_threads,
                                            host=host, outfolder=outfolder, indent=indent+1,
                                            hit_name_only=hit_name_only,
                                            user=user, driver=driver, write_intermediates=write_intermediates,
                                            fw_search_kwargs=fw_search_kwargs, rv_search_kwargs=rv_search_kwargs)
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
            # Initiate RecBlastRun Tasks
            for rec in rec_handle:
                if verbose > 2:
                    print('Sequence: ', rec.name)
                rb_queue.put(RecBlastRun(seq_record=rec, target_species=species))
    else:
        for rec in rec_handle:
            if verbose > 2:
                print('Sequence: ', rec.name)
            rb_queue.put(RecBlastRun(seq_record=rec, target_species=target_species))
    #########################################################################

    # Drop poison pills in queue
    for _ in rec_blast_instances:
        rb_queue.put(None)
    #########################################################################

    # Collect results
    recblast_out = list()
    write_args = {}
    if output_type == 'database':
        output_type = 'sqlite3'
    if output_type in ['sqlite3']:
        sql_file = outfolder.joinpath(''.join(run_name if run_name != 'default' else 'RecBlastOutput') + '.sqlite')
        write_args['table_name'] = ''.join(run_name if run_name != 'default'
                                           else 'RecBlastOutput').replace('-', '_')
        write_args['max_hit_col'] = 1
        write_args['col_to_make'] = 0
    write_args['row_number'] = 1
    while n_jobs:
        try:
            recblast_out.append(rb_results.get())
        except Exception as err:
            print(type(err), err)
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
                    print(RecBlastWriteError(err))
            if min_mem:
                recblast_out[-1] = 1
    #########################################################################

    # Kill any living threads at this point, for good measure
    for rcb in rec_blast_instances:
        if rcb.is_alive():
            rcb.join()
    #########################################################################

    # Combine all the rc_out records into one big one and return it:
    recblast_out = sum([rc for rc in recblast_out if rc != {}])
    progbar.done()
    return recblast_out
    #########################################################################


class RecSearch(object):
    def __init__(self, target_species, query_species, forward_search_type, reverse_search_type, verbose, **kwargs):
        # GLOBAL OPTIONS
        self.search_types = ("blast", "blat")
        self.start_servers = False
        self.max_processes = multiprocessing.cpu_count()/2
        self.records = []
        self.search_parameters = {"Forward": {}, "Reverse": {}}
        # Other parameters
        self.target_species = target_species
        self.query_species = query_species
        assert isinstance(target_species, list) or isinstance(target_species, str), "Target species must either be" \
                                                                                    "a list of target species names, " \
                                                                                    "or a single string indicating " \
                                                                                    "the target species name"
        self.verbose = verbose.lower().count('v') if isinstance(verbose, str) else verbose
        assert isinstance(verbose, int), 'Verbose was of type {}; must be either be an integer greater ' \
                                         'than or equal to zero, or a number of v\'s equal ' \
                                         'to the desired level of verbosity'.format(type(verbose))
        self.fw_search_type=forward_search_type
        self.rv_search_type=reverse_search_type
        assert any((search in self.fw_search_type for
                    search in self.search_types)), "Forward search type {0} is invalid, currently supported" \
                                                   "search types are: {1}.".format(self.fw_search_type,
                                                                                   ', '.join(self.verbose))
        assert any((search in self.rv_search_type for
                    search in self.search_types)), "Reverse search type {0} is invalid, currently supported" \
                                                   "search types are: {1}.".format(self.rv_search_type,
                                                                                   ', '.join(self.verbose))
        self.sequence_source = kwargs.pop("sequence_db", forward_search_type)

    def queries(self, *queries, infile_type):
        for query in queries:
            if isinstance(query, list):
                self.queries(*query)
                continue
            elif isinstance(query, str):
                if query.startswith(">"):
                    self.records += list(SeqIO.parse(StringIO(query), "fasta"))
                else:
                    self.queries(Path(query))
                continue
            elif isinstance(query, Path):
                assert query.exists(), "Sequence File {} does not exist!".format(str(query))
                assert query.is_file(), "Sequence File {} is not a file!".format(str(query))
                self.records += [i for i in SeqIO.parse(str(query), infile_type)]
                continue
            else:
                raise RecBlastException("Invalid Query type!")

    def
    def _calc_processes(self):
        # Calculation of processes to run
        if self.verbose > 1:
            print('Automatically calculating n_processes:')
        n_species = len(self.target_species) if isinstance(self.target_species, list) else 1
        n_rec = len(self.records)
        if n_rec < 1:
            raise RecBlastException("No query records have been set! Please use self.queries() to set query records.")
        n_jobs = n_rec * n_species
        if n_jobs > self.max_processes:
            self.n_processes = self.max_processes
        else:
            n_processes = n_jobs
            if self.verbose:
                print('Number of processes to be made: ', n_processes)


    def _blat_server_check_(self, target, direction, search_type, search_db, host, id_db_path, indent):
        assert isinstance(search_db, dict), "{0} search of type {1} must be a dictionary " \
                                            "of species-port key pairs".format(direction, search_type)
        server_activated = {}

        if search_db[target] == 'auto':
            rv_server_online = False
        else:
            rv_server_online = blat_server('auto', order='status', host=host, port=search_db[target],
                                           species=target, id_db_path=id_db_path, verbose=self.verbose,
                                           indent=1, type=search_type)
        if rv_server_online:
            return True
        else:
            print("BLAT server was offline!", indent=indent)
            if self.start_servers:
                print("Attempting to start server...", indent=indent)
                try:
                    search_db[target] = blat_server('auto', 'start', host=host, port=30000, species=target,
                                                    id_db_path=id_db_path, verbose=self.verbose,
                                                    indent=1, type=search_type)
                    return True
                except TimeoutError as err:
                    print("{0}: BLAT server start-up function timed out while starting server. "
                          "Full error message: \n".format(type(err)))
                    return False
                except BLATServerError as err:
                    print("{0}: BLAT server start-up failed! Full error message: \n".format(type(err)))
                    return False
                except Exception:
                    raise
            else:
                return False

    def __call__(self):
        assert len(self.records) > 0, "No query records have been set! Please use self.queries() to set query records."
        if self.verbose:
            print('RecBlast version: ', __version__)
            print('Using BioPython version: ', bp_version)
            print('Beginning RecBlastMP!')
        if self.verbose == 1:
            print('Basic verbose mode active. Will print only essential commentary.')
        elif self.verbose == 2:
            print('Verbose mode was set to 2. Will elaborate considerably about the script.')
        elif self.verbose == 3:
            print(
                'Debugging-level verbose mode set. You will be innunadated by text. '
                'Brace yourself, and hold on to your console.')
        elif self.verbose == 50:
            print("V FOR VERBOSE: \n"
                  "\"Voilà! In view, a humble vaudevillian veteran cast vicariously as both victim and villain by \n"
                  "the vicissitudes of Fate. This visage, no mere veneer of vanity, is a vestige of the vox populi, \n"
                  "now vacant, vanished. However, this valourous visitation of a bygone vexation stands vivified and \n"
                  "has vowed to vanquish these venal and virulent vermin vanguarding vice and vouchsafing the \n"
                  "violently vicious and voracious violation of volition! The only verdict is vengeance; a vendetta \n"
                  "held as a votive, not in vain, for the value and veracity of such shall one day vindicate the \n"
                  "vigilant and the virtuous. \n"
                  "Verily, this vichyssoise of verbiage veers most verbose, so let me simply add that it's my very \n"
                  "good honour to meet you and you may call me [Reciprocal-Best-Hit-BLAST Script].\" \n"
                  "\t - V \n"
                  "Moore, Alan, David Lloyd, Steve Whitaker, and Siobhan Dodds. V for Vendetta. New York: DC Comics, "
                  "2005.")
        if self.verbose > 3:
            multiprocessing.log_to_stderr(logging.DEBUG)

        # Converting perc_ident to integer because that's how these programs roll
        perc_ident = perc_ident * 100
        #########################################################################
        # Multiprocessing set-up
        self._calc_processes()
        if self.verbose > 1:
            print('Creating queues... ', end='')
        rb_queue = multiprocessing.JoinableQueue()
        rb_results = multiprocessing.Queue()
        if self.verbose > 1:
            print('Done!')
        #########################################################################
