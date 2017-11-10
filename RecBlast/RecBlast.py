import logging
from contextlib import redirect_stdout
from datetime import datetime
from operator import itemgetter
from pathlib import Path

import multiprocess as multiprocessing

from Bio import SearchIO, SeqIO
from Bio import __version__ as bp_version
from Bio.Blast import NCBIXML
from Bio.SeqRecord import SeqRecord
from BioSQL import BioSeqDatabase

from RecBlast import print, ProgressBar, merge_ranges, flatten, __version__
from RecBlast.FetchSeq import fetchseq, biosql_get_sub_db_names
from RecBlast.Search import Search, get_searchdb, blat_server
from RecBlast.RBC import RecBlastContainer
from RecBlast.WarningsExceptions import *
from RecBlast.Auxilliary import translate_annotation, percent_identity_searchio


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


class RecBlastMPThread(multiprocessing.Process):
    """
    RecBlast_MP_Thread_Handle is the first branch to be made. It will perform the actual RecBlast.
    """

    def __init__(self, proc_id, rb_queue, rb_results_queue, fw_search_db, infile_type, output_type, search_db_loc,
                 query_species, fw_search_type, rv_search_type, fw_search_local, rv_search_local, rv_search_db, expect,
                 perc_score, perc_span, outfolder, indent, reciprocal_method, hit_name_only, translate_hit_name,
                 perc_ident, perc_length, megablast, email, id_type, fw_source, id_db, fetch_batch_size, passwd,
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
        self.fw_source = fw_source
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
        self.search_db_loc = search_db_loc
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
                    output = rb_instance(fw_search_db=self.fw_search_db, search_db_loc=self.search_db_loc,
                                         infile_type=self.infile_type, output_type=self.output_type,
                                         query_species=self.query_species,
                                         fw_search_type=self.fw_search_type, rv_search_type=self.rv_search_type,
                                         fw_search_local=self.fw_search_local,
                                         rv_search_local=self.rv_search_local,
                                         rv_search_db=self.rv_search_db, expect=self.expect, perc_score=self.perc_score,
                                         perc_ident=self.perc_ident, perc_span=self.perc_span,
                                         perc_length=self.perc_length, megablast=self.megablast, email=self.email,
                                         id_type=self.id_type,
                                         fw_source=self.fw_source, id_db=self.id_db, 
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


class RecBlast(object):
    def __init__(self, seq_record, target_species):
        self.starttime = datetime.now()
        self.seq_record = seq_record
        transtab = str.maketrans('!@#$%^&*();:.,\'\"/\\?<>|[]{}-=+', '_____________________________')
        self.seq_record.id = self.seq_record.id.translate(transtab)
        self.target_species = target_species

    def __call__(self, fw_search_db, infile_type, output_type, search_db_loc, reciprocal_method,
                 query_species, fw_search_type, rv_search_type, fw_search_local, rv_search_local, rv_search_db, expect,
                 perc_score, indent, hit_name_only, perc_span, translate_hit_name,
                 perc_ident, perc_length, megablast, email, id_type, fw_source, id_db, fetch_batch_size, passwd,
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
                                                      search_db_path=search_db_loc,
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
                blat_2bit = get_searchdb(search_type=fw_search_type, species=target_species, db_loc=search_db_loc,
                                         verbose=verbose, indent=indent+1)
            elif isinstance(id_db, str):
                blat_2bit = id_db
            else:
                raise FileNotFoundError('Invalid 2bit file designation!')
            id_db = Path(search_db_loc, blat_2bit+'.2bit').absolute()
            id_db = str(id_db) if id_db.is_file() else None
            if id_db is None:
                raise FileNotFoundError('Invalid 2bit file!')
            if verbose > 1:
                print(id_db)
        if 'sql' in fw_source.lower():
            server = BioSeqDatabase.open_database(driver=driver, user=user, passwd=passwd,
                                                  host=host, db=id_db)
        try:
            if fw_source == 'strict':
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
                                                   id_type='brute', server=server, source=fw_source,
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
                                                          search_db_path=search_db_loc,
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
                 id_type='brute', fw_source='sql', id_db='bioseqdb', fetch_batch_size=50,
                 passwd='', hit_name_only=False, min_mem=False,
                 id_db_version='auto', search_db_loc='/usr/db/search_db_loc', indent=0, translate_hit_name=True,
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
                             id_type='brute', fw_source='2bit', id_db='auto', fetch_batch_size=50,
                             passwd='', hit_name_only=True, min_mem=True,
                             id_db_version='auto', search_db_loc='/usr/db/BLAT', indent=0,
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
    :param fw_source:
    :param id_db:
    :param fetch_batch_size:
    :param passwd:
    :param translate_hit_name:
    :param hit_name_only:
    :param min_mem:
    :param id_db_version:
    :param search_db_loc:
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
    elif isinstance(verbose, int) and verbose > 0:
        pass
    else:
        raise TypeError('Verbose must be either be an integer greater than or equal to zero, or a number of v\'s equal '
                        'to the desired level of verbosity')
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
                                           species=query_species, search_db_loc=search_db_loc, verbose=verbose,
                                           indent=1, type=fw_search_type)
        if not rv_server_online and 'blat' in rv_search_type.lower():
            rv_search_db[query_species] = blat_server('auto', 'start', host=host, port=30000, species=query_species,
                                                      search_db_loc=search_db_loc, verbose=verbose,
                                                      indent=1, type=rv_search_type)
            server_activated[query_species] = rv_search_db[query_species]
        elif not rv_server_online and 'tblat' in rv_search_type.lower():
            rv_search_db[query_species] = blat_server('auto', 'start', host=host, port=30000, species=query_species,
                                                      search_db_loc=search_db_loc, verbose=verbose,
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
                                                   species=species, search_db_loc=search_db_loc, verbose=verbose,
                                                   indent=1, type=fw_search_type)
                    print('Status of the forward server: ', fw_server_online)
                if not fw_server_online and 'blat' in fw_search_type.lower():
                    print('Forward server was not online, starting!')
                    fw_search_db[species] = blat_server('auto', 'start', host=host, port=20000, species=species,
                                                        search_db_loc=search_db_loc, verbose=verbose,
                                                        indent=1, type=fw_search_type)
                    server_activated[species] = fw_search_db[species]
                elif not fw_server_online and 'tblat' in fw_search_type.lower():
                    print('Forward server was not online, starting!')
                    fw_search_db[species] = blat_server('auto', 'start', host=host, port=20000, species=species,
                                                        search_db_loc=search_db_loc, verbose=verbose,
                                                        indent=1, type=fw_search_type)
                    server_activated[species] = fw_search_db[species]
            print(fw_search_db)
        else:
            if fw_search_db[target_species] == 'auto':
                print('fw_search_db was set to "auto!"')
                fw_server_online = False
            else:
                fw_server_online = blat_server('auto', 'status', host=host, port=fw_search_db[target_species],
                                               species=target_species, search_db_loc=search_db_loc,
                                               verbose=verbose, indent=1, type=fw_search_type)
                print('Status of the forward server: ', fw_server_online)
            if not fw_server_online and 'blat' in fw_search_type.lower():
                print('Forward server was not online, starting!')
                fw_search_db[target_species] = blat_server('auto', 'start', host=host, port=20000,
                                                           species=target_species, search_db_loc=search_db_loc,
                                                           verbose=verbose, indent=1, type=fw_search_type)
                server_activated[target_species] = fw_search_db[target_species]
            elif not fw_server_online and 'tblat' in fw_search_type.lower():
                print('Forward server was not online, starting!')
                fw_search_db[target_species] = blat_server('auto', 'start', host=host, port=20000,
                                                           species=target_species, search_db_loc=search_db_loc,
                                                           verbose=verbose, indent=1, type=fw_search_type)
                server_activated[target_species] = fw_search_db[target_species]
            print(fw_search_db)
    #########################################################################

    # RecBlast Thread init
    if verbose >= 1:
        print('Creating RecBlast Threads... ')
    rec_blast_instances = [RecBlastMPThread(proc_id=str(i + 1), rb_queue=rb_queue, rb_results_queue=rb_results,
                                            fw_search_db=fw_search_db, search_db_loc=search_db_loc,
                                            infile_type=infile_type, output_type=output_type,
                                            query_species=query_species,
                                            fw_search_type=fw_search_type, rv_search_type=rv_search_type,
                                            fw_search_local=fw_search_local, rv_search_local=rv_search_local,
                                            rv_search_db=rv_search_db, expect=expect, perc_score=perc_score,
                                            perc_ident=perc_ident, perc_span=perc_span,
                                            translate_hit_name=translate_hit_name,
                                            perc_length=perc_length, megablast=megablast, email=email, id_type=id_type,
                                            fw_source=fw_source, id_db=id_db, fetch_batch_size=fetch_batch_size,
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
