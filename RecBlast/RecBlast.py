import logging
from contextlib import redirect_stdout
from datetime import datetime
from pathlib import Path
from copy import deepcopy
from itertools import product

import multiprocess as multiprocessing

from Bio import SeqIO
from Bio import __version__ as bp_version
from io import StringIO
from Bio.SeqRecord import SeqRecord


from RecBlast import print, ProgressBar, __version__
from RecBlast.FetchSeq import fetchseq
from RecBlast.Search import Search, get_searchdb, blat_server, id_ranker
from RecBlast.RBC import RecBlastContainer
from RecBlast.WarningsExceptions import *
from RecBlast.Auxilliary import translate_annotation


class RecBlastMPThread(multiprocessing.Process):
    """
    RecBlastMPThread_Handle is the first branch to be made. It will perform the actual RecBlastRun.
    """

    def __init__(self, proc_id, rb_queue, rb_results_queue, query_species, forward_search_type, forward_search_criteria,
                 forward_search_settings, sequence_source_type, sequence_source_settings, reverse_search_type,
                 reverse_search_criteria, reverse_search_settings, reciprocal_method, output_type, outfolder,
                 translate_hit_name, verbose, memory_saver_level, start_rb, stop_rb, indent=0):
        multiprocessing.Process.__init__(self)
        self.name = proc_id
        self.rb_queue = rb_queue
        self.rb_results_queue = rb_results_queue
        self.query_species = query_species
        self.forward_search_type = forward_search_type
        self.forward_search_criteria = forward_search_criteria
        self.forward_search_settings = forward_search_settings
        self.sequence_source_type = sequence_source_type
        self.sequence_source_settings = sequence_source_settings
        self.reverse_search_type = reverse_search_type
        self.reverse_search_criteria = reverse_search_criteria
        self.reverse_search_settings = reverse_search_settings
        self.reciprocal_method = reciprocal_method
        self.output_type = output_type
        self.outfolder = outfolder
        self.translate_hit_name = translate_hit_name
        self.verbose = verbose
        self.memory_saver_level = memory_saver_level
        self.start_rb = start_rb
        self.stop_rb = stop_rb
        self.indent = indent

    def run(self):
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
                    output = rb_instance(proc_id=self.name,
                                         query_species=self.query_species,
                                         forward_search_type=self.forward_search_type,
                                         forward_search_criteria=self.forward_search_criteria,
                                         forward_search_settings=self.forward_search_settings,
                                         sequence_source_type=self.sequence_source_type,
                                         sequence_source_settings=self.sequence_source_settings,
                                         reverse_search_type=self.reverse_search_type,
                                         reverse_search_criteria=self.reverse_search_criteria,
                                         reverse_search_settings=self.reverse_search_settings,
                                         reciprocal_method=self.reciprocal_method,
                                         output_type=self.output_type,
                                         outfolder=self.outfolder,
                                         translate_hit_name=self.translate_hit_name,
                                         verbose=self.verbose,
                                         start_rb=self.start_rb,
                                         stop_rb=self.stop_rb,
                                         indent=self.indent,
                                         memory_saver_level=self.memory_saver_level)
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

    def _set_output_paths_(self, query_species, target_species, forward_search_type, output_type):
        return dict(
            forward_search_output=Path("{0}_recblast_out".format(target_species).replace(' ', '_'),
                                       "{0}_{1}_tmp".format(forward_search_type,
                                                            self.seq_record.id
                                                            ).replace(' ', '_')),
            forward_id_score_output=Path("{0}_recblast_out".format(target_species).replace(' ', '_'),
                                         "{0}_{1}_tmp".format(forward_search_type,
                                                              self.seq_record.id
                                                              ).replace(' ', '_'),
                                         "{0}_{1}_{2}_to_{3}.ID_Scores".format(forward_search_type,
                                                                               self.seq_record.id,
                                                                               query_species,
                                                                               target_species
                                                                               ).replace(' ', '_')),
            reverse_search_output=[],
            recblast_output_unanno=Path("{0}_recblast_out".format(target_species).replace(' ', '_'),
                                        "{0}_{1}_tmp".format(forward_search_type,
                                                             self.seq_record.id
                                                             ).replace(' ', '_'),
                                        "unannotated_{0}_{1}.fasta".format(forward_search_type,
                                                                           self.seq_record.id
                                                                           ).replace(' ', '_')),
            recblast_output=Path("{0}_recblast_out".format(target_species).replace(' ', '_'),
                                 "{0}_{1}.{2}".format(forward_search_type,
                                                      self.seq_record.id,
                                                      output_type
                                                      ).replace(' ', '_')),
            search_nohits=Path("{0}_recblast_out".format(target_species).replace(' ', '_'),
                               "{0}_{1}.no-hits".format(forward_search_type,
                                                        self.seq_record.id
                                                        ).replace(' ', '_'))
        )

    def _alt_start_(self, target_species, start_rb):
        if start_rb and start_rb[0].lower() in ["forward", "forward_search"]:
            try:
                output = (Search.load(Path(start_rb[1], "{0}_{1}.{2}".format(target_species.replace(' ', "_"),
                                                                             self.seq_record.name, start_rb[2]))),
                          None)
            except Exception as err:
                output = (None, err)
            return output
        elif start_rb[0].lower() in ["forward_id_ranker"]:
            id_ranked_path = Path(start_rb[1], "{0}_{1}.{2}".format(target_species.replace(' ', "_"),
                                                                 self.seq_record.name, start_rb[2]))
            f_id_ranked = []
            with id_ranked_path.open() as id_ranked_input:
                tmp_ids = (line.strip().split('\t') for line in id_ranked_input)
                for i, h_chr, h_start, h_end, h_name, h_score, h_strand, h_coverage in enumerate(tmp_ids):
                    try:
                        # Note: output is list of (seq_chr, seq_range, seq_name, seq_score, strand, seq_coverage)
                        f_id_ranked.append((h_chr, (h_start, h_end), h_name, h_score, h_strand, h_coverage))
                    except ValueError:
                        message = "In file {0}, line #{1}, too few items on line! Confirm that the input file is in " \
                                  "BED6 format and try again!".format(str(id_ranked_path), i)
                        print(message)
                        raise ValueError(message)
            return f_id_ranked
        elif start_rb[0].lower() in ['fetchseq']:
            seq_path = Path(start_rb[1], "{0}_{1}.{2}".format(target_species.replace(' ', "_"),
                                                           self.seq_record.name, start_rb[2]))
            if 'fa' in start_rb[2]:
                seq_type = 'fasta'
            elif start_rb[2] in ['gb', 'genbank']:
                seq_type = 'gb'
            else:
                seq_type = start_rb[2]
            seq_dict, missing_items = {SeqIO.to_dict(SeqIO.parse(str(seq_path), seq_type))}, []
            return seq_dict, missing_items

    def __call__(self, proc_id, query_species, forward_search_type, forward_search_criteria,
                 forward_search_settings, sequence_source_type, sequence_source_settings, reverse_search_type,
                 reverse_search_criteria, reverse_search_settings, reciprocal_method, output_type, outfolder,
                 translate_hit_name, verbose, memory_saver_level, indent, start_rb=None, stop_rb=None):
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
                print(kwarg, ':\t', value, indent=indent + 1)
        indent += 1
        if verbose > 1:
            print('Creating handles for intermediary outputs...', indent=indent)
        # Handle for output paths:
        rc_container['output_paths'] = self._set_output_paths_(query_species, target_species,
                                                               forward_search_type,
                                                               output_type)
        output_paths = rc_container['output_paths']
        forward_search = rc_container['forward_search']

        # FORWARD SEARCH
        try:
            if verbose:
                message = "Start was set to 'forward_search', fetching pre-generated search results from folder " \
                          "'{}'!".format(start_rb[1]) if start_rb and start_rb[0].lower() in ["forward", 
                                                                                              "forward_search"] else \
                    "Performing Forward search for {}... ".format(self.seq_record.id)
                print(message, indent=indent)
            fw_search = Search(search_type=forward_search_type)
            forward_search_params = {}
            forward_search_params.update(forward_search_settings)
            forward_search_params.update(forward_search_criteria)
            if start_rb and start_rb[0].lower() in ["forward", "forward_search"]:
                forwardsearchrecord, search_err = self._alt_start_(target_species, start_rb)
            else:
                forwardsearchrecord, search_err = fw_search(seq_record=self.seq_record, species=target_species,
                                                            search_type=forward_search_type,
                                                            blastoutput_custom=output_paths['forward_search_output'],
                                                            verbose=verbose, indent=indent + 1,
                                                            **forward_search_params)
        except ValueError:
            rc_container['recblast_unanno'] = [SeqRecord('')]
            return rc_container_full
        except Exception as err:
            print('WARNING! UNCATCHED EXCEPTION OCCURED!')
            print(type(err), err)
            return rc_container_full
        if search_err:
            forward_search['search_errors'] = search_err
            print('Forward search returned with an error!', search_err)
            return rc_container_full
        forward_search['search_results'] = forwardsearchrecord
        if not forwardsearchrecord:
            print('Forward search record was empty!')
            return rc_container_full
        if verbose:
            print('Forward search done!', indent=indent)
        if stop_rb and stop_rb.lower() in ["forward", "forward_search"]:
            if verbose:
                print("Stop was set to {}, ending run and returning results!".format(stop_rb))
            return rc_container_full

        # ID_RANKER
        if verbose:
            message = "Start was set to 'forward_search', fetching pre-generated search results from folder '{}'" \
                      "!".format(start_rb[1]) if start_rb and start_rb[0].lower() in ['forward_id_ranker'] else \
                      'Culling results based on given criteria...'
            print(message, indent=indent)
        try:
            if start_rb and start_rb[0].lower() in ['forward_id_ranker']:
                f_id_ranked = self._alt_start_(target_species, start_rb)
            else:
                # Note: output is list of (seq_chr, seq_range, seq_name, seq_score, strand, seq_coverage)
                f_id_ranked = id_ranker(forwardsearchrecord, verbose=verbose,
                                        indent=indent + 1, **forward_search_criteria)

        except Exception as err:
            print('WARNING! UNCATCHED ERROR IN ID_RANKER!')
            print(type(err), err)
            return rc_container_full
        if verbose:
            print('Done!', indent=indent)
        f_id_out_list = ['{chr}\t{s}\t{e}\t{n}\t{sc}\t{st}\n'.format(chr=id_i[0], s=id_i[1][0],
                                                                     e=id_i[1][1], n=id_i[2],
                                                                     sc=id_i[3], st=id_i[4]) for id_i in f_id_ranked]
        forward_ids = rc_container['forward_ids']
        forward_ids['ids'] = f_id_ranked
        forward_ids['pretty_ids'] = f_id_out_list
        if not f_id_ranked:
            print('Forward search for query {} yielded no hits!'.format(self.seq_record.id))
            return rc_container_full
        if stop_rb and stop_rb.lower() == "forward_id_ranker":
            if verbose:
                print("Stop was set to {}, ending run and returning results!".format(stop_rb))
            return rc_container_full

        # FETCHSEQ
        try:
            if verbose:
                message = "Start was set to 'fetchseq', reading in sequences from folder '{}'" \
                          "!".format(start_rb[1]) if start_rb and start_rb[0].lower() in ['forward_id_ranker'] else \
                    'Beginning Fetchseq!'
                print(message, indent=indent)

            if isinstance(sequence_source_settings['database'], dict):
                try:
                    sequence_source_settings['database'] = sequence_source_settings['database'][target_species]
                except KeyError:
                    raise DatabaseNotFoundError('No sequence source database for species {} '
                                                'was found in the provided dict!'.format(target_species))
            if start_rb and start_rb[0].lower() in ['fetchseq']:
                seq_dict, missing_items = self._alt_start_(target_species, start_rb)
            else:
                # Note: BioSQL is NOT thread-safe as implemented in fetchseq!
                # Throws tons of errors if executed with more than one thread!
                seq_dict, missing_items = fetchseq(ids=f_id_ranked, species=target_species, delim='\t',
                                                   source=sequence_source_type,
                                                   output_type=output_type, output_name='',
                                                   verbose=verbose, n_threads=1, indent=indent + 1,
                                                   n_subthreads=1, **sequence_source_settings)
            if verbose:
                print('Done with fetching!', indent=indent)
        except IndexError:
            print('WARNING! FETCHSEQ FAILED! ENDING THIS RUN!')
            return rc_container_full

        if not missing_items:
            forward_ids['missing_items'] = missing_items
            print('Items were missing!', indent=indent)
            for i in missing_items:
                print(i, indent=indent + 1)

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
        if stop_rb and stop_rb.lower() in ["fetchseq"]:
            if verbose:
                print("Stop was set to {}, ending run and returning results!".format(stop_rb))
            return rc_container_full

        # REVERSE SEARCH
        if verbose:
            print('Preparing for Reverse search...', indent=indent)
        for index, entry_record in enumerate(recblast_sequence):
            assert isinstance(entry_record, SeqRecord), 'Warning! entry_record is of type {} ' \
                                                        'and not a SeqRecord!'.format(str(type(entry_record)))
            if verbose:
                print("Entry {} in unannotated RecBlast Hits:\n".format(entry_record.id), indent=indent + 1)
                for item in [entry_record.id, entry_record.description,
                             entry_record.seq[0:10] + '...' + entry_record.seq[-1]]:
                    print(item, indent=indent + 2)
            output_paths['reverse_search_output'].append(
                Path("{0}_recblast_out".format(target_species).replace(' ', '_') +
                     '/' +
                     "{0}_{1}_tmp".format(reverse_search_type,
                                          self.seq_record.id
                                          ).replace(' ', '_') +
                     '/' +
                     "{0}_{1}_{3}_to_{2}_{4}.xml".format(reverse_search_type,
                                                         self.seq_record.id,
                                                         query_species,
                                                         target_species,
                                                         entry_record.name
                                                         ).replace(' ', '_')
                     ))
            if verbose:
                print("Performing Reverse search for {}... ".format(entry_record.name), indent=indent)
            reverse_search = rc_container['reverse_search']
            try:
                rv_search = Search(search_type=reverse_search_type)
                reverse_search_params = {}
                reverse_search_params.update(reverse_search_settings)
                reverse_search_params.update(reverse_search_criteria)
                reversesearchrecord, search_err = rv_search(seq_record=entry_record, species=query_species,
                                                            search_type=reverse_search_type,
                                                            blastoutput_custom=output_paths['reverse_search_output'],
                                                            verbose=verbose, indent=indent + 1,
                                                            **reverse_search_params)
            except ValueError:
                rc_container['recblast_unanno'] = [SeqRecord('')]
                return rc_container_full
            except Exception as err:
                print('WARNING! UNCATCHED EXCEPTION OCCURED!')
                print(type(err), err)
                return rc_container_full
            if search_err:
                print('Reverse search returned with errors!')
                reverse_search['search_errors'] = search_err
                return rc_container_full
            reverse_search['search_results'] = reversesearchrecord
            if verbose:
                print('Done with Reverse search!', indent=indent)
            if stop_rb and stop_rb.lower() in ["reverse", "reverse_search"]:
                if reciprocal_method in ['1:1', '1-to-1']:
                    if verbose:
                        print("Stop was set to {}, and reciprocal_method was set to '1:1', "
                              "ending run now!".format(stop_rb))
                    return rc_container_full
                else:
                    if verbose:
                        print("Stop was set to {}, skipping annotation and continuing to next hit!".format(stop_rb))
                    continue

            # Rank reverse hits
            if verbose:
                print('Culling results using reverse search criteria...', indent=indent)
            try:
                reverse_hits = id_ranker(reversesearchrecord, verbose=verbose,
                                         indent=indent + 1, method=reciprocal_method, **reverse_search_criteria)
            # Note: output is list of (seq_chr, seq_range, seq_name, seq_score, strand, seq_coverage)
            # old f_id_ranked: (seq_name, seq_range, strand, seq_score, seq_coverage)
            except Exception as err:
                print('No Reverse search hits were found for this hit!', indent=indent + 1)
                print('Continuing to next Sequence!', indent=indent + 1)
                continue
            print('Reverse search hits:', indent=indent + 1)
            print(reverse_hits, indent=indent + 2)
            reverse_search_annotations = []
            if stop_rb and stop_rb.lower() in ["reverse_id_ranker"]:
                if reciprocal_method in ['1:1', '1-to-1']:
                    if verbose:
                        print("Stop was set to {}, and reciprocal_method was set to '1:1', "
                              "ending run now!".format(stop_rb))
                    return rc_container_full
                else:
                    if verbose:
                        print("Stop was set to {}, skipping annotation and continuing to next hit!".format(stop_rb))
                    continue
            # ANNOTATION
            for anno in reverse_hits:
                try:
                    if translate_hit_name:
                        new_anno = translate_annotation(anno[0])  # Big note: this will fail if the "chr" category
                        # (ie the hit_id of the search, ie the header of the
                        # fasta file) is NOT a RefSeq ID!!
                        # TODO: add parameter to RecSearch to allow specification of translate_annotation kwargs!
                        # TODO: even better, allow the user to explicitly specify a function here for custom translation!
                except Exception as err:
                    new_anno = anno[0]
                finally:
                    reverse_search_annotations.append('\t |[ {0}:{1}-{2}{3} ]|'.format(new_anno, anno[1][0],
                                                                                       anno[1][1], anno[4]))
                if reciprocal_method in ['1-to-1', '1:1', 'best-hit', 'best hit']:
                    print('Best Hit Reciprocal BLAST was selected, will only use top hit for annotation!',
                          indent=indent)
                    continue
            if not reverse_search_annotations:
                print('No Reverse search hits were found for this hit!', indent=indent + 1)
                print('Continuing to next Sequence!', indent=indent + 1)
                continue
            else:
                if verbose > 1:
                    print('Done. Annotating RecBlast Hits:', indent=indent + 1)
            reverse_ids = rc_container['reverse_ids']
            reverse_ids['ids'] = reverse_search_annotations
            if reciprocal_method in ['1-to-1', '1:1']:
                print('1:1 Best Hit Reciprocal BLAST was selected, ending run after first hit!',
                      indent=indent)
                entry_record.description += '|-|' + reverse_search_annotations[0] if \
                    isinstance(reverse_search_annotations, list) \
                    else reverse_search_annotations
                break
            else:
                entry_record.description += '|-|'.join(reverse_search_annotations) if \
                    isinstance(reverse_search_annotations, list) \
                    else '|-|' + reverse_search_annotations
            if verbose > 3:
                print(entry_record, indent=indent + 2)

        if not isinstance(recblast_sequence, list):
            recblast_sequence = [recblast_sequence]
        if recblast_sequence == list():
            recblast_sequence.append(SeqRecord(''))
        rc_container['recblast_results'] = recblast_sequence
        if memory_saver_level > 1:
            rc_container_full = (target_species, self.seq_record.id,
                                 ["{h_chr}\t{h_start}\t{h_end}\t"
                                  "{h_query}\t{h_score}\t{h_strand}\t"
                                  "{h_anno}\t"
                                  "{h_coverage}".format(h_chr=rec.name,
                                                        h_start=rec.features[0].location.start,
                                                        h_end=rec.features[0].location.end,
                                                        h_query=self.seq_record.name,
                                                        h_score=str(rec.features[0].qualifiers['score']),
                                                        h_strand=rec.features[0].location.strand,
                                                        h_anno=rec.description.split('|[')[1].rstrip(']|').strip(),
                                                        h_coverage=rec.features[0].qualifiers['query_coverage']
                                                    ) for rec in recblast_sequence])
        run_time = str(datetime.now() - self.starttime)
        rc_container['run_time'] = run_time
        print('PID', 'SeqName', sep='\t')
        print(proc_id, self.seq_record.id, sep='\t')
        print('Run time: ', run_time)
        return rc_container_full
        # def run(self):


class RecSearch(object):
    """Performs a Reciprocal Search of sequences between a list of target species and a query species.

    The RecSearch object can be used to initiate and run various separate reciprocal searches. The

    Usage:
    >>> recblast = RecSearch(target_species="Species A", query_species="Homo sapiens",
    ...                      forward_search_type="tblat", reverse_search_type="blat-transcript",
    ...                      sequence_source="twobit", verbose=0)
    >>> recblast.set_queries("./Test/query/3protein.fasta", infile_type="fasta")
    >>> recblast.forward_search_settings['database'] = {"Loxodonta africana": "loxAfr3.2bit"}
    >>> recblast.sequence_source_settings['database'] = {"Loxodonta africana": "loxAfr3.2bit"}
    >>> recblast.forward_search_settings['database_port'] = {"Loxodonta africana":20008}
    >>> recblast.reverse_search_settings['database_port'] = {"Homo sapiens":20005}
    >>> rc = recblast("test", output_type=None, stop="forward_id_ranker")

    """
    def __init__(self, target_species, query_species, forward_search_type, reverse_search_type, sequence_source,
                 verbose):
        """

        :param target_species:
        :param query_species:
        :param forward_search_type:
        :param reverse_search_type:
        :param sequence_source:
        :param verbose:
        """
        # GLOBAL OPTIONS
        self.search_types = (
            'blastn',
            'blastn-transcript',
            'blastp',
            'blastx',
            'blat',
            'blat-transcript',
            'tblastn',
            'tblastn-transcript',
            'tblastx',
            'tblastx-transcript',
            'tblat',
            'tblat-transcript')
        self.sequence_sources = (
            'twobit',
            'sql',
            'ensemble',
            'fasta'
        )
        self.memory_saver_level = 0
        self.start_servers = False
        self.translate_hit_name = True
        self.max_processes = multiprocessing.cpu_count() / 2
        self.records = []
        self.search_parameters = {"Forward": {}, "Reverse": {}}
        # Other parameters
        self.target_species = target_species if isinstance(target_species, list) else [target_species]
        self.query_species = query_species
        self.verbose = verbose.lower().count('v') if isinstance(verbose, str) else verbose
        self.forward_search_type = forward_search_type.lower()
        self.reverse_search_type = reverse_search_type.lower()
        self._set_search_settings_()
        self.sequence_source = sequence_source
        self._set_sequence_source_settings_()

        ### Tests:
        assert isinstance(self.target_species, list), "Target species must either be a list of target species names, " \
                                                      "or a single string indicating the target species name"
        assert all((isinstance(species, str) for species in self.target_species)), "All target species must be strings!"

        assert isinstance(verbose, int), 'Verbose was of type {}; must be either be an integer greater ' \
                                         'than or equal to zero, or a number of v\'s equal ' \
                                         'to the desired level of verbosity'.format(type(verbose))

        assert self.forward_search_type in self.search_types, "Forward search type {0} is invalid, currently " \
                                                              "supported search types are: {1}.".format(
            self.forward_search_type,
            ', '.join(self.verbose))
        assert self.reverse_search_type in self.search_types, "Reverse search type {0} is invalid, currently " \
                                                              "supported search types are: {1}.".format(
            self.reverse_search_type,
            ', '.join(self.verbose))

    def set_queries(self, *queries, infile_type=None):
        infile_legal_types = ['fasta', 'fa', 'genbank', 'gb']
        if infile_type:
            assert infile_type.lower() in infile_legal_types, "infile_type must be one of the following" \
                                                              "filetypes: {}".format(", ".join(infile_legal_types))
        for query in queries:
            if isinstance(query, list):
                self.set_queries(*query, infile_type=infile_type.lower())
                continue
            elif isinstance(query, str):
                if query.startswith(">"):
                    self.records += list(SeqIO.parse(StringIO(query), "fasta"))
                else:
                    self.set_queries(Path(query), infile_type=infile_type.lower())
                continue
            elif isinstance(query, Path):
                assert query.exists(), "Sequence File {} does not exist!".format(str(query))
                assert query.is_file(), "Sequence File {} is not a file!".format(str(query))
                self.records += list(SeqIO.parse(str(query), infile_type.lower()))
                continue
            elif isinstance(query, SeqRecord):
                self.records.append(SeqRecord)
            else:
                raise RecBlastException("Invalid Query type!")

    def _calc_processes_(self):
        # Calculation of processes to run
        if self.verbose > 1:
            print('Automatically calculating n_processes:')
        n_species = len(self.target_species) if isinstance(self.target_species, list) else 1
        n_rec = len(self.records)
        if n_rec < 1:
            raise RecBlastException("No query records have been set! Please use self.queries() to set query records.")
        self.n_jobs = n_rec * n_species
        if self.n_jobs > self.max_processes:
            self.n_processes = self.max_processes
        else:
            self.n_processes = self.n_jobs
            if self.verbose:
                print('Number of processes to be made: ', self.n_processes)

    def _set_search_settings_(self):
        self.forward_search_criteria = dict(
            perc_score=0,
            perc_ident=0,
            perc_span=0,
            perc_length=0
        )
        if 'blat' in self.forward_search_type:
            self.forward_search_settings = dict(
                local="localhost",
                database_version="auto",
                database_path="/usr/db/BLAT",
                database="auto",
                database_port=dict()
            )
        elif 'blast' in self.forward_search_type:
            self.forward_search_settings = dict(
                local=True,
                email=None,
                megablast=True,
                database_version="auto",
                database_path="/usr/db/blastdb",
                database="auto"
            )
            self.forward_search_criteria["expect"] = 0
        else:
            raise NotImplementedError("Selection for forward search type {} "
                                      "is not yet implemented!".format(self.forward_search_type))
        self.reverse_search_criteria = dict(
            perc_score=0,
            perc_ident=0,
            perc_span=0,
            perc_length=0
        )
        if 'blat' in self.reverse_search_type:
            self.reverse_search_settings = dict(
                local="localhost",
                database_version="auto",
                database_path="/usr/db/BLAT",
                database="auto",
                database_port=dict()
            )
        elif 'blast' in self.reverse_search_type:
            self.reverse_search_settings = dict(
                local=True,
                email=None,
                megablast=True,
                database_version="auto",
                database_path="/usr/db/blastdb",
                database="auto"
            )
            self.reverse_search_criteria["expect"] = 0
        else:
            raise NotImplementedError("Selection for reverse search type {} "
                                      "is not yet implemented!".format(self.reverse_search_type))

    def _set_sequence_source_settings_(self):
        if self.sequence_source == 'sql':
            self.sequence_source_settings = dict(
                driver="psycopg2",
                host="localhost",
                database='bioseqdb',
                user='',
                password='',
                id_db_version="auto",
            )
        elif self.sequence_source == 'twobit':
            self.sequence_source_settings = dict(
                database_path="/usr/db/BLAT",
                database="auto"
            )
        elif self.sequence_source == 'fasta':
            self.sequence_source_settings = dict()
        elif self.sequence_source == 'ensemble':
            self.sequence_source_settings = dict(email="", batch_size=50)
        else:
            raise NotImplementedError("Selection for sequence source {} "
                                      "is not yet implemented!".format(self.sequence_source))
        self.sequence_source_settings.update({'id_type': "brute"})

    def _blat_server_check_(self, target, direction, search_type, database, database_port, host, database_path, indent):
        database_port = {target: database_port} if isinstance(database_port, int) else database_port
        assert isinstance(database_port, dict), "{0} search ports for search of type {1} must be a dictionary " \
                                                "of species-port key pairs. Provided object was of " \
                                                "type {2} with structure {3}".format(direction.capitalize(),
                                                                                     search_type,
                                                                                     type(database_port),
                                                                                     str(database_port))
        database = {target: database} if isinstance(database, str) else database
        assert isinstance(database, dict), "{0} search databases for search of type {1} must be a dictionary " \
                                           "of species-database key pairs. Provided object was of " \
                                           "type {2} with structure {3}".format(direction.capitalize(), search_type,
                                                                                type(database),
                                                                                str(database))
        server_activated = {}

        if database_port[target] == 'auto':
            server_online = False
        else:
            server_online = blat_server(database[target], order='status', host=host, port=database_port[target],
                                        species=target, database_path=database_path, verbose=self.verbose,
                                        indent=1, type=search_type)
        if server_online:
            return True
        else:
            print("BLAT server was offline!", indent=indent)
            if self.start_servers:
                print("Attempting to start server...", indent=indent)
                try:
                    database_port[target] = blat_server(database[target], 'start', host=host, port=30000,
                                                        species=target, database_path=database_path,
                                                        verbose=self.verbose, indent=1, type=search_type)
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

    def optimize_search_criteria(self, training_bed, smallest_step_size=0.1):
        # First, read the bed file and assemble the training set
        training_set = {}
        bedfile = Path(training_bed)
        assert bedfile.exists(), "Given bedfile path does not exist!"
        assert bedfile.is_file(), "Given bedfile path was not a file! Did you provide a directory?"
        with bedfile.open() as bed:
            hits = [line.strip().split('\t') for line in bed]
            for hit in hits:
                if hit[3] in training_set:
                    training_set[hit[3]].append(hit[0:3])
                else:
                    training_set[hit[3]] = hit[0:3]
        gene_records = {rec.name: rec for rec in self.records}
        for gene in training_set:
            assert gene in gene_records, "Training set contains the gene {}, but there is " \
                                         "no associated query with the same name!".format(gene)
        # Deep copy the current recblast for training
        test_recblast = deepcopy(self)
        test_recblast.set_queries((gene_records[gene] for gene in training_set))

        # Here's the magic: we're gonna run this in two phases: first, optimize parameters in the
        # forward direction, without performing the reciprocal search.
        # Next, we'll select the parameter values in the forward direct that maximizes the number of correct hits
        # Finally, we'll then feed in the results from the first search to the program to do the reverse,
        # where we'll optimize the reverse parameters.

        step_size_list = [round(1/ss, 3) for ss in range(2, int(1/smallest_step_size)+1)]
        step_size_start = [round(i * step_size_list[0], 3) for i in range(0, int(1 / step_size_list[0]) + 1)]

        fw_test_results = {}
        for params in product(step_size_start, repeat=4):
            test_recblast.forward_search_criteria.update({'perc_ident': params[0],
                                                          'perc_length': params[1],
                                                          'perc_score': params[2],
                                                          'perc_span': params[3]})

            r_name = "test" + "_".join(params)
            fw_test_results[params] = [rec[2]['forward_ids'] for rec in
                                       test_recblast(r_name, output_type=None, stop="forward_id_ranker")]

    def __call__(self, run_name, reciprocal_method="best-hit", output_location="./", output_type="bed-min",
                 start=None, stop=None):
        assert len(self.records) > 0, "No query records have been set! " \
                                      "Please use self.set_queries() to set query records."
        if start:
            assert isinstance(start, tuple) and len(start) == 3, "If not set to 'None', parameter 'start' must be a " \
                                                                 "tuple of length 2, with the first item indicating the" \
                                                                 " starting point of the program; the second item " \
                                                                 "indicating the folder of the intermediary files " \
                                                                 "necessary to proceed; and the third item indicating" \
                                                                 "the extension of the file. Items in folder must be " \
                                                                 "named in the format of {species}_{query_name} in " \
                                                                 "order for the program to correctly identify them."
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
        if "perc_ident" in self.forward_search_criteria:
            self.forward_search_criteria['perc_ident'] = int(self.forward_search_criteria['perc_ident'] * 100)
        if "perc_ident" in self.reverse_search_criteria:
            self.reverse_search_criteria['perc_ident'] = int(self.reverse_search_criteria['perc_ident'] * 100)
        #########################################################################
        # Multiprocessing set-up
        self._calc_processes_()
        if self.verbose > 1:
            print('Creating queues... ', end='')
        rb_queue = multiprocessing.JoinableQueue()
        rb_results = multiprocessing.Queue()
        if self.verbose > 1:
            print('Done!')
        #########################################################################
        # Make output folder
        date_str = datetime.now().strftime('%y-%m-%d_%I-%M-%p')
        outfolder = Path(output_location, '{0}'.format(run_name if run_name else date_str))
        try:
            outfolder.mkdir(parents=True)
        except FileExistsError:
            pass
        #########################################################################
        # Get Database Files if set to "auto"
        self.forward_search_settings['database'] = {
            species: get_searchdb(search_type=self.forward_search_type,
                                  species=species,
                                  db_loc=self.forward_search_settings['database_path'],
                                  verbose=self.verbose, indent=1).name
            for species in self.target_species}
        self.reverse_search_settings['database'] = {self.query_species:
                                                        get_searchdb(search_type=self.reverse_search_type,
                                                                     species=self.query_species,
                                                                     db_loc=self.reverse_search_settings['database_'
                                                                                                         'path'],
                                                                     verbose=self.verbose, indent=1).name}
        # If sequence_source= twobit and database=auto, set database equal to forward search db
        if self.sequence_source in ["twobit", "2bit"]:
            if self.sequence_source_settings['database'] == 'auto':
                self.sequence_source_settings['database'] = self.forward_search_settings['database']
        #########################################################################
        # Check BLAT server:
        # Forward:
        for species in self.target_species:
            self._blat_server_check_(target=species, direction="forward", search_type=self.forward_search_type,
                                     database=self.forward_search_settings['database'],
                                     database_port=self.forward_search_settings['database_port'],
                                     host=self.forward_search_settings['local'],
                                     database_path=self.forward_search_settings['database_path'],
                                     indent=1)
        # Reverse:
        self._blat_server_check_(target=self.query_species, direction="reverse", search_type=self.reverse_search_type,
                                 database=self.reverse_search_settings['database'],
                                 database_port=self.reverse_search_settings['database_port'],
                                 host=self.reverse_search_settings['local'],
                                 database_path=self.reverse_search_settings['database_path'], indent=1)
        #########################################################################
        # RecBlast Thread init
        if self.verbose >= 1:
            print('Creating RecBlast Threads... ')
        rec_blast_instances = [RecBlastMPThread(proc_id=str(i + 1), rb_queue=rb_queue, rb_results_queue=rb_results,
                                                query_species=self.query_species,
                                                forward_search_type=self.forward_search_type,
                                                forward_search_criteria=self.forward_search_criteria,
                                                forward_search_settings=self.forward_search_settings,
                                                sequence_source_type=self.sequence_source,
                                                sequence_source_settings=self.sequence_source_settings,
                                                reverse_search_type=self.reverse_search_type,
                                                reverse_search_criteria=self.reverse_search_criteria,
                                                reverse_search_settings=self.reverse_search_settings,
                                                reciprocal_method=reciprocal_method, output_type=output_type,
                                                outfolder=outfolder, start_rb=start, stop_rb=stop,
                                                translate_hit_name=self.translate_hit_name, verbose=self.verbose,
                                                memory_saver_level=self.memory_saver_level)
                               for i in range(self.n_processes)]
        for rcb in rec_blast_instances:
            rcb.start()
        #########################################################################

        # Progress bar
        progbar = ProgressBar(self.n_jobs, fmt=ProgressBar.FULL)
        if self.verbose:
            progbar()
        #########################################################################
        # Load species list and add RecBlasts to queue
        for species in self.target_species:
            # Initiate RecBlastRun Tasks
            for rec in self.records:
                if self.verbose > 2:
                    print('Sequence: ', rec.name)
                rb_queue.put(RecBlastRun(seq_record=rec, target_species=species))
        #########################################################################
        # Drop poison pills in queue
        for _ in rec_blast_instances:
            rb_queue.put(None)
        #########################################################################
        # Collect results
        recblast_out = list()
        while self.n_jobs:
            try:
                recblast_out.append(rb_results.get())
                self.n_jobs -= 1
            except Exception as err:
                print(type(err), err)
            else:
                progbar.current += 1
                if self.verbose:
                    progbar()
                # Get rid of any blank dicts that like to pop up inside the RBC every once in a while
                if isinstance(recblast_out[-1], dict):
                    _ = recblast_out[-1].pop('__dict__', None)
                try:
                    if recblast_out[-1] == {}:
                        del recblast_out[-1]
                    else:
                        if isinstance(recblast_out[-1], dict) and 'error' in recblast_out[-1].keys():
                            print('RBC with an error was returned!')
                            print('Proc_ID: ', recblast_out[-1]['proc_id'])
                            print(recblast_out[-1]['error'])
                        else:
                            if self.memory_saver_level > 1:
                                if output_location == '/dev/stdout':
                                    outfile = Path(output_location)
                                else:
                                    outfile = Path(outfolder, recblast_out[-1][0]+".out")
                                with outfile.open("a") as of:
                                    of.write("# "+recblast_out[-1][1] + '\n')
                                    of.write('\n'.join(recblast_out[-1][2])+'\n')
                            elif output_type:
                                recblast_out[-1].write(filetype=output_type, file_loc=outfolder, verbose = self.verbose)
                    if self.memory_saver_level:
                        recblast_out[-1] = 1

                except Exception as err:
                    print('WARNING! Could not write output of RecBlast #{}'.format(len(recblast_out)))
                    print(RecBlastWriteError(err))
                    continue
        #########################################################################
        # Kill any living threads at this point, for good measure
        for rcb in rec_blast_instances:
            if rcb.is_alive():
                rcb.join()
        #########################################################################
        # Combine all the rc_out records into one big one and return it:
        recblast_out = sum([rc for rc in recblast_out if rc != {}])
        if self.verbose:
            progbar.done()
        return recblast_out
        #########################################################################

    def __repr__(self):
        return ('RecSearch(target_species={0}, query_species={1}, '
                'forward_search_type={2}, reverse_search_type={3}, '
                'sequence_source={4}, n_records={5} verbose={6})'.format(self.target_species, self.query_species, 
                                                                         self.forward_search_type, 
                                                                         self.reverse_search_type, self.sequence_source,
                                                                         len(self.records), self.verbose)
                )
    def __str__(self):
        s = "RecSearch Version {version}\n" \
            "Query Species:\n" \
            "\t{query_species}\n" \
            "Target Species:\n" \
            "\t{target_species}\n" \
            "Queries:\n" \
            "\t{queries}\n" \
            "Forward Search Type: {f_search_type}\n" \
            "\t Parameters:\n" \
            "\t\t{f_search_param}\n" \
            "\t Criteria:\n" \
            "\t\t{f_search_crit}\n" \
            "Reverse Search Type: {r_search_type}\n" \
            "\t Parameters:\n" \
            "\t\t{r_search_param}\n" \
            "\t Criteria:\n" \
            "\t\t{r_search_crit}\n" \
            "Sequence Database: {seq_db}\n" \
            "\t{seq_set}" \
            "".format(version=__version__,
                      query_species=self.query_species,
                      target_species="\n\t".join(self.target_species),
                      queries="\n\t".join(("{0} ({1})".format(i.name, len(i)) for i in self.records)),
                      f_search_type=self.forward_search_type,
                      f_search_param="\n\t\t".join(
                          ("{0}: {1}".format(k, v) for k, v in self.forward_search_settings.items())
                      ),
                      f_search_crit="\n\t\t".join(
                          ("{0}: {1}".format(k, v) for k, v in self.forward_search_criteria.items())
                      ),
                      r_search_type=self.reverse_search_type,
                      r_search_param="\n\t\t".join(
                          ("{0}: {1}".format(k, v) for k, v in self.reverse_search_settings.items())
                      ),
                      r_search_crit="\n\t\t".join(
                          ("{0}: {1}".format(k, v) for k, v in self.reverse_search_criteria.items())
                      ),
                      seq_db=self.sequence_source,
                      seq_set="\n\t".join(
                          ("{0}: {1}".format(k, v) for k, v in self.sequence_source_settings.items())
                      )
                      )
        return s
        