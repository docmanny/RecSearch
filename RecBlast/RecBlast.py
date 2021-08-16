import logging
import warnings
from contextlib import redirect_stdout
from datetime import datetime
from pathlib import Path
from copy import deepcopy
from itertools import product
from operator import itemgetter

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
from RecBlast.Auxilliary import translate_ids


class RecBlastMPThread(multiprocessing.Process):
    """
    RecBlastMPThread_Handle is the first branch to be made. It will perform the actual RecBlastRun.
    """

    def __init__(self, proc_id, rb_queue, rb_results_queue, query_species, forward_search_type, forward_search_criteria,
                 forward_search_settings, sequence_source_type, sequence_source_settings, reverse_search_type,
                 reverse_search_criteria, reverse_search_settings, reciprocal_method, output_type, outfolder,
                 translate_annotation_params, verbose, memory_saver_level, indent=0):
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
        self.translate_annotation_params = translate_annotation_params
        self.verbose = verbose
        self.memory_saver_level = memory_saver_level
        self.indent = indent

    def run(self):
        if self.outfolder:
            master_out = self.outfolder.joinpath('Proc-{0}.log'.format(self.name)).absolute()
        else:
            master_out = Path("/dev/null")
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
                                         translate_annotation_params=self.translate_annotation_params,
                                         verbose=self.verbose,
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
    def __init__(self, seq_record, target_species, start_rb=None, stop_rb=None):
        self.starttime = datetime.now()
        self.seq_record = seq_record
        transtab = str.maketrans('!@#$%^&*();:.,\'\"/\\?<>|[]{}-=+', '_____________________________')
        self.seq_record.id = self.seq_record.id.translate(transtab)
        self.target_species = target_species
        self.start_rb = start_rb
        self.stop_rb = stop_rb

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

    @staticmethod
    def _translate_annotations_(hit_list, translate_annotation_params=None):
        """Translates list of hit ids to another standard based on parameters.

        :param hit_list:
        :param translate_annotation_params:
        :return:
        """

        if translate_annotation_params:
            if translate_annotation_params["method"].lower() == "mygene":
                return translate_ids(hit_list,
                                     orig=translate_annotation_params["anno_type_old"],
                                     to=translate_annotation_params["anno_type_new"],
                                     species=translate_annotation_params["species"])
            elif translate_annotation_params["method"].lower() == "table":
                trans_dict = translate_annotation_params["trans_dict"]
                return [trans_dict[rec] if rec in trans_dict else rec for rec in hit_list]
        else:
            return hit_list

    def _search_(self, record, target_species, search_type, settings, criteria, outpath, verbose, indent):
        search = Search(search_type=search_type)
        search_params = {}
        search_params.update(settings)
        search_params.update(criteria)
        if self.start_rb and "search" in self.start_rb[0].lower():
            searchrecord, search_err = (Search.load(Path(self.start_rb[1],
                                                         "{0}_{1}.{2}"
                                                         "".format(self.target_species.replace(' ', "_"),
                                                                   self.seq_record.name,
                                                                   self.start_rb[2]))),
                                        None)
        else:
            searchrecord, search_err = search(seq_record=record,
                                              species=target_species,
                                              search_type=search_type,
                                              blastoutput_custom=outpath,
                                              verbose=verbose, indent=indent + 1,
                                              **search_params)
        if search_err:
            raise SearchError(search_err)
        else:
            if searchrecord:
                return searchrecord
            else:
                raise SearchError("No Search Record returned for query {}!".format(self.seq_record.name))

    def _id_ranker_(self, search_record, verbose, return_only, search_criteria, indent):
        if self.start_rb and self.start_rb[0].lower() in ['forward_id_ranker']:
            id_ranked_path = Path(self.start_rb[1], "{0}_{1}.{2}".format(self.target_species.replace(' ', "_"),
                                                                         self.seq_record.name, self.start_rb[2]))
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
                        raise ValueError(message)
        else:
            # Note: output is list of (seq_chr, seq_range, seq_name, seq_score, strand, thickStart, thickEnd, rgb,
            #  blockcount, blockspans, blockstarts, query_span, seq_coverage)
            f_id_ranked = id_ranker(search_record, verbose=verbose,
                                    indent=indent + 1,
                                    return_only=return_only,
                                    **search_criteria)
        if f_id_ranked:
            return f_id_ranked
        else:
            message = "No hits meeting the forward search criteria were found for query {}!".format(self.seq_record.id)
            raise NoHitsError(message)

    def _get_sequence_(self, id_ranked, sequence_source_type, sequence_source_settings, output_type, verbose, indent):
        if isinstance(sequence_source_settings['database'], dict):
            try:
                sequence_source_settings['database'] = sequence_source_settings['database'][self.target_species]
            except KeyError:
                raise DatabaseNotFoundError('No sequence source database for species {} '
                                            'was found in the provided dict!'.format(self.target_species))
        if self.start_rb and self.start_rb[0].lower() in ['fetchseq']:
            seq_path = Path(self.start_rb[1], "{0}_{1}.{2}".format(self.target_species.replace(' ', "_"),
                                                                   self.seq_record.name, self.start_rb[2]))
            if 'fa' in self.start_rb[2]:
                seq_type = 'fasta'
            elif self.start_rb[2] in ['gb', 'genbank']:
                seq_type = 'gb'
            else:
                seq_type = self.start_rb[2]
            seq_list = list(SeqIO.parse(str(seq_path), seq_type))
            missing_items = []
        else:
            # Note: BioSQL is NOT thread-safe as implemented in fetchseq!
            # Throws tons of errors if executed with more than one thread!
            seq_list, missing_items = fetchseq(ids=id_ranked, species=self.target_species, delim='\t',
                                               source=sequence_source_type,
                                               output_type=output_type, output_name='',
                                               verbose=verbose, n_threads=1, indent=indent + 1,
                                               n_subthreads=1, **sequence_source_settings)
        if missing_items:
            print('Items were missing!', indent=indent)
            for i in missing_items:
                print(i, indent=indent + 1)
        if seq_list:
            return seq_list, missing_items
        else:
            message = 'No Sequence List was returned for species "{0}", query "{1}"!'.format(self.target_species,
                                                                                             self.seq_record.id)
            raise NoHitsError(message)

    @staticmethod
    def _annotate_forward_hits_(fw_hit_record, reverse_hit_list, verbose, indent):
        reverse_search_annotations = ",".join(('{0}:{1}-{2}({3})'.format(anno[0], anno[1][0],
                                                                         anno[1][1], anno[4])
                                               for anno in reverse_hit_list))
        fw_hit_record.features[0].qualifiers["annotations"] = reverse_search_annotations
        fw_hit_record.description += '\tAnnotations=' + reverse_search_annotations
        if verbose > 3:
            print("New Hit Record Description:", indent=indent)
            print(fw_hit_record.description, indent=indent + 1)
        return fw_hit_record

    def __call__(self, proc_id, query_species, forward_search_type, forward_search_criteria,
                 forward_search_settings, sequence_source_type, sequence_source_settings, reverse_search_type,
                 reverse_search_criteria, reverse_search_settings, reciprocal_method, output_type, outfolder,
                 translate_annotation_params, verbose, memory_saver_level, indent):
        # Simple shunt to minimize having to rewrite code.
        target_species = self.target_species
        self.query_species = query_species

        # Creating the RecBlast Container
        rc_container_full = RecBlastContainer(target_species=target_species, query_record=self.seq_record,
                                              query_species=query_species)
        # Shorthand reference for full container
        rc_container = rc_container_full[target_species][self.seq_record.id]
        rc_container.update(proc_id=proc_id)
        if verbose:
            message = '[Proc: {0}] [Seq.Name: {1}] [Target: {2}]'.format(proc_id, self.seq_record.id, target_species)
            print(message, indent=indent)
        indent += 1
        if verbose > 1:
            message = 'Creating handles for intermediary outputs...'
            print(message, indent=indent)
        # Handle for output paths:
        rc_container['output_paths'] = self._set_output_paths_(query_species, target_species,
                                                               forward_search_type,
                                                               output_type)
        output_paths = rc_container['output_paths']

        # FORWARD SEARCH
        try:
            if verbose:
                message = "Performing Forward Search for Record \"{}\": ".format(self.seq_record.id)
                print(message, indent=indent)
            forward_search_record = self._search_(record=self.seq_record,
                                                  target_species=self.target_species,
                                                  search_type=forward_search_type,
                                                  settings=forward_search_settings,
                                                  criteria=forward_search_criteria,
                                                  outpath=output_paths["forward_search_output"],
                                                  verbose=verbose,
                                                  indent=indent)
            rc_container['forward_search'] = None if memory_saver_level else forward_search_record
        except RecBlastException as err:
            if isinstance(err, SearchError):
                err = ForwardSearchError(err)
            print("{0}: {1}".format(str(type(err)), err))
            return rc_container_full

        # ID_RANKER
        if verbose:
            message = "Ranking Forward Hits: "
            print(message, indent=indent)
        try:
            f_id_ranked = self._id_ranker_(search_record=forward_search_record, verbose=verbose,
                                           return_only=1 if reciprocal_method in ['1:1', '1-to-1'] else None,
                                           search_criteria=forward_search_criteria, indent=indent)
            forward_ids = rc_container['forward_ids']
            forward_ids['ids'] = None if memory_saver_level else f_id_ranked
            f_id_out_list = ['{chr}\t{start}\t{end}\t{name}\t{score}\t{strand}\n'.format(chr=id_i[0], start=id_i[1][0],
                                                                                         end=id_i[1][1], name=id_i[2],
                                                                                         score=id_i[3], strand=id_i[4])
                             for id_i in f_id_ranked]
            forward_ids['pretty_ids'] = None if memory_saver_level else f_id_out_list
        except RecBlastException as err:
            print("{0}: {1}".format(str(type(err)), err))
            return rc_container_full

        # FETCHSEQ
        try:
            if verbose:
                message = "Fetching Sequences for Forward Hits:"
                print(message, indent=indent)
            recblast_sequence, missing_items = self._get_sequence_(id_ranked=f_id_ranked,
                                                                   sequence_source_type=sequence_source_type,
                                                                   sequence_source_settings=sequence_source_settings,
                                                                   output_type=output_type,
                                                                   verbose=verbose,
                                                                   indent=indent)
            forward_ids['missing_items'] = None if memory_saver_level else missing_items
            rc_container['recblast_unanno'] = None if memory_saver_level else recblast_sequence
        except RecBlastException as err:
            print("{0}: {1}".format(str(type(err)), err))
            return rc_container_full

        # RECIPROCAL SEARCH & ANNOTATION
        if verbose:
            message = 'Performing Reverse Searches:'
            print(message, indent=indent)
        if reciprocal_method in ['1-to-1', '1:1']:
            if verbose:
                message = '1:1 Reciprocal BLAST was selected, will only annotate top hit!'
                print(message, indent=indent)
            recblast_sequence = [recblast_sequence[0]]
        indent += 1
        old_indent = indent

        for index, fw_hit_record in enumerate(recblast_sequence):
            fw_hit_record.features[0].qualifiers["annotations"] = ""
            indent = old_indent
            assert isinstance(fw_hit_record, SeqRecord), ('Warning! entry_record is of type {} '
                                                          'rather than SeqRecord!').format(str(type(fw_hit_record)))
            if verbose:
                message = "Hit #{0}: {1}".format(index, fw_hit_record.name)
                print(message, indent=indent)
            indent += 1
            if verbose > 1:
                message = ">{}".format(fw_hit_record.description)
                print(message, indent=indent)
                message = fw_hit_record.seq[0:10] + '...' + fw_hit_record.seq[-10:]
                print(message, indent=indent)
            output_paths['reverse_search_output'].append(
                Path("{0}_recblast_out".format(self.target_species).replace(' ', '_'),
                     "{0}_{1}_tmp".format(type,
                                          self.seq_record.id
                                          ).replace(' ', '_'),
                     "{0}_{1}_{3}_to_{2}_{4}.xml".format(type,
                                                         self.seq_record.id,
                                                         self.query_species,
                                                         self.target_species,
                                                         fw_hit_record.name
                                                         ).replace(' ', '_')
                     )
            )

            # REVERSE SEARCH
            try:
                if verbose:
                    message = "Performing Reverse Search for Record \"{}\": ".format(fw_hit_record.id)
                    print(message, indent=indent)
                reverse_search_record = self._search_(record=fw_hit_record,
                                                      target_species=query_species,
                                                      search_type=reverse_search_type,
                                                      settings=reverse_search_settings,
                                                      criteria=reverse_search_criteria,
                                                      outpath=output_paths['reverse_search_output'][-1],
                                                      verbose=verbose,
                                                      indent=indent + 1)
            except NoHitsError as err:

                if isinstance(err, SearchError):
                    err = ReverseSearchError(err)
                print("{0}: {1}".format(str(type(err)), err), indent=indent)
                print("As nothing was found, continuing to next hit!", indent=indent)
                continue
            except RecBlastException as err:
                if isinstance(err, SearchError):
                    err = ReverseSearchError(err)
                print("{0}: {1}".format(str(type(err)), err))
                return rc_container_full

            # ID RANKER
            if verbose:
                message = "Ranking Reverse Hits: "
                print(message, indent=indent)
            try:
                r_id_ranked = self._id_ranker_(search_record=reverse_search_record, verbose=verbose,
                                               return_only=1 if reciprocal_method in ['1:1', '1-to-1'] else None,
                                               search_criteria=reverse_search_criteria, indent=indent+1)
                reverse_ids = rc_container['reverse_ids']
                reverse_ids['ids'] = None if memory_saver_level else r_id_ranked
                r_id_out_list = [
                    '{chr}\t{start}\t{end}\t{name}\t{score}\t{strand}\n'.format(chr=id_i[0], start=id_i[1][0],
                                                                                end=id_i[1][1], name=id_i[2],
                                                                                score=id_i[3], strand=id_i[4])
                    for id_i in r_id_ranked]
                reverse_ids['pretty_ids'] = None if memory_saver_level else r_id_out_list
            except NoHitsError as err:
                if isinstance(err, SearchError):
                    err = ReverseSearchError(err)
                print("{0}: {1}".format(str(type(err)), err), indent=indent)
                print('Continuing to next Sequence!', indent=indent)
                continue
            except RecBlastException as err:
                print("{0}: {1}".format(str(type(err)), err))
                return rc_container_full
            if reciprocal_method in ['best-hit', 'best hit']:
                if verbose:
                    message = 'Best Hit Reciprocal BLAST was selected, will only use top hit for annotation!'
                    print(message, indent=indent)
                r_id_ranked = [r_id_ranked[0]]

            # ANNOTATION
            if translate_annotation_params:
                try:
                    if verbose:
                        message = "Translating annotation using method {}".format(translate_annotation_params["method"])
                        print(message, indent=indent)
                    r_ids = [i[0] for i in r_id_ranked]
                    new_anno = self._translate_annotations_(r_ids,
                                                            translate_annotation_params=translate_annotation_params)
                    for i, new_id in enumerate(new_anno):
                        r_id_ranked[i][0] = new_id
                except Exception as err:
                    if verbose:
                        print("Failed to translate annotation!",
                              indent=indent + 1)
                        print("Error:", indent=indent + 1)
                        print("{0}: {1}".format(type(err), err), indent=indent + 2)
            try:
                if verbose:
                    print('Annotating RecBlast Hits:', indent=indent)
                recblast_sequence[index] = self._annotate_forward_hits_(fw_hit_record=fw_hit_record,
                                                                        reverse_hit_list=r_id_ranked,
                                                                        verbose=verbose,
                                                                        indent=indent+1)
            except RecBlastException as err:
                print("{0}: {1}".format(str(type(err)), err))
                return rc_container_full
        rc_container['recblast_results'] = recblast_sequence
        if memory_saver_level > 1:
            rc_container_full = (target_species,
                                 self.seq_record.id,
                                 "#{}".format(target_species) + "\n".join(
                                     ["{h_chr}\t{h_start}\t{h_end}\t"
                                      "{h_query}\t{h_score}\t{h_strand}\t"
                                      "{h_anno}\t"
                                      "{h_coverage}".format(h_chr=rec.name,
                                                            h_start=rec.features[0].location.start,
                                                            h_end=rec.features[0].location.end,
                                                            h_query=self.seq_record.name,
                                                            h_score=str(rec.features[0].qualifiers['score']),
                                                            h_strand=rec.features[0].location.strand,
                                                            h_anno=rec.features[0].qualifiers["annotations"],
                                                            h_coverage=rec.features[0].qualifiers['query_coverage']
                                                            ) for rec in recblast_sequence]
                                 ))
        run_time = str(datetime.now() - self.starttime)
        rc_container['run_time'] = run_time
        print('PID', 'SeqName', sep='\t')
        print(proc_id, self.seq_record.id, sep='\t')
        print('Run time: ', run_time)
        return rc_container_full


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
            'tblat-transcript'
        )
        self.sequence_sources = (
            'twobit',
            'sql',
            'ensemble',
            'fasta'
        )
        self.memory_saver_level = 0
        self.start_servers = False
        self.translate_annotation_params = None
        self.max_processes = max((multiprocessing.cpu_count() - 2, 1))
        self.records = []
        self.search_parameters = {"Forward": {}, "Reverse": {}}
        # Other parameters
        self.target_species = target_species if isinstance(target_species, (list, tuple)) else [target_species]
        if isinstance(self.target_species, tuple):
            self.target_species = list(self.target_species)
        self.query_species = query_species
        self.verbose = verbose.lower().count('v') if isinstance(verbose, str) else verbose
        self.forward_search_type = forward_search_type.lower()
        self.reverse_search_type = reverse_search_type.lower()
        self._set_search_settings_()
        self.sequence_source = sequence_source
        self._set_sequence_source_settings_()

        # Tests:
        assert isinstance(self.target_species, list), ("Target species must either be a list of target species names, "
                                                       "or a single string indicating the target species name")
        assert all((isinstance(species, str) for species in self.target_species)), ("All target species "
                                                                                    "must be strings! target_species "
                                                                                    "instead contains items of the "
                                                                                    "following types: "
                                                                                    ", ".join((str(type(spc))
                                                                                               for spc
                                                                                               in self.target_species)))

        assert isinstance(verbose, int), ('Verbose was of type {}; must be either be an integer greater '
                                          'than or equal to zero, or a number of v\'s equal '
                                          'to the desired level of verbosity').format(type(verbose))

        assert self.forward_search_type in self.search_types, ("Forward search type {0} is invalid, currently "
                                                               "supported search types are: {1}."
                                                               ).format(self.forward_search_type,
                                                                        ', '.join(self.verbose))
        assert self.reverse_search_type in self.search_types, ("Reverse search type {0} is invalid, currently "
                                                               "supported search types are: {1}."
                                                               ).format(self.reverse_search_type,
                                                                        ', '.join(self.verbose))

    def set_translation_annotation_parameters(self, method, **kwargs):
        self.translate_annotation_params = {"method": method} if method else None
        if not self.translate_annotation_params:
            return
        elif method.lower() == "mygene":
            self.translate_annotation_params.update(
                {
                    "species": kwargs.pop("species", "human"),
                    "anno_type_old": kwargs.pop("anno_type_old", "RefSeq"),
                    "anno_type_new": kwargs.pop("anno_type_new", "symbol")
                }
            )
        elif method.lower() == "table":
            self.translate_annotation_params["key_value_order"] = kwargs.pop("key_value_order", True)
            self.translate_annotation_params["tsv_location"] = kwargs.pop("tsv_location", "")
            if self.translate_annotation_params["tsv_location"]:
                table = Path(self.translate_annotation_params["tsv_location"])
                with table.open() as tbl:
                    lines = (line.strip().split("\t") for line in tbl.readlines())
                    if self.translate_annotation_params["key_value_order"]:
                        self.translate_annotation_params["trans_dict"] = {k: v for k, v in lines}
                    else:
                        self.translate_annotation_params["trans_dict"] = {k: v for v, k in lines}
        else:
            raise NotImplementedError()

    def set_queries(self, *queries, infile_type=None):
        infile_legal_types = ['fasta', 'fa', 'genbank', 'gb']
        if infile_type:
            assert infile_type.lower() in infile_legal_types, ("infile_type must be one of the following"
                                                               "filetypes: {}").format(", ".join(infile_legal_types))
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
        self.n_processes = int(self.n_processes)

    def _set_search_settings_(self):
        self.forward_search_criteria = dict(
            perc_score=0,
            perc_ident=0,
            perc_query_span=0
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
            perc_query_span=0
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
        assert isinstance(database_port, dict), ("{0} search ports for search of type {1} must be a dictionary "
                                                 "of species-port key pairs. Provided object was of "
                                                 "type {2} with structure {3}").format(direction.capitalize(),
                                                                                       search_type,
                                                                                       type(database_port),
                                                                                       str(database_port))
        database = {target: database} if isinstance(database, str) else database
        assert isinstance(database, dict), ("{0} search databases for search of type {1} must be a dictionary "
                                            "of species-database key pairs. Provided object was of "
                                            "type {2} with structure {3}").format(direction.capitalize(),
                                                                                  search_type,
                                                                                  type(database),
                                                                                  str(database))

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

    def optimize_parameters(self, training_bed, step_size):
        warnings.warn("Only does the first step of optimization, does not optimize reverse parameters!",
                      RecBlastWarning)
        # First, read the bed file and assemble the training set
        training_set = {}
        bedfile = Path(training_bed)
        assert bedfile.exists(), "Given bedfile path does not exist!"
        assert bedfile.is_file(), "Given bedfile path was not a file! Did you provide a directory?"
        with bedfile.open() as bed:
            hits = [line.strip().split('\t') for line in bed]
            for hit in hits:
                if hit[3] in training_set:
                    training_set[hit[3]].append((hit[0], int(hit[1]), int(hit[2])))
                else:
                    training_set[hit[3]] = [(hit[0], int(hit[1]), int(hit[2]))]
        gene_records = {rec.name: rec for rec in self.records}
        for gene in training_set:
            assert gene in gene_records, ("Training set contains the gene {}, but there is "
                                          "no associated query with the same name!").format(gene)
            training_set[gene] = sorted(sorted(sorted(training_set[gene],
                                                      key=itemgetter(2)),
                                               key=itemgetter(1)),
                                        key=itemgetter(0))
        # Deep copy the current recblast for training
        test_recblast = deepcopy(self)
        test_recblast.verbose = 0
        test_recblast.records = [gene_records[gene] for gene in training_set]

        # Here's the magic: we're gonna run this in two phases: first, optimize parameters in the
        # forward direction, without performing the reciprocal search.
        # Next, we'll select the parameter values in the forward direct that maximizes the number of correct hits
        # Finally, we'll then feed in the results from the first search to the program to do the reverse,
        # where we'll optimize the reverse parameters.
        param_tests = list(product(*[[round(i * step_size, 3) for i in range(0, int(1 / step_size) + 1)]
                                     for low, high in zip((0, 0, 0, 0), (1, 1, 1, 1))]))
        total_tests = len(param_tests)
        fw_test_results = {}

        for i, params in enumerate(param_tests):
            print("{0}/{1}".format(i+1, total_tests), end="... ")
            test_recblast.forward_search_criteria.update({'perc_ident': params[0],
                                                          'perc_query_span': params[1],
                                                          'perc_score': params[2]})

            r_name = "test" + "_fwd_" + "_".join(("{0}-{1}".format(k, v) for k, v in
                                                  test_recblast.forward_search_criteria.items()))
            fw_test_results[i] = [rec[2]['forward_ids']['ids'] for rec in
                                  test_recblast(r_name, output_type=None, stop="forward_id_ranker")]
            print("Done!")
        score_list = []
        for i, fwout in fw_test_results.items():

            score = {g: None for g in training_set}
            fwout = [sorted(sorted(sorted(a, key=lambda i: i[1][1]), key=lambda i: i[1][0]), key=itemgetter(0))
                     if len(a) > 1 else a for a in fwout]

            fw_results = {x[0][2]: [(y[0], y[1][0], y[1][1]) for y in x] for x in fwout if len(x) > 0}
            for gene in training_set:
                if gene not in fw_results:
                    fw_results[gene] = []
                max_score = len(set(training_set[gene]))
                hits = len(set(training_set[gene]) & set(fw_results[gene]))
                false_neg = len(set(training_set[gene]) - set(fw_results[gene]))
                false_pos = len(set(fw_results[gene]) - set(training_set[gene]))
                score[gene] = (hits/max_score, false_neg/max_score, false_pos/max_score)
            score_list.append([sum(i)/len(score) for i in zip(*score.values())]+[i])
        # First sort by ascending false positive scores (low->high), then by ascending false negative scores,
        # then by descending hit scores.
        # This ensures that the highest score maximizes the positive hits, followed by minimizing the false negatives,
        # and finally minimizing the false positives.
        highest_score = sorted(sorted(sorted(score_list,
                                             key=itemgetter(2)),
                                      key=itemgetter(1)),
                               key=itemgetter(0), reverse=True)[0]
        best_params_ind = [i for i, j in enumerate(score_list) if j[0:3] == highest_score[0:3]]
        best_params = [param_tests[i] for i in best_params_ind]

        return [sorted(set(i)) for i in zip(*best_params)]

    def dump_paramfile(self, location):
        location = Path(location) if not isinstance(location, Path) else location
        with location.open("w") as outf:
            outf.write(str(self))

    def _makeoutputfolders(self, run_name, location):
        rootfolder = Path(location, '{0}'.format(run_name))
        logfolder = Path(rootfolder, "logs/")
        outfolder = Path(rootfolder, "output/")
        for folder in (rootfolder, logfolder, outfolder):
            try:
                folder.mkdir(parents=True)
            except FileExistsError:
                pass
        return (rootfolder, logfolder, outfolder)

    def __call__(self, run_name=None, reciprocal_method="best-hit", output_location="./", output_type="bed-min",
                 start=None, stop=None):
        date_str = datetime.now().strftime('%y-%m-%d_%I-%M-%p')
        run_name = run_name if run_name else date_str
        assert len(self.records) > 0, ("No query records have been set! "
                                       "Please use self.set_queries() to set query records.")
        if start:
            assert isinstance(start, tuple) and len(start) == 3, ("If not set to 'None', parameter 'start' must be a "
                                                                  "tuple of length 2, with the first item indicating "
                                                                  "the starting point of the program; the second item "
                                                                  "indicating the folder of the intermediary files "
                                                                  "necessary to proceed; and the third item indicating"
                                                                  "the extension of the file. Items in folder must be "
                                                                  "named in the format of {species}_{query_name} in "
                                                                  "order for the program to correctly identify them.")
        if self.verbose:
            print('RecBlast version: ', __version__)
            print('Using BioPython version: ', bp_version)
            print('Beginning RecSearchMP!')
        if self.verbose == 1:
            print('Basic verbose mode active. Will print only essential commentary.')
        elif self.verbose == 2:
            print('Verbose mode was set to 2. Will elaborate considerably about the script.')
        elif self.verbose == 3:
            print('Debugging-level verbose mode set. You will be innunadated by text. '
                  'Brace yourself, and hold on to your console.')
        elif self.verbose >= 49:
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
        if "perc_ident" in self.forward_search_criteria and not isinstance(self.forward_search_criteria['perc_ident'],
                                                                           float):
            self.forward_search_criteria['perc_ident'] = int(self.forward_search_criteria['perc_ident'] * 100)
        if "perc_ident" in self.reverse_search_criteria and not isinstance(self.reverse_search_criteria['perc_ident'],
                                                                           float):
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

        (rootfolder, logfolder, outfolder) = self._makeoutputfolders(run_name, output_location) if output_type else None

        # Dump parameters for future reference
        self.dump_paramfile(Path(logfolder,run_name+".parameters"))
        #########################################################################
        # Get Database Files if set to "auto"
        if self.forward_search_settings['database'] == "auto":
            self.forward_search_settings['database'] = {
                species: get_searchdb(search_type=self.forward_search_type,
                                      species=species,
                                      db_loc=self.forward_search_settings['database_path'],
                                      verbose=self.verbose, indent=1).name
                for species in self.target_species}
        if self.reverse_search_settings['database'] == "auto":
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
            print('Creating RecSearch Threads... ')
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
                                                outfolder=logfolder,
                                                translate_annotation_params=self.translate_annotation_params,
                                                verbose=self.verbose,
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
                rb_queue.put(RecBlastRun(seq_record=rec, target_species=species, start_rb=start, stop_rb=stop))
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
                                    of.write(recblast_out[-1][2])
                            elif output_type:
                                recblast_out[-1].write(filetype=output_type, file_loc=outfolder, verbose=self.verbose)
                    if self.memory_saver_level:
                        recblast_out[-1] = 1

                except Exception as err:
                    print('WARNING! Could not write output of RecSearch #{}'.format(len(recblast_out)))
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
                'sequence_source={4}, n_records={5}, verbose={6})'.format(self.target_species,
                                                                          self.query_species,
                                                                          self.forward_search_type,
                                                                          self.reverse_search_type,
                                                                          self.sequence_source,
                                                                          len(self.records),
                                                                          self.verbose)
                )

    def __str__(self):
        # 2021-08-15: doesn't seem like t is used, rather, only tap is
        #t = ("{0}: {1}".format(k, v) for k, v in self.translate_annotation_params.items())
        tap="\n\t".join(i if len(i) <= 200 else i[0:200] + " ... " + i[-10:-1] for i in t) if \
            isinstance(self.translate_annotation_params, dict) else "None"
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
            "\t{seq_set}\n" \
            "Annotation Translation Parameters:\n" \
            "\t{tap}\n" \
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
                      ),
                      tap=tap
                      )
        return s