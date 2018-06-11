import bz2
import hashlib
import sqlite3
from functools import partial, reduce
from copy import deepcopy
from pathlib import Path
from datetime import datetime

import dill as pickle
from pickle import UnpicklingError
from Bio import SeqIO, SeqFeature
from Bio.SeqRecord import SeqRecord
from Bio.Blast.Record import Blast as BioBlastRecord

from RecBlast import print, __version__
from RecBlast.WarningsExceptions import *


class RecBlastRecord(dict):
    def __init__(self, proc_id=None, version=None, run_time=None, query_record=None, query_species=None,
                 target_species=None, forward_search=None, forward_ids=None, recblast_unanno=None,
                 reverse_search=None, reverse_ids=None, recblast_results=None, output_paths=None):
        super(dict, self).__init__()
        self.proc_id = proc_id if proc_id else ""
        self.version = version if version else __version__
        self.run_time = run_time if run_time else 0
        self.query_record = query_record if query_record else SeqRecord(seq="")
        self.query_species = query_species if query_species else ""
        self.target_species = target_species if target_species else ""
        self.forward_search = forward_search if forward_search else dict(search_results=BioBlastRecord,
                                                                         search_errors='')
        self.forward_ids = forward_ids if forward_ids else dict(ids=list(),
                                                                missing_ids=list(),
                                                                pretty_ids=list())
        self.recblast_unanno = recblast_unanno if recblast_unanno else []
        self.reverse_search = reverse_search if reverse_search else dict(search_results=BioBlastRecord,
                                                                         search_errors='')
        self.reverse_ids = reverse_ids if reverse_ids else dict(ids=list(), missing_ids=list(), pretty_ids=list())
        self.recblast_results = recblast_results if recblast_results else [],
        self.output_paths = output_paths if output_paths else dict(forward_search_output=Path(),
                                                                   forward_id_score_output=Path(),
                                                                   recblast_output_unanno=Path(),
                                                                   reverse_search_output=list(),
                                                                   recblast_output=Path(),
                                                                   search_nohits=Path()
                                                                   )


class RecBlastContainer(dict):
    """RecBlastContainer class containing all intermediary and final RecBlast outputs."""

    def __init__(self, target_species, query_record, **kwargs):
        super(dict, self).__init__()
        if isinstance(query_record, SeqIO.SeqRecord):
            self[target_species] = {query_record.id: dict(
                proc_id=kwargs.pop('proc_id', str()),
                version=kwargs.pop('version', __version__),
                run_time=kwargs.pop('run_time', 0),
                query_record=kwargs.pop('query_record', query_record),
                target_species=kwargs.pop('target_species', target_species),
                query_species=kwargs.pop('query_species', str()),
                forward_search=kwargs.pop('forward_search',
                                          dict(search_results=BioBlastRecord,
                                               search_errors='')),
                forward_ids=kwargs.pop('forward_ids',
                                       dict(ids=list(),
                                            missing_ids=list(),
                                            pretty_ids=list())),
                recblast_unanno=kwargs.pop('recblast_unanno', list()),
                reverse_search=kwargs.pop('reverse_search',
                                          dict(search_results=BioBlastRecord,
                                               search_errors='')),
                reverse_ids=kwargs.pop('reverse_ids',
                                       dict(ids=list(),
                                            missing_ids=list(),
                                            pretty_ids=list())),
                recblast_results=kwargs.pop('recblast_results', list()),
                output_paths=kwargs.pop('output_paths',
                                        dict(forward_search_output=Path(),
                                             forward_id_score_output=Path(),
                                             recblast_output_unanno=Path(),
                                             reverse_search_output=list(),
                                             recblast_output=Path(),
                                             search_nohits=Path()
                                             ))
            )}
        elif isinstance(query_record, dict):
            self[target_species] = query_record
        else:
            raise AssertionError('Got Query Record of type {}, must either be a query dict '
                                 'from another RecBlastContainer or a SeqRecord Object!'.format(type(query_record)))

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

    def __iter__(self):
        for species in self.keys():
            for query in self[species].keys():
                yield (species, query, self[species][query])

    def species(self):
        return list(self.keys())

    def queries(self):
        q = []
        _ = [q.extend(self[i].keys()) for i in self.keys()]
        q = list(set(q))
        return q

    def result_filter(self, func, *args, species=None, query=None, **kwargs):
        replace_internal = kwargs.pop('replace_internal', False)
        if '__dict__' in self.keys():
            del self['__dict__']

        if species is None:
            if replace_internal:
                for spec in self.keys():
                    self.result_filter(func=func, *args, species=spec, query=query,
                                       replace_internal=replace_internal, **kwargs)
            else:
                rc = (self.result_filter(func=func, *args, species=spec, query=query,
                                         replace_internal=replace_internal, **kwargs) for spec in self.keys())
                return sum(rc)
        elif query is None:
            if replace_internal:
                for q in self[species].keys():
                    self.result_filter(func=func, *args, species=species, query=q, replace_internal=replace_internal,
                                       **kwargs)
            else:
                rc = (self.result_filter(func=func, *args, species=species,
                                         query=q, replace_internal=replace_internal,
                                         **kwargs) for q in self[species].keys())
                return sum(rc)
        else:
            if replace_internal:
                self[species][query]['recblast_results'] = list(filter(partial(func, *args, **kwargs),
                                                                       self[species][query]['recblast_results']))
                return self
            else:
                query_record = dict()
                query_record[query] = deepcopy(self[species][query])
                summary_statistic = kwargs.pop('summary_statistic', None)
                recblast_object = kwargs.pop('recblast_object', 'recblast_results')
                assert isinstance(recblast_object, str), ('recblast_object was of type {}, '
                                                          'must be a str!').format(type(recblast_object))
                if summary_statistic:
                    try:
                        stat = summary_statistic(self[species][query][recblast_object], *args, **kwargs)
                    except KeyError:
                        raise KeyError('Record {0}, {1} has no key {2}'.format(species, query, recblast_object))
                    query_record[query]['recblast_results'] = list(filter(partial(func, *args, stat=stat, **kwargs),
                                                                          self[species][query]['recblast_results']))
                else:
                    query_record[query]['recblast_results'] = list(filter(partial(func, *args, **kwargs),
                                                                          self[species][query]['recblast_results']))
                return RecBlastContainer(target_species=species, query_record=query_record)

    def result_map(self, func, *args, species=None, query=None, **kwargs):
        replace_internal = kwargs.pop('replace_internal', False)
        if '__dict__' in self.keys():
            del self['__dict__']

        if species is None:
            if replace_internal:
                for spec in self.keys():
                    self.result_map(func=func, *args, species=spec, query=query,
                                    replace_internal=replace_internal, **kwargs)
            else:
                rc = (self.result_map(func=func, *args, species=spec, query=query,
                                      replace_internal=replace_internal, **kwargs) for spec in self.keys())
                return sum(rc)
        elif query is None:
            if replace_internal:
                for q in self[species].keys():
                    self.result_map(func=func, *args, species=species, query=q, replace_internal=replace_internal,
                                    **kwargs)
            else:
                rc = (self.result_map(func=func, *args, species=species,
                                      query=q, replace_internal=replace_internal,
                                      **kwargs) for q in self[species].keys())
                return sum(rc)
        else:
            if replace_internal:
                self[species][query]['recblast_results'] = list(map(partial(func, *args, **kwargs),
                                                                    self[species][query]['recblast_results']))
                return self
            else:
                query_record = dict()
                query_record[query] = deepcopy(self[species][query])
                summary_statistic = kwargs.pop('summary_statistic', None)
                recblast_object = kwargs.pop('recblast_object', 'recblast_results')
                assert isinstance(recblast_object, str), 'recblast_object must be a str!'
                if summary_statistic:
                    try:
                        stat = summary_statistic(self[species][query][recblast_object], *args, **kwargs)
                    except KeyError:
                        raise KeyError('Record {0}, {1} has no key {2}'.format(species, query, recblast_object))
                    query_record[query]['recblast_results'] = list(map(partial(func, *args, stat=stat, **kwargs),
                                                                       self[species][query]['recblast_results']))
                else:
                    query_record[query]['recblast_results'] = list(map(partial(func, *args, **kwargs),
                                                                       self[species][query]['recblast_results']))
                return RecBlastContainer(target_species=species, query_record=query_record)

    def result_reduce(self, func, *args, species=None, query=None, **kwargs):
        replace_internal = kwargs.pop('replace_internal', False)
        if '__dict__' in self.keys():
            del self['__dict__']

        if species is None:
            if replace_internal:
                for spec in self.keys():
                    self.result_reduce(func=func, *args, species=spec, query=query,
                                       replace_internal=replace_internal, **kwargs)
            else:
                rc = (self.result_reduce(func=func, *args, species=spec, query=query,
                                         replace_internal=replace_internal, **kwargs) for spec in self.keys())
                return sum(rc)
        elif query is None:
            if replace_internal:
                for q in self[species].keys():
                    self.result_reduce(func=func, *args, species=species, query=q, replace_internal=replace_internal,
                                       **kwargs)
            else:
                rc = (self.result_reduce(func=func, *args, species=species,
                                         query=q, replace_internal=replace_internal,
                                         **kwargs) for q in self[species].keys())
                return sum(rc)
        else:
            if replace_internal:
                self[species][query]['recblast_results'] = list(reduce(partial(func, *args, **kwargs),
                                                                       self[species][query]['recblast_results']))
                return self
            else:
                query_record = dict()
                query_record[query] = deepcopy(self[species][query])
                summary_statistic = kwargs.pop('summary_statistic', None)
                recblast_object = kwargs.pop('recblast_object', 'recblast_results')
                assert isinstance(recblast_object, str), 'recblast_object must be a str!'
                if summary_statistic:
                    try:
                        stat = summary_statistic(self[species][query][recblast_object], *args, **kwargs)
                    except KeyError:
                        raise KeyError('Record {0}, {1} has no key {2}'.format(species, query, recblast_object))
                    query_record[query]['recblast_results'] = list(reduce(partial(func, *args, stat=stat, **kwargs),
                                                                          self[species][query]['recblast_results']))
                else:
                    query_record[query]['recblast_results'] = list(reduce(partial(func, *args, **kwargs),
                                                                          self[species][query]['recblast_results']))
                return RecBlastContainer(target_species=species, query_record=query_record)

    def write(self, file_loc=None, filetype='fasta', verbose=1, **kwargs):
        if filetype is None:
            return 0
        if file_loc is None:
            date_str = datetime.now().strftime('%y-%m-%d_%I-%M-%p')
            file_loc = Path('./RecBlast_output/{0}/'.format(date_str)).absolute()
            try:
                file_loc.mkdir(parents=True)
            except FileExistsError:
                pass
        elif isinstance(file_loc, Path):
            pass
        elif isinstance(file_loc, str):
            file_loc = Path(file_loc)
        else:
            raise TypeError('file_loc was of type {}, must be either a string or a Path object!'.format(file_loc))
        if self == dict():
            if verbose:
                print('rc_container is empty!')
            return 0
        else:
            if filetype.lower() in 'sqlite3':
                tbn = kwargs.pop('table_name', 'RecBlastOutput')
                odb = kwargs.pop('outdb', None)
                nwrite = self._write_sqlite(outdb=odb, sqlfile=file_loc, table_name=tbn, verbose=verbose, **kwargs)
            elif 'bed' in filetype.lower():
                col = kwargs.pop('col', 12)
                custom = kwargs.pop('custom', None)
                filename = kwargs.pop('filename', 'RecBlastOutput.bed').replace(' ', '_')
                filename += '' if filename.endswith('.bed') else '.bed'
                if filetype.lower() == 'bed-min':
                    nwrite = self._write_bed(file_loc=file_loc, filename=filename, col=4, custom=custom,
                                             verbose=verbose)
                elif filetype.lower() == 'bed-complete':
                    nwrite = self._write_bed(file_loc=file_loc, filename=filename, col=col, custom='complete',
                                             verbose=verbose)
                else:
                    nwrite = self._write_bed(file_loc=file_loc, filename=filename, col=col, custom=custom,
                                             verbose=verbose)
            elif filetype.lower() in 'gff3':
                nwrite = self._write_gff3(file_loc=file_loc, verbose=verbose, **kwargs)
            else:
                nwrite = self._write_files(file_loc=file_loc, filetype=filetype, verbose=verbose, **kwargs)
        return nwrite

    def _write_bed(self, file_loc, filename, custom, col, verbose):
        nwrite = 0
        for species in self.keys():
            bed = []
            recblast_output = file_loc.joinpath(species.replace(' ', '_') + '_' + str(filename))
            try:
                recblast_output.parent.mkdir(parents=True)
            except FileExistsError:
                pass
            if verbose:
                print('Output Location for bed file of {0}:\t{1}'.format(species, str(recblast_output)))
            for query, record in self[species].items():
                if not record['recblast_results']:
                    continue
                for index, hit in enumerate(record['recblast_results']):
                    try:
                        feat = hit.features[0]
                    except IndexError:
                        feat = SeqFeature.SeqFeature()
                    except AttributeError:
                        print(type(hit))
                        print(hit.__dict__)
                        print(hit)
                        raise
                    try:
                        loc = feat.location
                        try:
                            start = str(loc.start) if loc.start is not None else '0'
                        except AttributeError:
                            start = '0'
                        try:
                            end = str(loc.end) if loc.end is not None else '0'
                        except AttributeError:
                            end = '0'
                        try:
                            strand = str(loc.strand) if loc.strand is not None else '.'

                            if strand == "1":
                                strand = "+"
                            elif strand == "-1":
                                strand = "-"
                            else:
                                strand = '.'
                        except AttributeError:
                            strand = '.'
                        loc = (str(start), str(end), str(strand))
                    except AttributeError:
                        loc = ('0', '0', '.')
                    try:
                        score = str(feat.qualifiers['score'])
                    except KeyError:
                        score = '.'
                    try:
                        thick_start = str(feat.qualifiers['thickStart'])
                    except KeyError:
                        thick_start = '.'
                    try:
                        thick_end = str(feat.qualifiers['thickEnd'])
                    except KeyError:
                        thick_end = '.'
                    try:
                        item_rgb = str(feat.qualifiers['itemRGB'])
                    except KeyError:
                        item_rgb = '.'
                    try:
                        block_count = str(feat.qualifiers['blockCount'])
                    except KeyError:
                        block_count = '.'
                    try:
                        block_sizes = str(feat.qualifiers['blockSizes'])
                    except KeyError:
                        block_sizes = '.'
                    try:
                        block_starts = str(feat.qualifiers['blockStarts'])
                    except KeyError:
                        block_starts = '.'
                    extra = []
                    if custom:
                        new_cols = []
                        if custom == 'complete':
                            for val in feat.qualifiers.keys():
                                if val in ['score', 'thickStart', 'thickEnd', 'itemRGB', 'blockCount', 'blockSizes',
                                           'blockStarts']:
                                    continue
                                else:
                                    new_cols.append('{0}={1}'.format(str(val), str(feat.qualifiers[val])))
                        for val in custom:
                            if val in feat.qualifiers.keys():
                                new_cols.append('{0}={1}'.format(str(val), str(feat.qualifiers[val])))
                        extra.append(';'.join(new_cols))
                    items = [hit.name,
                             loc[0],
                             loc[1],
                             str(query + '_' + str(index)),
                             score,
                             loc[2],
                             thick_start,
                             thick_end,
                             item_rgb,
                             block_count,
                             block_sizes,
                             block_starts]
                    if extra:
                        items += extra
                        items = [str(i) for i in items]
                    else:
                        items = [str(i) for i in items][0:col]
                    line = '\t'.join(items) + '\n'
                    bed.append(line)
                nwrite += 1
            with recblast_output.open('a+') as rc_out:
                rc_out.writelines(bed)
        return nwrite

    def _write_sqlite(self, outdb, sqlfile, table_name, verbose, **kwargs):
        nwrite = 0
        max_hit_col = kwargs.pop('max_hit_col', 1)
        # col_to_make = kwargs.pop('col_to_make', 0)
        row_number = kwargs.pop('row_number', 1)
        if outdb is None:
            outdb = sqlite3.connect(str(sqlfile))
        cursor = outdb.cursor()
        # Create the table if it doesn't exist
        command = ('CREATE TABLE IF NOT EXISTS {tn} (id INT PRIMARY KEY , target_species TEXT, query_species TEXT, '
                   'query_name TEXT, query_record TEXT, hit_1 TEXT)').format(tn=table_name)
        command = command.replace('-', '_')
        cursor.execute(command)
        outdb.commit()
        species = [i for i in self.keys() if i != '__dict__']
        for spc in species:
            targets = list(self[spc].keys())
            for target in targets:
                rc = self[spc][target]
                ts = spc
                qn = target
                qs = rc['query_species']
                qr = rc['query_record'].format('fasta')
                hits = rc['recblast_results']
                if hits == list():
                    if verbose:
                        print('No RecBlast hits in species {0} for sequence {1}!'.format(spc, rc['query_record'].id))
                    continue
                else:
                    try:
                        n_hits = len(hits)
                        if n_hits == 0:
                            hits = ['NA']
                        # Check if there's more hits than columns, and if so, add more columns to the table
                        elif n_hits > max_hit_col:
                            col_to_make = n_hits - max_hit_col
                            while col_to_make > 0:
                                col_to_make -= 1
                                max_hit_col += 1
                                hc_name = 'hit_{}'.format(max_hit_col)
                                try:
                                    cursor.execute("ALTER TABLE {tn} ADD COLUMN '{cn}' TEXT".format(tn=table_name,
                                                                                                    cn=hc_name))
                                except sqlite3.OperationalError as err:
                                    if 'duplicate column name' in str(err).lower():
                                        continue
                                    else:
                                        raise
                        # Construct value string to pass to db
                        colvals = "{id}, '{ts}', '{qs}', '{qn}', '{qr}'".format(id=row_number,
                                                                                ts=ts,
                                                                                qs=qs,
                                                                                qn=qn,
                                                                                qr=qr)
                        colvals += ', ' + ', '.join(["'{}'".format(hit.format('fasta')) for hit in hits])
                        colnames = 'id, target_species, query_species, query_name, query_record'
                        if n_hits > 0:
                            colnames += ', ' + ', '.join(['hit_{}'.format(n + 1) for n in range(n_hits)])
                        else:
                            colnames += ', hit_1'
                        # Add row
                        cursor.execute("INSERT OR IGNORE INTO {tn} ({cols}) VALUES ({vals})".format(tn=table_name,
                                                                                                    cols=colnames,
                                                                                                    vals=colvals))
                        nwrite += 1
                        row_number += 1
                        outdb.commit()
                    except AttributeError as err:
                        print('WARNING! Could not write output of RecBlastContainer[{0}, {1}]'.format(ts, qn))
                        print('Error:\n', type(err), ": ", err, indent=1)
                        continue
        outdb.commit()
        outdb.close()
        return nwrite

    def _write_gff3(self, file_loc, verbose, **kwargs):
        # Example row
        # seqid\tRecBlast\tduplicated_pseudogene,supported_by_sequence_similarity\tstart\tend\t0\t+\t0\t
        nwrite = 0
        source = 'RecBlast_MP: ' + __version__

        raise NotImplementedError('RBC write_gff3 method has not yet been implemented!')

    def _write_files(self, file_loc, filetype, verbose, **kwargs):
        """

        :param file_loc:
        :param filetype:
        :return:
        """
        nwrite = 0
        species = [i for i in self.keys() if i != '__dict__']
        for spc in species:
            targets = list(self[spc].keys())
            for target in targets:
                rc_local = self[spc][target]
                recblast_sequence = rc_local['recblast_results']
                if recblast_sequence == list():
                    if verbose:
                        print('No RecBlast hits in species {0} for sequence {1}!'.format(spc,
                                                                                         rc_local['query_record'].id))
                    continue
                else:
                    recblast_output = file_loc.joinpath(rc_local['output_paths']['recblast_output'])
                    if verbose:
                        print('Output Location:\t', str(recblast_output))
                    try:
                        recblast_output.parent.mkdir(parents=True)
                    except FileExistsError:
                        pass
                with recblast_output.open('a+') as rc_out:
                    if isinstance(recblast_sequence, list):
                        if isinstance(recblast_sequence[0], SeqRecord):
                            nwrite += SeqIO.write(recblast_sequence, rc_out, filetype)
                        else:
                            nwrite += rc_out.write('\n'.join(['>{}'.format(i) for i in recblast_sequence]))
                    elif isinstance(recblast_sequence, SeqRecord):
                        nwrite += SeqIO.write(recblast_sequence, rc_out, filetype)
                    else:
                        nwrite += rc_out.write('\n'.join(['>{}'.format(i) for i in recblast_sequence]))
        return nwrite

    def __str__(self):
        super(RecBlastContainer, self).__str__()
        strobj = ''
        i = '\t'
        for species, species_dict in self.items():
            strobj += species + ':\n'
            for query, query_dict in species_dict.items():
                strobj += 1 * i + query + ":\n"
                for key, value in query_dict.items():
                    strobj += 2 * i + key + ":\n"
                    if isinstance(value, dict):
                        for subkey, subvalue in value.items():
                            strobj += 3 * i + subkey + ":\n"
                            if isinstance(subvalue, dict):
                                for subsubkey, subsubvalue in value.items():
                                    strobj += 4 * i + subsubkey + ":\n"
                                    val = str(subsubvalue).replace('\n', '\n' + 5 * i)
                                    strobj += 5 * i + val + "\n"
                            else:
                                val = str(subvalue).replace('\n', '\n' + 4 * i)
                                strobj += 4 * i + val + "\n"
                    else:
                        val = str(value).replace('\n', '\n' + 3 * i)
                        strobj += 3 * i + val + "\n"
        return strobj

    def __add__(self, other):
        assert isinstance(other, RecBlastContainer), "Cannot add a non-RecBlastContainer object to a RecBlastContainer!"
        for species in other.keys():
            if species in self.keys():
                other_sub = other[species]
                self_sub = self[species]
                for query in other_sub.keys():
                    if query in self_sub.keys():
                        RecBlastWarning('The record corresponding to Species {0} with '
                                        'Query {1} is duplicated and will be dropped'.format(species, query))
                    else:
                        self[species].update({query: other_sub[query]})
            else:
                self.update({species: other[species]})
        return self

    def __radd__(self, other):
        if isinstance(other, int):
            pass
        else:
            for species in other.keys():
                if species in self.keys():
                    other_sub = other[species]
                    self_sub = self[species]
                    for query in other_sub.keys():
                        if query in self_sub.keys():
                            RecBlastWarning('The record corresponding to Species {0} with '
                                            'Query {1} is duplicated and will be dropped'.format(species, query))
                        else:
                            self[species].update({query: other_sub[query]})
                else:
                    self.update({species: other[species]})
        return self


def RBC_load(rbc_file, sha256_checksum, compressed=True):
    def _load(_rbc_file, _compressed, _sha256_checksum):
        fmethod = bz2.BZ2File if _compressed else open
        with fmethod(_rbc_file, 'rb', buffering=0) as f:
            f.seek(-32, 2)
            csum = f.read()
            f.seek(0, 0)
            if _sha256_checksum != 'ignore' and csum != _sha256_checksum:
                raise CheckSumError('SHA256 checksums do not match!')
            while True:
                try:
                    yield pickle.load(f)
                except EOFError:
                    break
                except UnpicklingError:
                    yield f.read()

    return sum((RecBlastContainer(i[0], **i[2]) for i in _load(rbc_file,
                                                               compressed,
                                                               sha256_checksum) if all([isinstance(i[0], str),
                                                                                        isinstance(i[2], dict)]
                                                                                       )))


def RBC_dump(rbc_object, filename='RBC', compressed=True):
    sha256 = hashlib.sha256()
    fmethod = bz2.BZ2File if compressed else open
    suffix = '.pbz2' if compressed else 'p'
    with fmethod(filename + suffix, 'wb', buffering=0) as f:
        for c in rbc_object:
            s = pickle.dumps(c, protocol=pickle.HIGHEST_PROTOCOL)
            f.write(s)
            sha256.update(s)
        rbc_checksum = sha256.digest()
        f.write(rbc_checksum)
    return rbc_checksum


class RecBlastStats(object):
    def __init__(self, recblastrecord):
        self.stats = {}
        for species, query_record in recblastrecord.items():
            for query, record in query_record.items():
                self.stats[(species, query)] = {}
                self.stats[(species, query)]['query_len'] = len(record['query_record'])
        self.n_searches = sum([len(recblastrecord[species]) for species in recblastrecord.keys()])
