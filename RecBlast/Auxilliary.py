import re
import sys
from collections import namedtuple, Counter, OrderedDict
from operator import itemgetter
from math import log
from Bio import SeqIO
from RecBlast import print, merge_ranges
from RecBlast.Search import id_search
from itertools import chain, islice
import mygene
from pathlib import Path
from RecBlast.RBC import RecBlastContainer
from io import StringIO
import pandas as pd
import sqlite3


def cleanup_fasta_input(handle, filetype='fasta', write=True):

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


def massively_translate_fasta(SeqIter):

    mg = mygene.MyGeneInfo()
    all_genes = []
    def chunks(iterable, size=1000):
        iterator = iter(iterable)
        for first in iterator:
            yield list(chain([first], islice(iterator, size - 1)))

    for x in chunks(SeqIter):
        out = mg.querymany([a.id for a in x], scopes='refseq', fields='symbol', species='Homo sapiens', returnall=True)
        tdict={}
        for a in out['out']:
            try:
                tdict[a['query']] = a['symbol']
            except KeyError:
                continue
        for i in x:
            try:
                i.id = tdict[i.id]
            except KeyError:
                continue
        all_genes += x
    return all_genes


def translate_annotation(annotation, orig='refseq', to='symbol', species='human'):
    """
    Converts a name from one type to another using mygene.
    :param annotation: 
    :param orig: 
    :param to: 
    :return: 
    """
    mg = mygene.MyGeneInfo()
    out = mg.querymany(annotation, scopes=orig, fields=to, species=species)
    try:
        trans = out[0]
        return trans['symbol']
    except KeyError:
        raise Exception('No symbol found for {}'.format(annotation))


def nr_by_longest(handle, filetype='fasta', write=True):

    oldlist = SeqIO.parse(handle, filetype)
    seqdict = {}

    for seq in oldlist:
        if seq.seq == 'Sequenceunavailable':
            print('Seq Unavailable:\t', seq.name)
            continue
        try:
            seq.id, seq.description = seq.id.split('|')[0], seq.id.split('|')[1]
        except IndexError:
            seq.id, seq.description = seq.id.split(' ')[0], ''.join(seq.id.split('|')[1:len(seq.id.split('|'))])
        assert seq.id != 'gi' or seq.id != 'emb' or seq.id != 'acc'
        if seq.id in seqdict:
            if len(seq)>len(seqdict[seq.id]):
                seqdict[seq.id] = seq
            else:
                continue
        else:
            seqdict[seq.id] = seq
    newlist = (seq for _, seq in seqdict.items())
    if write:

        outhandle = 'nr_' + str(Path(handle).name)
        with Path(outhandle).open('w') as outf:
            SeqIO.write(newlist, outf, filetype)
    return newlist


def cull_reciprocal_best_hit(recblast_out):
    """
    returns a recblast_out container that only has the reciprocal best hits.
    :param recblast_out:
    :return:
    """
    pat = re.compile('\|\[(.*?)\]\|')  # regex for items in annotation
    if isinstance(recblast_out, list):
        rc_out_list = []
        for index, rc in enumerate(recblast_out):
            rc_out_list.append(cull_reciprocal_best_hit(rc))
        return rc_out_list
    else:
        #assert isinstance(recblast_out, RecBlastContainer), "Items must be RecBlastContainer Objects!"
        for species, rc_spec_rec in recblast_out.items():
            #print('Species:\t', species, indent=0)
            for query, rc_rec in rc_spec_rec.items():
                #print('Query:\t', query, indent=1)
                try:
                    rc_out = rc_rec['recblast_results']
                except KeyError:
                    print('No entries in recblast_results for query {0} in species {1}'.format(query, species))
                    continue
                tmprecord = []
                for record in rc_out:
                    try:
                        #print(record.description, indent=3)
                        target_id, annotations = record.description.split('|-|')
                        #print('Target ID:\t', target_id, indent=4)
                        #print('Annotations:', annotations.lstrip('\t'), indent=4)
                    except ValueError:
                        print(record.description, indent=2)
                        print('Could not unpack annotations!', indent=2)
                        continue
                    id_lst = pat.findall(annotations)
                    #print('id_list:\t', id_lst, indent=4)
                    if id_lst:
                        if query in id_lst[0]:
                            tmprecord.append(record)
                        else:
                            print("For query {0}, target {1} was not a reciprocal best hit!".format(query,
                                                                                                    target_id))
                            continue
                    else:
                        print('No annotations found for record {0} in species {1}, query {2}'.format(record.name,
                                                                                                     species,
                                                                                                     query))
                        continue
                recblast_out[species][query]['recblast_results'] = tmprecord
        return recblast_out


def simple_struct(recblast_out, verbose=True):
    """Returns a nice diagram of queries, targets, and annotations"""
    master_dict = {}
    pat = re.compile('\|\[(.*?)\]\|')  # regex for items in annotation
    if isinstance(recblast_out, list):
        # Prepare a list of dictionaries of length recblast_out, along with a list of respective species
        master_count = [dict] * len(recblast_out)

        for index, rc in enumerate(recblast_out):
            try:
                master_count[index] = simple_struct(rc)
            except AttributeError:
                master_count[index] = rc
        for subdict in master_count:
            for species, species_dict in subdict.items():
                if isinstance(species_dict, Exception):
                    continue
                try:
                    comb_spec_dict = master_dict[species]
                except KeyError:
                    master_dict[species] = dict()
                    comb_spec_dict = master_dict[species]
                for query, query_dict in species_dict.items():
                    try:
                        comb_query_dict = comb_spec_dict[query]
                    except KeyError:
                        comb_spec_dict[query] = dict()
                        comb_query_dict = comb_spec_dict[query]
                    for target_id, annotation_list in query_dict.items():
                        try:
                            comb_anno_list = comb_query_dict[target_id]
                        except KeyError:
                            comb_query_dict[target_id] = list()
                            comb_anno_list = comb_query_dict[target_id]
                        comb_anno_list += annotation_list if isinstance(annotation_list, list) else [annotation_list]
        return master_dict

    else:
        """
        Structure:
            master_dict:
                Species|    species_dict:
                                Query|  query_dict:
                                            target_id|  annotations_list
        """
        #assert isinstance(recblast_out, RecBlastContainer), 'Item in recblast_out was not a RecBlastContainer object!'
        try:
            recblast_out.__delitem__('__dict__')
        except KeyError:
            pass
        for species, rc_spec_rec in recblast_out.items():
            #print('Species:\t', species, indent=0)
            try:
                species_dict = master_dict[species]
            except KeyError:
                master_dict[species] = dict()
                species_dict = master_dict[species]
            for query, rc_rec in rc_spec_rec.items():
                #print('Query:\t', query, indent=1)
                try:
                    query_dict = species_dict[query]
                except KeyError:
                    species_dict[query] = dict()
                    query_dict = species_dict[query]
                try:
                    rc_out = rc_rec['recblast_results']
                except KeyError:
                    print('No entries in recblast_results for query {0} in species {1}'.format(query, species))
                    continue
                for record in rc_out:
                    try:
                        #print(record.description, indent=3)
                        target_id, annotations = record.description.split('|-|')
                        #print('Target ID:\t', target_id, indent=4)
                        #print('Annotations:', annotations.lstrip('\t'), indent=4)
                    except ValueError:
                        print(record.description, indent=2)
                        #print('Could not unpack annotations!', indent=2)
                        continue
                    try:
                        target_list = query_dict[target_id]
                    except KeyError:
                        query_dict[target_id] = list()
                        target_list = query_dict[target_id]
                    id_lst = pat.findall(annotations)
                    #print('id_list:\t', id_lst, indent=4)
                    if id_lst:
                        target_list += id_lst
                    else:
                        print('No annotations found for record {0} in species {1}, query {2}'.format(record.name,
                                                                                                     species,
                                                                                                     query))
        if verbose:
            print('*******************************************')
            for species, species_dict in master_dict.items():
                print(species, indent=0)
                for query, query_dict in species_dict.items():
                    print(query, indent=1)
                    for target_id, annotation_list in query_dict.items():
                        print(target_id, indent=2)
                        tmp = []
                        for annotation in annotation_list:
                            p, item, seq_range, id_type = id_search(annotation, id_type='brute', verbose=0)
                            if id_type == 'symbol':
                                tmp.append(item)
                            else:
                                tmp.append(item)
                        query_dict[target_id] = tmp
                        for annotation in query_dict[target_id]:
                            print(annotation, indent=3)
            print('*******************************************')
        return master_dict


def rc_out_stats(rc_out):
    # Todo: use 'from Collections import Counter' to rapidly count duplicates
    if isinstance(rc_out, list):
        holder =[]
        for rc in rc_out:
            holder.append(rc_out_stats(rc))
        c_hit_list, c_multihit_list = zip(holder)
        hit_perc = sum(c_hit_list)/len(c_hit_list)
        multihit_perc = sum(c_multihit_list)/len(c_multihit_list)

    # Percentage of searches with reciprocal hits, regardless of number:
    # Percentage of searches with more than one hit:
    elif isinstance(rc_out, RecBlastContainer):
        c_hit = 0
        c_multihit = 0
        for species, queries_dict in rc_out.items():
            for query, results in rc_out.items():
                try:
                    record_list = results['recblast_results']
                except KeyError:
                    return (0, 0)
                has_run = 0
                for record in record_list:
                    if not has_run:
                        c_hit += 1
                        has_run = 0
                    c_multihit += 1

    else:
        return None



def count_dups(recblast_out):
    """ Inverts target-annotation dictionary to find out, for every best-hit annotation, how many targets there are"""
    species_anno_target_dict = {}
    species_anno_count_dict = {}
    master_dict = simple_struct(recblast_out, verbose=False)

    for species, species_dict in master_dict.items():
        try:
            anno_target_dict = species_anno_target_dict[species]
        except KeyError:
            species_anno_target_dict[species] = {}
            anno_target_dict = species_anno_target_dict[species]
        print(species_dict, indent=0)
        for query, query_dict in species_dict.items():
            # ignoring query
            print(query_dict, indent=1)
            for target_id, annotation_list in query_dict.items():
                print(annotation_list, indent=2)
                tophit = annotation_list[0]
                print(tophit, indent=2)
                try:
                    anno_target_dict[tophit] += [target_id]
                except KeyError:
                    anno_target_dict[tophit] = list()
                    anno_target_dict[tophit].append(target_id)
                print(anno_target_dict[tophit], indent=3)
    for species, anno_dict in species_anno_target_dict.items():
        print(species, indent=0)
        try:
            anno_count_dict = species_anno_count_dict[species]
        except KeyError:
            species_anno_count_dict[species] = {}
            anno_count_dict = species_anno_count_dict[species]
        for annotation, target_list in anno_dict.items():
            print(annotation, '\t\t\t', len(target_list))
            anno_count_dict[annotation] = len(target_list)
    return species_anno_target_dict, species_anno_count_dict


class FilterRBHs(object):
    def __init__(self, **kwargs):
        """Convenience class for use with RecBlastContainer.result_filter(). Removes non-Reciprocal Best Hits from RBC.
        """
        self._recblast_object = 'query_record'
        self.args = {'func': self.fun, 'summary_statistic': self._stat, 'recblast_object': self._recblast_object}
        for k,v in kwargs.items():
            self.args[k] = v

    def _stat(self, query_record):
        """ Summary Statistic function for filter_RBH. Requires setting recblast_object='query_record'

        :param query_record:
        :return:
        """
        return query_record.name
    def fun(self, hit, stat, verbose=False):


        pat = re.compile('\|\[(.*?):.*\]\|')  # regex for items in annotation
        try:
            hit_split = hit.description.split('|-|')
            top_anno = hit_split[1]
        except ValueError:
            print(hit.description, indent=2)
            print('Could not unpack annotations!', indent=2)
            return False
        except IndexError:
            print(hit.description, indent=2)
            print('Could not unpack annotations!', indent=2)
            return False
        id_lst = pat.findall(top_anno)[0].strip()
        if id_lst:
            _, hit_symbol, _, _ = id_search(id_lst, id_type='symbol', verbose=verbose)

            if stat == hit_symbol:
                return True
        else:
            return False


def map_ranges(hit):
    """ Convenience function for RBC.results_map(). Replaces results with a tup of result descriptions and loci."""
    _, h_id, h_range, _ = id_search(hit.description, verbose=False)
    h_start = h_range[0]
    h_end = h_range[1]
    h_strand = h_range[2]
    h_d = (hit.description, h_id, h_start, h_end, h_strand)
    return h_d


def RBC_drop_many_to_one_hits(RBC):
    loci_dict_RBC = {}
    for species, query, rec in RBC.result_map(map_ranges):
        r = rec['recblast_results']
        for index, hit in enumerate(r):
            loci_dict_RBC[(hit[1], hit[2], hit[3], ''.join((query, str(index))))] = (species, query, index)
    filtered_loci_dict_RBC = drop_overlaps_bed(loci_dict_RBC)
    filter_dict = {}
    for hit_loc in filtered_loci_dict_RBC.values():
        species, query, index = hit_loc
        if (species, query) in filter_dict.keys():
            filter_dict[(species, query)].append(index)
        else:
            filter_dict[(species, query)] = [index]
    for (species, query), indexes  in filter_dict.items():
        for hit_index, hit in enumerate(RBC[species][query]['recblast_results']):
            if hit_index in indexes:
                continue
            else:
                del RBC[species][query]['recblast_results'][hit_index]


def count_reciprocal_best_hits(recblast_out):
    pat = re.compile('\|\[(.*?)\]\|')  # regex for items in annotation
    species_counters = {}
    for species, species_dict in recblast_out.items():
        species_counters[species] = Counter()
        for query, query_dict in species_dict.items():
            try:
                rc_out = query_dict['recblast_results']
            except KeyError:
                print('No entries in recblast_results for query {0} in species {1}'.format(query, species))
                continue
            for hit in rc_out:
                try:
                    hit_split = hit.description.split('|-|')
                    target_id = hit_split[0]
                    annotations = hit_split[1]
                except ValueError:
                    print(hit.description, indent=2)
                    print('Could not unpack annotations!', indent=2)
                    continue
                except IndexError:
                    print(hit.description, indent=2)
                    print('Could not unpack annotations!', indent=2)
                    continue
                id_lst = ''.join(pat.findall(annotations))
                if id_lst:
                    _, hit_symbol, _, _ = id_search(id_lst, id_type='symbol', verbose=0)

                else:
                    print('No annotations found for record {0} in species {1}, query {2}'.format(hit.name,
                                                                                                 species,
                                                                                                 query))
                    continue

                if query == hit_symbol:
                    species_counters[species].update({query:1})
    return species_counters


def export_count_as_csv(rec_hit_counter_dict, filename='RecBlastCount'):
    # First get a list of all the genes, period.
    allgenes = []
    for species, species_counter in rec_hit_counter_dict.items():
        for key, value in species_counter.items():
            if key in allgenes:
                continue
            else:
                allgenes.append(key)
    # Next, make a dict with a tuple of counts per species

    genedict = {gene: tuple((rec_hit_counter_dict[species][gene]
                             for species in rec_hit_counter_dict.keys())) for
                gene in allgenes}
    all_lines = ['Gene\t'+'\t'.join([species for species in rec_hit_counter_dict.keys()])+'\n']
    all_lines += ['{Gene}\t{counts_str}\n'.format(Gene=key,counts_str='\t'.join([str(i) for i in value]))
                  for key, value in genedict.items()]
    with open(filename+'.tsv', 'w') as outf:
        outf.writelines(all_lines)

def count_reciprocal_best_hits_from_pandas(pandas_df):

    pat = re.compile('\|\[(.*?)\]\|')  # regex for items in annotation
    spec_list = list(pandas_df.target_species.unique())
    species_counters = {}
    for species in spec_list:
        species_counters[species] = Counter()
        species_results = pandas_df.loc[pandas_df['target_species']==species]
        query_list = list(species_results.query_name.unique())
        for query in query_list:
            print(query)
            query_results = species_results.loc[species_results['query_name']==query].ix[:,5:-1]
            rc_out=[]
            for i, d in query_results.iterrows():
                rc_out += d.tolist()
            # Annoying shunt
            rc_out_asfasta = '\n'.join(['>'+i for i in rc_out if i is not None])
            tmp = StringIO(rc_out_asfasta)
            rc_out = SeqIO.parse(tmp, 'fasta')
            for hit in rc_out:
                try:
                    hit_split = hit.description.split('|-|')
                    id_lst = ''.join(pat.findall(hit_split[1]))
                except ValueError:
                    print(hit.description, indent=2)
                    print('Could not unpack annotations!', indent=2)
                    continue
                if id_lst:
                    _, hit_symbol, _, _ = id_search(id_lst, id_type='symbol', verbose=0)

                else:
                    print('No annotations found for record {0} in species {1}, query {2}'.format(hit.name,
                                                                                                 species,
                                                                                                 query))
                    continue
                if query == hit_symbol:
                    species_counters[species].update({query:1})
    return species_counters


def sqlite_to_pandas(sql_file, table_name):


    conn = sqlite3.connect(sql_file)
    df = pd.read_sql_query("select * from {0};".format(table_name), conn)
    return df


def filter_hits_pandas(pandas_df):
    def filter_func(row):
        qrec = row.query_record
        qrec = SeqIO.read(StringIO(qrec), 'fasta')
        min_len = 0.25 * len(qrec)
        intro = row.iloc[0:6].tolist()
        hits = row.iloc[5:-1].tolist()
        new_hits = []
        for hit in hits:

            if hit == 'NA':
                new_hits.append(None)
                continue
            elif hit is not None:
                tmp = '>' + hit
            else:
                new_hits.append(None)
                continue
            hit = SeqIO.read(StringIO(tmp), 'fasta')
            id_lst = hit.id
            _, hit_symbol, seq_range, _ = id_search(id_lst, id_type='brute', verbose=0)
            try:
                seq_range = seq_range[hit_symbol]
            except KeyError:
                new_hits.append(None)
                continue
            seq_len = abs(int(seq_range[1]) - int(seq_range[0]))
            new_hits.append(hit.description if seq_len >= min_len else None)
        full = intro + new_hits
        return full
    return pandas_df.apply(filter_func, axis=1)


class DataIntegratorParser(object):
    def __init__(self, file):
        transtab = str.maketrans('!@#$%^&*();:.,\'\"/\\?<>|[]{}-=+', '_____________________________')
        if isinstance(file, str):
            self.file = Path(file)
            assert self.file.exists(), file + ' is an invalid file path or does not exist!'
            assert self.file.is_file(), file + ' is not a valid file!'
        elif isinstance(file, Path):
            assert self.file.exists(), str(file) + ' is an invalid file path or does not exist!'
            assert self.file.is_file(), str(file) + ' is not a valid file!'
        else:
            raise TypeError('File must be either a str or Path object!')
        self.regions = []
        with self.file.open() as f:
            for index, line in enumerate(f):
                line = line.strip()
                if index == 0:
                    self.header = line.lstrip('# ')
                    continue
                elif line.startswith('# region='):
                    region = line.lstrip('# region=').translate(transtab)
                    if getattr(self, region, None) is None:
                        self.regions.append(region)
                        setattr(self, region, [])
                    continue
                elif line.startswith('#') and not line.startswith('# '):
                    cnames = line.lstrip('#').translate(transtab)
                    ColNames = namedtuple('ColNames', cnames.split('\t'))
                    self.colnames = ColNames._fields
                    continue
                elif line.startswith ('# No data'):
                    newitem = getattr(self, region, []) + [ColNames(*[None]*len(self.colnames))]
                    setattr(self, region, newitem)
                    continue
                else:
                    try:
                        newitem = getattr(self, region, []) + [ColNames(*line.split('\t'))]
                        setattr(self, region, newitem)
                    except NameError as err:
                        raise NameError(str(err) + '\nParser encountered a line of data before either the column names '
                                                   'or the genomic region was declared in the file!')
                    except TypeError:
                        print(line, file=sys.stderr)
                        raise
                    continue

    def rename_regions_via_bedfile(self, bedfile):
        transtab = str.maketrans('!@#$%^&*();:.,\'\"/\\?<>|[]{}-=+', '_____________________________')
        if isinstance(bedfile, str):
            self.bedfile = Path(bedfile)
            assert self.bedfile.exists(), bedfile + ' is an invalid file path or does not exist!'
            assert self.bedfile.is_file(), bedfile + ' is not a valid file!'
        elif isinstance(bedfile, Path):
            assert self.bedfile.exists(), str(bedfile) + ' is an invalid file path or does not exist!'
            assert self.bedfile.is_file(), str(bedfile) + ' is not a valid file!'
        else:
            raise TypeError('File must be either a str or Path object!')
        bed_trans = {}
        with self.bedfile.open() as f:
            for line in f:
                line = line.strip().split('\t')
                bed_trans['{0}_{1}_{2}'.format(line[0], str(int(line[1])+1), line[2])] = line[3].translate(transtab)
        self.regions = []
        for oldid in bed_trans:
            self.regions.append(bed_trans[oldid])
            setattr(self, bed_trans[oldid], getattr(self, oldid, []))
            delattr(self, oldid)

    def count_stats_per_record(self, attr_name):
        counts = OrderedDict()
        for region in sorted(self.regions):
            rec = getattr(self, region)
            c = Counter([getattr(r, attr_name) for r in rec])
            counts[region] = c
        return counts

    def __iter__(self):
        for region in self.regions:
            yield getattr(self, region)

    def __str__(self):
        string = ''
        for region in self.regions:
            content = getattr(self, region)
            string += "{0}:\t {1} ... {2} ({3})\n".format(region,
                                                          content[0][0],
                                                          content[-1][0],
                                                          len(content))
        return string


def read_bed(bedfile, key_col = 3):
    """Returns a dict using the given 0-indexed key_column"""
    d = {}
    bedfile = Path(bedfile)
    assert bedfile.exists(), "Given bedfile path does not exist!"
    assert bedfile.is_file(), "Given bedfile path was not a file! Did you provide a directory?"
    with bedfile.open() as bed:
        for line in bed:
            items = line.strip().split('\t')
            for i, j in enumerate(items):
                try:
                    new_j = int(j)
                    items[i] = new_j
                except ValueError:
                    try:
                        new_j = float(j)
                        items[i] = new_j
                    except ValueError:
                        continue
            if isinstance(key_col, slice):
                key = tuple(items[key_col])
                if key in d.keys():
                    raise KeyError('Duplicate keys in dictionary!')
                else:
                    d[key] = items
            else:
                if items[key_col] in d.keys():
                    raise KeyError('Duplicate keys in dictionary!')
                else:
                    d[items[key_col]] = items
    return d

def drop_overlaps_bed(bedfile):
    d = bedfile if isinstance(bedfile, dict) else read_bed(bedfile, key_col=slice(0,3))
    d_new = []
    dlocs = {}
    for loc in d.keys():
        if loc[0] in dlocs.keys():
            dlocs[loc[0]].append([int(loc[1]), int(loc[2]), loc[3]])
        else:
            dlocs[loc[0]] = [[int(loc[1]), int(loc[2]), loc[3]]]
    for k, v in dlocs.items():
        if len(v) > 1:
            v = [sorted(i[0:2])+[i[2]] for i in v]
            # comparison matrix
            t = [[max(v[i][0], j[0]) <= min(v[i][1], j[1]) for j in v] for i in range(0, len(v))]
            # set diagonal identities to False
            for index in range(0,len(t)):
                t[index][index] = False
            # sum per column of matrix
            t_sums = [sum(i) for i in zip(*t)]
            # Select only items which have a zero in the t_sums index
            filtered_v = [v[i] for i in range(0,len(t_sums)) if t_sums[i] == 0]
            d_new += [(k, i[0], i[1], i[2]) for i in filtered_v]
        else:
            try:
                v = v[0]
                d_new.append((k, v[0], v[1], v[2]))
            except Exception:
                print(k, v)
                raise
    filtered_d = {}
    for item in d_new:
        if item in d.keys():
            filtered_d[item] = d[item]
        elif (item[0], item[2], item[1]) in d.keys():
            filtered_d[(item[0], item[2], item[1])] = d[(item[0], item[2], item[1])]
        else:
            print(item)
            raise Exception
    return filtered_d


def calc_effective_copy_number_by_coverage(query_record):
    # get list of ranges
    if len(query_record['recblast_results']) == 0:
        return None
    else:
        raw_ranges = (hit.features[0].qualifiers['query_coverage'] for hit in query_record['recblast_results'])
        ranges = []
        for r in raw_ranges:
            try:
                rng = (int(r[0]), int(r[1]))
                ranges.append(sorted(rng))
            except IndexError:
                continue
        coverage = list(merge_ranges(ranges))
        sum_coverage = sum([i[1] - i[0] for i in coverage])
        if sum_coverage == 0:
            return 0
        else:
            sum_nuc = sum([sum([sum([s in range(r[0], r[1]) for s in range(i[0], i[1])]) for i in ranges]) for r in coverage])
            return round(sum_nuc/sum_coverage, 2)


def bed_get_flanking_regions(bedfile, left_range, right_range, genome_file=None):
    """Returns two new bedfiles with ranges left-and-right of each item of the original file, respectively.

    :param str bedfile:
    :param left_range: Either a single positive integer indicating the left-most number of bases in range;
                        or a tuple of two integers indicating the left-and-right bound of the range.
    :param right_range: Either a single positive integer indicating the right-most number of bases in range;
                        or a tuple of two integers indicating the right-and-left bound of the range.
    :return:
    """
    if isinstance(left_range, int):
        left_range = (left_range, 0)
    if isinstance(right_range, int):
        right_range = (0, right_range)

    assert isinstance(left_range, tuple), "Parameter 'left_range' must either be an integer or a tuple!"
    assert len(left_range) == 2, "Parameter 'left_range' must be a tuple of length 2!"
    assert left_range[0] > left_range[1] or left_range == (0,0), "The left-side range modifier of left_range must be " \
                                                                 "less than the right-side!"
    assert isinstance(right_range, tuple), "Parameter 'right_range' must either be an integer or a tuple!"
    assert len(right_range) == 2, "Parameter 'right_range' must be a tuple of length 2!"
    assert right_range[0] < right_range[1] or right_range == (0,0), "The right-side range modifier of left_range must" \
                                                                    " be greater than the left-side!"
    bedfile = Path(bedfile)
    assert bedfile.exists(), "Given bedfile path does not exist!"
    assert bedfile.is_file(), "Given bedfile path was not a file! Did you provide a directory?"
    leftbed = bedfile.with_name(bedfile.stem +
                                "_left_Offset{0}_Size{1}".format(left_range[1], 
                                                                 left_range[0] - left_range[1]) + 
                                bedfile.suffix)
    rightbed = bedfile.with_name(bedfile.stem +
                                "_right_Offset{0}_Size{1}".format(right_range[1], 
                                                                 right_range[0] - right_range[1]) + 
                                bedfile.suffix)
    granges = {chrm:int(size) for chrm, size
               in [line.strip().split("\t") for line in open(genome_file)]} if genome_file else None
    with bedfile.open() as bf, leftbed.open("w") as lbf, rightbed.open("w") as rbf:
        records = (line.strip().split('\t')[0:4] for line in bf)
        for (chr, s, e, id) in records:
            if left_range != (0,0):
                left = [chr,
                        int(s) - left_range[0],
                        int(s) - left_range[1],
                        id + "_left"]
                ldiff = 0
                if left[2]>left[1]>0:
                    left[3] += "_offset-{0}_size-{1}".format(left_range[1],
                                                             left[2] - left[1])
                else:
                    if left[1]<0:
                        ldiff = -left[1]  # note its '-' because left[1] is negative
                        left[2] += ldiff
                        left[2] = left[2] if left[2] <= int(s) else int(s)
                        left[1] = 0
                        if left[1] == left[2]:
                            left[2] += 1
                            ldiff -= 1
                        left[3] += "_offset-{0}_size-{1}".format(left_range[1] - ldiff,
                                                                 left[2]-left[1])
                    else:
                        left[3] += "_offset-{0}_size-{1}".format(left_range[1],
                                                                 left[2] - left[1])
                left = (str(i) for i in left)
                lbf.write('\t'.join(left) + "\n")
            if right_range != (0,0):
                right = [chr,
                         int(e) + right_range[0],
                         int(e) + right_range[1],
                         id + "_right"]
                if granges:
                    if granges[chr] <= right[2] or granges[chr] <= right[1]:
                        rdiff = granges[chr] - right[2]
                        right[2] = granges[chr]
                        right[1] += rdiff
                        right[1] = right[1] if right[1]>=int(e) else int(e)
                        if right[2] == right[1]:
                            right[1] -= 1
                            rdiff -= 1
                        right[3] += "_offset-{0}_size-{1}".format(right_range[0] + rdiff,
                                                                  right[2] - right[1])
                    else:
                        right[3] += "_offset-{0}_size-{1}".format(right_range[0],
                                                                  right[2] - right[1])
                else:
                    right[3] += "_offset-{0}_size-{1}".format(right_range[0],
                                                              right[2] - right[1])
                right = (str(i) for i in right)
                rbf.write('\t'.join(right) + "\n")
    return


def bed_rename_list_post_merge(bedfile):
    pass



def bed_extract_duplicates(bedfile, outfile="", verbose = False):
    bedfile = Path(bedfile)
    assert bedfile.exists(), "Given bedfile path does not exist!"
    assert bedfile.is_file(), "Given bedfile path was not a file! Did you provide a directory?"
    bed_dict = read_bed(bedfile)
    hits = sorted(bed_dict.keys())
    counts = Counter((''.join(hit.split("_")[:-1]) for hit in hits))
    duphits = (hit for hit in hits if counts[hit.split("_")[0]] > 1)
    outfile = Path(outfile) if outfile else bedfile.with_suffix(".bed.dups")
    try:
        first = next(duphits)
        if verbose:
            print(first, "\t", counts[first.split("_")[0]])
        with outfile.open("w") as of:
            of.write("\t".join((str(i) for i in bed_dict[first])) + "\n")
            for hit in duphits:
                if verbose:
                    print(hit, "\t", counts[hit.split("_")[0]])
                of.write("\t".join((str(i) for i in bed_dict[hit]))+"\n")
    except StopIteration:
        if verbose:
            print("No duplicates found in file!")

def merge_ids(fasta):
    outfasta = Path(fasta)
    with outfasta.with_name(outfasta.name+"_joined").open('w') as outfile:
        from Bio import SeqIO
        bla = SeqIO.parse(fasta, "fasta")
        newrec = {}
        for rec in bla:
            rec.id = rec.id.split("_left")[0].split("_right")[0]
            if rec.id in newrec:
                newrec[rec.id].seq += rec.seq
                newrec[rec.id].description += "\t" + rec.description
            else:
                newrec[rec.id] = rec
        SeqIO.write((v for v in newrec.values()), outfile, "fasta")


class BLASTSearchParameters(object):
    def __init__(self, blast_type, blastdb_path, blast_db="auto", expect=10, perc_score=0.009, perc_span=0.1,
                 ncbi_search=False, perc_ident=0.69, perc_length=0.001, megablast=True, blastdb_version='auto',
                 email = '', **kwargs):
        self.search_type = blast_type
        self.search_local = not ncbi_search
        self.email = email
        self.expect = expect
        self.perc_score = perc_score
        self.perc_ident = perc_ident
        self.perc_span = perc_span
        self.perc_length = perc_length
        self.megablast = megablast
        self.id_db_version = blastdb_version
        self.id_db_path = blastdb_path
        self.search_db = blast_db if isinstance(blast_db, dict) or isinstance(blast_db, str) else "auto"
        for k, v in kwargs:
            setattr(self, k, v)
        if ncbi_search:
            assert "@" in self.email, "If using NCBI for remote BLAST searching, a valid email must be set!"

class BLATSearchParameters(object):
    def __init__(self, blat_type, twobit_path, twobit_port_dict, gfserver_host="localhost",
                 expect=10, perc_score=0.009, perc_span=0.1, perc_ident=0.69,
                 perc_length=0.001, twobit_file_dict="auto", twobit_version='auto'):
        self.search_type = blat_type
        self.expect = expect
        self.perc_score = perc_score
        self.perc_ident = perc_ident
        self.perc_span = perc_span
        self.perc_length = perc_length
        self.search_local = gfserver_host
        self.id_db_version = twobit_version
        self.id_db_path = twobit_path
        self.id_db = twobit_file_dict if (isinstance(twobit_file_dict, dict) or
                                          isinstance(twobit_file_dict, str)) else "auto"
        self.search_db = twobit_port_dict
        self.id_source = "twobit"


class SQLServerParameters(object):
    def __init__(self, host='localhost', id_db='bioseqdb', user='postgres', driver='psycopg2',
                 password='', id_db_version='auto'):
        self.id_source = 'sql'
        self.driver = driver
        self.host = host
        self.id_db = id_db
        self.user = user
        self.password = password
        self.id_db = id_db
        self.id_db_version = id_db_version
