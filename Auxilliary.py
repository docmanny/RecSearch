import re
import sys
from builtins import print as _print



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


def print(*objects, indent=0, markup='', **print_kwargs):
    _print('\t'*indent, markup, *objects, **print_kwargs)


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


def nr_by_longest(handle, filetype='fasta', write=True):
    from Bio import SeqIO
    oldlist = SeqIO.parse(handle, filetype)
    seqdict = {}
    for seq in oldlist:
        if seq.seq == 'Sequenceunavailable':
            print('Seq Unavailable:\t', seq.name)
            continue
        seq.id, seq.description = seq.id.split('|')[0], seq.id.split('|')[1]
        assert seq.id != 'gi' or seq.id != 'emb' or seq.id != 'acc'
        if seq.id in seqdict:
            if len(seq)>len(seqdict[seq.id]):
                seqdict[seq.id] = seq
            else:
                continue
        else:
            seqdict[seq.id] = seq
    newlist = []
    for key, seq in seqdict.items():
        newlist.append(seq)
    if write:
        from pathlib import Path
        outhandle = 'nr_' + str(Path(handle).name)
        with Path(outhandle).open('w') as outf:
            SeqIO.write(newlist, outf, filetype)
    return newlist


def simple_struct(recblast_out, verbose=True):
    """Returns a nice diagram of queries, targets, and annotations"""

    from recblast_MP import id_search, RecBlastContainer
    master_dict = {}
    pat = re.compile('\|\[(.*?)\]\|')  # regex for items in annotation
    if isinstance(recblast_out, list):
        # Prepare a list of dictionaries of length recblast_out, along with a list of respective species
        master_count = [dict] * len(recblast_out)

        for index, rc in enumerate(recblast_out):
            master_count[index] = simple_struct(rc)
        for subdict in master_count:
            for species, species_dict in subdict.items():
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
                        comb_anno_list += annotation_list
        return master_dict

    else:
        """
        Structure:
            master_dict:
                Species|    species_dict:
                                Query|  query_dict:
                                            target_id|  annotations_list
        """
        assert isinstance(recblast_out, RecBlastContainer), 'Item in recblast_out was not a RecBlastContainer object!'
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
                        query_dict[target_id] = [''.join(id_search(annotation, id_type='brute', verbose=0)[1][0])
                                                 for annotation in annotation_list]
                        for annotation in query_dict[target_id]:
                            print(annotation, indent=3)
            print('*******************************************')
        return master_dict


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






