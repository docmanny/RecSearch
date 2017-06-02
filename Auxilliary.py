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


def count_dups(recblast_out):
    from recblast_MP import id_search
    master_dict = {}
    pat = re.compile('\|\[(.*)\]\|')  # regex for items in annotation

    for rc in recblast_out:
        try:
            rc.__delitem__('__dict__')
        except KeyError:
            pass
        for species, rc_spec_rec in rc.items():
            print('Species:\t', species, indent=0)
            try:
                species_dict = master_dict[species]
            except KeyError:
                master_dict[species] = dict()
                species_dict = master_dict[species]
            for gene, rc_rec in rc_spec_rec.items():
                print('Gene:\t' ,gene, indent=1)
                try:
                    gene_dict = species_dict[gene]
                except KeyError:
                    species_dict[gene] = dict()
                    gene_dict = species_dict[gene]
                try:
                    rc_out = rc_rec['recblast_results']
                except KeyError:
                    print('No entries in recblast_results for query {0} in species {1}'.format(gene, species))
                    continue
                for record in rc_out:
                    try:
                        print(record.description, indent=3)
                        target_id, annotations = record.description.split('|-|')
                        print(target_id, indent=4)
                        print(annotations.lstrip('\t'), indent=4)
                    except ValueError:
                        print(record.description, indent=2)
                        print('Could not unpack annotations!', indent=2)
                        continue
                    try:
                        target_list = gene_dict[target_id]
                    except KeyError:
                        gene_dict[target_id] = list()
                        target_list = gene_dict[target_id]
                    id_lst = pat.findall(annotations)
                    print(id_lst, indent=4)
                    if id_lst:
                        target_list += id_lst
                    else:
                        print('No annotations found for record {0} in species {1}, gene {2}'.format(record.name,
                                                                                                    species,
                                                                                                    gene))
    for species, species_dict in master_dict.items():
        for gene, gene_dict in species_dict.items():
            for target_id, annotation_list in gene_dict.items():
                for annotation in annotation_list:
                    _, id_list_ids, seq_range, _ = id_search(annotation, id_type='brute', verbose=0)
                    print(id_list_ids)




