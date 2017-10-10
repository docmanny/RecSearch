# Copyright 2017 by Juan Manuel Vazquez. All Rights Reserved.
# This package is governed by its license as mentioned in the
# LICENSE file that should have been included as a part of this
# package.
""" A collection of tools to simplify and analyze Reciprocal Best-Hit Searches.

https://www.vazquez.bio/RecBlast
"""

__version__ = 'v1.3dev'

# Functions that are essential yet unrelated to recblast_MP.py

class ProgressBar(object):
    """Adapted from Romuald Brunet at StackExchange"""
    DEFAULT = 'Progress: %(bar)s %(percent)3d%%'
    FULL = '%(bar)s %(current)d/%(total)d (%(percent)3d%%) %(remaining)d to go'

    def __init__(self, total, width=100, fmt=DEFAULT, symbol='=',
                 output=sys.stdout):
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

