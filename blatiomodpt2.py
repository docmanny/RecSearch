def get_processor(format, mapping):
    """Returns the object to process the given format according to the mapping.
    :param format: format name
    :type format: string, lower case
    :param mapping: mapping of format name and its processor object
    :type mapping: dictionary {string: object}
    """
    # map file format to iterator name
    try:
        obj_info = mapping[format]
    except KeyError:
        # handle the errors with helpful messages
        if format is None:
            raise ValueError("Format required (lower case string)")
        elif not isinstance(format, basestring):
            raise TypeError("Need a string for the file format (lower case)")
        elif format != format.lower():
            raise ValueError("Format string %r should be lower case" %
                             format)
        else:
            raise ValueError("Unknown format %r. Supported formats are "
                             "%r" % (format, "', '".join(mapping)))

    mod_name, obj_name = obj_info
    mod = __import__('RecBlast.etc.blatiomod', fromlist=[''])

    return getattr(mod, obj_name)


from __future__ import print_function

import sys
import warnings

from Bio import BiopythonExperimentalWarning
from Bio.File import as_handle
from Bio._py3k import basestring

warnings.warn('Bio.SearchIO is an experimental submodule which may undergo '
              'significant changes prior to its future official release.',
              BiopythonExperimentalWarning)

__all__ = ('read', 'parse', 'to_dict', 'index', 'index_db', 'write', 'convert')

# dictionary of supported formats for parse() and read()
_ITERATOR_MAP = {
    'blast-tab': ('BlastIO', 'BlastTabParser'),
    'blast-text': ('BlastIO', 'BlastTextParser'),
    'blast-xml': ('BlastIO', 'BlastXmlParser'),
    'blat-psl': ('RecBlast.etc.blatiomod.py', 'BlatPslParser'),
    'exonerate-cigar': ('ExonerateIO', 'ExonerateCigarParser'),
    'exonerate-text': ('ExonerateIO', 'ExonerateTextParser'),
    'exonerate-vulgar': ('ExonerateIO', 'ExonerateVulgarParser'),
    'fasta-m10': ('FastaIO', 'FastaM10Parser'),
    'hmmer2-text': ('HmmerIO', 'Hmmer2TextParser'),
    'hmmer3-text': ('HmmerIO', 'Hmmer3TextParser'),
    'hmmer3-tab': ('HmmerIO', 'Hmmer3TabParser'),
    # for hmmer3-domtab, the specific program is part of the format name
    # as we need it distinguish hit / target coordinates
    'hmmscan3-domtab': ('HmmerIO', 'Hmmer3DomtabHmmhitParser'),
    'hmmsearch3-domtab': ('HmmerIO', 'Hmmer3DomtabHmmqueryParser'),
    'phmmer3-domtab': ('HmmerIO', 'Hmmer3DomtabHmmqueryParser'),
}

# dictionary of supported formats for index()
_INDEXER_MAP = {
    'blast-tab': ('BlastIO', 'BlastTabIndexer'),
    'blast-xml': ('BlastIO', 'BlastXmlIndexer'),
    'blat-psl': ('BlatIO', 'BlatPslIndexer'),
    'exonerate-cigar': ('ExonerateIO', 'ExonerateCigarIndexer'),
    'exonerate-text': ('ExonerateIO', 'ExonerateTextIndexer'),
    'exonerate-vulgar': ('ExonerateIO', 'ExonerateVulgarIndexer'),
    'fasta-m10': ('FastaIO', 'FastaM10Indexer'),
    'hmmer2-text': ('HmmerIO', 'Hmmer2TextIndexer'),
    'hmmer3-text': ('HmmerIO', 'Hmmer3TextIndexer'),
    'hmmer3-tab': ('HmmerIO', 'Hmmer3TabIndexer'),
    'hmmscan3-domtab': ('HmmerIO', 'Hmmer3DomtabHmmhitIndexer'),
    'hmmsearch3-domtab': ('HmmerIO', 'Hmmer3DomtabHmmqueryIndexer'),
    'phmmer3-domtab': ('HmmerIO', 'Hmmer3DomtabHmmqueryIndexer'),
}

# dictionary of supported formats for write()
_WRITER_MAP = {
    'blast-tab': ('BlastIO', 'BlastTabWriter'),
    'blast-xml': ('BlastIO', 'BlastXmlWriter'),
    'blat-psl': ('BlatIO', 'BlatPslWriter'),
    'hmmer3-tab': ('HmmerIO', 'Hmmer3TabWriter'),
    'hmmscan3-domtab': ('HmmerIO', 'Hmmer3DomtabHmmhitWriter'),
    'hmmsearch3-domtab': ('HmmerIO', 'Hmmer3DomtabHmmqueryWriter'),
    'phmmer3-domtab': ('HmmerIO', 'Hmmer3DomtabHmmqueryWriter'),
}


def parse(handle, format=None, **kwargs):
    """Turns a search output file into a generator that yields QueryResult
    objects.
     - handle - Handle to the file, or the filename as a string.
     - format - Lower case string denoting one of the supported formats.
     - kwargs - Format-specific keyword arguments.
    This function is used to iterate over each query in a given search output
    file:
    >>> from Bio import SearchIO
    >>> qresults = SearchIO.parse('Blast/mirna.xml', 'blast-xml')
    >>> qresults
    <generator object ...>
    >>> for qresult in qresults:
    ...     print("Search %s has %i hits" % (qresult.id, len(qresult)))
    ...
    Search 33211 has 100 hits
    Search 33212 has 44 hits
    Search 33213 has 95 hits
    Depending on the file format, `parse` may also accept additional keyword
    argument(s) that modifies the behavior of the format parser. Here is a
    simple example, where the keyword argument enables parsing of a commented
    BLAST tabular output file:
    >>> from Bio import SearchIO
    >>> for qresult in SearchIO.parse('Blast/mirna.tab', 'blast-tab', comments=True):
    ...     print("Search %s has %i hits" % (qresult.id, len(qresult)))
    ...
    Search 33211 has 100 hits
    Search 33212 has 44 hits
    Search 33213 has 95 hits
    """
    # get the iterator object and do error checking
    iterator = get_processor(format, _ITERATOR_MAP)

    # HACK: force BLAST XML decoding to use utf-8
    handle_kwargs = {}
    if format == 'blast-xml' and sys.version_info[0] > 2:
        handle_kwargs['encoding'] = 'utf-8'

    # and start iterating
    with as_handle(handle, 'rU', **handle_kwargs) as source_file:
        generator = iterator(source_file, **kwargs)

        for qresult in generator:
            yield qresult


qresults = list(parse('/home/manny/Downloads/test.pslx', 'blat-psl', pslx=True))
