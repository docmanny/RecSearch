import os

import pytest
from Bio.SeqRecord import SeqRecord
# from pytest_mock import mocker
from BioSQL.BioSeq import DBSeqRecord

import recblast_MP as rb

try:
    BLASTDB=os.environ['BLASTDB']
    BLATDB=os.environ['BLATDB']
except KeyError:
    print([key for key, _ in os.environ.items()])
    BLASTDB='~/db/blastdb'
    BLATDB='~/db/blatdb'
root_dir = os.getcwd()
print('Root Dir:', root_dir)


"""
def pytest_generate_tests(metafunc):
    # called once per each test function

    funcarglist = metafunc.cls.params[metafunc.function.__name__]
    argnames = sorted(funcarglist[0])
    metafunc.parametrize(argnames, [[funcargs[name] for name in argnames]
            for funcargs in funcarglist])
"""


def pytest_generate_tests(metafunc):
    idlist = []
    argvalues = []
    for scenario in metafunc.cls.scenarios:
        idlist.append(scenario[0])
        items = scenario[1].items()
        argnames = [x[0] for x in items]
        argvalues.append(([x[1] for x in items]))
    metafunc.parametrize(argnames, argvalues, ids=idlist, scope="class")


class Test_percent_identity_searchio(object):
    scenarios = [('',{})]
    pass


class Test_getsearchdb(object):
    scenarios = [('', {})]
    pass


class Test_blatserver(object):
    scenarios = [('',{})]
    pass


class Test_id_ranker(object):
    scenarios = [('',{})]
    pass


class Test_biosql_get_sub_db_names(object):
    scenarios = [('',{})]
    pass

class Test_biosql_DBSeqRecord_to_SeqRecord(object):
    """Class containing all tests related to biosql_DBSeqRecord_to_SeqRecord"""
    scenarios = [('InRecord = SeqRecord',{'inrecord':SeqRecord(seq='AATTGGCC', id='test'),
                                          'off':False}),
                 ('InRecord = DBSeqRecord; no conversion',{'inrecord':DBSeqRecord,
                                                           'off':True}),
                 ('InRecord = DBSeqRecord; yes conversion', {'inrecord':DBSeqRecord,
                                                             'off':False}),
                 ('InRecord = str',{'inrecord':str,
                                    'off':True})]

    def test_check_inrecord_is_DBSeqRecord(self, inrecord, off):
        """Tests if raised AssertionError if inrecord is not DBSeqRecord"""
        if isinstance(inrecord, DBSeqRecord):
            if off:
                assert isinstance(rb.biosql_DBSeqRecord_to_SeqRecord(DBSeqRecord_=inrecord,
                                                                     off=off),
                                  DBSeqRecord)
            else:
                assert isinstance(rb.biosql_DBSeqRecord_to_SeqRecord(DBSeqRecord_=inrecord,
                                                                     off=off),
                                  SeqRecord)
        else:
            with pytest.raises(AssertionError):
                rb.biosql_DBSeqRecord_to_SeqRecord(DBSeqRecord_=inrecord,
                                                   off=off)
class Test_RecBlastContainer(object):
    scenarios = [('',{})]
    pass


class Test_biosql_seq_lookup_cascade(object):
    scenarios = [('',{})]
    pass


class Test_GetSeqMP(object):
    scenarios = [('',{})]
    pass


class Test_BioSeqLookupCascade(object):
    scenarios = [('',{})]
    pass


class Test_biosql_get_record_mp(object):
    scenarios = [('',{})]
    pass


class Test_id_search(object):
    scenarios = [('',{})]

    def test_id(self):
        pass
    def test_gi(self):
        pass
    def test_accession(self):
        pass
    def test_chr(self):
        pass
    def test_scaffold(self):
        pass
    def test_seqrange(self):
        pass


class Test_FetchSeqMP(object):
    scenarios = [('',{})]
    pass

class Test_FetchSeq(object):
    scenarios = [('',{})]
    pass


class Test_fetchseqMP(object):
    scenarios = [('',{})]
    pass


class Test_RecBlastMP_Thread(object):
    scenarios = [('',{})]
    pass


class Test_RecBlast(object):
    scenarios = [('',{})]
    pass


class Test_recblastMP(object):
    scenarios = [('',{})]
    pass


class Test_RecBlastControl(object):
    scenarios = [('',{})]
    pass
