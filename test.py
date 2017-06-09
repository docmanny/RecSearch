import os

import pytest
from Bio.Blast.Record import Blast as BioBlastRecord
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
    scenarios = [('', {})]
    pass


class Test_getsearchdb(object):
    scenarios = [('', {})]
    pass


class Test_blatserver(object):
    scenarios = [('', {})]
    pass


class Test_id_ranker(object):
    scenarios = [('', {})]
    pass


class Test_biosql_get_sub_db_names(object):
    scenarios = [('', {})]
    pass


class Test_biosql_DBSeqRecord_to_SeqRecord(object):
    """Class containing all tests related to biosql_DBSeqRecord_to_SeqRecord"""
    scenarios = [('InRecord = SeqRecord', {'inrecord': SeqRecord(seq='AATTGGCC', id='test'),
                                          'off': False}),
                 ('InRecord = DBSeqRecord; no conversion', {'inrecord': DBSeqRecord,
                                                           'off' : True}),
                 ('InRecord = DBSeqRecord; yes conversion', {'inrecord': DBSeqRecord,
                                                             'off': False}),
                 ('InRecord = str', {'inrecord': str,
                                     'off': True})]

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
    scenarios = [('', {})]
    pass


class Test_Blast(object):
    """
    full_list = list(itertools.product(['CDN1B.faprt', 'INS.faprt', 'UCHL3.faprt'],
                                       ['Myotis lucifugus', 'Pteropus vampyrus', 'Homo sapiens'],
                                       ['blastp'], [False], ['refseq_protein'],['Homo sapiens']))
    full_list += list(itertools.product(['CDN1B.faprt', 'INS.faprt', 'UCHL3.faprt'],
                                        ['Myotis lucifugus', 'Pteropus vampyrus', 'Homo sapiens'],
                                        ['tblastn'], [False], ['refseq_genome'],['Homo sapiens']))
    full_list += list(itertools.product(['CDKN1B.fa', 'INS.fa', 'UCHL3.fa'],
                                        ['Myotis lucifugus', 'Pteropus vampyrus', 'Homo sapiens'],
                                        ['blastx'], [False], ['refseq_protein'],['Homo sapiens']))
    full_list += list(itertools.product(['CDKN1B.fa', 'INS.fa', 'UCHL3.fa'],
                                        ['Myotis lucifugus', 'Pteropus vampyrus', 'Homo sapiens'],
                                        ['tblastx'], [False], ['refseq_genome'],['Homo sapiens']))
    full_list += list(itertools.product(['CDKN1B.fa', 'INS.fa', 'UCHL3.fa'],
                                        ['Myotis lucifugus', 'Pteropus vampyrus', 'Homo sapiens'],
                                        ['blastn'], [False], ['refseq_genome'],['Homo sapiens']))

    scenarios = [('{blast_type}-{seq_record}-{target_species}-{local}-'.format(seq_record=i[0],
                                                                               target_species=i[1],
                                                                               blast_type=i[2],
                                                                               local='local' if i[3] else
                                                                                     'remote').replace(' ', '_') +
                  '{database}-{query_species}'.format(query_species=i[5],
                                                      database=i[4]).replace(' ', '_'),
                  dict(seq_record='Test/query/{}'.format(i[0]),
                       target_species=i[1],
                       blast_type=i[2],
                       database=i[4],
                       local_blast=i[3],
                       query_species=i[5])
                  ) for i in full_list]"""
    scenarios =[('blastp-CDN1B.faprt-Myotis_lucifugus-remote-refseq_protein-Homo_sapiens',
                 {'blast_type': 'blastp',
                  'database': 'refseq_protein',
                  'local_blast': False,
                  'query_species': 'Homo sapiens',
                  'seq_record': 'Test/query/CDN1B.faprt',
                  'target_species': 'Myotis lucifugus'}),
                ('blastp-CDN1B.faprt-Myotis_lucifugus-local-Myotis_lucifugus_protein_v2.0.fa-Homo_sapiens',
                 {'blast_type': 'blastp',
                  'database': 'Myotis_lucifugus_protein_v2.0.fa',
                  'local_blast': True,
                  'query_species': 'Homo sapiens',
                  'seq_record': 'Test/query/CDN1B.faprt',
                  'target_species': 'Myotis lucifugus'})
                ]
    @pytest.mark.long
    def test_search_remote_real(self, seq_record, target_species, database, query_species, blast_type, local_blast,
                                BLASTDB=globals()['BLASTDB']):
        results, err = rb.blast(seq_record=seq_record, target_species=target_species, database=database,
                                query_species=query_species, blast_type=blast_type, local_blast=local_blast,
                                BLASTDB=BLASTDB)
        assert isinstance(results, BioBlastRecord)
        assert err == '' or err == [] or err == None

    def test_search_local_mocked(self, seq_record, target_species, database, query_species, blast_type, local_blast,
                               BLASTDB=globals()['BLASTDB']):
        pass


class Test_biosql_seq_lookup_cascade(object):
    scenarios = [('', {})]
    pass


class Test_GetSeqMP(object):
    scenarios = [('', {})]
    pass


class Test_BioSeqLookupCascade(object):
    scenarios = [('', {})]
    pass


class Test_biosql_get_record_mp(object):
    scenarios = [('', {})]
    pass


class Test_id_search(object):
    scenarios = [('', {})]

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
    scenarios = [('', {})]
    pass


class Test_FetchSeq(object):
    scenarios = [('', {})]
    pass


class Test_fetchseqMP(object):
    scenarios = [('', {})]
    pass


class Test_RecBlastMP_Thread(object):
    scenarios = [('', {})]
    pass


class Test_RecBlast(object):
    scenarios = [('', {})]
    pass


class Test_recblastMP(object):
    scenarios = [('', {})]
    pass


class Test_RecBlastControl(object):
    scenarios = [('', {})]
    pass
