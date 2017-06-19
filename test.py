import os
from io import StringIO

import pytest
from Bio.Blast.Record import Blast as BioBlastRecord
from Bio.SeqRecord import SeqRecord
# from pytest_mock import mocker
from BioSQL.BioSeq import DBSeqRecord

import Auxilliary as aux
import recblast_MP as rb

try:
    BLASTDB = os.environ['BLASTDB']
    BLATDB = os.environ['BLATDB']
except KeyError:
    print([key for key, _ in os.environ.items()])
    BLASTDB = '~/db/blastdb'
    BLATDB = '~/db/blatdb'
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
                                                            'off': True}),
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
    """('blastp-CDN1B.faprt-Myotis_lucifugus-remote-refseq_protein-Homo_sapiens',
                 {'blast_type': 'blastp',
                  'database': 'refseq_protein',
                  'local_blast': False,
                  'query_species': 'Homo sapiens',
                  'seq_record': 'Test/query/CDN1B.faprt',
                  'target_species': 'Myotis lucifugus'}),"""
    scenarios = [('tblastn-CDN1B.faprt-Myotis_lucifugus-local-Myotis_lucifugus_protein_v2.0.fa-Homo_sapiens',
                 {'blast_type': 'tblastn',
                  'database': 'Myotis_lucifugus_genome_v2.0.fa',
                  'local_blast': True,
                  'query_species': 'Homo sapiens',
                  'seq_record': 'CDN1B.faprt',
                  'target_species': 'Myotis lucifugus'})
                 ]

    @pytest.mark.long
    def test_search_remote_real(self, seq_record, target_species, database, query_species, blast_type, local_blast,
                                BLASTDB=globals()['BLASTDB']):
        results, err = rb.blast(seq_record='Test/query/{}'.format(seq_record), target_species=target_species,
                                database=database,
                                query_species=query_species, blast_type=blast_type, local_blast=local_blast,
                                BLASTDB=BLASTDB)
        assert isinstance(results, BioBlastRecord)
        assert err == '' or err == [] or err is None

    def test_search_local_mocked(self, mocker, seq_record, target_species, database, query_species, blast_type,
                                 local_blast, BLASTDB=globals()['BLASTDB']):

        file_handle = 'Test/BLAST/'
        file_handle += '{blast_type}_{local_blast}_'.format(blast_type=blast_type,
                                                            local_blast='Local' if local_blast else 'Remote')
        file_handle += '{query_species}_{target_species}_{seq_record}.xml'.format(query_species=query_species,
                                                                                  target_species=target_species,
                                                                                  seq_record=seq_record.split('.fa')[0])
        file_handle = file_handle.replace(' ', '_')
        with open(file_handle) as fin:
            mockingbird = fin.readlines()
        if not local_blast:
            mockingbird = StringIO(mockingbird)

        mocker.patch('recblast_MP.subprocess.Popen', return_value=(mockingbird, None))
        mocker.patch('recblast_MP.subprocess.check_output', return_value=(mockingbird, None))
        mocker.patch('recblast_MP.NCBIWWW.qblast', return_value=mockingbird)
        results, err = rb.blast(seq_record='Test/query/{}'.format(seq_record), target_species=target_species,
                                database=database, query_species=query_species, blast_type=blast_type,
                                local_blast=local_blast, BLASTDB=BLASTDB)
        assert isinstance(results, BioBlastRecord)
        assert err == '' or err == [] or err is None


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

    def test_id(self, id, id_type):
        if id_type == 'id':
            pass
        else:
            pass

    def test_gi(self, id, id_type):
        if id_type == 'gi':
            pass
        else:
            pass

    def test_accession(self, id, id_type):
        if id_type == 'accession':
            pass
        else:
            pass

    def test_chr(self, id, id_type):
        if id_type == 'chr':
            pass
        else:
            pass

    def test_scaffold(self, id, id_type):
        if id_type == 'scaffold':
            pass
        else:
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
    full_list = []
    scenarios = [('',
                  dict(seqfile=i[0],
                       target_species=i[1],
                       fw_blast_db=i[2],
                       rv_blast_db=i[3],
                       infile_type=i[4],
                       output_type=i[5],
                       host=i[6],
                       user=i[7],
                       driver=i[8],
                       query_species=i[9],
                       blast_type_1=i[10],
                       blast_type_2=i[11],
                       local_blast_1=i[12],
                       local_blast_2=i[13],
                       expect=i[14],
                       perc_score=i[15],
                       perc_ident=i[16],
                       perc_length=i[17],
                       megablast=i[18],
                       email=i[19],
                       id_type=i[20],
                       fw_source=i[21],
                       fw_id_db=i[22],
                       fetch_batch_size=i[23],
                       passwd=i[24],
                       fw_id_db_version=i[25],
                       BLASTDB=i[26],
                       indent=i[27],
                       verbose=i[28],
                       max_n_processes=i[29],
                       n_threads=i[30],
                       write_intermediates=i[31],
                       write_final=i[32],
                       fw_blast_kwargs=i[33],
                       rv_blast_kwargs=i[34])
                  ) for i in full_list]
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


class Test_count_dups(object):
    scenarios = [('Not a RecBlastContainer', {'rc_out': [dict()],
                                              'empirical_count': None}),
                 ('Empty RBC', {'rc_out': rb.RecBlastContainer(target_species=None,
                                                               query_record=SeqRecord(name='', seq='')),
                                'empirical_count': {'': 0}}),
                 ('', {}),
                 ('', {}),
                 ('', {})]
    pass

    def test_input(self, rc_out, empirical_count):
        if isinstance(rc_out, rb.RecBlastContainer):
            pass
        else:
            pytest.raises(AssertionError, message='Expecting AssertionError: '
                                                  'Item in recblast_out was not a RecBlastContainer object!')

    def test_count_dups(self, rc_out, empirical_count):
        if empirical_count is None:
            pass
        else:
            count_dict = aux.count_dups(rc_out)
        pass
