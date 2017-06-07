import itertools
import os

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




def test_get_searchdb():
    """Checks get_searchdb()"""
    print('Todo: write test_get_searchdb')


def test_blast(mocker, db_loc=BLASTDB):
    """Checks BLAST searches done using the blast() function"""

    # Internal test variables
    target_list = ['Myotis lucifugus', 'Pteropus vampyrus', 'Homo sapiens']
    prot_query_list = ['CDN1B.faprt', 'INS.faprt', 'UCHL3.faprt']
    nuc_query_list = ['CDKN1B.fa', 'INS.fa', 'UCHL3.fa']
    blast_list = ['blastn', 'blastx', 'tblastn', 'tblastx', 'blastp']
    remote_list = [True, False]
    db = [rb.get_searchdb(search_type, species, db_loc=db_loc, verbose=0) for search_type, species in
          list(itertools.product(blast_list, target_list))]
    outlist = []
    
    #print('Sequence', 'Target_species', 'Blast_type', 'Local?', 'Pass/Fail', sep='\t')
    for r in remote_list:
        if not r:
            mocker.patch('recblast_MP.NCBIWWW.qblast', return_value=(1,0))
            d = 'refseq_genomic'
            for target in target_list:
                for bt in blast_list:
                    if bt in ['blastn', 'blastx', 'tblastx']:
                        for nuc in nuc_query_list:
                            arglist = [nuc, target, bt, r]
                            out, err = rb.blast(seq_record=root_dir+'/Test/query/'+nuc,
                                                target_species=target, blast_type=bt,
                                                    local_blast=r, database=d)
                            assert out
                            print(['{}\t'.format(i) for i in arglist], bool(out))
                    else:
                        for prot in prot_query_list:
                            arglist = [prot, target, bt, r]
                            out, err = rb.blast(seq_record=root_dir+'/Test/query/'+prot,
                                                target_species=target, blast_type=bt,
                                                    local_blast=r, database=d)
        else:
            mocker.patch('recblast_MP.subprocess.Popen', return_value=(1,0))
            mocker.patch('recblast_MP.subprocess.check_output', return_value=(1,0))
            for d in db:
                for target in target_list:
                    for bt in blast_list:
                        if bt in ['blastn', 'blastx', 'tblastx']:
                            for nuc in nuc_query_list:
                                arglist = [nuc, target, bt, r]
                                out, err = rb.blast(seq_record=root_dir+'/Test/query/'+nuc,
                                                    target_species=target, blast_type=bt,
                                                        local_blast=r, database=d)
                        else:
                            for prot in prot_query_list:
                                arglist = [prot, target, bt, r]
                                out, err = rb.blast(seq_record=root_dir+'/Test/query/'+prot,
                                                    target_species=target, blast_type=bt,
                                                        local_blast=r, database=d)


def pytest_generate_tests(metafunc):
    idlist = []
    argvalues = []
    for scenario in metafunc.cls.scenarios:
        idlist.append(scenario[0])
        items = scenario[1].items()
        argnames = [x[0] for x in items]
        argvalues.append(([x[1] for x in items]))
    metafunc.parametrize(argnames, argvalues, ids=idlist, scope="class")

target_list = ['Myotis lucifugus', 'Pteropus vampyrus', 'Homo sapiens']
prot_query_list = ['CDN1B.faprt', 'INS.faprt', 'UCHL3.faprt']
nuc_query_list = ['CDKN1B.fa', 'INS.fa', 'UCHL3.fa']
blast_list = ['blastn', 'blastx', 'tblastn', 'tblastx', 'blastp']
remote_list = [True, False]
db = [rb.get_searchdb(search_type, species, db_loc=db_loc, verbose=0) for search_type, species in
      list(itertools.product(blast_list, target_list))]

scenario_kwarglist = [dict(seqfile = scenario[0],
                           target_species = scenario[1],
                           fw_blast_db = scenario[2],
                           rv_blast_db = scenario[3],
                           infile_type = scenario[4],
                           output_type = scenario[5],
                           host = scenario[6],
                           user = scenario[7],
                           driver = scenario[8],
                           query_species = scenario[9],
                           blast_type_1 = scenario[10],
                           blast_type_2 = scenario[11],
                           local_blast_1 = scenario[12],
                           local_blast_2 = scenario[13],
                           expect = scenario[14],
                           perc_score = scenario[15],
                           perc_ident = scenario[16],
                           perc_length = scenario[17],
                           megablast = scenario[18],
                           email = scenario[19],
                           id_type = scenario[20],
                           fw_source = scenario[21],
                           fw_id_db = scenario[22],
                           fetch_batch_size = scenario[23],
                           passwd = scenario[24],
                           fw_id_db_version = scenario[25],
                           BLASTDB = scenario[26],
                           indent = scenario[27],
                           verbose = scenario[28],
                           max_n_processes = scenario[29],
                           n_threads = scenario[30],
                           write_intermediates = scenario[31],
                           write_final = scenario[32],
                           fw_blast_kwargs = scenario[33],
                           rv_blast_kwargs = scenario[34]
                           )
                      for scenario in scenario_list]


class TestBlast(object):
    """def __init__(self):
        self.target_list = ['Myotis lucifugus', 'Pteropus vampyrus', 'Homo sapiens']
        self.prot_query_list = ['CDN1B.faprt', 'INS.faprt', 'UCHL3.faprt']
        self.nuc_query_list = ['CDKN1B.fa', 'INS.fa', 'UCHL3.fa']
        self.blast_list = ['blastn', 'blastx', 'tblastn', 'tblastx', 'blastp']
        #self.remote_list = [True, False]"""
    def test_local_blast(self):
        print('local')
    def test_remote_blast(self):
        print('remote')
    def __run__(self):
        self.test_local_blast()
        self.test_remote_blast()

####################################################################################3
import os
try:
    BLASTDB=os.environ['BLASTDB']
    BLATDB=os.environ['BLATDB']
except KeyError:
    print([key for key, _ in os.environ.items()])
    BLASTDB='~/db/blastdb'
    BLATDB='~/db/blatdb'

seqfile_gene = ['CDN1B.fa', ['CDN1B.fa'], 'INS.fa', 'UCHL3.fa', '3gene.fasta', ['CDN1B.fa']]
seqfile_prot = ['CDN1B.faprt', 'INS.faprt', 'UCHL3.faprt', '3protein.fasta', ['CDN1B.faprt']]
target_species = ['Myotis lucifugus', 'Pteropus vampyrus', 'Homo sapiens', ['Myotis lucifugus']],
fw_blast_db = ['auto']
rv_blast_db = ['auto']
fw_blat_db = [{'Myotis lucifugus': 20001, 'Pteropus vampyrus': 20003, 'Homo sapiens': 20004},
              {'Myotis lucifugus': 30001, 'Pteropus vampyrus': 30003, 'Homo sapiens': 30004},
              {'Myotis lucifugus': 'auto', 'Pteropus vampyrus': 'auto', 'Homo sapiens': 'auto'}]
rv_blat_db = [{'Myotis lucifugus': 20001, 'Pteropus vampyrus': 20003, 'Homo sapiens': 20004},
              {'Myotis lucifugus': 30001, 'Pteropus vampyrus': 30003, 'Homo sapiens': 30004},
              {'Myotis lucifugus': 'auto', 'Pteropus vampyrus': 'auto', 'Homo sapiens': 'auto'}]
infile_type = ['fasta'],
output_type = ['fasta'],
host = ['localhost', '127.0.0.1'],
user = ['postgres'], # Todo: config .travis.yml & add db: https://docs.travis-ci.com/user/database-setup/#PostgreSQL
driver = ['psycopg2'], #Todo: config requirements.txt with pscyopg2
query_species = ['Homo sapiens', ['Homo sapiens']]
blast_type_1 = ['blastn', 'blastx', 'tblastn', 'tblastx', 'blastp', 'blat', 'tblat']
blast_type_2 = ['blastn', 'blastx', 'tblastn', 'tblastx', 'blastp', 'blat', 'tblat']
local_blast_1 = [True, False],
local_blast_2 = [True, False],
expect = [10],
perc_score = [0.5],
perc_ident = [50],
perc_length = [0.5],
megablast = [True, False],
email = ['johndoe@example.com'],
id_type = ['gi', 'accession', 'scaffold', 'id', 'chr', 'brute'],
fw_source = ['sql', 'fasta', '2bit'],
fw_id_db = ['auto'],
fetch_batch_size = [50],
passwd = [''],
fw_id_db_version = ['auto'],
BLASTDB = [BLASTDB, BLATDB],
indent = [0],
verbose = ['vvv'],
max_n_processes = [1,5],
n_threads = [1],
write_intermediates = [True, False],
write_final = [True, False],
fw_blast_kwargs = [dict()],
rv_blast_kwargs = [dict()]