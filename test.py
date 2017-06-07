import recblast_MP as rb
import Auxilliary
from Bio.Blast import Record
import itertools
from pytest_mock import mocker
import os

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


def test_blast_1(mocker, db_loc=BLASTDB):
    """Checks BLAST searches done using the blast() function"""

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
