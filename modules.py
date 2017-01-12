"""
Author: Juan Manuel Vazquez
Date Created: 11/01/2016
Date Last Modified: 11/02/2016
"""

"""
crosscheck is a function that looks through a series of .fasta files, and checks which IDs in the headers are common
between them
"""
#def crosscheck():
    # TODO: Write Crosscheck


"""
RecBlast is a function that will perform a reciprocal blast of a given sequence between two species. It takes a file
containing sequences as its input, in addition to the type of sequence and the names of species. First, it blasts
the sequence of choice from Species 1 to Species 2; it then takes all hits above the user-selected threshold,
and stores them in Fasta format. Then, for each hit, it will blast the FASTA sequence of Species 2 against Species 1;
it will then obtain the names of the top hits and append them in order to the header of the Species 2 FASTA sequence.
"""


def RecBlast(seqfile, species2, database1, filetype="fasta", species1="Homo Sapiens", blasttype='blastp',
             localblast1=True, localblast2=False, database2="refseq_proteins", evalue=0.001, identitiesperc=0.75,
             scoreperc=0.75, lengthperc=0.75, idthres='Score'):
    import os
    from Bio import SeqIO
    from Bio.Blast import NCBIWWW, NCBIStandalone
    from Bio.Blast.Applications import NcbiblastxCommandline,NcbiblastnCommandline,NcbiblastpCommandline, \
        NcbitblastxCommandline,NcbitblastnCommandline
    print("Now starting...")
    # First loop will iterate over each sequence in a file, preferably FASTA but also allows for GenBank
    for index, seq_record in enumerate(SeqIO.parse(seqfile, filetype)):
        print("{}: {}".format(index,seq_record.name))

        """
        recblast_out is the file where the final Blast output will go to for this particular sequence. This is the
        annotated FASTA, which will contain the top hits of the first BLAST, annotated with the name of the
        second BLAST
        """
        recblast_out = os.path.join(os.getcwd(),
                                    "{0}_recblast_out".format(species2),
                                    "{0}_{1}.fasta".format(blasttype, seq_record.name))

        # Begin by opening recblast_out, and then start with the primary BLAST
        with open(recblast_out, "w+") as blastout:
            length = len(seq_record)
            if localblast1:
                if blasttype == "blastn":
                    NcbiblastnCommandline(query=seq_record.seq, db=database1, evalue=evalue, outfmt=5,
                                          out="{0}to{1}blast_{2}.xml".format(species1, species2,
                                                                             seq_record.name))
                elif blasttype == "blastp":
                    NcbiblastpCommandline(query=seq_record.seq, db=database1, evalue=evalue, outfmt=5,
                                          out="{0}to{1}blast_{2}.xml".format(species1, species2,
                                                                             seq_record.name))
                elif blasttype == "blastx":
                    NcbiblastxCommandline(query=seq_record.seq, db=database1, evalue=evalue, outfmt=5,
                                          out="{0}to{1}blast_{2}.xml".format(species1, species2,
                                                                             seq_record.name))
                elif blasttype == "tblastx":
                    NcbitblastxCommandline(query=seq_record.seq, db=database1, evalue=evalue, outfmt=5,
                                           out="{0}to{1}blast_{2}.xml".format(species1, species2,
                                                                              seq_record.name))
                elif blasttype == "tblastn":
                    NcbitblastnCommandline(query=seq_record.seq, db=database1, evalue=evalue, outfmt=5,
                                           out="{0}to{1}blast_{2}.xml".format(species1, species2,
                                                                              seq_record.name))
                else:
                    raise Exception("Invalid blast choice!")
            else:
                with open(os.path.join(os.getcwd(),
                                       "recblast_out",
                                       "{0}to{1}blast_{2}.xml".format(species1, species2, seq_record.name)),
                          "w+") as fxml:
                    fxml.write(NCBIWWW.qblast(program=blasttype, database=database1, sequence=seq_record.seq,
                                              entrez_query=(species2 + "[ORGN]"), expect=evalue))

            # Start Reciprocal Blast
            # First load the primary BLAST hits to a handle.
            blasthits = open("{0}to{1}blast_{2}.xml".format(species1, species2, seq_record.name),"r")

            """
            Next, load up the error parsers:
                errorfile is the error log for the first BLAST parser, which will be used to get top hits of the
                        forward direction BLAST, for use in the reverse BLAST
                errorfile2 is the error log for the second BLAST parser, which will be used to get the top hits of the
                        reverse BLAST, which then will be used to annotate the hits from the forward BLAST
            """
            errorfile = open(os.path.join(os.getcwd(), "recblast_out", "{}_recblasterror1.log".format(seqfile)), "w+")
            errorfile2 = open(os.path.join(os.getcwd(), "recblast_out", "{}_recblasterror2.log".format(seqfile)), "w+")
            blast_error_parser = NCBIStandalone.BlastErrorParser(errorfile)
            blast_error_parser2 = NCBIStandalone.BlastErrorParser(errorfile2)

            # Iterate through all the blast hits
            for blastrecord in NCBIStandalone.Iterator(blasthits, blast_error_parser):
                for alignment in blastrecord.alignments:
                    scorelist = []
                    topscore = 0
                    # We will loop through each alignment twice, first to get
                    for hsp in alignment.hsps:
                        scorelist.append(hsp.score)
                    scorelist.sort()
                    topscore = scorelist[0]
                    for hsp in alignment.hsps:
                        if (hsp.score >= (scoreperc * topscore) & hsp.expect <= evalue &
                            hsp.align_length >= (length * lengthperc) & hsp.identities >=identitiesperc):

                            if localblast2:
                                if blasttype == "blastn":
                                    NcbiblastnCommandline(query=hsp.query,db=database2,evalue=evalue,outfmt=5,
                                                          out="{0}to{1}blast_{2}.xml".format(species2,species1,
                                                                                          seq_record.name))
                                elif blasttype == "blastp":
                                    NcbiblastpCommandline(query=hsp.query,db=database2,evalue=evalue,outfmt=5,
                                                          out="{0}to{1}blast_{2}.xml".format(species2,species1,
                                                                                          seq_record.name))
                                elif blasttype == "blastx":
                                    NcbiblastxCommandline(query=hsp.query,db=database2,evalue=evalue,outfmt=5,
                                                          out="{0}to{1}blast_{2}.xml".format(species2,species1,
                                                                                          seq_record.name))
                                elif blasttype == "tblastx":
                                    NcbitblastxCommandline(query=hsp.query,db=database2,evalue=evalue,outfmt=5,
                                                          out="{0}to{1}blast_{2}.xml".format(species2,species1,
                                                                                          seq_record.name))
                                elif blasttype == "tblastn":
                                    NcbitblastnCommandline(query=hsp.query,db=database2,evalue=evalue,outfmt=5,
                                                          out="{0}to{1}blast_{2}.xml".format(species2,species1,
                                                                                          seq_record.name))
                                else:
                                    raise Exception("Invalid blast choice!")
                            else:
                                with open(os.path.join(os.getcwd(), "recblast_out", "{0}to{1}blast_{2}.xml".format(
                                        species2,species1,seq_record.name), "w+")) as fxml:
                                    fxml.write(NCBIWWW.qblast(program=blasttype,database=database2,sequence= hsp.query,
                                               entrez_query=(species1 + "[ORGN]"),expect=evalue))
                            recblasthits = open("{1}to{0}blast_{2}.xml".format(species1, species2, seq_record.name))
                            for recblastrecord in NCBIStandalone.Iterator(recblasthits, blast_error_parser2):
                                hit_names = []
                                if str(idthres).isnumeric():
                                    it_idthres=0
                                    while idthres>it_idthres:
                                        hit_names.append(recblastrecord.alignments[idthres].title)
                                        idthres += 1
                                else:
                                    for alignment2 in recblastrecord.alignments:
                                        hit_names = []
                                        scorelist2 = []
                                        topscore2 = 0
                                        for hsp2 in alignment2.hsps:
                                            scorelist2.append(hsp2.score)
                                        scorelist2.sort()
                                        topscore2 = scorelist2[0]
                                        for hsp2 in alignment2.hsps:
                                            if (hsp2.score >= scoreperc * topscore2 & hsp2.expect <= evalue &
                                                hsp2.align_length >= length * lengthperc & 
                                                hsp2.identities >= identitiesperc):
                                                hit_names.append(alignment2.title)
                            recblasthits.close()
                            hit_title = ''.join(alignment.title,hit_names)
                            blastout.write(''.join('****Alignment****','\n'))
                            blastout.write(''.join('sequence:', hit_title,'\n'))
                            blastout.write(''.join('length:', alignment.length,'\n'))
                            blastout.write(''.join('e value:', hsp.expect,'\n'))
                            blastout.write(''.join(hsp.query,'\n'))
                            blastout.write(''.join(hsp.match,'\n'))
                            blastout.write(''.join(hsp.sbjct,'\n'))
            blasthits.close()
            errorfile.close()
            errorfile2.close()


"""
RecBlast2 uses the experimental Biopython SearchIO module. Honestly it looks so much cleaner from the docs so I kinda
jumped ship to fix RecBlast without the whole "for x in y / for y in z / for z in a / ..." hell.
"""


def RecBlast2(seqfile, species2, database1, filetype="fasta", species1="Homo Sapiens", blasttype='blastp',
              localblast1=True, localblast2=False, database2="refseq_proteins", evalue=0.001, identitiesperc=0.75,
              scoreperc=0.75, lengthperc=0.75, idthres='Score'):
    import os
    from Bio import SeqIO
    from Bio.Blast import NCBIWWW, NCBIStandalone
    from Bio.Blast.Applications import NcbiblastxCommandline,NcbiblastnCommandline,NcbiblastpCommandline,NcbitblastxCommandline,NcbitblastnCommandline
    print("Now starting...")
    # First loop will iterate over each sequence in a file, preferably FASTA but also allows for GenBank
    for index, seq_record in enumerate(SeqIO.parse(seqfile, filetype)):
        print("{}: {}".format(index,seq_record.name))

        # recblast_out is the file where the final Blast output will go to for this particular sequence. This is the annotated
        #  FASTA, which will contain the top hits of the first BLAST, annotated with the name of the second BLAST
        # errorfile is the error log for the first BLAST parser, which will be used to get top hits for the reciprocal
        #  BLAST
        # errorfile2 is the error log for the second BLAST parser, which will be used to get the top hits of the
        #  reciprocal BLAST, which then will be used to annotate the hits from the primary BLAST
        recblast_out = os.path.join(os.getcwd(), "{}_recblast_out".format(species2), "{0}_{1}.fasta".format(blasttype,seq_record.name))
        errorfile = open(os.path.join(os.getcwd(), "recblast_out", "{}_recblasterror1.log".format(seqfile)),"w+")
        errorfile2 = open(os.path.join(os.getcwd(), "recblast_out", "{}_recblasterror2.log".format(seqfile)), "w+")

        # Begin by opening recblast_out, and then start with the primary BLAST
        with open(recblast_out, "w+") as blastout:
            length = len(seq_record)
            if localblast1:
                if blasttype == "blastn":
                    NcbiblastnCommandline(query=seq_record.seq,db=database1,evalue=evalue,outfmt=5,
                                          out="{0}to{1}blast_{2}.xml".format(species1, species2,
                                                                             seq_record.name))
                elif blasttype == "blastp":
                    NcbiblastpCommandline(query=seq_record.seq,db=database1,evalue=evalue,outfmt=5,
                                          out="{0}to{1}blast_{2}.xml".format(species1, species2,
                                                                             seq_record.name))
                elif blasttype == "blastx":
                    NcbiblastxCommandline(query=seq_record.seq,db=database1,evalue=evalue,outfmt=5,
                                          out="{0}to{1}blast_{2}.xml".format(species1, species2,
                                                                             seq_record.name))
                elif blasttype == "tblastx":
                    NcbitblastxCommandline(query=seq_record.seq,db=database1,evalue=evalue,outfmt=5,
                                           out="{0}to{1}blast_{2}.xml".format(species1, species2,
                                                                              seq_record.name))
                elif blasttype == "tblastn":
                    NcbitblastnCommandline(query=seq_record.seq,db=database1,evalue=evalue,outfmt=5,
                                           out="{0}to{1}blast_{2}.xml".format(species1, species2,
                                                                              seq_record.name))
                else:
                    raise Exception("Invalid blast choice!")
            else:
                with open(os.path.join(os.getcwd(), "recblast_out", "{0}to{1}blast_{2}.xml".format(species1, species2,
                                        seq_record.name)), "w+") as fxml:
                    fxml.write(NCBIWWW.qblast(program=blasttype,database=database1,sequence=seq_record.seq,
                                              entrez_query=(species2 + "[ORGN]"),expect=evalue))

            # Start Reciprocal Blast
            # First load the primary BLAST hits to a handle.
            blasthits = open("{0}to{1}blast_{2}.xml".format(species1, species2, seq_record.name),"r")
            # Next, load up the error parsers
            blast_error_parser = NCBIStandalone.BlastErrorParser(errorfile)
            blast_error_parser2 = NCBIStandalone.BlastErrorParser(errorfile2)

            # Iterate through all the blast hits


""""
SirBlastALot() is a function that, for a given sequence, BLASTS it against a list of organisms of one's choice; the
default list has an extensive vertebrate coverage. It also has an option to use BLAT instead. By default it just does
does a unidirectional blast of your sequence to each organism; you can set it to do a Reciprocal Blast as well.
"""
