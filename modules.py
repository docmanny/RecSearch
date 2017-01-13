"""
Author: Juan Manuel Vazquez
Date Created: 11/01/2016
Date Last Modified: 11/02/2016
"""

"""
crosscheck is a function that looks through a series of .fasta files, and checks which IDs in the headers are common
between them
"""


def crosscheck():
    # TODO: Write Crosscheck
    pass


"""
blast just does as it implies - runs a blast from Species 1 to Species 2 for a series of things in the sequence.
Importantly, it has an option, "recblast", that serves to
"""


def blast(seqfile, target_species, database, query_species="Homo sapiens", filetype="fasta", blasttype='blastn',
             localblast=True, evalue=0.005, recblast_on=False, megablast=True, blastoutput_custom=""):
    import os
    from Bio import SeqIO
    from Bio.Blast import NCBIWWW
    from Bio.Blast.Applications import NcbiblastxCommandline, NcbiblastnCommandline, NcbiblastpCommandline, \
        NcbitblastxCommandline, NcbitblastnCommandline

    if recblast_on:
        length_list = []
    else:
        print("Now starting...")
    if blasttype != 'blastn':
        megablast = False
    # First loop will iterate over each sequence in a file, preferably FASTA but also allows for GenBank
    for index, seq_record in enumerate(SeqIO.parse(seqfile, filetype)):
        if recblast_on:
            # recblastout = open(os.path.join(os.getcwd(),
            #                            "{0}_recblast".format(target_species),
            #                            "{0}_{1}.list".format(blasttype, seq_record.name)), "a+")
            # length_list.append(len(seq_record))
            pass
        else:
            print("{}: {}".format(index, seq_record.name))

        # Begin by opening recblast_out, and then start with the primary BLAST
        if blastoutput_custom == '':
            blastoutput_custom = os.path.join(os.getcwd(),
                                              "{0}_blast".format(target_species),
                                              "{0}_{1}_{2}_to_{3}.xml".format(blasttype, seq_record.name,
                                                                              query_species, target_species))
        if localblast:
            if blasttype == "blastn":
                NcbiblastnCommandline(query=seq_record.seq, db=database, evalue=evalue, outfmt=5, megablast=megablast,
                                      out=blastoutput_custom)
            elif blasttype == "blastp":
                NcbiblastpCommandline(query=seq_record.seq, db=database, evalue=evalue, outfmt=5, megablast=megablast,
                                      out=blastoutput_custom)
            elif blasttype == "blastx":
                NcbiblastxCommandline(query=seq_record.seq, db=database, evalue=evalue, outfmt=5, megablast=megablast,
                                      out=blastoutput_custom)
            elif blasttype == "tblastx":
                NcbitblastxCommandline(query=seq_record.seq, db=database, evalue=evalue, outfmt=5, megablast=megablast,
                                       out=blastoutput_custom)
            elif blasttype == "tblastn":
                NcbitblastnCommandline(query=seq_record.seq, db=database, evalue=evalue, outfmt=5, megablast=megablast,
                                       out=blastoutput_custom)
            else:
                raise Exception("Invalid blast choice!")
        else:
            with open(blastoutput_custom, "w+") as fxml:
                fxml.write(NCBIWWW.qblast(program=blasttype, database=database, sequence=seq_record.seq,
                                          entrez_query=(target_species + "[ORGN]"), expect=evalue, megablast=megablast))
        # if recblast_on:
            # recblastout.close()
            # return length_list


def recblast(seqfile, target_species, database1, filetype="fasta", query_species="Homo sapiens", blasttype='blastn',
             localblast1=True, localblast2=False, database2="RefSeq_Genes", evalue=0.001, identitiesperc=0.75,
             scoreperc=0.75, lengthperc=0.75, idthres='Score', megablast=True):
    """By Reciprocal BLAST, finds orthologs in Species 2 of a list of genes from Species 1 and annotates them.

    Reciprocal BLAST involves using a primary BLAST to identify putative orthologs in the "target_species" using
     sequences from the "query_species", which by default is "Homo sapiens".
    Input is a list of genes ("seqfile") saved as a specified "filetype" (defaults to FASTA), to be searched against an
     indicated database. Other options include:
    blasttype -- BLAST program to be used. ("blastn", "blastp", "blastx", "tblastx", "tblastn")
    localblast1 -- Should the Forward BLAST be done locally or at NCBI? (default True)
    localblast2 -- Should the Reverse BLAST be done locally or at NCBI? (default False)
    database2 -- Database to be queried for the Reverse Blast to ID putative orthologs. (default "RefSeq_Genes")
    evalue -- Maximum E-Value accepted from HSPs
    identitiesperc -- Minimum percent identity accepted from HSPs
    scoreperc -- Minimum percentage from the top score that will be used as a cut-off for putative orthologs.
    lenghtperc -- Minimum fraction of the total length of the alignment that will be accepted.
    idthres -- TODO: clarify
    """

    import os
    from Bio import SeqIO
    from Bio.Blast import NCBIXML

    print("Now starting...")

    # First loop will iterate over each sequence in a file, preferably FASTA but also allows for GenBank
    for index, seq_record in enumerate(SeqIO.parse(seqfile, filetype)):
        print("{}: {}".format(index, seq_record.name))
        length = len(seq_record)    # For use in calculating the length percentage of every HSP

        forward_blast_output = os.path.join(os.getcwd(),
                                            "{0}_recblast_out".format(target_species),
                                            "{0}_{1}_{2}_to_{3}.xml".format(blasttype, seq_record.name, query_species,
                                                                            target_species))
        forward_id_score_output = os.path.join(os.getcwd(),
                                       "{0}_recblast_out".format(target_species),
                                       "{0}_{1}_{2}_to_{3}.ID_Scores.tmp".format(blasttype, seq_record.name,
                                                                                 query_species, target_species))
        recblast_output = os.path.join(os.getcwd(),
                                    "{0}_recblast_out".format(target_species),
                                    "{0}_{1}.fasta".format(blasttype, seq_record.name))
        # Forward Blast:
        blast(seqfile=seqfile, target_species=target_species, database=database1, query_species="Homo sapiens",
              filetype=filetype, blasttype=blasttype, localblast=localblast1, evalue=evalue, recblast_on=True,
              megablast=megablast, blastoutput_custom=forward_blast_output)

        # Easy part's over - now we need to get the top hits from the forward BLAST, ID them, then compile a new
        # FASTA file with sequences from Species 2 that will be annotated via the Reverse BLAST against Species 1.

        # First we load the primary BLAST XML results to a handle, read the file, then loop over all alignments
        # to get the top scoring HSPs for each (I don't trust NCBI to always give me a pre-sorted list beforehand).
        with open(forward_blast_output, "r") as forward_blasthits:
            blastrecord = NCBIXML.read(forward_blasthits)
            align_scorelist = []
            for alignment in blastrecord.alignments:
                hsp_scorelist = sorted([hsp.score for hsp in alignment.hsps])
                align_scorelist.append(hsp_scorelist[0])
        # Two parts to this next loop: first we loop for each alignment. Next, we look though the HSPs in each
        # alignment file. If the HSP being considered has a score above the thresholds, we note down the ID and
        # score of that HSP and corresponding alignment; once we do that for one HSP in the series, we update the
        # "blast_has_run" variable and proceed to skip to the next alignment result. This goes on until all
        # alignments have been considered, and so we now have a complete list of putative orthologs.
            with open(forward_id_score_output, "w") as f_id_out:
                for align_index, alignment in enumerate(blastrecord.alignments):
                    blast_has_run = False  # Every time we consider a new alignment we
                    for hsp in alignment.hsps:
                        if blast_has_run:
                            break
                        if (hsp.score >= (scoreperc * align_scorelist[align_index]) and hsp.expect <= evalue and
                                hsp.align_length >= (length * lengthperc)
                                and (int(hsp.identities)/int(hsp.align_length)) >= identitiesperc):
                            f_id_out.write('{0}\nSCORE:{1};\n'.format(alignment.title,hsp.score))
                            blast_has_run = True
                        else:
                            continue
        # Now, equiped with the list of hits, we need to look these up on a database and get their sequences as a
        # FASTA file.
        with open(recblast_output, "w+") as blastout: # TODO: WRITE CODE TO GET SEQUENCES FROM LIST OF ID'S AND SAVE AS FASTA.
            pass # TODO: Remove this

        # Now that we have the sequences we can do the Reverse BLAST:
        # Big caveat though: we need to do each target individually... FIXME: DO SOMETHING ABOUT THIS.
        seq_number = 0 # simple counter to figure out how many sequences I have
        for entry_index, entry_record in enumerate(SeqIO.parse(recblast_output,"fasta")):
            reverse_blast_output = os.path.join(os.getcwd(),
                                                "{0}_recblast_out".format(target_species),
                                                "{0}_{1}_{3}_to_{2}_{4}.xml".format(blasttype, seq_record.name,
                                                                                   query_species, target_species,
                                                                                   entry_index)) # FIXME: entry_index is indexing every entry in recblast_output, rather than indexing only important hits. does this matter?
            reverse_id_score_output = os.path.join(os.getcwd(),
                                                   "{0}_recblast_out".format(target_species),
                                                   "{0}_{1}_{3}_to_{2}.ID_Scores.tmp".format(blasttype,
                                                                                                 seq_record.name,
                                                                                                 query_species,
                                                                                                 target_species)) # FIXME: this is definitely broken, need to figure out how to keep track of id_score without overwriting every loop.
            blast(seqfile=entry_record, target_species=query_species, database=database2,
                  query_species=target_species, filetype=filetype, blasttype=blasttype, localblast=localblast1,
                  evalue=evalue, recblast_on=True, megablast=megablast, blastoutput_custom=reverse_blast_output)
            with open(reverse_blast_output, "r") as reverse_blast_hits:
                blastrecord2 = NCBIXML.read(reverse_blast_hits)
                align_scorelist = []
                for alignment in blastrecord2.alignments:
                    hsp_scorelist = sorted([hsp.score for hsp in alignment.hsps])
                    align_scorelist.append(hsp_scorelist[0])
                with open(reverse_id_score_output, "w") as r_id_out:
                    for align_index, alignment in enumerate(blastrecord2.alignments):
                        blast_has_run = False
                        for hsp in alignment.hsps:
                            if blast_has_run:
                                break
                            if (hsp.score >= (scoreperc * align_scorelist[align_index]) and hsp.expect <= evalue and
                                    hsp.align_length >= (length * lengthperc)
                                    and (int(hsp.identities)/int(hsp.align_length)) >= identitiesperc):
                                r_id_out.write('{0}\nSCORE:{1};\n'.format(alignment.title,hsp.score))
                                blast_has_run = True
                            else:
                                continue

            # Now we can finally use the

        """ Old stuff, just ignore this stuff down here for now.

            for recblastrecord in NCBIStandalone.Iterator(reverse_blast_hits, blast_error_parser2):
                hit_names = []
                if str(idthres).isnumeric():
                    it_idthres = 0
                    while idthres > it_idthres:
                        hit_names.append(recblastrecord.alignments[idthres].title)
                        idthres -= 1
                else:
                    for alignment2 in recblastrecord.alignments:
                        hit_names = []
                        scorelist2 = []
                        for hsp2 in alignment2.hsps:
                            scorelist2.append(hsp2.score)
                        scorelist2.sort()
                        topscore2 = scorelist2[0]
                        for hsp2 in alignment2.hsps:
                            if (hsp2.score >= scoreperc * topscore2 and hsp2.expect <= evalue and
                                    hsp2.align_length >= length * lengthperc and
                                    hsp2.identities >= identitiesperc):
                                hit_names.append(alignment2.title)
            reverse_blast_hits.close()
            hit_title = '{0} {1}'.format(alignment.title, hit_names)
            blastout.write('****Alignment****\n')
            blastout.write('sequence: {}\n'.format(hit_title))
            blastout.write('length: {}\n'.format(alignment.length))
            blastout.write('e value: {}\n'.format(hsp.expect))
            blastout.write(hsp.query)
            blastout.write(hsp.match)
            blastout.write(hsp.sbjct)
            """

"""
recblast2 uses the experimental Biopython SearchIO module. Honestly it looks so much cleaner in the docs so I kinda
jumped ship to fix recblast without the whole "for x in y / for y in z / for z in a / ..." hell.
"""


def recblast2(seqfile, target_species, database1, filetype="fasta", query_species="Homo sapiens", blasttype='blastp',
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
        recblast_out = os.path.join(os.getcwd(), "{}_recblast_out".format(target_species), "{0}_{1}.fasta".format(blasttype,seq_record.name))
        errorfile = open(os.path.join(os.getcwd(), "recblast_out", "{}_recblasterror1.log".format(seqfile)),"w+")
        errorfile2 = open(os.path.join(os.getcwd(), "recblast_out", "{}_recblasterror2.log".format(seqfile)), "w+")

        # Begin by opening recblast_out, and then start with the primary BLAST
        with open(recblast_out, "w+") as blastout:
            length = len(seq_record)
            if localblast1:
                if blasttype == "blastn":
                    NcbiblastnCommandline(query=seq_record.seq,db=database1,evalue=evalue,outfmt=5,
                                          out="{0}to{1}blast_{2}.xml".format(query_species, target_species,
                                                                             seq_record.name))
                elif blasttype == "blastp":
                    NcbiblastpCommandline(query=seq_record.seq,db=database1,evalue=evalue,outfmt=5,
                                          out="{0}to{1}blast_{2}.xml".format(query_species, target_species,
                                                                             seq_record.name))
                elif blasttype == "blastx":
                    NcbiblastxCommandline(query=seq_record.seq,db=database1,evalue=evalue,outfmt=5,
                                          out="{0}to{1}blast_{2}.xml".format(query_species, target_species,
                                                                             seq_record.name))
                elif blasttype == "tblastx":
                    NcbitblastxCommandline(query=seq_record.seq,db=database1,evalue=evalue,outfmt=5,
                                           out="{0}to{1}blast_{2}.xml".format(query_species, target_species,
                                                                              seq_record.name))
                elif blasttype == "tblastn":
                    NcbitblastnCommandline(query=seq_record.seq,db=database1,evalue=evalue,outfmt=5,
                                           out="{0}to{1}blast_{2}.xml".format(query_species, target_species,
                                                                              seq_record.name))
                else:
                    raise Exception("Invalid blast choice!")
            else:
                with open(os.path.join(os.getcwd(), "recblast_out", "{0}to{1}blast_{2}.xml".format(query_species, target_species,
                                        seq_record.name)), "w+") as fxml:
                    fxml.write(NCBIWWW.qblast(program=blasttype,database=database1,sequence=seq_record.seq,
                                              entrez_query=(target_species + "[ORGN]"),expect=evalue))

            # Start Reciprocal Blast
            # First load the primary BLAST hits to a handle.
            blasthits = open("{0}to{1}blast_{2}.xml".format(query_species, target_species, seq_record.name),"r")
            # Next, load up the error parsers
            blast_error_parser = NCBIStandalone.BlastErrorParser(errorfile)
            blast_error_parser2 = NCBIStandalone.BlastErrorParser(errorfile2)

            # Iterate through all the blast hits


""""
SirBlastALot() is a function that, for a given sequence, BLASTS it against a list of organisms of one's choice; the
default list has an extensive vertebrate coverage. It also has an option to use BLAT instead. By default it just does
does a unidirectional blast of your sequence to each organism; you can set it to do a Reciprocal Blast as well.
"""
