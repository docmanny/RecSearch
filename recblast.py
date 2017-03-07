from misc_code import *  # Not good practice, I know, but necessary as a stand-in until I sort things out


def recblast(seqfile, target_species, fw_blast_db='chromosome', infile_type="fasta", output_type="fasta",
             query_species="Homo sapiens", blast_type='blastn', local_blast_1=False, local_blast_2=False,
             rv_blast_db="nt", expect=10, perc_score=0.50, perc_ident=50, perc_length=0.5,
             megablast=True, email='', id_type='brute', fw_source="psql", fw_id_db="", batch_size=50,
             passwd='', fw_id_db_version='1.0', verbose=True, BLASTDB='/usr/db/blastdb', **kwargs):
    """By Reciprocal BLAST, finds orthologs in Species 2 of a list of genes from Species 1 and annotates them.
    Reciprocal BLAST involves using a primary BLAST to identify putative orthologs in the "target_species" using
     sequences from the "query_species", which by default is "Homo sapiens".
    Input is a list of genes ("seqfile") saved as a specified "infile_type" (defaults to FASTA), to be searched against
    an indicated database. Other options include:
    blast_type -- BLAST program to be used. ("blastn", "blastp", "blastx", "tblastx", "tblastn")
    local_blast_1 -- Should the Forward BLAST be done locally or at NCBI? (default True)
    local_blast_2 -- Should the Reverse BLAST be done locally or at NCBI? (default False)
    rv_blast_db -- Database to be queried for the Reverse Blast to ID putative orthologs. (default "RefSeq_Genes")
    expect -- Maximum E-Value accepted from HSPs
    identitiesperc -- Minimum percent identity accepted from HSPs
    perc_score -- Minimum percentage from the top score that will be used as a cut-off for putative orthologs.
    lenghtperc -- Minimum fraction of the total length of the alignment that will be accepted.
    idthres -- TODO: clarify
    """

    from pathlib import Path
    from Bio import SeqIO, __version__
    from Bio.Blast import NCBIXML
    from operator import itemgetter

    if verbose:
        print("Now starting RecBlast...")
        print('BioPython Version: ', __version__)
    if isinstance(seqfile, list):
        seq_gen = ((index, seq_record) for (index, seq_record) in enumerate(seqfile))
    else:
        seqfile_path = Path(seqfile)
        if seqfile_path.exists() and seqfile_path.is_file():
            seq_gen = ((index, seq_record) for (index, seq_record) in enumerate(SeqIO.parse(str(seqfile_path),
                                                                                            infile_type)))
        else:
            raise FileNotFoundError

    # First loop will iterate over each sequence in a file, preferably FASTA but also allows for GenBank
    for index, seq_record in seq_gen:
        if verbose:
            print("Forward BLAST - {}: {}".format(index + 1, seq_record.name))
        forward_blast_output = Path("{0}_recblast_out".format(target_species).replace(' ', '_') + '/' +
                                    "{0}_{1}_tmp".format(blast_type, seq_record.name).replace(' ', '_') + '/' +
                                    "{0}_{1}_{2}_to_{3}.xml".format(blast_type, seq_record.name, query_species,
                                                                    target_species).replace(' ', '_'))

        forward_id_score_output = Path("{0}_recblast_out".format(target_species).replace(' ', '_') + '/' +
                                       "{0}_{1}_tmp".format(blast_type, seq_record.name).replace(' ', '_') + '/' +
                                       "{0}_{1}_{2}_to_{3}.ID_Scores.tmp".format(blast_type, seq_record.name,
                                                                                 query_species,
                                                                                 target_species).replace(' ', '_'))

        recblast_output_unanno = Path("{0}_recblast_out".format(target_species).replace(' ', '_') + '/' +
                                      "{0}_{1}_tmp".format(blast_type, seq_record.name).replace(' ', '_') + '/' +
                                      "unannotated_{0}_{1}.tmp".format(blast_type, seq_record.name).replace(' ', '_'))

        try:
            forward_blast_output.absolute().parent.mkdir(parents=True)
        except FileExistsError:
            pass
        try:
            forward_id_score_output.absolute().parent.mkdir(parents=True)
        except FileExistsError:
            pass
        try:
            recblast_output_unanno.absolute().parent.mkdir(parents=True)
        except FileExistsError:
            pass

        # Forward Blast:
        if fw_blast_db == 'skip':
            if verbose:
                print("Skipping Forward Blast!")
            pass
        else:
            blast(seq_record=seq_record, target_species=target_species, database=fw_blast_db,
                  query_species=query_species, filetype=infile_type, blast_type=blast_type, local_blast=local_blast_1,
                  expect=expect, megablast=megablast, blastoutput_custom=str(forward_blast_output), write=True,
                  perc_ident=perc_ident, verbose=verbose, BLASTDB=BLASTDB, **kwargs)
            if verbose:
                print('Forward blast done!')
        # Easy part's over - now we need to get the top hits from the forward BLAST, ID them, then compile a new
        # FASTA file with sequences from Species 2 that will be annotated via the Reverse BLAST against Species 1.

        # First we load the primary BLAST XML results to a handle, read the file, then loop over all alignments
        # to get the top scoring HSPs for each (I don't trust NCBI to always give me a pre-sorted list beforehand).
        # In addition, to really get to the crux of what this script should be doing, I also need to get the query
        # start and end points for each HSP, to tile them over the query, in order to get the true query coverage.
        # Furthermore I need to do the same for subject start and end so I can get a specific subrange for the sequence.
        with forward_blast_output.open("r") as forward_blasthits:
            if verbose:
                print('Opening Forward blast output: ', str(forward_blast_output.absolute()))
            blastrecord = NCBIXML.read(forward_blasthits)
        align_scorelist = []
        hsp_scorelist = []
        subject_range = []
        query_start_end = []
        for alignment in blastrecord.alignments:
            if verbose:
                print('Sorting through alignment\'s HSPs to get top scores of all alignments...')
            subject_range_hsp = []
            query_start_end_hsp = []
            for hsp in alignment.hsps:
                hsp_scorelist.append(hsp.score)
                subject_range_hsp.append(hsp.sbjct_start)
                subject_range_hsp.append(hsp.sbjct_end)
                query_start_end_hsp.append((hsp.query_start, hsp.query_end))
            hsp_scorelist.sort(reverse=True)
            query_start_end.append(i for i in merge_ranges(query_start_end_hsp))
            subject_range.append((subject_range_hsp[0], subject_range_hsp[-1]))
            if verbose:
                print("HSP Score List: \n\t", hsp_scorelist)
            align_scorelist.append(hsp_scorelist[0])
            if verbose:
                print("Alignment Score List: \n\t", align_scorelist)
        if verbose:
            print('Done with sorting!')
        # Two parts to this next loop: first we loop for each alignment. Next, we look though the HSPs in each
        # alignment file. If the HSP being considered has a score above the thresholds, we note down the ID and
        # score of that HSP and corresponding alignment; once we do that for one HSP in the series, we update the
        # "blast_got_hit" variable and proceed to skip to the next alignment result. This goes on until all
        # alignments have been considered, and so we now have a complete list of putative orthologs.
        with forward_id_score_output.open("w") as f_id_out:
            if verbose:
                print('Searching through alignments to get top-scoring hit IDs')
            has_written = False
            for align_index, alignment in enumerate(blastrecord.alignments):
                blast_got_hit = False  # Every time we consider a new alignment
                for hsp in alignment.hsps:
                    if blast_got_hit:
                        break
                    if ((hsp.score >= (perc_score * align_scorelist[align_index])) and (hsp.expect <= expect) and
                            (sum([i[-1] - i[0] for i in query_start_end[align_index]]) / blastrecord.query_length
                                 >= perc_length)):
                        if verbose:
                            print('Found annotation above threshold!')
                        f_id_out.write('{0}\t{1}\t{2}\n'.format(alignment.title.replace('/t', ' '),
                                                                ':{0}-{1}'.format(subject_range[align_index][0],
                                                                                  subject_range[align_index][-1]),
                                                                hsp.score))
                        has_written = True
                        blast_got_hit = True
                    else:
                        continue
                if not blast_got_hit:
                    print('NOTE: FOR ALIGNMENT {}, NO HITS WERE FOUND!'.format(alignment.title))
            if not has_written:
                print('WARNING! FOR THIS RUN, NO HITS WERE WRITTEN TO FILE, CONTINUING TO NEXT SEQUENCE IN LIST!')
                continue
        # Now, equiped with the list of hits, we need to look these up on a database and get their sequences as a
        # FASTA file.
        if verbose:
            print('Fetching sequences for ID\'ed hits...')
        try:
            fetchseq(id_file=str(forward_id_score_output), species=target_species, email=email, source=fw_source,
                     output_type=output_type, output_name=str(recblast_output_unanno), db=fw_id_db, delim='\t',
                     id_type=id_type, batch_size=batch_size, passwd=passwd, version=fw_id_db_version, verbose=verbose)
            if verbose:
                print('Done with fetching!')
        except IndexError:
            print('WARNING! FETCHSEQ FAILED! SKIPPING THIS SEQUENCE!')
            continue
        # Little caveat: fetchseq by design appends a .[output_type] to the end of the file so we need to add that on:
        recblast_output_unanno = str(recblast_output_unanno) + '.{}'.format(output_type)
        # Now that we have the sequences we can do the Reverse BLAST:
        # Big caveat though: we need to do each target individually.
        if verbose:
            print('Preparing for Reverse BLAST...')
        recblast_output = Path("{0}_recblast_out".format(target_species).replace(' ', '_') + '/' +
                               "{0}_{1}.{2}".format(blast_type, seq_record.name, output_type).replace(' ', '_'))
        try:
            recblast_output.absolute().parent.mkdir(parents=True)
        except FileExistsError:
            pass
        for entry_index, entry_record in enumerate(SeqIO.parse(str(recblast_output_unanno), "fasta")):
            if entry_record.seq:
                pass
            else:
                print(Warning('Entry {0} in unnanotated recblast file {1} came '
                              'back empty'.format(entry_record.name,
                                                  str(recblast_output_unanno))))
                continue
            if verbose:
                print("Entry #{} in unannotated RecBlast Hits:\n".format(entry_index))
                for item in [entry_record.id, entry_record.description, entry_record.seq]:
                    print('\t', item)
            reverse_blast_output = Path("{0}_recblast_out".format(target_species).replace(' ', '_') + '/' +
                                        "{0}_{1}_tmp".format(blast_type, seq_record.name).replace(' ', '_') + '/' +
                                        "{0}_{1}_{3}_to_{2}_{4}.xml".format(blast_type, seq_record.name,
                                                                            query_species, target_species,
                                                                            entry_index).replace(' ', '_'))
            try:
                reverse_blast_output.absolute().parent.mkdir(parents=True)
            except FileExistsError:
                pass
            if verbose:
                print('Performing Reverse Blast:')
            if rv_blast_db == 'skip':
                pass
            elif rv_blast_db == 'stop':
                print('Not performing reverse blast!')
                continue
            else:
                blast(seq_record=entry_record, target_species=query_species, database=rv_blast_db,
                      query_species=target_species, filetype=infile_type, blast_type=blast_type,
                      local_blast=local_blast_2, write=True,
                      expect=expect, megablast=megablast, blastoutput_custom=str(reverse_blast_output),
                      perc_ident=perc_ident, BLASTDB=BLASTDB, **kwargs)
            if verbose:
                print('Done with Reverse Blast!')
            with reverse_blast_output.open("r") as reverse_blast_hits:
                if verbose:
                    print('Getting top scores for each alignment...')
                blastrecord2 = NCBIXML.read(reverse_blast_hits)
                align_scorelist2 = []
                hsp_scorelist2 = []
                subject_range2 = []
                query_start_end2 = []
                for alignment2 in blastrecord2.alignments:
                    if verbose > 4:
                        print('Sorting through alignment\'s HSPs to get top scores of all alignments...')
                    subject_range_hsp2 = []
                    query_start_end_hsp2 = []
                    for hsp2 in alignment2.hsps:
                        hsp_scorelist2.append(hsp2.score)
                        subject_range_hsp2.append(hsp2.sbjct_start)
                        subject_range_hsp2.append(hsp2.sbjct_end)
                        query_start_end_hsp2.append((hsp2.query_start, hsp2.query_end))
                    hsp_scorelist2.sort(reverse=True)
                    query_start_end2.append([i for i in merge_ranges(query_start_end_hsp2)])
                    subject_range2.append((subject_range_hsp2[0], subject_range_hsp2[-1]))
                    if verbose > 4:
                        print("HSP Score List: \n\t", hsp_scorelist2)
                    align_scorelist2.append(hsp_scorelist2[0])
                    if verbose > 4:
                        print("Alignment Score List: \n\t", align_scorelist2)
                        print("Query_start_end: \n\t", query_start_end2)
                        print("Subject Range: \n\t", subject_range2)
                if verbose:
                    print('Done with sorting!')
                # Now we have a list of the top score of each alignment for the current entry_record.
            with recblast_output.open("w+") as rb_out:
                if verbose:
                    print('Annotating BLAST results')
                has_written2 = False
                reverse_blast_annotations = list()
                for align_index2, alignment2 in enumerate(blastrecord2.alignments):
                    blast_got_hit2 = False
                    for hsp2 in alignment2.hsps:
                        if (hsp2.score >= (perc_score * align_scorelist2[align_index2])):
                            if verbose > 4:
                                print('hsp score above threshold')
                            if (hsp2.expect <= expect):
                                if verbose > 4:
                                    print('hsp expect below threshold')
                                if verbose > 4:
                                    print('HSP Length: ', query_start_end2[align_index2])
                                length_alignment = sum([i[-1] - i[0] for i in query_start_end2[align_index2]])
                                align_len_threshold = blastrecord2.query_length * perc_length
                                if verbose > 4:
                                    print(length_alignment)
                                    print(align_len_threshold)
                                if length_alignment >= align_len_threshold:
                                    print('hsp perc length above threshold')
                                    if verbose:
                                        print('Found hit!')
                                    reverse_blast_annotations.append((alignment2.title, '[:{0}-{1}]'.format(
                                        subject_range2[align_index2][0],
                                        subject_range2[align_index2][1]),
                                                                      hsp2.score))
                                    has_written2 = True
                                    blast_got_hit2 = True
                                else:
                                    print('WARNING HSP LENGTH BELOW THRESHOLD')
                                    print(length_alignment,
                                          ' not greater than ', align_len_threshold)
                            else:
                                print('WARNING HSP EXPECT ABOVE THRESHOLD')
                                print(hsp2.expect, 'not less than', expect)
                        else:
                            print('WARNING HSP SCORE BELOW THRESHOLD')
                            print(hsp2.score, ' not greater than ', (perc_score * align_scorelist2[align_index2]))
                            # else:
                            # continue
                    if not blast_got_hit2:
                        print('NOTE: Alignment {} was not used to annotate.'.format(alignment2.title))
                if reverse_blast_annotations:
                    sorted(reverse_blast_annotations, reverse=True, key=itemgetter(2))
                    annotations = ['\t||{0} {1} ({2})'.format(anno[0], anno[1], anno[2]) for anno
                                   in reverse_blast_annotations]
                    if verbose:
                        print('********************************************')
                        print(annotations)
                        print('********************************************')
                    entry_record.description += ''.join(annotations)
                    if verbose > 3:
                        print(entry_record)
                    SeqIO.write(entry_record, rb_out, output_type)
                if not has_written2:
                    print(Warning('NONE OF THE REVERSE BLAST HITS FOR THIS RUN MET ANNOTATION CRITERIA!'))
                    continue
        if verbose:
            print('DONE!!!!')
            # Done!