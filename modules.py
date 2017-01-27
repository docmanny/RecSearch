"""
Author: Juan Manuel Vazquez
Date Created: 11/01/2016
Date Last Modified: 01/19/2017
"""


def biosql_addrecord(sub_db_name, description, file, passwd, filetype='fasta', driver="psycopg2", user="postgres",
                     host="localhost", db="bioseqdb", verbose=True, pretend=False):  # TODO: FILL OUT DOCSTRING
    """Wrapper for adding records to a BioSQL database.

    :param sub_db_name: testing.
    :param description:
    :param file:
    :param passwd:
    :param filetype:
    :param driver:
    :param user:
    :param host:
    :param db:
    :param verbose:
    :param pretend:
    :return:
    """
    from Bio import SeqIO
    from BioSQL import BioSeqDatabase
    from pathlib import Path
    from sys import exc_info
    if verbose:
        print("Beginning addition of {0} to main db {1}".format(filetype, db))
        print("Opening BioSeqDB server...")
    try:
        server = BioSeqDatabase.open_database(driver=driver, user=user, passwd=passwd, host=host, db=db)
    except ImportError:
        if verbose:
            print("Import Error! The driver you selected isn't correct")
        raise
    except:
        if verbose:
            print("Oops! Something went wrong with opening the server! Are you use all the right statements?")
        raise
    else:
        if verbose:
            print("Database opened!")

    if verbose:
        print("Creating new sub-database for file...")
    try:
        if pretend:
            print('Pretend is active, this is where I would have committed the info to the server!')
        else:
            db = server.new_database(sub_db_name, description=description)
    except:
        if verbose:
            print('Failed to create new server!')
        raise
    else:
        if verbose:
            print("Successfully generated new sub-database {0}!".format(sub_db_name))

    if verbose:
        print("Committing sub-database to server...")
    try:
        if pretend:
            print('Pretend is active, this is where I would have committed the new sub-database to the server!')
        else:
            server.commit()
    except:
        if verbose:
            print('Couldn\'t commit new database!')
        raise
    else:
        if verbose:
            print("Sub-database successfully committed!")

    if verbose:
        print("Parsing file now for entry into {}... (this takes a while)".format(sub_db_name))
    infile = Path(file)
    try:
        if infile.exists() and infile.is_file():
            try:
                if pretend:
                    print('Pretend is active, this is where I would have tried to load the data!')
                else:
                    count = db.load(SeqIO.parse(str(infile), filetype))
            except:
                if verbose:
                    print("Problem loading data!")
                raise
            else:
                if pretend:
                    print('Pretend is active, this is where I would have said that records were loaded!')
                else:
                    if verbose:
                        print("Loaded {} records".format(count))
            if verbose:
                print("Commiting new data to db {}".format(sub_db_name))
            try:
                if pretend:
                    print('Pretend is active, this is where I would have committed the info to the server!')
                else:
                    server.commit()
            except:
                if verbose:
                    print('Couldn\'t commit new database!')
                raise
            else:
                if verbose:
                    print("Sub-database successfully committed!")
        else:
            print('Sorry, file {} does not seem to exist...'.format(infile))
    except:
        print('Whoops! Something happened trying to open file {}:'.format(infile), exc_info())
    # End of Function


def biosql_recordids(sub_db_name, passwd, dumpfile=True, driver="psycopg2", user="postgres", host="localhost",
                     db="bioseqdb"):  # TODO: FILL OUT DOCSTRING
    """

    :param sub_db_name:
    :param passwd:
    :param dumpfile:
    :param driver:
    :param user:
    :param host:
    :param db:
    :return:
    """
    from BioSQL import BioSeqDatabase

    server = BioSeqDatabase.open_database(driver=driver, user=user, passwd=passwd, host=host, db=db)
    db = server[sub_db_name]
    if dumpfile:
        import sys
        og_stdout = sys.stdout
        outfile = open('./'+sub_db_name+'.index', "w")
        sys.stdout = outfile
    print("This database contains {} records".format(len(db)))
    for key, record in db.items():
        if key == 1:
            print("As an example of what items are in the database, here is the first entry:")
            print(key, record)
        print("Key {0} maps to a sequence record with id {0}".format(key, record.id))
    if dumpfile:
        sys.stdout = og_stdout
        outfile.close()
    # End of Function


def biosql_getrecord(sub_db_name, passwd, id_list=[], id_type='accession', driver="psycopg2", user="postgres",
                     host="localhost", db="bioseqdb", verbose=False):  # TODO: FILL OUT DOCSTRING
    """

    :param sub_db_name:
    :param passwd:
    :param id_list:
    :param id_type:
    :param driver:
    :param user:
    :param host:
    :param db:
    :param verbose:
    :return:
    """
    from BioSQL import BioSeqDatabase
    server = BioSeqDatabase.open_database(driver=driver, user=user, passwd=passwd, host=host, db=db)
    db = server[sub_db_name]
    seqdict = {}
    for identifier in id_list:
        try:
            seqdict[identifier] = db.lookup(**{id_type: identifier})
        except IndexError as err:
            if verbose:
                print('Attempting to search using Primary ID instead of declared type:')
            try:
                seqdict[identifier] = db.lookup(primary_id=identifier)
            except IndexError:
                if verbose:
                    print('Primary ID search didn\'t work! Troubleshooting...')
                try:
                    if verbose:
                        print('Attempting to search using name instead of declared type:')
                    seqdict[identifier] = db.lookup(name=identifier)
                except IndexError:
                    if verbose:
                        print('Name search didn\'t work! Troubleshooting...')
                    try:
                        id_type = input('Last shot, chose an ID type: '
                                        '[accession, primary_id, gi, version, display_id, name]')
                        seqdict[identifier] = db.lookup(**{id_type: identifier})
                    except IndexError:
                        print("WARNING: couldn't find {0}: \n full error: {1}".format(identifier, err))

                finally:
                    if verbose:
                        print('Error resolved! Got sequences')
            finally:
                if verbose:
                    print('Error resolved! Got sequences')
        finally:
            if verbose:
                print('Error resolved! Got sequences')
    return seqdict
    # End of Function


def biosql_addmultirecord(base_dir, passwd='', description_base='Record imported from the file', filetype='fasta',
                          driver="psycopg2", user="postgres", host="localhost", db="bioseqdb", verbose=True,
                          pretend=False):
    """
    Wrapper for adding FASTA files in a directory hierarchy to a BioSQL database in a automated fashion.
    :param base_dir: The base directory from which the function will begin its search. It assumes that this location
                        contains a series of sub-directories, and that immediately in those directories lie the files
                        to be added as records into the database.
    :param description_base: A basic description of the data to be added. This descriptor will be applied to all entries
                                to be added, with the addition of "[filename] located at [file location]" added on.
                                (Default: 'Record imported from the file'
    :param passwd: The password of the database to be accessed. (Default: '')
    :param filetype: The format of the files to be imported. (Default: 'fasta')
    :param driver: The python driver to be used to interact with the database. By default, it uses the psycopg2 module
                    to communicate with a PostgreSQL database, but in principle this could also be a MySQL or some other
                    SQL server type.
    :param user:
    :param host:
    :param db:
    :param verbose:
    :param pretend:
    :return:
    """
    from pathlib import Path
    from sys import exc_info
    if verbose:
        print('Changing directory to {}'.format(base_dir))
    p_base = Path(base_dir)
    for p_sub in [subpath for subpath in p_base.iterdir() if subpath.is_dir()]:
        if verbose:
            print('Found path {}, checking for files'.format(str(p_sub)))
        p_file = [p for p in sorted(p_sub.glob('*.{}*'.format(filetype[0:2]))) if p.is_file()]
        if len(p_file) > 1:
            checkyoself = True
            while checkyoself:
                for a, b in enumerate(p_file):
                    print('{0}: {1}'.format(a+1, b))
                p_file_index = input("Multiple files found for {}, please select one:".format(str(p_file[0].parent)))
                if p_file_index in [str(i) for i in range(1, len(p_file)+1)]:
                    checkyoself = False
                elif p_file_index in ['quit', 'q', 'e', 'exit']:
                    escape = input("Are you sure you want to stop? (y/n)")
                    if escape == 'y':
                        raise Exception("Program ended by User")
                else:
                    print('Invalid selection!')
        elif (len(p_file) == 0) and verbose:
            print('No file found, continuing...')
            continue
        else:
            p_file = p_file[0]
            if verbose:
                print("File found: {}".format(str(p_file.name)))
            p_file_lst = p_file.name.split('_')
            sub_db_name = ''.join([x[0:3].title() for x in p_file_lst[0:2]] +
                                  [p_file_lst[len(p_file_lst)-1].lstrip('v').rstrip(p_file.suffix)])
            description_full = description_base + str(p_file.name) + 'located at' + str(p_file.parent)
            try:
                biosql_addrecord(sub_db_name=sub_db_name, description=description_full, file=str(p_file),
                                 passwd=passwd, filetype=filetype, driver=driver, user=user,
                                 host=host, db=db, verbose=verbose, pretend=pretend)
            except:
                print("Unexpected error:", exc_info()[0])
                continue
    if verbose:
        print("Done!")


def fetchseq(id_file, species, email='', source="psql", output_type="fasta", output_name="outfile",
             db="nucleotide", delim='\t', id_type='accession', batch_size=50, passwd='', version='1.0', verbose=True):
    # Todo: Fill out docstring
    """

    :param id_file:
    :param species:
    :param email:
    :param source:
    :param output_type:
    :param output_name:
    :param db:
    :param delim:
    :param id_type:
    :param batch_size:
    :param passwd:
    :param version:
    :param verbose:
    :return:
    """
    import re
    from os import strerror
    from errno import ENOENT
    from Bio import SeqIO
    from pathlib import Path

    success_status = 1  # Lets downstream functions know if this worked 100% even though most errors will be caught
    in_file = Path(id_file)
    out_file = Path(output_name + '.' + output_type)
    if verbose:
        print("Loading ID File...")
    if in_file.exists():
        if verbose:
            print('ID File found successfully: ', str(in_file.absolute()))
    else:
        raise FileNotFoundError(ENOENT, strerror(ENOENT), str(in_file.name))

    # Read ID file to compile lists:
    with in_file.open('r') as infile_handle:
        if verbose:
            print('Reading ID File...')
        id_prelist = [line.strip() for line in infile_handle]     # list of each line in the file
    id_list = [str(item.split(delim)) for item in id_prelist]    # Breaks the tab sep in the lines into strings

    # Define the regex functions
    p = [re.compile('(gi)([| :_]+)(\d\d+\.?\d*)(.*)'),      # regex for gi
         re.compile('(NW|ref)([| _:]+)(\d\d+\.?\d*)(.*)'),  # regex for accession
         re.compile('(id)([| :_]+)(\d\d+\.?\d*)(.*)'),      # regex for generic ID
         re.compile(':(\d+)-(\d+)'),                             # regex for sequence range
         ]
    id_list_ids = []    # Initialized list of IDs
    seq_range = {}      # Initialized dict of sequence ranges

    # Begin search:
    if verbose:
        print('ID File Loaded, performing regex search for identifiers...')
        print('ID Specified as: ', id_type)
    if id_type == 'brute':
        if bool(p[0].findall(id_list[0])):
            id_type = 'gi'
        if bool(p[1].findall(id_list[0])):
            id_type = 'accession'
        if bool(p[2].findall(id_list[0])):
            id_type = 'id'
        else:
            id_type = 'other'
        if verbose:
            print('Brute Force was set, tested strings for all pre-registered IDs. ID was selected as type ', id_type)
    if id_type == 'gi':
        if bool(p[0].findall(id_list[0])):
            found_id = True
            if verbose:
                print('Successfully found GI numbers, list has been compiled!')
            for item in id_list:
                id_list_ids.append(''.join(p[0].findall(item)[0][0:3]))
                if bool(p[3].findall(id_list[0])):
                    # Seq_range will be a list of tuples where the second element is the range, and the first
                    # is the ID. This way, the function accommodates sequences with a subrange and sequences without a
                    # subrange.
                    seq_range[''.join(p[0].findall(item)[0][0:3])] = p[3].findall(item)[0]
                    if verbose:
                        print('Found sequence delimiters in IDs!')
        else:
            found_id = False
    elif id_type == 'accession':
        if bool(p[1].findall(id_list[0])):
            found_id = True
            if verbose:
                print('Successfully found accession numbers, list has been compiled!')
            for item in id_list:
                id_list_ids.append(''.join(p[1].findall(item)[0][0:3]))
                if bool(p[3].findall(id_list[0])):
                    seq_range[''.join(p[1].findall(item)[0][0:3])] = p[3].findall(item)[0]
                    if verbose:
                        print('Found sequence delimiters in IDs!')
        else:
            found_id = False
    elif id_type == 'id':
        if bool(p[2].findall(id_list[0])):
            found_id = True
            if verbose:
                print('Successfully found ID numbers, list has been compiled!')
            for item in id_list:
                id_list_ids.append(''.join(p[2].findall(item)[0][0:3]))
                if bool(p[3].findall(id_list[0])):
                    seq_range[''.join(p[2].findall(item)[0][0:3])] = p[3].findall(item)[0]
                    if verbose:
                        print('Found sequence delimiters in IDs!')
        else:
            found_id = False
    else:
        found_id = False
    while not found_id:
        print('Header identified for first sequence ID:', id_list[0], sep='\n')
        custom_regex = input('Couldn\'t find ID using preset patterns... Please enter ID pattern for regex search:')
        if custom_regex[0].lower() == 'q':
            exit()
        print('Will try again...')
        p.append(re.compile(custom_regex))
        if bool(p[4].findall(id_list[0])):
            id_type = input('ID name:')
            found_id = True
            if verbose:
                print('Successfully found custom ID numbers, list has been compiled!')
            for item in id_list:
                id_list_ids.append(''.join(p[4].findall(item)[0][0:3]))
                if bool(p[3].findall(str(item))):
                    seq_range[''.join(p[4].findall(item)[0][0:3])] = p[3].findall(item)[0]
                    if verbose:
                        print('Found sequence delimiters in IDs!')
        else:
            print('Sorry, still can\'t find it...')
    if verbose:
        print('ID list: ')
        for index, ID_item in enumerate(id_list_ids):
            print(index+1, ': ', ID_item)

    # Armed with the ID list, we fetch the sequences from the appropriate source
    if source == "Entrez":  # Todo: Make sure this will actually output the correct sequence range...
        if verbose:
            print('Source selected was Entrez. Beginning search now:')
        from Bio import Entrez
        from urllib.error import HTTPError
        from time import sleep
        Entrez.email = email
        if verbose:
            print('Entrez email set as: ', email)
        id_str = ",".join(id_list_ids)
        search_results = Entrez.read(Entrez.epost(db, id=id_str))
        if verbose:
            print('EPost with IDs for database {} submitted to Entrez'.format(db))
        webenv = search_results["WebEnv"]
        query_key = search_results["QueryKey"]
        with out_file.open("a+") as out_handle:
            if verbose:
                print('Opened outfile ', str(out_file.name))
                print('Commencing download:')
            for start in range(0, len(id_list_ids), batch_size):
                if verbose:
                    print('Fetching sequences {0}-{1}'.format(start, start + batch_size))
                attempt = 0
                while attempt < 3:
                    if verbose:
                        print('Attempt #', str(attempt+1))
                    attempt += 1
                    try:
                        fetch_handle = Entrez.efetch(db=db, rettype="fasta", retmode="text", retstart=start,
                                                     retmax=batch_size, webenv=webenv, query_key=query_key)
                    except HTTPError as err:
                        if 500 <= err.code <= 599:
                            print("Received error from server ", err)
                            print("Attempt {} of 3".format(attempt))
                            print('Will wait before next attempt...')
                            sleep(15)
                        else:
                            print('could\'t get sequences, omitting', id_list[start:start+batch_size])
                            success_status = 0
                            continue
                data = fetch_handle.read()
                fetch_handle.close()
                out_handle.write(data)
    elif source.lower() == "psql":
        sub_db_name = ''.join([i[0:3] for i in species.title().split(' ')]) + version
        seqdict = biosql_getrecord(sub_db_name=sub_db_name, id_list=id_list_ids, id_type=id_type, passwd=passwd,
                                   driver="psycopg2", user="postgres", host="localhost", db="bioseqdb")
        itemsnotfound = [x for x in id_list_ids if x not in seqdict.keys()]
        if itemsnotfound:
            if verbose:
                print('Some items were not found. List of items will be saved to the file items_not_found.output')
                for item in itemsnotfound:
                    print(item)
            with open(str(out_file.cwd()) + 'items_not_found.output', 'w') as missingitems:
                missingitems.writelines(itemsnotfound)
        keys = [k for k in seqdict.keys()]
        if verbose:
            print("Sequence Dictionary keys:")
            print(keys)
        if bool(seq_range):
            seqrange_ids = [ids for ids in seq_range.keys()]
            if verbose:
                print('Sequence Range IDs:')
                print(seqrange_ids)
            for k in keys:
                if seqdict[k].id in seqrange_ids:
                    if verbose:
                        print('For sequence {}, found a sequence range!'.format(str(seqdict[k].id)))
                        print('Full length of sequence: {}'.format(len(seqdict[k])))
                    if id_type == 'gi':
                        seq_description_full = p[0].findall(seqdict[k].description)[0]
                    elif id_type == 'accession':
                        seq_description_full = p[1].findall(seqdict[k].description)[0]
                    elif id_type == 'id':
                        seq_description_full = p[2].findall(seqdict[k].description)[0]
                    else:
                        seq_description_full = p[4].findall(seqdict[k].description)[0]
                else:
                    if verbose:
                        print('No sequence range found, continuing...')
                    continue
                id_range = ':' + '-'.join(seq_range[k])
                seqdict[k] = seqdict[k][int(seq_range[k][0]):int(seq_range[k][1])]
                if verbose:
                    print('Seq_description_full: ', seq_description_full)
                    print('id_range: ', id_range[1:])
                seqdict[k].description = ''.join(seq_description_full[0:3]) + id_range + str(seq_description_full[3])
                if verbose:
                    print(seqdict[k].description)
                seqdict[k].id += id_range
                if verbose:
                    print(seqdict[k].id)
        if verbose:
            print('Sequence Record post-processing, to be saved:')
            print(seqdict)
            print('Length of subsequence with range{0}: {1}'.format(id_range, len(seqdict[k])))
        SeqIO.write([seqdict[key] for key in seqdict.keys()], str(out_file), output_type)
    elif source == "fasta":  # Note: anecdotally, this doesn't run terribly fast - try to avoid.
        seqdict = SeqIO.index(db, source,
                              key_function=lambda identifier: p[0].search(p[2].search(identifier).group()).group())
        itemsnotfound = [x for x in id_list_ids if x not in seqdict.keys()]
        if itemsnotfound:
            if verbose:
                print('Some items were not found. List of items will be saved to the file items_not_found.output')
                for item in itemsnotfound:
                    print(item)
            with open(str(out_file.cwd()) + 'items_not_found.output', 'w') as missingitems:
                missingitems.writelines(itemsnotfound)
        keys = [k for k in seqdict.keys()]
        if verbose:
            print("Sequence Dictionary keys:")
            print(keys)
        if bool(seq_range):
            seqrange_ids = [ids for ids in seq_range.keys()]
            if verbose:
                print('Sequence Range IDs:')
                print(seqrange_ids)
            for k in keys:
                if seqdict[k].id in seqrange_ids:
                    if verbose:
                        print('For sequence {}, found a sequence range!'.format(str(seqdict[k].id)))
                        print('Full length of sequence: {}'.format(len(seqdict[k])))
                    if id_type == 'gi':
                        seq_description_full = p[0].findall(seqdict[k].description)[0]
                    elif id_type == 'accession':
                        seq_description_full = p[1].findall(seqdict[k].description)[0]
                    elif id_type == 'id':
                        seq_description_full = p[2].findall(seqdict[k].description)[0]
                    else:
                        seq_description_full = p[4].findall(seqdict[k].description)[0]
                else:
                    if verbose:
                        print('No sequence range found, continuing...')
                    continue
                id_range = ':' + '-'.join(seq_range[k])
                seqdict[k] = seqdict[k][int(seq_range[k][0]):int(seq_range[k][1])]
                seqdict[k].description = ''.join(seq_description_full[0:3]) + id_range + str(seq_description_full[3])
                seqdict[k].id += id_range
                if verbose:
                    print('Length of subsequence with range{0}: {1}'.format(id_range, len(seqdict[k])))
        SeqIO.write([seqdict[key] for key in seqdict.keys()], str(out_file), output_type)
    else:
        print('Not a valid database source!')
    if verbose:
        print('Done!')
    return success_status

def fetchseq_multi():
    # Function that either:
    # 1) Collects sequences for a single protein from various species, given a list of species;
    # 2) Collects various sequences from one species, given a list of sequence IDs;
    # 3) Collects various sequences, given sequence IDs.
    pass  # TODO: Write fetchseq_multi


def crosscheck():
    """Looks through a series of .fasta files, and checks which IDs in the headers are common between them.
    """
    # TODO: Write Crosscheck
    pass


def blast(seqfile, target_species, database, query_species="Homo sapiens", filetype="fasta", blasttype='blastn',
          localblast=True, evalue=0.005, recblast_on=False, megablast=True, blastoutput_custom="", verbose=True):
    # TODO: TEST IT 
    # TODO: WRITE DOCSTRING
    # TODO: AFTER RECBLAST WORKS, CLEAN UP COMMENTED-OUT CODE
    # TODO: Check if you can get query coverage info from blast file
    from pathlib import Path
    from Bio import SeqIO
    from Bio.Blast import NCBIWWW
    from Bio.Blast.Applications import NcbiblastxCommandline, NcbiblastnCommandline, NcbiblastpCommandline, \
        NcbitblastxCommandline, NcbitblastnCommandline

    if recblast_on:
        length_list = []
    if verbose:
        print("Now starting BLAST...")
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
        if verbose:
            print("Current BLAST: \n {}: {}".format(index+1, seq_record.name))

        # Begin by opening recblast_out, and then start with the primary BLAST
        if blastoutput_custom == '':
            blastoutput_custom = Path("{0}_blast".format(target_species),
                                      "{0}_{1}_{2}_to_{3}.xml".format(blasttype, seq_record.name,
                                                                      query_species, target_species))
        else:
            blastoutput_custom = Path(blastoutput_custom)
        recstring = str(seq_record.seq)
        if localblast:
            if blasttype == "blastn":
                NcbiblastnCommandline(query=recstring, db=database, evalue=evalue, outfmt=5, megablast=megablast,
                                      out=blastoutput_custom)
            elif blasttype == "blastp":
                NcbiblastpCommandline(query=recstring, db=database, evalue=evalue, outfmt=5, megablast=megablast,
                                      out=str(blastoutput_custom))
            elif blasttype == "blastx":
                NcbiblastxCommandline(query=recstring, db=database, evalue=evalue, outfmt=5, megablast=megablast,
                                      out=str(blastoutput_custom))
            elif blasttype == "tblastx":
                NcbitblastxCommandline(query=recstring, db=database, evalue=evalue, outfmt=5, megablast=megablast,
                                       out=str(blastoutput_custom))
            elif blasttype == "tblastn":
                NcbitblastnCommandline(query=recstring, db=database, evalue=evalue, outfmt=5, megablast=megablast,
                                       out=str(blastoutput_custom))
            else:
                raise Exception("Invalid blast choice!")
        else:
            blast_result = NCBIWWW.qblast(program=str(blasttype), database=str(database), sequence=recstring,
                                          entrez_query='"{}"[ORGN]'.format(target_species), expect=evalue,
                                          megablast=megablast)
            if verbose:
                print('Done, saving to outfile....')
            with blastoutput_custom.open("w") as fxml:
                fxml.write(blast_result.read())
        print("Done!")


        # if recblast_on:
            # recblastout.close()
            # return length_list


def recblast(seqfile, target_species, fw_blast_db, infile_type="fasta", output_type="fasta",
             query_species="Homo sapiens", blasttype='blastn', localblast1=True, localblast2=False,
             rv_blast_db="RefSeq_Genes", evalue=0.001, identitiesperc=0.75, scoreperc=0.75, lengthperc=0.75,
             megablast=True, email='', id_type='bruteforce', fw_source="psql", fw_id_db="", batch_size=50,
             passwd='', fw_id_db_version='1.0', rv_id_db="", rv_id_db_version='1.0', verbose=True):
    """By Reciprocal BLAST, finds orthologs in Species 2 of a list of genes from Species 1 and annotates them.

    Reciprocal BLAST involves using a primary BLAST to identify putative orthologs in the "target_species" using
     sequences from the "query_species", which by default is "Homo sapiens".
    Input is a list of genes ("seqfile") saved as a specified "infile_type" (defaults to FASTA), to be searched against an
     indicated database. Other options include:
    blasttype -- BLAST program to be used. ("blastn", "blastp", "blastx", "tblastx", "tblastn")
    localblast1 -- Should the Forward BLAST be done locally or at NCBI? (default True)
    localblast2 -- Should the Reverse BLAST be done locally or at NCBI? (default False)
    rv_blast_db -- Database to be queried for the Reverse Blast to ID putative orthologs. (default "RefSeq_Genes")
    evalue -- Maximum E-Value accepted from HSPs
    identitiesperc -- Minimum percent identity accepted from HSPs
    scoreperc -- Minimum percentage from the top score that will be used as a cut-off for putative orthologs.
    lenghtperc -- Minimum fraction of the total length of the alignment that will be accepted.
    idthres -- TODO: clarify
    """

    from pathlib import Path
    from Bio import SeqIO
    from Bio.Blast import NCBIXML

    print("Now starting...")

    # First loop will iterate over each sequence in a file, preferably FASTA but also allows for GenBank
    for index, seq_record in enumerate(SeqIO.parse(seqfile, infile_type)):
        print("{}: {}".format(index, seq_record.name))
        length = len(seq_record)    # For use in calculating the length percentage of every HSP

        forward_blast_output = Path("{0}_recblast_out".format(target_species) + '/' +
                                    "{0}_{1}_{2}_to_{3}.xml".format(blasttype, seq_record.name, query_species,
                                                                    target_species))
        forward_id_score_output = Path("{0}_recblast_out".format(target_species) + '/' +
                                       "{0}_{1}_{2}_to_{3}.ID_Scores.tmp".format(blasttype, seq_record.name,
                                                                                 query_species, target_species))
        recblast_output = Path("{0}_recblast_out".format(target_species) + '/' +
                               "{0}_{1}.fasta".format(blasttype, seq_record.name))
        # Forward Blast:
        blast(seqfile=seqfile, target_species=target_species, database=fw_blast_db, query_species="Homo sapiens",
              filetype=infile_type, blasttype=blasttype, localblast=localblast1, evalue=evalue, recblast_on=True,
              megablast=megablast, blastoutput_custom=str(forward_blast_output))

        # Easy part's over - now we need to get the top hits from the forward BLAST, ID them, then compile a new
        # FASTA file with sequences from Species 2 that will be annotated via the Reverse BLAST against Species 1.

        # First we load the primary BLAST XML results to a handle, read the file, then loop over all alignments
        # to get the top scoring HSPs for each (I don't trust NCBI to always give me a pre-sorted list beforehand).
        with forward_blast_output.open("r") as forward_blasthits:
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
            with forward_id_score_output.open("w") as f_id_out:
                for align_index, alignment in enumerate(blastrecord.alignments):
                    blast_has_run = False  # Every time we consider a new alignment we
                    for hsp in alignment.hsps:
                        if blast_has_run:
                            break
                        if (hsp.score >= (scoreperc * align_scorelist[align_index]) and hsp.expect <= evalue and
                                hsp.align_length >= (length * lengthperc)
                                and (int(hsp.identities)/int(hsp.align_length)) >= identitiesperc):
                            f_id_out.write('{0}\t{1}\n'.format(alignment.title.replace('/t',' '), hsp.score))
                            blast_has_run = True
                        else:
                            continue
        # Now, equiped with the list of hits, we need to look these up on a database and get their sequences as a
        # FASTA file.
        fetchseq(id_file=str(forward_id_score_output), species=target_species, email=email, source=fw_source,
                 output_type=output_type, output_name=str(recblast_output), db=fw_id_db, delim='\t', id_type=id_type,
                 batch_size=batch_size, passwd=passwd, version=fw_id_db_version, verbose=verbose)

        # Now that we have the sequences we can do the Reverse BLAST:
        # Big caveat though: we need to do each target individually...  FIXME: DO SOMETHING ABOUT THIS.

        seq_number = 0 # simple counter to figure out how many sequences I have
        for entry_index, entry_record in enumerate(SeqIO.parse(str(recblast_output),"fasta")):
            reverse_blast_output =Path("{0}_recblast_out".format(target_species),
                                       "{0}_{1}_{3}_to_{2}_{4}.xml".format(blasttype, seq_record.name,
                                                                           query_species, target_species, entry_index))  # FIXME: entry_index is indexing every entry in recblast_output, rather than indexing only important hits. does this matter?
            reverse_id_score_output = Path("{0}_recblast_out".format(target_species),
                                           "{0}_{1}_{3}_to_{2}.ID_Scores.tmp".format(blasttype, seq_record.name,
                                                                                     query_species, target_species))  # FIXME: this is definitely broken, need to figure out how to keep track of id_score without overwriting every loop.
            blast(seqfile=entry_record, target_species=query_species, database=rv_blast_db,
                  query_species=target_species, filetype=infile_type, blasttype=blasttype, localblast=localblast2,
                  evalue=evalue, recblast_on=True, megablast=megablast, blastoutput_custom=str(reverse_blast_output))
            with reverse_blast_output.open("r") as reverse_blast_hits:
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
                                r_id_out.write('{0}\t{1}\n'.format(alignment.title, hsp.score))
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
def recblast2(seqfile, target_species, fw_blast_db, infile_type="fasta", query_species="Homo sapiens", blasttype='blastp',
              localblast1=True, localblast2=False, rv_blast_db="refseq_proteins", evalue=0.001, identitiesperc=0.75,
              scoreperc=0.75, lengthperc=0.75, idthres='Score'):
    # recblast2 uses the experimental Biopython SearchIO module. Honestly it looks so much cleaner in the docs so I
    # kinda
    # jumped ship to fix recblast without the whole "for x in y / for y in z / for z in a / ..." hell.
    import os
    from Bio import SeqIO
    from Bio.Blast import NCBIWWW, NCBIStandalone
    from Bio.Blast.Applications import NcbiblastxCommandline,NcbiblastnCommandline,NcbiblastpCommandline, /
                                        NcbitblastxCommandline,NcbitblastnCommandline
    print("Now starting...")
    # First loop will iterate over each sequence in a file, preferably FASTA but also allows for GenBank
    for index, seq_record in enumerate(SeqIO.parse(seqfile, infile_type)):
        print("{}: {}".format(index,seq_record.name))

        # recblast_out is the file where the final Blast output will go to for this particular sequence. This is the
        # annotated
        #  FASTA, which will contain the top hits of the first BLAST, annotated with the name of the second BLAST
        # errorfile is the error log for the first BLAST parser, which will be used to get top hits for the reciprocal
        #  BLAST
        # errorfile2 is the error log for the second BLAST parser, which will be used to get the top hits of the
        #  reciprocal BLAST, which then will be used to annotate the hits from the primary BLAST
        recblast_out = os.path.join(os.getcwd(), "{}_recblast_out".format(target_species),
                                    "{0}_{1}.fasta".format(blasttype,seq_record.name))
        errorfile = open(os.path.join(os.getcwd(), "recblast_out", "{}_recblasterror1.log".format(seqfile)),"w+")
        errorfile2 = open(os.path.join(os.getcwd(), "recblast_out", "{}_recblasterror2.log".format(seqfile)), "w+")

        # Begin by opening recblast_out, and then start with the primary BLAST
        with open(recblast_out, "w+") as blastout:
            length = len(seq_record)
            if localblast1:
                if blasttype == "blastn":
                    NcbiblastnCommandline(query=recstring,db=fw_blast_db,evalue=evalue,outfmt=5,
                                          out="{0}to{1}blast_{2}.xml".format(query_species, target_species,
                                                                             seq_record.name))
                elif blasttype == "blastp":
                    NcbiblastpCommandline(query=recstring,db=fw_blast_db,evalue=evalue,outfmt=5,
                                          out="{0}to{1}blast_{2}.xml".format(query_species, target_species,
                                                                             seq_record.name))
                elif blasttype == "blastx":
                    NcbiblastxCommandline(query=recstring,db=fw_blast_db,evalue=evalue,outfmt=5,
                                          out="{0}to{1}blast_{2}.xml".format(query_species, target_species,
                                                                             seq_record.name))
                elif blasttype == "tblastx":
                    NcbitblastxCommandline(query=recstring,db=fw_blast_db,evalue=evalue,outfmt=5,
                                           out="{0}to{1}blast_{2}.xml".format(query_species, target_species,
                                                                              seq_record.name))
                elif blasttype == "tblastn":
                    NcbitblastnCommandline(query=recstring,db=fw_blast_db,evalue=evalue,outfmt=5,
                                           out="{0}to{1}blast_{2}.xml".format(query_species, target_species,
                                                                              seq_record.name))
                else:
                    raise Exception("Invalid blast choice!")
            else:
                with open(os.path.join(os.getcwd(), "recblast_out", "{0}to{1}blast_{2}.xml".format(query_species,
                                                                                                    target_species,
                                                                                                    seq_record.name)),
                          "w+") as fxml:
                    fxml.write(NCBIWWW.qblast(program=blasttype,database=fw_blast_db,sequence=seq_record.seq,
                                              entrez_query=(target_species + "[ORGN]"),expect=evalue))

            # Start Reciprocal Blast
            # First load the primary BLAST hits to a handle.
            blasthits = open("{0}to{1}blast_{2}.xml".format(query_species, target_species, seq_record.name),"r")
            # Next, load up the error parsers
            blast_error_parser = NCBIStandalone.BlastErrorParser(errorfile)
            blast_error_parser2 = NCBIStandalone.BlastErrorParser(errorfile2)

            # Iterate through all the blast hits
"""


def sirblastalot():
    # is a function that, for a list of sequences, BLASTS them against a list of organisms of one's choice; the
    # default list has an extensive vertebrate coverage. It also has an option to use BLAT instead. By default it just
    # does a unidirectional blast of your sequence to each organism; you can set it to do a Reciprocal Blast as well.
    pass
    # TODO: WRITE SIRBLASTALOT
