"""
Author: Juan Manuel Vazquez
Date Created: 11/01/2016
Date Last Modified: 01/19/2017
"""


def merge_ranges(ranges):
    """
    Merge overlapping and adjacent ranges and yield the merged ranges in order.
    The argument must be an iterable of pairs (start, stop).
    (Source: Gareth Rees, StackExchange)

    >>> list(merge_ranges([(5,7), (3,5), (-1,3)]))
    [(-1, 7)]
    >>> list(merge_ranges([(5,6), (3,4), (1,2)]))
    [(1, 2), (3, 4), (5, 6)]
    >>> list(merge_ranges([]))
    []
    """
    ranges = iter(sorted(ranges))
    try:
        current_start, current_stop = next(ranges)
    except StopIteration:  # ranges is empty
        return
    for start, stop in ranges:
        if start > current_stop:
            # Gap between segments: output current segment and start a new one.
            yield current_start, current_stop
            current_start, current_stop = start, stop
        else:
            # Segments adjacent or overlapping: merge.
            current_stop = max(current_stop, stop)
    yield current_start, current_stop


def biosql_DBSeqRecord_to_SeqRecord(DBSeqRecord_, off=False):
    """
    As I wrote this script I realized very quickly that there are a great many issues with DBSeqRecords throwing bugs
    left and right, so I got tired of trying to debug a blackbox and instead decided to convert it to something that
    works.
    :param DBSeqRecord_:
    :return:
    """
    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import Seq
    if off:
        return DBSeqRecord_
    else:
        return SeqRecord(seq=Seq(str(DBSeqRecord_.seq)), id=DBSeqRecord_.id, name=DBSeqRecord_.name,
                         description=DBSeqRecord_.description, dbxrefs=DBSeqRecord_.dbxrefs,
                         features=DBSeqRecord_.features, annotations=DBSeqRecord_.annotations,
                         letter_annotations=DBSeqRecord_.letter_annotations)


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

    count = 0

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
    from pathlib import Path

    server = BioSeqDatabase.open_database(driver=driver, user=user, passwd=passwd, host=host, db=db)
    db = server[sub_db_name]
    if dumpfile:
        import sys
        og_stdout = sys.stdout
        outfile = Path('./'+sub_db_name+'.index').open("w")
        sys.stdout = outfile
    print("This database contains {} records".format(len(db)))
    for key, record in db.items():
        if key == 1:
            print("As an example of what items are in the database, here is the first entry:")
            print(key, record)
        print("Key {0} maps to a sequence record with id {1} and name {2}".format(key, record.id, record.name))
    if dumpfile:
        sys.stdout = og_stdout
        outfile.close()
    # End of Function


def biosql_getrecord(sub_db_name, passwd, id_list=list(), id_type='accession', driver="psycopg2", user="postgres",
                     host="localhost", db="bioseqdb", verbose=True, parallel=False):  # TODO: FILL OUT DOCSTRING
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
    # Note: DBSeqRecord has many bugs that seem to center on the common thread:
    #   "'DBSeq' object has no attribute '_data'"
    # So I'm wrapping all fetching functions in biosql_DBSeqRecord_to_SeqRecord() to avoid the issue.

    from BioSQL import BioSeqDatabase
    import time
    if id_list:
        pass
    else:
        raise Exception('Received List was empty')
    try_counter = 0
    if parallel:
        from random import randrange
        time_to_wait = randrange(1,10,step=0.5)
    else:
        time_to_wait = 0
    while try_counter <=3:
        time.sleep(time_to_wait)
        try_counter += 1
        try:
            if verbose:
                print("Opening database server ", db)
            server = BioSeqDatabase.open_database(driver=driver, user=user, passwd=passwd, host=host, db=db)
            if verbose:
                print("Server opened successfully!")
            try_counter = 4
        except:
            print('Sorry, couldn\'t open the database! Retrying in {} seconds...'.format(str(time_to_wait**
                                                                                             try_counter)))
    try:
        server
    except NameError:
        raise Exception('Server could not be opened!')
    try_counter = 0
    if parallel:
        from random import randrange
        time_to_wait = randrange(1, 10, step=0.5)
    else:
        time_to_wait = 0
    while try_counter <= 3:
        time.sleep(time_to_wait)
        try_counter += 1
        try:
            if verbose:
                print("Opening sub-database ", sub_db_name)
            dtbse = server[sub_db_name]
            if verbose:
                print("Successfully opened sub-database!")
            try_counter = 4
        except:
            print('Sorry, couldn\'t open the sub-database! Retrying in {} seconds...'.format(str(time_to_wait**
                                                                                             try_counter)))
    try:
        dtbse
    except NameError:
        raise Exception('sub-database could not be opened!')
    seqdict = {}
    for identifier in id_list:
        try_get_id = True
        if try_get_id:
            try:
                if verbose:
                    print("Now searching database {0} for {1}: {2}".format(sub_db_name, id_type, identifier))
                seqdict[identifier] = biosql_DBSeqRecord_to_SeqRecord(dtbse.lookup(**{id_type: identifier}))
                if verbose:
                    print('Got sequence for {}!'.format(identifier))
                try_get_id = True
            except IndexError as err:
                print("WARNING: couldn't find {0} using given ID type... \n Full error: {1}".format(identifier, err))
        if try_get_id:
            identifier_sans_subnumber = identifier.split('.')[0]
            if verbose:
                print('Seeing if removing any sub-numbers (acc: xxxxxx.1 for example) helps...')
                print('Identifier: ', identifier_sans_subnumber)
            try:
                if verbose:
                    print("Now searching database {0} for {1}: {2}".format(sub_db_name, id_type,
                                                                           identifier_sans_subnumber))
                seqdict[identifier] = biosql_DBSeqRecord_to_SeqRecord(dtbse.lookup(**{id_type: identifier_sans_subnumber}))
                if verbose:
                    print('Got sequence for {}!'.format(identifier))
                try_get_id = True
            except IndexError as err1:
                print("WARNING: couldn't find {0} using abbreviated ID... \n Full error: {1}"
                      .format(identifier_sans_subnumber,
                              err1))
        if try_get_id:
            try:
                if verbose:
                    print('Attempting to search using Primary ID instead of declared type:')
                seqdict[identifier] = biosql_DBSeqRecord_to_SeqRecord(dtbse.lookup(primary_id=identifier))
                if verbose:
                    print('Got sequence for {}!'.format(identifier))
                try_get_id = False
            except IndexError as err2:
                print("WARNING: couldn't find {0} using Primary ID... \n full error: {1}".format(identifier, err2))
                if verbose:
                    print('Primary ID search didn\'t work! Troubleshooting...')
        if try_get_id:
            try:
                if verbose:
                    print('Attempting to search using name instead of declared type:')
                seqdict[identifier] = biosql_DBSeqRecord_to_SeqRecord(dtbse.lookup(name=identifier))
                if verbose:
                    print('Got sequence for {}!'.format(identifier))
                try_get_id = False
            except IndexError as err3:
                print("WARNING: Still couldn't find {0}: \n full error: {1}".format(identifier, err3))
                if verbose:
                    print('Name search didn\'t work! Troubleshooting...')
        if try_get_id:
            try:
                id_type = input('Last shot, chose an ID type: '
                                '[accession, primary_id, gi, version, display_id, name]')
                if id_type == 'exit':
                    exit(exit(), 'Script ended!')
                seqdict[identifier] = biosql_DBSeqRecord_to_SeqRecord(dtbse.lookup(**{id_type: identifier}))
                if verbose:
                    print('Got sequence for {}!'.format(identifier))
                try_get_id = False
            except IndexError as err5:
                print("WARNING: COULD NOT FIND SEQUENCES FOR ID:{0}: \n full error: {1}".format(identifier, err5))
        if try_get_id:
            continue
    server.close()
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
             db="nucleotide", delim='\t', id_type='accession', batch_size=50, passwd='', version='1.0', verbose=True,
             parallel=False):
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
        id_prelist = list(filter(None, id_prelist))
        if verbose:
            print('Full header for Entry 1:')
            try:
                print(id_prelist[0])
            except IndexError:
                print('No items found!')
                raise
    # Check to make sure list is not empty
    if verbose and (not id_prelist or id_prelist is None):
        print('id_prelist is empty!')

    id_list = [str(item.split(delim)) for item in id_prelist]    # Breaks the tab sep in the lines into strings

    # Define the regex functions
    p = [re.compile('(gi)([| :_]+)(\d\d+\.?\d*)(.*)'),                      # regex for gi
         re.compile('([AXNYZ][MWRCPGTZ]|ref)([| _:]+)(\d\d+\.?\d*)(.*)'),   # regex for accession
         re.compile('(id)([| :_]+)(\d\d+\.?\d*)(.*)'),                      # regex for generic ID
         re.compile(':(\d+)-(\d+)'),                                        # regex for sequence range
         ]
    id_list_ids = []    # Initialized list of IDs
    seq_range = {}      # Initialized dict of sequence ranges

    # Begin search:
    if verbose:
        print('ID File Loaded, performing regex search for identifiers...')
        print('ID Specified as: ', id_type)
    if id_type == 'brute':
        if bool(p[1].findall(id_list[0])):
            id_type = 'accession'
            if verbose:
                print(p[1].findall(id_list[0]))
        elif bool(p[0].findall(id_list[0])):
            id_type = 'gi'
            if verbose:
                print(p[1].findall(id_list[0]))
        elif bool(p[2].findall(id_list[0])):
            id_type = 'id'
            if verbose:
                print(p[1].findall(id_list[0]))
        else:
            id_type = 'other'
        if verbose:
            print('Brute Force was set, tested strings for all pre-registered IDs. ID was selected as type ', id_type)
    if id_type == 'gi':
        if bool(p[0].findall(id_list[0])):
            found_id = True
            if verbose:
                print('Successfully found GI numbers, compiling list!')
            for item in id_list:
                item_parts = p[0].findall(item)
                if verbose:
                    print('Item:\t', item_parts)
                id_list_ids.append(item_parts[0][0:3])
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
                print('Successfully found accession numbers, compiling list!')
            for item in id_list:
                item_parts = p[1].findall(item)
                if verbose:
                    print('Item:\t', item_parts)
                id_list_ids.append(item_parts[0][0:3])
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
                print('Successfully found ID numbers, compiling list!')
            for item in id_list:
                item_parts = p[2].findall(item)
                if verbose:
                    print('Item:\t', item_parts)
                id_list_ids.append(item_parts[0][0:3])
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
                print('Successfully found custom ID numbers, compiling list!')
            for item in id_list:
                item_parts = p[4].findall(item)
                if verbose:
                    print('Item:\t', item_parts)
                id_list_ids.append(item_parts[0][0:3])
                if bool(p[3].findall(str(item))):
                    seq_range[''.join(p[4].findall(item)[0][0:3])] = p[3].findall(item)[0]
                    if verbose:
                        print('Found sequence delimiters in IDs!')
        else:
            print('Sorry, still can\'t find it...')
    if verbose:
        print('ID list: ')
        for index, ID_item in enumerate(id_list_ids):
            print(index+1, ': ', ''.join(ID_item))

    # Armed with the ID list, we fetch the sequences from the appropriate source
    if source.lower() == "entrez":  # Todo: Make sure this will actually output the correct sequence range...
        if verbose:
            print('Source selected was Entrez. Beginning search now:')
        from Bio import Entrez
        from urllib.error import HTTPError
        from time import sleep
        Entrez.email = email
        if verbose:
            print('Entrez email set as: ', email)
        id_str = ",".join([i[2] for i in id_list_ids])
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
        if verbose:
            print('Searching for sequences in local PostgreSQL db...')
        sub_db_name = ''.join([i[0:3] for i in species.title().split(' ')]) + version
        id_list_search = [''.join(i[0:3]) for i in id_list_ids]
        seqdict = biosql_getrecord(sub_db_name=sub_db_name, id_list=id_list_search, id_type=id_type,
                                   passwd=passwd, driver="psycopg2", user="postgres", host="localhost", db="bioseqdb",
                                   parallel=parallel, verbose=verbose)
        itemsnotfound = [''.join(x) for x in id_list_ids if ''.join(x) not in seqdict.keys()]
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
                if int(seq_range[k][0]) > int(seq_range[k][1]):
                    tmp_id = seqdict[k].id
                    tmp_name = seqdict[k].name
                    tmp_desc = seqdict[k].description
                    tmp_dbxrefs = seqdict[k].dbxrefs
                    tmp_feat = seqdict[k].features
                    tmp_annotations = seqdict[k].annotations
                    tmp_let_anno = seqdict[k].letter_annotations
                    seqdict[k].seq = seqdict[k][int(seq_range[k][1]):int(seq_range[k][0])].seq.reverse_complement()

                else:
                    seqdict[k] = seqdict[k][int(seq_range[k][0]):int(seq_range[k][1])]
                if verbose:
                    print('Seq_description_full: ', seq_description_full)
                    print('id_range: ', id_range[1:])
                if int(seq_range[k][0]) > int(seq_range[k][1]):
                    seqdict[k].description = ''.join(seq_description_full[0:3]) + id_range + '(-)' + \
                                                     str(seq_description_full[3])
                else:
                    seqdict[k].description = ''.join(seq_description_full[0:3]) + id_range + '(+)' + \
                                             str(seq_description_full[3])
                if verbose:
                    print('Sequence Description: \n\t', seqdict[k].description)
                seqdict[k].id += id_range
                if verbose:
                    print('Sequence ID: \n\t', seqdict[k].id)
                    if id_range:
                        print('Length of subsequence with range {0}: {1}'.format(id_range, len(seqdict[k])))
        if verbose:
            print('Sequence Record post-processing, to be saved:')
            print(seqdict)

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
                        print('\tFull length of sequence: {}'.format(len(seqdict[k])))
                    if id_type == 'gi':
                        seq_description_full = p[0].findall(seqdict[k].description)[0]
                    elif id_type == 'accession':
                        seq_description_full = p[1].findall(seqdict[k].description)[0]
                    elif id_type == 'id':
                        seq_description_full = p[2].findall(seqdict[k].description)[0]
                    else:
                        seq_description_full = p[4].findall(seqdict[k].description)[0]
                if verbose:
                    print(int(seq_range[k][0]))
                    print(int(seq_range[k][1]))
                id_range = ':' + '-'.join(seq_range[k])
                seqdict[k] = seqdict[k][int(seq_range[k][0]):int(seq_range[k][1])]
                seqdict[k].description = ''.join(seq_description_full[0:3]) + id_range + str(seq_description_full[3])
                seqdict[k].id += id_range
                if verbose:
                    print('\tLength of subsequence with range{0}: {1}'.format(id_range, len(seqdict[k])))
                else:
                    if verbose:
                        print('No sequence range found, continuing...')
                    continue

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


def blast(seq_record, target_species, database, query_species="Homo sapiens", filetype="fasta", blasttype='blastn',
          localblast=False, expect=0.005, megablast=True, blastoutput_custom="", perc_ident=75,
          verbose=True):

    from pathlib import Path
    from Bio import SeqIO
    from Bio.Blast import NCBIWWW
    from Bio.Blast.Applications import NcbiblastxCommandline, NcbiblastnCommandline, NcbiblastpCommandline, \
        NcbitblastxCommandline, NcbitblastnCommandline
    if isinstance(seq_record, SeqIO.SeqRecord):
        pass
    else:
        seq_record = SeqIO.read(seq_record, filetype)
    if verbose:
        print("Now starting BLAST...")
    if blasttype != 'blastn':
        megablast = False
    # Begin by opening recblast_out, and then start with the primary BLAST
    if blastoutput_custom == '':
        blastoutput_custom = Path("{0}_blast".format(target_species),
                                  "{0}_{1}_{2}_to_{3}.xml".format(blasttype, seq_record.name,
                                                                  query_species, target_species))
    else:
        blastoutput_custom = Path(blastoutput_custom)
    recstring = str(seq_record.seq)
    if localblast:
        if verbose:
            print('Running Local Blast...')
        if blasttype == "blastn":
            NcbiblastnCommandline(query=recstring, db=database, expect=expect, outfmt=5, megablast=megablast,
                                  out=str(blastoutput_custom), perc_ident=perc_ident)
        elif blasttype == "blastp":
            NcbiblastpCommandline(query=recstring, db=database, expect=expect, outfmt=5, megablast=megablast,
                                  out=str(blastoutput_custom), perc_ident=perc_ident)
        elif blasttype == "blastx":
            NcbiblastxCommandline(query=recstring, db=database, expect=expect, outfmt=5, megablast=megablast,
                                  out=str(blastoutput_custom), perc_ident=perc_ident)
        elif blasttype == "tblastx":
            NcbitblastxCommandline(query=recstring, db=database, expect=expect, outfmt=5, megablast=megablast,
                                   out=str(blastoutput_custom), perc_ident=perc_ident)
        elif blasttype == "tblastn":
            NcbitblastnCommandline(query=recstring, db=database, expect=expect, outfmt=5, megablast=megablast,
                                   out=str(blastoutput_custom), perc_ident=perc_ident)
        else:
            raise Exception("Invalid blast choice!")
    else:
        blast_result = NCBIWWW.qblast(program=str(blasttype), database=str(database), sequence=recstring,
                                      entrez_query='"{}"[ORGN]'.format(target_species), expect=expect,
                                      megablast=megablast, perc_ident=perc_ident)
        if verbose:
            print('Done, saving to outfile....')
        with blastoutput_custom.open("w") as fxml:
            fxml.write(blast_result.read())


def blast_many(seqfile, target_species, database, query_species="Homo sapiens", filetype="fasta", blasttype='blastn',
               localblast=False, expect=0.005, megablast=True, blastoutput_custom="", perc_ident=75, verbose=True):
    # TODO: TEST IT 
    # TODO: WRITE DOCSTRING
    # TODO: AFTER RECBLAST WORKS, CLEAN UP COMMENTED-OUT CODE
    from Bio import SeqIO

    if verbose:
        print("Now starting BLAST...")

    # First loop will iterate over each sequence in a file, preferably FASTA but also allows for GenBank
    for index, seq_record in enumerate(SeqIO.parse(seqfile, filetype)):
        blast(seq_record=seq_record, target_species=target_species, database=database, query_species=query_species,
              filetype=filetype, blasttype=blasttype, localblast=localblast, expect=expect,
              megablast=megablast, blastoutput_custom=blastoutput_custom, perc_ident=perc_ident, verbose=verbose)
    # Done


def recblast(seqfile, target_species, fw_blast_db='chromosome', infile_type="fasta", output_type="fasta",
             query_species="Homo sapiens", blasttype='blastn', localblast1=False, localblast2=False,
             rv_blast_db="nt", expect=10, scoreperc=0.50, perc_ident=50, perc_length=0.5,
             megablast=True, email='', id_type='brute', fw_source="psql", fw_id_db="", batch_size=50,
             passwd='', fw_id_db_version='1.0', verbose=True, parallel=False):
    """By Reciprocal BLAST, finds orthologs in Species 2 of a list of genes from Species 1 and annotates them.

    Reciprocal BLAST involves using a primary BLAST to identify putative orthologs in the "target_species" using
     sequences from the "query_species", which by default is "Homo sapiens".
    Input is a list of genes ("seqfile") saved as a specified "infile_type" (defaults to FASTA), to be searched against
    an indicated database. Other options include:
    blasttype -- BLAST program to be used. ("blastn", "blastp", "blastx", "tblastx", "tblastn")
    localblast1 -- Should the Forward BLAST be done locally or at NCBI? (default True)
    localblast2 -- Should the Reverse BLAST be done locally or at NCBI? (default False)
    rv_blast_db -- Database to be queried for the Reverse Blast to ID putative orthologs. (default "RefSeq_Genes")
    expect -- Maximum E-Value accepted from HSPs
    identitiesperc -- Minimum percent identity accepted from HSPs
    scoreperc -- Minimum percentage from the top score that will be used as a cut-off for putative orthologs.
    lenghtperc -- Minimum fraction of the total length of the alignment that will be accepted.
    idthres -- TODO: clarify
    """

    from pathlib import Path
    from Bio import SeqIO, __version__
    from Bio.Blast import NCBIXML
    if parallel:
        pass

    if verbose:
        print("Now starting RecBlast...")
        print('BioPython Version: ', __version__)
    if isinstance(seqfile,list):
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
                                    "{0}_{1}_tmp".format(blasttype, seq_record.name).replace(' ', '_') + '/' +
                                    "{0}_{1}_{2}_to_{3}.xml".format(blasttype, seq_record.name, query_species,
                                                                    target_species).replace(' ', '_'))

        forward_id_score_output = Path("{0}_recblast_out".format(target_species).replace(' ', '_') + '/' +
                                       "{0}_{1}_tmp".format(blasttype, seq_record.name).replace(' ', '_') + '/' +
                                       "{0}_{1}_{2}_to_{3}.ID_Scores.tmp".format(blasttype, seq_record.name,
                                                                                 query_species,
                                                                                 target_species).replace(' ', '_'))

        recblast_output_unanno = Path("{0}_recblast_out".format(target_species).replace(' ', '_') + '/' +
                                      "{0}_{1}_tmp".format(blasttype, seq_record.name).replace(' ', '_') + '/' +
                                      "unannotated_{0}_{1}.tmp".format(blasttype, seq_record.name).replace(' ', '_'))

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
                  query_species=query_species, filetype=infile_type, blasttype=blasttype, localblast=localblast1,
                  expect=expect, megablast=megablast, blastoutput_custom=str(forward_blast_output),
                  perc_ident=perc_ident)
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
                    if ((hsp.score >= (scoreperc * align_scorelist[align_index])) and (hsp.expect <= expect) and
                            (sum([i[-1]-i[0] for i in query_start_end[align_index]])/blastrecord.query_length
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
                     id_type=id_type, batch_size=batch_size, passwd=passwd, version=fw_id_db_version, verbose=verbose,
                     parallel=parallel)
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
                               "{0}_{1}.{2}".format(blasttype, seq_record.name, output_type).replace(' ', '_'))
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
                                        "{0}_{1}_tmp".format(blasttype, seq_record.name).replace(' ', '_') + '/' +
                                        "{0}_{1}_{3}_to_{2}_{4}.xml".format(blasttype, seq_record.name,
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
                      query_species=target_species, filetype=infile_type, blasttype=blasttype, localblast=localblast2,
                      expect=expect, megablast=megablast, blastoutput_custom=str(reverse_blast_output),
                      perc_ident=perc_ident)
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
                for alignment in blastrecord2.alignments:
                    if verbose:
                        print('Sorting through alignment\'s HSPs to get top scores of all alignments...')
                    subject_range_hsp2 = []
                    query_start_end_hsp2 = []
                    for hsp2 in alignment.hsps:
                        hsp_scorelist2.append(hsp2.score)
                        subject_range_hsp2.append(hsp2.sbjct_start)
                        subject_range_hsp2.append(hsp2.sbjct_end)
                        query_start_end_hsp2.append((hsp2.query_start, hsp2.query_end))
                    hsp_scorelist2.sort(reverse=True)
                    query_start_end2.append(i for i in merge_ranges(query_start_end_hsp2))
                    subject_range2.append((subject_range_hsp2[0], subject_range_hsp2[-1]))
                    if verbose:
                        print("HSP Score List: \n\t", hsp_scorelist2)
                    align_scorelist2.append(hsp_scorelist2[0])
                    if verbose:
                        print("Alignment Score List: \n\t", align_scorelist2)
                if verbose:
                    print('Done with sorting!')
                # Now we have a list of the top score of each alignment for the current entry_record.
                with recblast_output.open("w+") as rb_out:
                    if verbose:
                        print('Annotating BLAST results')
                    has_written2 = False
                    for align_index2, alignment in enumerate(blastrecord2.alignments):
                        blast_got_hit2 = False
                        for hsp2 in alignment.hsps:
                            if (hsp2.score >= (scoreperc * align_scorelist2[align_index2])):
                                print('hsp score above threshold')
                                if (hsp2.expect <= expect):
                                    print('hsp expect below threshold')
                                    if (sum([i[-1] - i[0] for i in query_start_end2[align_index2]])/blastrecord2.query_length >= perc_length):
                                        print('hsp perc length above threshold')
                                        if verbose:
                                            print('Found hit!')
                                        entry_record.id += '\t||{0} ({1})'.format(alignment.title +
                                                                          ' [:{0}-{1}]'.format(
                                                                              subject_range2[align_index2][0],
                                                                              subject_range2[align_index2][0]
                                                                                              ), hsp.score)
                                        SeqIO.write(entry_record, rb_out, output_type)
                                        has_written2 = True
                                        blast_got_hit2 = True
                                    else:
                                        print('WARNING HSP LENGTH BELOW THRESHOLD')
                                        print(sum([i[-1] - i[0] for i in
                                                   query_start_end2[align_index2]]) / blastrecord2.query_length,' not greater than ',perc_length)
                                else:
                                    print('WARNING HSP EXPECT ABOVE THRESHOLD')
                                    print(hsp2.expect, 'not less than', expect)
                            else:
                                print('WARNING HSP SCORE BELOW THRESHOLD')
                                print(hsp2.score, ' not greater than ', (scoreperc * align_scorelist2[align_index2]))

                            #else:
                                #continue
                        if not blast_got_hit2:
                            print('NOTE: Alignment {} was not used to annotate.'.format(alignment.title))
                    if not has_written2:
                        print(Warning('NONE OF THE REVERSE BLAST HITS FOR THIS RUN MET ANNOTATION CRITERIA!'))
                        continue
        if verbose:
            print('DONE!!!!')
        #Done!


def sirblastalot():
    # is a function that, for a list of sequences, BLASTS them against a list of organisms of one's choice; the
    # default list has an extensive vertebrate coverage. It also has an option to use BLAT instead. By default it just
    # does a unidirectional blast of your sequence to each organism; you can set it to do a Reciprocal Blast as well.
    pass
    # TODO: WRITE SIRBLASTALOT
