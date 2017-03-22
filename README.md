RecBlast
========

Recblast is a module containing functions that can do Reciprocal-Best-Hit-BLAST
for a given list of sequences between a query species, and a list of target species.
Unlike other implementations of Reciprocal-Best-Hit-BLAST, RecBlast is designed to
be flexible in its functionality. Rather than simply identifying a lone reciprocal
best-hit, RecBlast compiles all hits above a set of given criteria in order of
descending scores, and then annotates them with a list of reciprocal BLAST hits
that were compiled and sorted in the same fashion.

Additionally, RecBlast is designed to be flexible in how it goes about performing
the searches; in its final form, RecBlast will be able to not only use BLAST,
but will also be able to use BLAT and other custom programs, either locally or
on remote servers.

  (NOTE: As of the latest release, V1.0, RecBlastMP is the main module containing
the function recblastmp(), that can perform the Reciprocal-BLASTs in a
parallel fashion. Recblast.py should not be used, and instead,
RecBlastMP should be used with n_processes = 1 if one desires a
single-process version of the script.)

  <br>

###  ___Features:___
--------
  - Reciprocal BLASTs of any number of sequences between a query
species and any number of target species;
  - Multiprocessing-enabled using Python's multiprocessing module;
  - Options for specifying or automatically selecting local databases to use at each step;
  - Local or remote BLAST with customizable options.
  - Sequence retrival from either FASTA or SQL databases (local-only as of V1.0).

<br>

###  ___How it works:___
------------

  The script begins by initializing various RecBlast processes, one for each
RecBlast to be run, up to a mutable maximum limit. Once initialized, each RecBlast
process will then do the following:

  1. __Forward BLAST__ of the given sequence against the target species BLAST database;
  2. __Filter hits__ based on given parameters;
  3. __Get the full-length sequences__ of all filtered hits;
  4. __Reverse BLAST__ each filtered hit against the query species BLAST database;
  5. __Filter reciprocal hits__ based on given parameters;
  6. __Annotate hits__ using the filtered, sorted reciprocal hits.

  The user can choose to either write all intermediate files to disk; write
only the final, annotated file to record; or to write nothing at all. The
program returns a RecBlastContainer object with all the intermediate outputs
stored within so that the user can inspect the results in memory after
running the program.

<br>

###  ___Parameters and Settings:___
----------

##### Forward and Reverse BLAST

  The module has its own function, simply called blast(), that handles
the search functionality and calls to BLAST. In addition to being able
to specify whether to use a local BLAST installation or NCBI's servers,
both the Forward and Reverse BLAST calls take a dictionary parameter,
fw_blast_kwargs and rv_blast_kwargs, respectively, that can be used to
modify the two searches independently.

TODO: add more detail here?

<br>

### ___Roadmap:___
----------------------

  - [V1.5] Entrez implementation for remote sequence searches
  - [V2.0] BLAT compatibility
  - [Debated] Use hub-and-spoke model to coordinate all searches in a
single SearchMaster Process
  - [VX.x] Tools to use to analyze output

