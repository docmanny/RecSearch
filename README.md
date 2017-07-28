[![DOI](https://zenodo.org/badge/82135779.svg)](https://zenodo.org/badge/latestdoi/82135779)
[![GitHub release](https://img.shields.io/github/release/docmanny/RecBlast.svg)]()
[![Build Status](https://travis-ci.org/docmanny/RecBlast.svg?branch=master)](https://travis-ci.org/docmanny/RecBlast)

RecBlast
========

Recblast is a module containing functions that can do Reciprocal-Best-Hit-BLAST
for a given list of sequences between a query species, and a list of target species.
Unlike other implementations of Reciprocal-Best-Hit-BLAST, RecBlast is designed to
be flexible in its functionality; the function provides a variety of
changable parameters for filtering and annotating hits, allowing one
to either perform a exclusive best-hit annotation of all forward hits,
or else annotating all forward hits with a list of ranked reciprocal hits
that meet the specified criteria.

RecBlast is also designed to use either BLAST or BLAT as the search algorithm, and
allows mixing and matching of various different search types. The query sequences
can be either protein or nucleotide queries, and one can specify forward searches
against either protein, nucleotide, translated protein, or translated nucleotide
databases by specifying what kind of BLAST or BLAT search one wishes to employ;
in the reverse direction, one can also perform the reciprocal search with either
of the two algorithms, and against any kind of database one wishes given the output
type (either protein or nucelotide) of the first search.

The program is designed to automate much of the argument selection, with functions
to automatically select the correct database files for search algorithms and sequence
acquisition based on the selected search types and sequence selection. Finally, the
module was designed using the 'multiprocess' Python framework, which greatly speeds up
the Reciprocal Best-Hit BLAST/BLAT process by running different query-species searches
in parallel.

  <br>

###  ___Features:___
--------
  - Reciprocal BLAST or BLAT for any number of sequences between a query
species,  and any number of target species;
  - Parallelized searches using the Multiprocess framework;
  - Either manual or automatic specification of databases and/or various parameters;
  - Choice of either BLAST or BLAT searches, and of any combination of nucleotide
  or protein search queries;
  - Searches can be run either locally or on remote databases;
  - Sequence retrival from 2bit, FASTA, or SQL databases;
  - Additional functions for quick analysis of output.
  - A convenience function, RecBlastControl, allows the user to provide a control
  file that then dictates what searches to run with what criteria.
  - Integration of the MyGene package for convenient renaming of reciprocal hit
  annotations.

<br>

###  ___How it works:___
------------

  The script begins by initializing various RecBlast processes, one for each
RecBlast to be run, up to a mutable maximum limit. Once initialized, each RecBlast
process will then do the following:

  1. __Forward BLAST__ of the given sequence against the target species BLAST database;
  2. __Filters hits__ based on given parameters;
  3. __Fetches the full-length sequences__ of all filtered hits;
  4. __Reverse BLAST__ of each filtered hit against the query species BLAST database;
  5. __Filters reciprocal hits__ based on given parameters;
  6. __Annotates hits__ using the filtered, sorted reciprocal hits.

<br>

###  ___Parameters and Settings:___
----------

##### Sequence Input:
The main recblastMP() function accepts either a string or a Path object
pointing to a file with the input sequences; this file can either in either the
FASTA or GenBank file formats, and can have any number of included sequences. The
species to be queried can be provided as a list or tuple of strings, or as a
single string for a single species.

##### Forward and Reverse BLAST/BLAT:

The module has its own function, simply called blast(), that handles
the search functionality and calls to BLAST. In addition to being able
to specify whether to use a local BLAST installation or NCBI's servers,
both the Forward and Reverse BLAST calls take a dictionary parameter,
fw_blast_kwargs and rv_blast_kwargs, respectively, that can be used to
modify the two searches independently.

##### Hit Filtering:

The module provides various criteria for filtering hits, including: percent identity;
percent of query spanned by the hit; percent length spanned by the hit compared to
the longest hit; percent of score earned by the hit versus the top-scoring hit;
and the percent length spanned by a hit's HSP compared to the longest HSP in the
same hit. Using recblastMP(), one can set a unified set of criteria for all searches;
or if one wants to run different searches using different criteria, one can use
RecBlastControl with a control file to do so.

##### Sequence Fetching:

RecBlast uses the coordinates of the forward hits to obtain a full sequence
for the hit to use for the reverse search. To do this, the user can specify which
of various methods the fetchseq() function should employ to get the full sequence.
Currently, the module supports fetching sequences from FASTA files, SQL databases,
or from .2bit files when BLAT is used.


##### Annotations:
There are currently two methods that can be used for annotating forward hits:
the default "kitchen-sink" method, and the "best hit" method. The default method
will gather all the reverse hits associated with a given Reverse BLAST/BLAT, filter
them based on the forward search criteria, and then sort them based on their score.
The forward hit will then be annotated with said ordered list.
The "best hit" method will also filter and sort the results, but it will only annotate
the forward hit with the top-scoring hit.

##### Output:
When run, RecBlast will first create a folder named RecBlast_output in the active
directory, with a subfolder with the date and time of execution; within this
folder, it will create various log files for each process opened by the script
that will document all of its generated output (controlled by the verbosity level
assigned to the function). Inside the date-time folder it will also create
subfolders with the naming structure {species}\_recblast_out, where 'species'
represents the full species name given to the program for the RecBlast run. Each
species folder will then contain FASTA files with the naming structure
{search type}_{sequence name}.fasta, where 'search type' is the forward search
type used in the RecBlast run, and 'sequence name' is the name associated with the
sequence used as the query in the run. Additionally, the function itself will return
a RecBlastContainer object that contains all the intermediate and final results from
the executed code.
For very large runs, where storing various RecBlastContainer objects in memory would
result in excessive loads, one can set the parameter 'min_mem' to 'True', which then
purges RecBlastContainer objects from memory after the output files are written to
disk. Additionally, if one does not care for the full sequence information of the hits
and only wants the hit location information with annotations, setting 'hit_name_only'
to 'True' will result in the program only writing the hit name, location, and
annotations to disk in FASTA format as headers sans sequence.


<br>

### ___Roadmap:___
----------------------

  - Entrez implementation for remote sequence searches
  - Additional tools to use to analyze output
