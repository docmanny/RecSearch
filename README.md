---
---

RecBlast
========
[![DOI](https://zenodo.org/badge/82135779.svg)](https://zenodo.org/badge/latestdoi/82135779)

RecBlast is a module that provides a Reciprocal searching framework for finding homologues
of sequences between a query species, and a list of target species. RecBlast is designed
to be flexible by allowing the user to customize each step of the pipeline, from
the forward and reverse search algorithms, to the sequence database for hits. The output
of the search is also highly customizable, and can be set to either provide 1:1 Best-Hits or
1:Many hits. The package is written in a modular, parallelized fashion to maximize search
efficiency while also allowing for easy additions and modifications to the pipeline.

  <br>

###  ___Features:___
--------
  - Reciprocal searching of any number of sequences between a query
species and any number of target species
  - Runs on multiple cores
  - Uses either BLAT or BLAST for the searches, with room for other options such as HMMER
  - Local or remote search host options
  - Can use FASTA files, BLAT or BLAST databases, SQL databases, or Ensembl for sequence retreval


<br>

###  ___How it works:___
------------

The core of the module is the RecSearch object, which organizes the various settings for the
reciprocal search. After initializing it with the target and query species, the search types,
and the sequence source, one can then change the criteria for the searches or change advanced
settings.

``` python
from RecBlast.RecBlast import RecSearch
# Initialize the RecSearch Object
# We are using protein sequences from Humans to query the African Elephant genome,
# using BLAT as the algorithm, and obtaining results from the twobit file.
recblast = RecSearch(target_species="Loxodonta africana", query_species="Homo sapiens",
                     forward_search_type="tblat", reverse_search_type="blat-transcript",
                     sequence_source="twobit", verbose=0)
# Read query sequences from a FASTA file
recblast.set_queries("Test/query/3protein.fasta", infile_type="fasta")
# For the forward BLAT search, set 'database' to the loxAfr 2bit file
# RecSearch will by default automatically search for a 2bit file otherwise.
recblast.forward_search_settings['database'] = {"Loxodonta africana": "loxAfr3.2bit"}
# For the forward and reverse BLAT searches, set the ports of the gfServer for each species
# By default, it looks for the server at localhost, but it can also be remote if one wishes
recblast.forward_search_settings['database_port'] = {"Loxodonta africana":20008}
recblast.reverse_search_settings['database_port'] = {"Homo sapiens":20005}
# Obtain the sequence for the forward hit from the twobit file.
recblast.sequence_source_settings['database'] = {"Loxodonta africana": "loxAfr3.2bit"}
```

 The search is then initialized by running the object, and will return a
RecBlastContainer object with all the results and intermediary outputs.

``` python
rc = recblast("Example", output_type="BED)
```

  The script begins by initializing various RecSearch worker processes, one for each
species-sequence pair, up to a mutable maximum limit (by default, it is half the number of
CPUs detected by `multiprocess`). Once initialized, each RecSearch process executes the
following protocol:

  1. __Forward Search__ of the given sequence against the target species database;
  2. __Filter hits__ based on given parameters;
  3. __Get the full-length sequences__ of all filtered hits;
  4. __Reverse Search__ each filtered hit against the query species database;
  5. __Filter reciprocal hits__ based on given parameters;
  6. __Annotate hits__ using the filtered, sorted reciprocal hits.



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


<br>

### ___Roadmap:___
----------------------

  - Entrez implementation for remote sequence searches
  - Additional tools for analyzing output
---
---

RecBlast
========
[![DOI](https://zenodo.org/badge/82135779.svg)](https://zenodo.org/badge/latestdoi/82135779)

RecBlast is a module that provides a Reciprocal searching framework for finding homologues
of sequences between a query species, and a list of target species. RecBlast is designed
to be flexible by allowing the user to customize each step of the pipeline, from
the forward and reverse search algorithms, to the sequence database for hits. The output
of the search is also highly customizable, and can be set to either provide 1:1 Best-Hits or
1:Many hits. The package is written in a modular, parallelized fashion to maximize search
efficiency while also allowing for easy additions and modifications to the pipeline.

  <br>

###  ___Features:___
--------
  - Reciprocal searching of any number of sequences between a query
species and any number of target species
  - Runs on multiple cores
  - Uses either BLAT or BLAST for the searches, with room for other options such as HMMER
  - Local or remote search host options
  - Can use FASTA files, BLAT or BLAST databases, SQL databases, or Ensembl for sequence retreval


<br>

###  ___How it works:___
------------

The core of the module is the RecSearch object, which organizes the various settings for the
reciprocal search. After initializing it with the target and query species, the search types,
and the sequence source, one can then change the criteria for the searches or change advanced
settings.

``` python
from RecBlast.RecBlast import RecSearch
# Initialize the RecSearch Object
# We are using protein sequences from Humans to query the African Elephant genome,
# using BLAT as the algorithm, and obtaining results from the twobit file.
recblast = RecSearch(target_species="Loxodonta africana", query_species="Homo sapiens",
                     forward_search_type="tblat", reverse_search_type="blat-transcript",
                     sequence_source="twobit", verbose=0)
# Read query sequences from a FASTA file
recblast.set_queries("Test/query/3protein.fasta", infile_type="fasta")
# For the forward BLAT search, set 'database' to the loxAfr 2bit file
# RecSearch will by default automatically search for a 2bit file otherwise.
recblast.forward_search_settings['database'] = {"Loxodonta africana": "loxAfr3.2bit"}
# For the forward and reverse BLAT searches, set the ports of the gfServer for each species
# By default, it looks for the server at localhost, but it can also be remote if one wishes
recblast.forward_search_settings['database_port'] = {"Loxodonta africana":20008}
recblast.reverse_search_settings['database_port'] = {"Homo sapiens":20005}
# Obtain the sequence for the forward hit from the twobit file.
recblast.sequence_source_settings['database'] = {"Loxodonta africana": "loxAfr3.2bit"}
```

 The search is then initialized by running the object, and will return a
RecBlastContainer object with all the results and intermediary outputs.

``` python
rc = recblast("Example", output_type="BED)
```

  The script begins by initializing various RecSearch worker processes, one for each
species-sequence pair, up to a mutable maximum limit (by default, it is half the number of
CPUs detected by `multiprocess`). Once initialized, each RecSearch process executes the
following protocol:

  1. __Forward Search__ of the given sequence against the target species database;
  2. __Filter hits__ based on given parameters;
  3. __Get the full-length sequences__ of all filtered hits;
  4. __Reverse Search__ each filtered hit against the query species database;
  5. __Filter reciprocal hits__ based on given parameters;
  6. __Annotate hits__ using the filtered, sorted reciprocal hits.



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


<br>

### ___Roadmap:___
----------------------

  - Entrez implementation for remote sequence searches
  - Additional tools for analyzing output
