RecBlast
========

Recblast is a module containing functions that can do Reciprocal-Best-Hit-BLAST for a given list of sequences between a query species, and a list of target species.

  As of the latest release, V1.0, RecBlastMP is the main module containing the function recblastmp(), that can perform the Reciprocal-BLASTs in a parallelized fashion.
Recblast.py should not be used, and instead, RecBlastMP should be used with n_processes = 1 if one desires a single-process version of the script.

  <br>