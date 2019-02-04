

************************************
Processing variant-association files
************************************

This is an introduction to processing variant-association files using cimr.

===========
Input files
===========

* genome-wide association study results
* expression quantitative loci mapping results
* splicing quantitative loci mapping results
* protein quantitative loci mapping results
* topologically associating domain annotations


============
Output files
============

* uniformly formatted association results for gene tests
* summary of missing or mis-formatted data
* (recommended) integration into a public repository (cimr-d)
* (optional) annotated results


=======================
Parameters and defaults
=======================

Following are the list of arguements available when using cimr as a standalone 
tool::

    file-name:
    
    data-type: {gwas, eqtl}

    genome-build:




=============
Example usage
=============

Association files can be checked 

>>> cimr processor --data-type eqtl --file-name eqtl.txt


