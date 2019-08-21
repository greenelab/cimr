

.. doc for processing variant-association files, created by YoSon Park


************************************
Processing variant-association files
************************************

This is an introduction to processing variant-association files using cimr. See
`guide to install <https://github.com/greenelab/cimr/doc/install.rst>`_.


===================
Input file examples
===================

* genome-wide association study results
* expression quantitative loci mapping results
* splicing quantitative loci mapping results
* topologically associating domain annotations


======================
Output (file) examples
======================

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
        Name of the input file. cimr processor can process variant-based association 
        summary statistis. Additionally annotation files such as genomic coordinates 
        generated from Hi-C experiments can be integrated into cimr-d for public
        sharing and research use.
    
    data-type: {gwas, eqtl, gene, tad}

    genome-build: {b37, b38}
        Currently, cimr supports two versions of the human genome references:
        b37 (GRCh37/HG19) and b38 (GRCh38/HG20/HG38). 
        `cimr-d <https://github.com/greenelab/cimr-d>`_ contains data annotated to 
        b38 and all data processed through 
        `cimr <https://github.com/greenelab/cimr>`_ can be converted from
        b37 to b38 using the command --update-map.




=============
Example usage
=============

Association files can be checked for common formatting errors.

-----------------------------
Example local file processing
-----------------------------



>>> 


