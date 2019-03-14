

****
cimr
****


***************************************
cimr is not yet released for public use
***************************************


================================================================================
continuous integration and analysis using variant association summary statistics
================================================================================

==========
YoSon Park
==========

**Useful links**:
`Source repository <https://github.com/greenelab/cimr>`_ |
`Issues & Ideas <https://github.com/greenelab/cimr/issues>`_ | 
`Documentation <https://cimr.readthedocs.io>`_ | 
`cimr-d <https://github.com/greenelab/cimr-d>`_


cimr is a convenience tool for continuous analyses of variant-based 
association results from GWAS (genome-wide association studies), eQTL 
(expression-quantitative trait loci mapping) or other association studies. 
cimr began as a python module to run large-scale Mendelian randomization 
analysis (hence the name). As the project developed, it became more 
evident that there are many parts preceding the analyses that require 
pipelining. So the current incarnation of cimr aims to streamline the 
pre-analysis processing steps, provide standardized input files and write 
example scripts to run various downstream methods seamlessly.



============
Installation
============

-----------------
Installing python
-----------------

cimr requires python :math: `\ge` 3.6. Installation of data analysis bundles 
such as `miniconda <https://conda.io/miniconda.html>`_ or 
`anaconda <https://www.anaconda.com/download/>`_ are recommended and will 
install all python packages cimr depends on. However, all required python 
packages can be downloaded and installed with setup.py or requirements.txt 
provided here.


------------------
Installing git lfs
------------------

cimr-d and some functionalities in cimr may use 
`git large file storage (LFS) <https://git-lfs.github.com/>`_ . 
See how to install `git <https://www.atlassian.com/git/tutorials/install-git>`_ .
git-lfs is not required for using cimr as a standalone tool without cimr-d.


To install git-lfs on Ubuntu, run:


>>> curl -s https://packagecloud.io/install/repositories/github/git-lfs/script.deb.sh | sudo bash
>>> sudo apt-get install -y git git-lfs
>>> git-lfs install


Alternatively, you can install git-lfs through conda:

>>> conda install -c conda-forge git-lfs && git lfs install


---------------
Installing cimr
---------------

You can use pip to install the latest stable release of cimr.

>>> pip3 install cimr


If you want to try out the nightly build of cimr at your own risk, 
clone the repository from git.


>>> git clone https://github.com/greenelab/cimr.git
>>> cd cimr
>>> python3 setup.py build
>>> python3 setup.py install


=================
Analysis examples
=================

------------------------------------------------------------------------
Quality assurance and processing of association summary statistics files
------------------------------------------------------------------------

cimr contains various functionalities in 
`processor <https://cimr.readthedocs.io/cimr/processor>`_ 
for processing summary statistics files for downstream analysis.


====================
Contributing to cimr
====================

-----------------
Contributing data
-----------------

You may contribute summary statistics from GWAS, eQTL and other similar studies. 
cimr currently expects hg20/GRCh38 reference for genomic position mapping.
However, variants mapped to hg19/GRCh37 may be used if updated using the
following command:


>>> cimr processor --datatype {datatype} --filename {filename} --update-map


Following columns are expected for association summary statistics files::

  gene : gene id in ensembl format
  rsnum : rs id of the variant
  constant_id : chromosome\_position\_referenceallele\_alternateallele\_genomebuild 
  e.g. chr2_128747549_G_T_hg19
  inc_allele : allele with respect to which variant's effect sizes are estimated
  inc_afrq : allele frequency of inc_allele
  beta : beta coefficient estimate for the association effect of the variant 
  se : standard error of the beta
  pval : p-value of the beta estimate



Here is an eQTL input file example::

  gene_id rsnum constant_id inc_allele inc_afrq beta se pval  
  GPR17 rs17262104 chr2_128747549_G_T G 0.06457 0.73698 0.11432743 5.5415e-10



----------------------------------
Contributing to cimr python module
----------------------------------


Contribute to cimr-d code or resources `here <https://github.com/greenelab/cimr-d>`_ .
Guidelines are provided `here <https://github.com/greenelab/cimr-d/CONTRIBUTING.md>`_ .

Contribute to cimr code or resources `here <https://github.com/greenelab/cimr>`_ .
Guidelines are provided `here <https://github.com/greenelab/cimr/CONTRIBUTING.md>`_ .

