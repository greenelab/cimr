

****
cimr
****


***************************************
cimr is not yet released for public use
***************************************

=====================================================
continuous integration and analysis of complex traits
=====================================================

==========
YoSon Park
==========

**Useful links**:
`Source repository <https://github.com/greenelab/cimr>`_ |
`Issues & Ideas <https://github.com/greenelab/cimr/issues>`_ |
`Documentation <https://cimr.readthedocs.io>`_ |
`cimr-d <https://github.com/greenelab/cimr-d>`_


*cimr* (continuously integrated meta-resource) is a convenience tool
for continuous analyses of variant-based association results from
GWAS (genome-wide association studies), eQTL (expression-quantitative
trait loci mapping) or other association studies. cimr aims to
streamline the pre-analysis processing steps, provide standardized
input files and automate scripting for standard downstream analyses.



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
See how to install
`git <https://www.atlassian.com/git/tutorials/install-git>`_ .


To install git-lfs on Ubuntu, run::

    curl -s https://packagecloud.io/install/repositories/github/git-lfs/script.deb.sh | sudo bash
    sudo apt-get install -y git git-lfs
    git-lfs install


Alternatively, you can install git-lfs through conda::

    conda install -c conda-forge git-lfs && git lfs install


---------------
Installing cimr
---------------

You can use pip to install the latest stable release of cimr::

    pip3 install cimr


If you want to try out the nightly build of cimr at your own risk,
clone the repository from git::

    git clone https://github.com/greenelab/cimr.git
    cd cimr
    python3 setup.py build
    python3 setup.py install



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


Contribute to cimr-d code or resources `here <https://github.com/greenelab/cimr-d>`_ .
Guidelines are provided `here <https://github.com/greenelab/cimr-d/CONTRIBUTING.md>`_ .

Contribute to cimr code or resources `here <https://github.com/greenelab/cimr>`_ .
Guidelines are provided `here <https://github.com/greenelab/cimr/CONTRIBUTING.md>`_ .

