

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

You can use cimr to standardize public datasets using a yaml file, e.g.::

    # example.yaml

    data_file:
        description: >-
            Global Lipid Genetics Consortium GWAS results for high-density
            cholesterol levels
        location:
            url: https://zenodo.org/record/3338180/files/HDL_Cholesterol.txt.gz
            md5: 2b28816a0a363db1a09ad9a6ba1a6620
        columns:
            variant_id: panel_variant_id
            variant_chrom: chromosome
            variant_pos: position
            rsnum: variant_id

    data_info:
        citation: 10.1038/ng.2797
        data_source: http://lipidgenetics.org/
        data_type: gwas
        context: hdl cholesterol
        build: b38
        sample_size: 187167
        n_cases: na
        can_be_public: true

    method:
        name: linear regression
        tool: PLINK;SNPTEST;EMMAX;Merlin;GENABEL;MMAP
        website: >-
            http://zzz.bwh.harvard.edu/plink/download.shtml;
            https://mathgen.stats.ox.ac.uk/genetics_software/snptest/snptest.html;
            https://genome.sph.umich.edu/wiki/EMMAX;
            https://csg.sph.umich.edu/abecasis/Merlin/tour/assoc.html;
            http://www.genabel.org/sites/default/files/html_for_import/GenABEL_tutorial_html/GenABEL-tutorial.html;
            https://mmap.github.io/

    contributor:
        name: Contributor Name
        github: contributorgithub
        email: contributoremail@emaildomain.emailextension



Details can be found in the
`cimr-d contributions.md <https://github.com/greenelab/cimr-d/blob/master/doc/contributing.md>`_.


Once the yaml file is prepared, you can run cimr locally::

    cimr processor -process -yaml-file example.yaml


