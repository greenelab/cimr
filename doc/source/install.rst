
************
Installation
************

===============================
Installing Python and git(-lfs)
===============================

----------------------
Python version support
----------------------

Officially Python 3.6 and 3.7.

------------------
Installing git lfs
------------------

cimd-d and cimr integrator functionalities may use git large file
storage (LFS).
Install `git <https://www.atlassian.com/git/tutorials/install-git>`_
and
`git-lfs <https://git-lfs.github.com/>`_ to use these functionalities.

To install git-lfs on Ubuntu, run::

    curl -s https://packagecloud.io/install/repositories/github/git-lfs/script.deb.sh
    sudo apt-get install -y git git-lfs
    git-lfs install

Alternatively, conda can be used to install git-lfs::

    conda install -c conda-forge git-lfs && git lfs install



=================
Installing cimr
=================

cimr requires Python :math:`\ge` 3.6.
`Miniconda <https://conda.io/miniconda.html>`_ or
`Anaconda <https://www.anaconda.com/download/>`_ are convenient
options for scientific computing with Python. However, if a basic
Python version has been installed, requirements.txt and setup.py
files can be used to download specific packages required for cimr
functions.

cimr can be installed using pip3::

    pip3 install cimr


Or by cloning the git repository for the nightly builds::

    git clone https://github.com/greenelab/cimr.git
    cd cimr
    python3 setup.py build
    python3 setup.py install

