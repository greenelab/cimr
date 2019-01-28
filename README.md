
# cimr
## continuous integration for mendelian randomization
## YoSon Park


cimr is a convenience tool for continuous analyses of variant-based 
association results from GWAS (genome-wide association studies), eQTL 
(expression-quantitative trait loci mapping) or other association studies. 
cimr began as a python module to run large-scale Mendelian randomization 
analysis (hence the name). As the project developed, it became more 
evident that there are many parts preceding the analyses that require 
pipelining. So the current incarnation of cimr aims to streamline the 
pre-analysis processing steps convenient, providing standardized input 
files and example scripts to run various downstream methods seamlessly.


*cimr is currently not released for public use.*



<!--ts-->

* [Installation](#installation)
  * [Installing python](#installing-python)
  * [Installing git lfs](#installing-git-lfs)
  * [Installing cimr](#installing-cimr)

* [Analyses](#Analyses)

* [cimr python module](#cimr-python-module)

* [Contributing to cimr](#contributing-to-cimr)
  * [Contributing data](#contributing-data)
  * [Contributing to cimr python module](#contributing-to-cimr-python-module)

<!--te-->


cimr is a convenience tool to process, integrate and maintain 
data in [cimr-d](https://github.com/greenelab/cimr-d).


## Installation

Python >= 3.6 is used for cimr. Additionally, in order to use cimr integrator 
functions or cimr-d, you may need git lfs installed.


### Installing python

cimr requires python >= 3.6. I recommend you have either miniconda ([download page](https://conda.io/miniconda.html)) 
or anaconda ([download page](https://www.anaconda.com/download/)) installed. However, all required python 
packages will be downloaded and installed as you build cimr using setup.py 
or requirements.txt provided here.


### Installing git lfs

cimr-d and some functionalities in cimr may use [git large file storage (LFS)](https://git-lfs.github.com/). 
See how to install git [here](https://www.atlassian.com/git/tutorials/install-git). 

To install git-lfs on Ubuntu, run:


```
curl -s https://packagecloud.io/install/repositories/github/git-lfs/script.deb.sh | sudo bash
sudo apt-get install -y git git-lfs
git-lfs install
```

Alternatively, you can install git-lfs through conda:

```
conda install -c conda-forge git-lfs && git lfs install
```

### Installing cimr

You can use pip to install the latest stable release of cimr.

```
pip3 install cimr
```

If you want to try out the nightly build of cimr at your own risk, 
clone the repository from git.

```
git clone https://github.com/greenelab/cimr.git
cd cimr
python3 setup.py build
python3 setup.py install
```


## Analysis examples

### Quality assurance and processing of association summary statistics files

cimr contains various functionalities in [processor](cimr/processor/README_processor.md) for
processing summary statistics files for downstream analysis.


## Contributing to cimr

### Contributing data

You may contribute summary statistics from GWAS, eQTL and other similar studies. 
cimr currently expects hg20/GRCh38 reference for genomic position mapping.
However, variants mapped to hg19/GRCh37 may be used if updated using the
following command:

```
cimr processor --datatype {datatype} --filename {filename} --update-map
```

Following columns are expected for association summary statistics files

```
gene : gene id in ensembl format
rsnum : rs id of the variant
constant_id : chromosome\_position\_referenceallele\_alternateallele\_genomebuild
e.g. chr2_128747549_G_T_hg19
inc_allele : allele with respect to which variant's effect sizes are estimated
inc_afrq : allele frequency of inc_allele
beta : beta coefficient estimate for the association effect of the variant on the gene 
se : standard error of the beta
pval : p-value of the beta estimate
```

Here is an eQTL input file example   

```
gene_id rsnum constant_id inc_allele inc_afrq beta se pval  
GPR17 rs17262104 chr2_128747549_G_T G 0.06457 0.73698 0.11432743 5.5415e-10
```

GWAS input files are not required to have gene_id column. Here is a GWAS input file example  

```
rsnum constant_id inc_allele inc_afrq beta se pval  
rs1172982 chr1_100230111_T_C T 0.3219 0.0043 0.0055 0.4689  
```

### Contributing to cimr python module

Read this [documentation regarding contributions](./CONTRIBUTING.md) for more information on contributing to cimr.
Read this [documentation regarding contributions](https://github.com/greenelab/cimr-d/CONTRIBUTING.md) for information
regarding contributions to cimr-d.


