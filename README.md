
# cimr
## continuous integration for mendelian randomization
## YoSon Park


cimr is a convenience tool for continuous analyses of variant-based 
association results from GWAS (genome-wide association studies), eQTL 
(expression-quantitative trait loci mapping) or other related studies. 
cimr began as a python module to run large-scale Mendelian randomization 
analysis (hence the name). As the project developed, it became more 
evident that there are many parts preceding the analyses that require 
pipelining. So the current incarnation of cimr aims to streamline the 
pre-analysis processing steps convenient, providing standardized input 
files and example scripts to run various downstream methods seamlessly.


\* cimr is currently not released for public use.
\* cimr is a convenience tool to process, integrate and maintain 
data in [cimr-d](https://github.com/greenelab/cimr-d).


<!--ts-->

* [Installation](#installation)
  * [Installing python via miniconda/anaconda](#installing-python-via-miniconda-or-anaconda)
  * [Installing git lfs](#installing-git-lfs)
  * [Installing cimr](#installing-cimr)

* [Analyses](#Analyses)

* [cimr python module](#cimr-python-module)

* [Contributing to cimr](#contributing-to-cimr)
  * [Contributing data](#contributing-data)
  * [Contributing to cimr python module](#contributing-to-cimr-python-module)

<!--te-->


## Installation

Python >= 3.6 is used for cimr. Additionally, in order to use cimr integrator 
functions or cimr-d, you may need git lfs installed.


### Installing python via miniconda or anaconda

cimr requires python >= 3.6. I recommend you have either miniconda ([download page](https://conda.io/miniconda.html)) 
or anaconda ([download page](https://www.anaconda.com/download/)) installed. However, all required python 
packages will be downloaded and install as you build cimr using setup.py 
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

If you want to try out the nightly build of cimr at your own risk, you may clone the repository from git.

```
git clone https://github.com/greenelab/cimr.git
cd cimr
python3 setup.py build
python3 setup.py install
```


## analyses



### 

[processor](cimr/processor/README_processor.md)


## Contributing to cimr

### Contributing data

You may contribute summary statistics from gwas, eqtl and other association studies. cimr currently expects hg20/grch38 reference for genomic position mapping.


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
gene rsnum constant_id inc_allele inc_afrq beta se pval  
GPR17 rs17262104 chr2_128747549_G_T G 0.06456953 0.736983583560685 0.11432743 5.541546921310881e-10  
HAX1 rs12749691 chr1_154251259_T_A T 0.27152318 0.280817746771117 0.05622703 1.08387876813775e-06  
```
Here is a GWAS input file example  

```
rsnum constant_id inc_allele inc_afrq beta se pval  
rs1172982 chr1_100230111_T_C T 0.3219 0.0043 0.0055 0.4689  
rs1172981 chr1_100230197_T_C T 0.06069 0.0057 0.0103 0.7688  
```

### Contributing to cimr python module

Read this [documentation regarding contributions](./CONTRIBUTIONG.md) for more information on contributing to cimr.
Read this [documentation regarding contributions](https://github.com/greenelab/cimr-d/CONTRIBUTING.md) for information
regarding contributions to cimr-d.


