
# cimr
## YoSon Park


continuous analyses of gwas, eqtl or other summary statistics

<!--ts-->

* [installation](#installation)
  * [install miniconda/anaconda](#install-miniconda-or-anaconda)
* [install git lfs](#install-git-lfs)
* [install cimr](#install-cimr)
* [contributing to cimr](#contributing-to-cimr)
  * [summary statistics](#summary-statistics)
* [analyses](#analyses)
* [cimr python module](#cimr-python-module)

<!--te-->

## installation

### install miniconda or anaconda
cimr requires python3 and conda to manage dependencies. we recommend you have either miniconda ([download page](https://conda.io/miniconda.html)) or anaconda ([download page](https://www.anaconda.com/download/)) installed. 

## install git lfs

cimr uses git and [git large file storage (LFS)](https://git-lfs.github.com/). See how to install git [here](https://www.atlassian.com/git/tutorials/install-git). To install git-lfs on Ubuntu, run:

```
curl -s https://packagecloud.io/install/repositories/github/git-lfs/script.deb.sh | sudo bash
sudo apt-get install -y git git-lfs
git-lfs install
```

alternatively, you can install git-lfs through conda:

```
conda install -c conda-forge git-lfs && git lfs install
```

## install cimr


### cimr as a python package

you may clone the repository from git and using as a python module:

```
git clone https://github.com/ypar/cimr.git
cd cimr
python3 setup.py build
python3 setup.py install
```


## contributing to cimr

### contributing summary statistics

you may contribute summary statistics from gwas, eqtl and other association studies. cimr currently expects hg20/grch38 reference for genomic position mapping.


following columns are expected for association summary statistics files

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

eqtl input file example   

```
gene rsnum constant_id inc_allele inc_afrq beta se pval  
GPR17 rs17262104 chr2_128747549_G_T G 0.06456953 0.736983583560685 0.11432743 5.541546921310881e-10  
HAX1 rs12749691 chr1_154251259_T_A T 0.27152318 0.280817746771117 0.05622703 1.08387876813775e-06  
```
> gwas input file example  

```
rsnum constant_id inc_allele inc_afrq beta se pval  
rs1172982 chr1_100230111_T_C T 0.3219 0.0043 0.0055 0.4689  
rs1172981 chr1_100230197_T_C T 0.06069 0.0057 0.0103 0.7688  
```

### contributing to twas methods

you may contribute to downstream analyses using summary statistics. currently incorporated methods are as follows:


### cimr python module




