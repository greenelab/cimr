

"""Utility functions related to annotating the input files or checking 
user-provided annotations such as gene IDs, gene names, etc.

(c) YoSon Park
"""


import sys
import json
import logging
import requests


class Querier:
    """The base class for annotation checks and remapping functions.

    Parameter
    ---------

    headers : headers used for the requests. default provided by mygene is
              {'content-type':'application/x-www-form-urlencoded'}

    url : the current default is v3 and url is set as below
          'https://mygene.info/v3/query'
          v2 is provided below (outdated, but not deprecated as of Jan 2019)
          'https://mygene.info/v2/query'

    species : the default is 9606 or human for cimr analyses.
              for general mygene usage, the doc is available in the below site
              docs.mygene.info/en/latest/doc/data.html#species

    scopes : the search space of the queired terms can be listed here.
             default setting includes that the query involves a gene with
             identifiers in the below described scopes
             symbol - official gene symbol (e.g. APOE)
             entrezgene - entrez gene ID
             ensembl.gene - Ensembl gene ID starting with ENSG
             ensembl.transcript - Ensembl transcript ID starting with ENST

    fields : the below table is from the mygene documentation indicated below
             docs.mygene.info/en/latest/doc/query_service.html#available-fields
             I have reordered the table alphabetically and shortedned some
             descriptions.
             The fields and scopes parameters share the same variable names.

    field             | description         | examples
    ------------------------------------------------------------------------------
    accession         | NCBI GeneBank number| accession:AA810989
    alias             | Gene alias          | alias:p33
    ensembl.gene      | Ensembl gene ID     | ensembl.gene:ENSG00000123374
    ensembl.protein   | Ensembl protein ID  | ensembl.protein:ENSP00000243067
    ensembl.transcript| Ensembl transcriptID| ensembl.transcript:ENST00000266970
    entrezgene        | Entrez gene ID      | entrezgene:1017   
    flybase           | Drosophila gen(om)es| flybase:FBgn0004107&species=fruitfly
    go                | Gene Ontology ID    | go:0000307
    hgnc              | HUGO                | hgnc:1771
    homologene        | NCBI HomoloGene ID  | homologene:74409
    hprd              | Human Protein Ref DB| hprd:00310
    interpro          | InterPro ID         | interpro:IPR008351
    mgi               | Mouse Genome Inf    | mgi:MGI\\:88339
    mim               | OMIM ID             | mim:116953
    mirbase           | miRNA DB            | mirbase:MI0017267
    name              | Gene name           | name:cyclin-dependent
    pdb               | PDB ID              | pdb:1AQ1
    pfam              | PFam ID             | pfam:PF00069
    pharmkb           | PharmGKB ID         | pharmgkb:PA101
    prosite           | Prosite ID          | prosite:PS50011
    reagent           | GNF reagent ID      | reagent:GNF282834
    refseq            | NCBI RefSeq ID      | refseq:NP_439892
    reporter          | Affymetrix probeset | reporter:204252_at
    retired           | Retired EntrezGeneID| retired:84999
    rgd               | Rat Genome Database | rgd:620620
    summary           | Gene summary text   | summary:insulin
    symbol            | Official gene symbol| symbol:cdk2
    tair              | Arabidopsis Info Res| tair:AT3G48750&species=thale-cress
    unigene           | NCBI UniGene ID     | unigene:Hs.19192
    uniprot           | UniProt ID          | uniprot:P24941
    wormbase          | C.elgans+nematode DB| wormbase:WBGene00057218&species=31234
    xenbase           | Xenopus l. + t. DB  | xenbase:XB-GENE-1001990&species=frog
    zfin              | Zebrafish DB   | zfin:ZDB-GENE-980526-104&species=zebrafish


    Notes
    -----

    mygene.info website documentation can be found here:
    docs.mygene.info/en/latest/index.html

    mygene.info is a website for gene annotation queries using the biothings
    backend. biothings documentation is here
    biothingsapi.readthedocs.io/en/latest

    cimr currently only utilizes the batch query function using direct queries
    to the website api

    """

    FIELDS = {'accession', 'alias', 'ensembl.gene', 'ensembl.protein', 
              'ensembl.transcript', 'entrezgene', 'flybase', 'go', 'hgnc',
              'homologene', 'hprd', 'interpro', 'mgi', 'mim', 'mirbase', 
              'name', 'pdb', 'pfam', 'pharmgkb', 'prosite', 'reagent',
              'refseq', 'reporter', 'retired', 'rgd', 'summary', 'symbol',
              'tair', 'unigene', 'uniprot', 'wormbase', 'xenbase', 'zfin'
              }


    def __init__(self, 
                 genes,
                 headers = {'content-type':'application/x-www-form-urlencoded'}, 
                 url = 'https://mygene.info/v3/query', 
                 species = 9606,
                 scopes = 'symbol+entrezgene+ensembl.gene+ensembl.transcript',
                 fields = 'name+symbol+taxid+entrezgene+ensembl+alias+refseq'
                 ):
        """Initialize the Annotator"""
        self.genes = genes
        self.headers = headers
        self.url = url
        self.species = species
        self.scopes = scopes
        self.fields = fields


    @staticmethod
    def subtract_plus(joined_fields):
        """Split the feature list by the delimiter to run per-element checks"""
        return joined_fields.split('+')
    

    @staticmethod
    def add_plus(split_fields):
        """Insert delimiters to convert the list (or pandas Series, etc.)
        into a suitable format to be used
        Accepted delimiters include a comma (,), a space ( ) or a plus (+).
        The default set by cimr is a plus (+).
        """
        return '+'.join(map(str, split_fields))

    
    def check_fields(self):
        """Check if all items indicated in the fields parameter are legit.
        The list has been updated based on mygene.info documentation on
        January 2019.
        """
        splitted = self.subtract_plus(self.fields)
        for field in splitted:
            if field not in self.FIELDS:
                raise ValueError(' %s is not a valid field supported in cimr.' % field)
            else:
                splitted.remove(field)
        self.fields = self.add_plus(splitted)


    @classmethod
    def from_csv(cls, features):
        if len(features) == 1:
            feature_list = features.split(',')
            queried = cls(feature_list).genes_for_query(feature_list)
            return queried
        else:
            raise ValueError(' features appear to be already split')        


    def form_query(self):
        """Form a query based on parameters. cimr uses the batch query of a gene
        list to the website api of the database.
        By default, only gene_ids in the format of symbol, entrez id, or
        ensembl (gene/transcript) are accepted but other options may be indicated
        by scopes option if cimr AnnotateGene is used without cimr-d integration.
        """
        params = 'q=' + str(self.genes) +\
                 '&scopes=' + str(self.scopes) +\
                 '&fields=' + str(self.fields) + \
                 '&species=' + str(self.species)
        print(self.url)
        print(self.headers)
        print(params)
        queried = requests.post(self.url, headers = self.headers, data = params)
        return queried


        


