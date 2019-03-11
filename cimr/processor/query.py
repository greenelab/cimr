

"""Utility functions related to annotating the input files or checking 
user-provided annotations such as gene IDs, gene names, etc.

(c) YoSon Park
"""


import sys
import json
import pandas
import logging
import requests


def subtract_plus(joined_fields):
    """Split the feature list by the delimiter to run per-element checks"""
    return joined_fields.split('+')


def add_plus(split_fields):
    """Insert delimiters to convert the list (or pandas Series, etc.)
    into a suitable format to be used
    Accepted delimiters include a comma (,), a space ( ) or a plus (+).
    The default set by cimr is a plus (+).
    """
    return '+'.join(map(str, split_fields))


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
              https://docs.mygene.info/en/latest/doc/data.html#species

    scopes : the search space of the queried terms can be listed here.
             default setting includes that the query involves a gene with
             identifiers in the below described scopes
             symbol - official gene symbol (e.g. APOE)
             entrezgene - entrez gene ID
             ensembl.gene - Ensembl gene ID starting with ENSG
             ensembl.transcript - Ensembl transcript ID starting with ENST

    fields : the below table is from the mygene documentation indicated below
             https://docs.mygene.info/en/latest/doc/query_service.html#available-fields
             I have reordered the table alphabetically and shortened some
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
    zfin              | Zebrafish DB        | zfin:ZDB-GENE-980526-104&species=zebrafish

    Methods
    -------

    

    Notes
    -----

    https://mygene.info website documentation can be found here:
    https://docs.mygene.info/en/latest/index.html

    https://mygene.info is a website for gene annotation queries using the biothings
    backend. biothings documentation is here
    https://biothingsapi.readthedocs.io/en/latest

    cimr currently only utilizes the batch query function using direct queries
    to the website api

    """

    FIELDS = {
        'accession', 'alias', 'ensembl.gene', 'ensembl.protein', 
        'ensembl.transcript', 'entrezgene', 'flybase', 'go', 'hgnc',
        'homologene', 'hprd', 'interpro', 'mgi', 'mim', 'mirbase', 
        'name', 'pdb', 'pfam', 'pharmgkb', 'prosite', 'reagent',
        'refseq', 'reporter', 'retired', 'rgd', 'summary', 'symbol',
        'tair', 'unigene', 'uniprot', 'wormbase', 'xenbase', 'zfin'
    }


    def __init__(self, 
                 genes,
                 headers={'content-type':'application/x-www-form-urlencoded'}, 
                 url='https://mygene.info/v3/query', 
                 species=9606,
                 scopes='alias+symbol+entrezgene+ensembl.gene+ensembl.transcript',
                 fields='name+symbol+taxid+entrezgene+ensembl+alias+refseq'
                 ):
        """Initialize the Querier"""
        self.genes = genes
        self.headers = headers
        self.url = url
        self.species = species
        self.scopes = scopes
        self.fields = fields


    def make_string(self):
        """Take the gene_id column input from the data table and turn it
        into a string to be used in form_query
        """
        if type(self.genes) is list:
            logging.info(f' converting a gene list into a queriable string.')
            self.genestring = add_plus(self.genes)
        else:
            logging.error(f' is the input a list with gene_id\'s?')

    
    def check_fields(self):
        """Check if all items indicated in the fields parameter are legit.
        The list has been updated based on mygene.info documentation on
        January 2019.
        """
        splitted = subtract_plus(self.fields)
        for field in splitted:
            if field not in self.FIELDS:
                raise ValueError(' %s is not a valid field supported in cimr.' % field)
            else:
                splitted.remove(field)
        self.fields = add_plus(splitted)
    

    def check_ensembl(self):
        """Commonly used Ensembl ID may include a version identifier. Check for any
        Ensembl IDs in the feature list and remove the version identifier.
        """
        if any('ENSG' in g for g in self.genes):
            self.gene_db = 'ensembl.gene'
        elif any('ENST' in g for g in self.genes):
            self.gene_db = 'ensembl.transcript'
        elif any('ENSP' in g for g in self.genes):
            self.gene_db = 'ensembl.protein'
        else:
            logging.info(f' no human ensembl ID found in the gene_id column.')

        if not self.gene_db in subtract_plus(self.scopes):
            self.scopes = str(self.scopes) + '+' + str(self.gene_db)
        

    @classmethod
    def from_csv(cls, features):
        if len(features) == 1:
            feature_list = features.split(',')
            queried = cls(feature_list)
            return queried
        else:
            raise ValueError(' features appear to be already split.')        


    def form_query(self):
        """Form a query based on parameters. cimr uses the batch query of a gene
        list to the website api of the database.
        By default, only gene_ids in the format of symbol, entrez id, or
        ensembl (gene/transcript) are accepted but other options may be indicated
        by scopes option if cimr Querier is used without cimr-d integration.
        """
        try:
            self.make_string()
        except:
            raise ValueError(' An error occurred while converting the gene list.')

        params = ('q=' + str(self.genestring) +
                 '&scopes=' + str(self.scopes) +
                 '&fields=' + str(self.fields) + 
                 '&species=' + str(self.species)
        )

        self.queried = requests.post(self.url, headers=self.headers, data=params)
        self.jsoned = json.loads(self.queried.text)


    def write_json(self, annot_file_out):
        """Write the results of the query into a file.
        The json dump indent default is set as 4.
        """
        try:
            with open(annot_file_out, 'w') as outfile:
                json.dump(self.jsoned, outfile, indent=4)
        except:
            raise ValueError(' an error occurred while writing the queried gene list.')


    def write_gene(self, annot_file_out):
        """Write official gene symbol, entrez IDs and ensembl gene IDs 
        for each gene_id.
        """
        with open(annot_file_out, 'w') as outfile:

            try:    
                for gene in range(0, len(self.jsoned)):

                    try:
                        symbol = self.jsoned[gene]['symbol']
                    except KeyError:
                        logging.info(f' %s does not have an official gene symbol' % self.jsoned[gene]['query'])
                        symbol = 'NA'

                    try:
                        entrez = self.jsoned[gene]['entrezgene']
                    except KeyError:
                        logging.info(f' %s does not have an entrez ID.' % self.jsoned[gene]['query'])
                        entrez = 'NA'

                    try:
                        ensembl = self.jsoned[gene]['ensembl']['gene']
                    except KeyError:
                        logging.info(f' %s does not have an ensembl ID.' % self.jsoned[gene]['query'])
                        ensembl = 'NA'

                    except TypeError:
                        ensembl = self.jsoned[gene]['ensembl'][0]['gene']
                    except:
                        logging.info(f' %s has a type error' % self.jsoned[gene]['query'])
                        ensembl = 'NA'
                    
                    try: 
                        alias = self.jsoned[gene]['alias']
                    except KeyError:
                        logging.info(f' %s does not have an alias.' % self.jsoned[gene]['query'])
                        alias = 'NA'

                    print(symbol, entrez, ensembl, alias, file=outfile)

            except:
                print(self.jsoned[gene]["query"])
                raise ValueError(' an error occurred while writing the parsed gene list')


class Snpper:
    """Query refsnp variation database directly to update RS IDs of SNPs.
    
    Returns
    -------
    Either the original (if newest) or updated RS ID of the SNP along with
    its position in the GRCh38.p12 reference genome.

    Notes
    -----
    Primarily queries the refsnp variation database for the latest RS IDs of
    GRCh38 reference genome.

    """

    REFSNP_URL = 'https://api.ncbi.nlm.nih.gov/variation/v0/beta/refsnp/'

    def __init__(self, rsid):
        """Initializing Snpper class"""
        self.rsid = rsid


    def update_rsid(self):
        """Find updated mapping of a SNP RS ID"""
        rsid = str(self.rsid).replace('rs','')
        response = requests.get(self.REFSNP_URL + rsid)
        json_data = json.loads(response.text)

        if 'merged_snapshot_data' in json_data.keys():
            updated_rsid = json_data['merged_snapshot_data']['merged_into'][0]
        else:
            updated_rsid = json_data['refsnp_id']

        genome_build = json_data['primary_snapshot_data']['placements_with_allele'][0]['placement_annot']['seq_id_traits_by_assembly'][0]['assembly_name']

        if genome_build == 'GRCh38.p12':            
            refseq_chrom = json_data['primary_snapshot_data']['placements_with_allele'][0]['alleles'][0]['allele']['spdi']['seq_id']
            pos = json_data['primary_snapshot_data']['placements_with_allele'][0]['alleles'][0]['allele']['spdi']['position']
            refseq_mapping_file = '/work/drug/gwas/gwas_catalog/chromosome_to_refseq_mapping.txt'
            refseq_mapping = pandas.read_csv(refseq_mapping_file, sep='\t', header=0, na_values=['NA'])
            chrom = refseq_mapping[refseq_mapping['RefSeq_sequence'] == refseq_chrom]['Molecule_name'].values[0]
        else:
            logging.error('')

        return 'rs' + str(updated_rsid), str(chrom) + ':' + str(pos)



