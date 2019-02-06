

"""Utilities for checking input file tagged as tad.
Topologically associated domain (tad) annotations are obtained from
publicly available resources and are used for various downstream
analyses.

(c) YoSon Park
"""

import pandas
import logging

from .utils import check_numeric


class Tadpole:
    """Process TAD annotations.

    Parameters
    ==========

    file_name: annotation file name

    Expects the following information in the data file

        chrom: chromosome id in the format of chr+number
        start: starting position of the tad block
        end: ending position of the tad block
        (optional) id: unique id of the tad block
        (optional) comment: optional column containing additional information

        e.g.
        chr1	1960001	2400000	tad0|hg19|chr1:1960001-2400000	1000

    cell_type: cell or tissue type wherein TAD information was measured
    species: species. default is 9606 (homo sapiens; human)

    Notes
    =====

    """

    HEADERS = ['chrom', 'start', 'end', 'id', 'comment', 'species', 'cell_type', 'pmid']

    def __init__(self, file_name, species, cell_type):
        self.file_name = file_name
        self.template = pandas.DataFrame(columns=self.HEADERS)
        self.species = species
        self.cell_type = cell_type


    def read_file(self):
        """Read the annotation file to integrate into cimr-d."""
        input_annot = pandas.read_csv(self.file_name, sep='\t', header=0)

        input_annot.drop_duplicates(inplace=True)
        input_annot['species'] = str(self.species)
        input_annot['cell_type'] = str(self.cell_type)
        check_numeric(input_annot, 'start')
        check_numeric(input_annot, 'end')

        null_id = pandas.isnull(input_annot['id'])
        chroms = input_annot['chrom'].astype(str)
        starts = input_annot['start'].astype(str)
        ends = input_annot['end'].astype(str)
        input_annot.loc[(null_id), 'id'] = chroms + ':' + starts + '-' + ends
                                                        
        self.annot = self.template.append(input_annot, ignore_index=True)[self.template.columns]

        




