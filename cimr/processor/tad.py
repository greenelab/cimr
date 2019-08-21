

"""Utilities for checking input file tagged as tad.
Topologically associated domain (tad) annotations are obtained from
publicly available resources and are used for various downstream
analyses.

(c) YoSon Park
"""

import os
import csv
import sys
import pandas
import logging
import datetime

from .utils import check_numeric
from .utils import set_chrom_dict


class Tadpole:
    """Process TAD annotations.

    Parameters
    ==========

    file_name: annotation file name

    Expects the following information in the data file

        chrom: chromosome id in the format of chr+number
        start: starting position of the tad block
        end: ending position of the tad block

        e.g.
        chr1	1960001	2400000

    Following columns can be provided with appropriate cimr arguements

        study_id: unique ID of the study, preferably a GEO unique identifier
        pub_id: DOI of the article published / preprinted using the data
                this ID will be listed as a part of the list of recommended
                citations when cimr-d is used for any other studies
        species: species. Default is homo_sapiens (9606, human)
        cell_type: cell or tissue type wherein TAD information was measured
        data_type: data type. Default is tad.

    Notes
    =====

    """

    HEADERS = [
        'chrom',
        'start',
        'end',
        'feature_id'
    ]

    METAHEADERS = [
        'file_name',
        'study_id',
        'pub_id',
        'species',
        'cell_type',
        'data_type',
        'protocol',
        'entry_count',
        'processed_date'
    ]

    METADATA_URL = 'https://raw.githubusercontent.com/greenelab/cimr-d/master/metadata/'


    def __init__(self,
                 file_name,
                 study_id,
                 pub_id,
                 species,
                 cell_type,
                 data_type,
                 protocol
                 ):
        self.file_name = file_name
        self.template = pandas.DataFrame(columns=self.HEADERS)
        self.study_id = study_id
        self.pub_id = pub_id
        self.species = species
        self.cell_type = cell_type
        self.data_type = data_type
        self.protocol = protocol
        self.processed_date = datetime.date.today().strftime('%Y-%m-%d')
        self.metafile = self.METADATA_URL + data_type + '.txt'


    def open_metadata(self):
        """Load or create a metadata for the new submission."""
        try:
            self.metadata = pandas.read_csv(self.metafile, sep='\t', header=0)
            logging.info(f' metadata is loaded.')
        except:
            logging.info(f' metadata is not available for the current data_type.')
            logging.info(f' new metadata is created.')
            self.metadata = pandas.DataFrame(columns=self.METAHEADERS)


    def write_metadata(self, tempdir='.'):
        """Write updated metadata file for a new PR.
        Default action assumes a function call within .travis.yml and a cloned
        cim-r working directory.
        """
        try:
            self.metadata.to_csv(
                tempdir+'/'+self.data_type+'.txt',
                sep='\t',
                header=True,
                index=False,
                na_reps='NA'
            )
        except:
            logging.error(f' not able to write the metadata.')
            sys.exit()


    def update_metadata(self):
        """Update metadata for the data_type with the new submission."""
        self.open_metadata()
        submitted = {
            'file_name' : self.file_name,
            'study_id' : self.study_id,
            'pub_id' : self.pub_id,
            'species' : self.species,
            'cell_type' : self.cell_type,
            'protocol' : self.protocol,
            'entry_count' : self.entry_count,
            'processed_date' : self.processed_date
        }
        self.metadata.append(submitted, ignore_index=True)


    def read_file(self):
        """Read the annotation file to integrate into cimr-d."""
        submitted_data = pandas.read_csv(
            self.file_name,
            sep='\t',
            header=None,
            names=['chrom', 'start', 'end'],
            dtype={'chrom' : str}
        )

        submitted_data.drop_duplicates(inplace=True)
        self.entry_count = len(submitted_data)

        check_numeric(submitted_data, 'start')
        check_numeric(submitted_data, 'end')

        submitted_data['feature_id'] = submitted_data.apply(
            lambda row: str(row.chrom) + ":" + str(row.start) + "-" + str(row.end),
            axis=1
        )

        chrom_dict, _ = set_chrom_dict()
        submitted_data['chrom'] = submitted_data['chrom'].map(
            chrom_dict, na_action='ignore'
        ).fillna(submitted_data['chrom'])

        self.processed_data = self.template.append(
            submitted_data,
            ignore_index=True
        )[self.template.columns]

        self.update_metadata()

    def write_file(self):
        """If checks pass, write the processed file into its new dir"""
        if str(self.file_name).startswith('cimr-d/submitted_data'):
            processed_file = str(self.file_name).replace('submitted_data', 'processed_data')
        else:
            processed_file = str(self.file_name) + '.cimr.txt'

        self.processed_data.to_csv(processed_file, sep='\t', header=False, index=False)


