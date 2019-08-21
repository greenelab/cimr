"""Reading and parsing through the contributor's yaml file(s).

By default within cimr-d ci processing, 'submitted' dir will be
scanned for file names ending with 'yaml' and 'yml'.

This is the default uploading skim for single and bulk files
using zenodo. For PR-based test file uploader,
check .circleci/deploy.sh
and .circleci/process_submitted_data.py

(c) YoSon Park
"""


import os
import sys
import yaml
import pandas
import pathlib
import logging

from pandas.api.types import is_numeric_dtype

from ..defaults import DATA_TYPES
from ..defaults import CONFIG_FILE_EXTENSION
from ..defaults import BULK_EXTENSION
from ..defaults import FILE_EXTENSION


def check_yaml_before_commit():
    """A git-status-dependent function used when locally applying
    parse_yaml.py. It searches for a new or modified yml/yaml file and
    returns its pathlib path.
    """
    import subprocess

    status_check = 'git status --porcelain'
    jobsplit = subprocess.check_output(
        status_check,
        stderr=subprocess.STDOUT,
        shell=True,
        universal_newlines=True
    ).replace('\n', '').split('?? ')

    for job in jobsplit:
        if job.endswith(CONFIG_FILE_EXTENSION):
            yaml_file = pathlib.Path('submitted/' + job.split('/')[-1])

    return yaml_file


def check_yaml_in_ci():
    """A git-status-dependent function used during ci processing.
    It searches for a new or modified yml/yaml file from a new pr.
    User-defined yaml files can be stored in the following dir:
    submitted/
    """
    import subprocess

    status_check = 'git diff origin --name-only'
    jobsplit = subprocess.check_output(
        status_check,
        stderr=subprocess.STDOUT,
        shell=True,
        universal_newlines=True
    ).split('\n')

    for job in jobsplit:
        if job.endswith(CONFIG_FILE_EXTENSION):
            yaml_file = pathlib.Path(
                'submitted/' + job.split('/')[-1]
            )

    return yaml_file


def find_yaml_in_dir():
    """Considering multiple yaml files in a submitted/ dir."""
    yaml_files = []
    yaml_dir = os.path.join(os.getcwd(), 'submitted/')

    for yaml_file in os.listdir(yaml_dir):
        yaml_file = os.path.join(yaml_dir, yaml_file)
        if yaml_file.endswith('.place_holder'):
            continue
        if yaml_file.endswith(CONFIG_FILE_EXTENSION):
            yaml_files.append(yaml_file)
        else:
            raise Exception(f' {yaml_file} is not an acceptible yaml file.')

    return yaml_files


def load_yaml(yaml_file):
    """Read the found yaml file and return the read object."""
    with open(yaml_file, 'r') as args:
        try:
            return yaml.safe_load(args)
        except yaml.YAMLError as exc:
            logging.error(f' {exc}')
            sys.exit(1)


def validate_data_type(data_type):
    """Validate data_type variable for cimr compatibility."""
    if data_type in DATA_TYPES:
        return True
    else:
        logging.error(f' check your data_type.')
        sys.exit(1)


def verify_weblink(path):
    """Verify the provided link to the contributed file."""
    import requests
    try:
        response = requests.get(path)
        response.raise_for_status()
        return True
    except Exception as e:
        print(e)
        return False


def trim_zenodo_link(path):
    """Trim a zenodo download link to extract file name.

    e.g. https://zenodo.org/record/3369410/files/gwas.txt.gz?download=1
    -> https://zenodo.org/record/3369410/files/gwas.txt.gz
    """
    path = path.replace('?download=1', '')
    return path


def download_gdrive_file(
    path,
    outdir,
    filename
    ):
    """Given a link starting with https://drive.google.com,
    initialize a google drive file download.

    need two arguments:
      * file id from the share link
      * file name (may be different from the original file name)
    """
    import os
    path = path.replace('https://drive.google.com/file/d/', '')
    path = path.replace('/view?usp=sharing','')

    with open('download_gdrive_file.sh', 'w') as downloader:
        downloader.write(r'''
        #!/bin/bash

        FILEID=$1
        FILENAME=$2

        CONFIRM=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate "https://docs.google.com/uc?export=download&id=$FILEID" -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1\n/p')

        wget --load-cookies /tmp/cookies.txt "https://docs.google.com/uc?export=download&confirm=$CONFIRM&id=$FILEID" -O $FILENAME

        rm -rf /tmp/cookies.txt

        '''
        )

    run_cmd = 'bash download_gdrive_file.sh ' + path + ' ' + outdir + filename
    os.system(run_cmd)
    os.system('rm download_gdrive_file.sh')


def download_file(path, outdir='./'):
    """Download data based on the provided link.

    Note;
    Progress bars added based on the following reference:
    https://stackoverflow.com/questions/37573483/progress-bar-while-down
    load-file-over-http-with-requests/37573701
    """
    from tqdm import tqdm
    import requests
    import math

    r = requests.get(path, stream=True)
    total_size = int(r.headers.get('content-length', 0))
    block_size = 1024
    wrote = 0
    file_name = path.split('/')[-1]
    file_path = outdir + file_name

    with open(file_path, 'wb') as f:
        for data in tqdm(r.iter_content(block_size),
                         total=math.ceil(total_size//block_size),
                         unit='KB',
                         leave=True,
                         ncols=72,
                         unit_scale=True,
                         unit_divisor=1024):
            wrote = wrote + len(data)
            f.write(data)

    if total_size != 0 and wrote != total_size:
        logging.error(f' check the file link and try again.')
        sys.exit(1)


def validate_hash(path, hash):
    """Validate a file against an MD5 hash value."""
    import hashlib

    md5 = hashlib.md5()

    with open(path, 'rb') as f:
        while True:
            chunk = f.read(10000000)
            if not chunk:
                break
            md5.update(chunk)

    return md5.hexdigest() == hash


def verify_dir(tarred_data):
    """Check directory tree of tarball containing multiple files"""
    for member in tarred_data.getmembers():
        # skip sub-directories
        if member.isdir():
            continue

        # No any non-file types (links, devices, FIFOs, etc) allowed
        if not member.isfile():
            logging.error(f' illegal file type: {member.name}')
            sys.exit(1)

        file_path = member.name
        # Leading '/' in file path NOT allowed
        if file_path.startswith('/'):
            logging.error(f' Leading / found in archived file.')
            sys.exit(1)

        # Remove leading './' from file path
        if file_path.startswith('./'):
            file_path = file_path.replace('./', '')

        # Ensure that a regular file's name in archive is always in the
        # format of "<data_type>/<filename>"
        tokens = file_path.split('/')
        if len(tokens) != 2:
            logging.error(f' Illegal file system hierarchy in archive: {member.name}')
            sys.exit(1)

        if tokens[0] not in DATA_TYPES:
            logging.error(f' Unknown data_type in archive: {tokens[0]}')
            sys.exit(1)


def convert_yaml(yaml_files):
    """Convert yaml parameters to cimr arguments"""
    fileset = []
    for yaml_file in yaml_files:
        yaml_file = pathlib.Path(yaml_file)
        if yaml_file.is_file():
            yaml_file_path = yaml_file.resolve(strict=True)
            logging.info(f' processing {yaml_file_path}')
        else:
            logging.error(f' {yaml_file} is not accessible.')
            sys.exit(1)

        yaml_data = load_yaml(yaml_file)
        y = Yamler(yaml_data)
        y.check_data_file()

        if hasattr(y, 'columnset'):
            columnset = y.columnset
        else:
            columnset = {}

        if hasattr(y, 'genome_build'):
            genome_build = y.genome_build

        fileset = [*fileset, *y.fileset]

    return genome_build, fileset, columnset


def standardize_context(context):
    """Standardizing the context description"""
    context = str(context).lower().replace(' ', '_')
    # context = context.split(';')
    return context


class Yamler:
    """A collection of utilities to parse the yaml file, check metadata
    and trigger cimr processing of the contributed file
    """
    def __init__(self, yaml_data):
        self.yaml_data = yaml_data
        self.data_type = None
        self.genome_build = None
        self.keys = None
        self.hash = None
        self.sub_datatype_dir = None
        self.downloaded_file = None


    def pick_keys(self):
        """List keys for the dictionarized yaml data."""
        self.keys = self.yaml_data.keys()


    def set_data_type(self):
        """Pull out data_type variable value from yaml"""
        try:
            data_type = self.yaml_data['data_info']['data_type']
            if validate_data_type(data_type):
                self.data_type = data_type
        except ValueError:
            logging.error(f' there is no data_type indicated.')
            sys.exit(1)


    def set_genome_build(self):
        """Pull out genome build variable value from yaml"""
        if 'build' in self.yaml_data['data_info']:
            self.genome_build = self.yaml_data['data_info']['build']
        else:
            logging.info(f' genome build is not provided in yaml')
            sys.exit(1)


    def check_hash(self):
        """Compare md5 of the downloaded file to the provided value"""
        if validate_hash(self.downloaded_file, self.hash):
            logging.info(f' md5sum verified.')
        else:
            raise ValueError(' provided md5sum didn\'t match.')


    def download(self):
        """Check if provided weblink to the file exists.
        Download if verified.
        """
        self.file_link = self.yaml_data['data_file']['location']['url']
        if 'zenodo.org' in str(self.file_link):
            self.file_link = trim_zenodo_link(self.file_link)

        if 'input_name' in self.yaml_data['data_file'].keys():
            self.infile = self.yaml_data['data_file']['input_name']
        else:
            self.infile = self.file_link.split('/')[-1]

        self.submitted_dir = 'submitted_data/'
        pathlib.Path(self.submitted_dir).mkdir(exist_ok=True)
        self.sub_datatype_dir = self.submitted_dir + str(self.data_type) + '/'
        pathlib.Path(self.sub_datatype_dir).mkdir(exist_ok=True)

        self.hash = self.yaml_data['data_file']['location']['md5']
        self.downloaded_file = self.sub_datatype_dir + self.infile
        self.fileset = [self.downloaded_file,]

        if os.path.isfile(self.sub_datatype_dir + self.infile):
            logging.info(f' file exists in {self.sub_datatype_dir}')
            return

        if verify_weblink(self.file_link):
            if not os.path.isfile(self.sub_datatype_dir + self.infile):
                logging.info(f' starting download')
                if 'drive.google.com' in self.file_link:
                    download_gdrive_file(
                        self.file_link,
                        self.sub_datatype_dir,
                        self.infile
                    )
                else:
                    download_file(
                        self.file_link,
                        self.sub_datatype_dir,
                    )
            else:
                logging.info(f' file found in {self.sub_datatype_dir}')

            self.hash = self.yaml_data['data_file']['location']['md5']
            self.downloaded_file = self.sub_datatype_dir + self.infile
            self.fileset = [self.downloaded_file,]
        else:
            logging.error(f' file unavailable')
            sys.exit(1)


    def make_archive_dir(self):
        """Make a temp dir for downloaded files to be archived"""
        self.archive_dir = 'submitted_data/downloaded_archive/'
        pathlib.Path(self.archive_dir).mkdir(exist_ok=True)


    def extract_bulk(self):
        """Extract downloaded bulk file"""
        import tarfile

        if tarfile.is_tarfile(self.downloaded_file):
            tarred_data = tarfile.open(
                self.downloaded_file,
                mode='r:*'
            )
            # 'multiple' option is added for edge cases where multiple
            # data type files in a tarfile share the same yaml
            if self.data_type == 'multiple':
                verify_dir(tarred_data)
                tarred_data.extractall(path=self.submitted_dir)
                self.fileset = [self.submitted_dir + x for x in tarred_data.getnames()]
            else:
                self.fileset = []
                for member in tarred_data.getmembers():
                    if member.isreg():
                        member.name = os.path.basename(member.name)
                        tarred_data.extract(member, path=self.sub_datatype_dir)
                        self.fileset.append(self.sub_datatype_dir + member.name)

        else:
            raise Exception(' invalid archive file for upload_bulk')

        # moving the downloaded archive file to subdir
        # 'downloaded_archive' to avoid being processed later
        # ref: https://github.com/greenelab/cimr-d/pull/12
        self.make_archive_dir()
        new_path = self.archive_dir + self.infile
        os.rename(
            self.downloaded_file,
            new_path
        )
        self.downloaded_file = new_path


    def get_colnames(self):
        """Initializing submitted column names to cimr variables"""
        if 'columns' in self.yaml_data['data_file'].keys():
            columnset = self.yaml_data['data_file']['columns']
            self.columnset = {
                v: k for k, v in columnset.items() if v != 'na'
            }
            logging.info(f' following columns will be renamed: ')
            logging.info(f' {self.columnset.keys()}')
        else:
            logging.info(f' column header change is not indicated.')


    def make_metatable(self):
        """Collect and store information for downstream analyses
        This function is meant to create a catalog for processed data
        in cimr-d. The updated table will be committed back to repo
        """
        from ..defaults import META_HEADER
        from ..defaults import CATALOG

        metadata_file = 'cimr-d_catalog.txt'

        if os.path.isfile(metadata_file):
            logging.info(f' reading a local cimr-d catalog file.')
            metadata = pandas.read_csv(
                metadata_file,
                header=0,
                index_col=None,
                sep='\t'
            )
        elif verify_weblink(CATALOG):
            logging.info(f' reading cimr-d catalog from cimr-d repo.')
            metadata = pandas.read_csv(
                CATALOG,
                header=0,
                index_col=None,
                sep='\t'
            )
        else:
            logging.info(f' creating a new cimr-d catalog.')
            metadata = pandas.DataFrame(columns=META_HEADER)

        for file_name in self.fileset:
            new_row = {
                'file_name': file_name.split('/')[-1],
                'submitted_data_url': self.file_link,
                'submitted_data_md5': self.hash,
                'build': self.genome_build
            }

            new_row['data_type'] = self.data_type

            if 'context' in self.yaml_data['data_info'].keys():
                context = self.yaml_data['data_info']['context']
                context = standardize_context(context)
                new_row['context'] = context
            else:
                logging.error(f' context description is required.')
                sys.exit(1)

            if 'description' in self.yaml_data['data_file'].keys():
                new_row['description'] = self.yaml_data['data_file']['description']
            else:
                logging.info(f' data description is not provided.')

            if 'sample_size' in self.yaml_data['data_info'].keys():
                new_row['sample_size'] = self.yaml_data['data_info']['sample_size']
            else:
                logging.info(f' sample_size is not provided.')

            if 'n_cases' in self.yaml_data['data_info'].keys():
                new_row['n_cases'] = self.yaml_data['data_info']['n_cases']
            else:
                logging.info(f' n_cases is not provided.')

            if 'citation' in self.yaml_data['data_info'].keys():
                new_row['citation'] = self.yaml_data['data_info']['citation']
            else:
                logging.info(f' citation is not provided.')

            if 'data_source' in self.yaml_data['data_info'].keys():
                new_row['data_source'] = self.yaml_data['data_info']['data_source']
            else:
                logging.info(f' data_source is not provided.')

            if 'name' in self.yaml_data['method'].keys():
                new_row['method_name'] = self.yaml_data['method']['name']
            else:
                logging.info(f' method name is not provided.')

            if 'tool' in self.yaml_data['method'].keys():
                new_row['method_tool'] = self.yaml_data['method']['tool']
            else:
                logging.info(f' method tool is not provided.')

            logging.info(f' updating cimr-d catalog.txt for {file_name}.')
            metadata = metadata.append(new_row, ignore_index=True)
            metadata.reset_index(inplace=True, drop=True)

            metadata = metadata[META_HEADER]

            metadata.to_csv(
                metadata_file,
                header=True,
                index=False,
                sep='\t',
                na_rep='NA',
                mode='w'
            )


    def check_data_file(self):
        """Standard set of Yamler functions to check information on the
        contributed data file for ci cimr processing.
        """
        self.set_data_type()
        self.set_genome_build()
        self.download()
        self.check_hash()

        if self.infile.endswith(BULK_EXTENSION):
            self.extract_bulk()

        self.get_colnames()
        self.make_metatable()


if __name__ == '__main__':

    logging.basicConfig(level='DEBUG')

    if len(sys.argv) == 1:
        yaml_files = find_yaml_in_dir()
    else:
        yaml_files = [pathlib.Path(sys.argv[1]),]

    convert_yaml(yaml_files)
