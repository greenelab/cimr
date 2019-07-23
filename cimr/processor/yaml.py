"""Reading and parsing through the contributor's yaml file.
(c) YoSon Park

This is the default uploading skim for single and bulk files 
using zenodo. For PR-based file uploader, check .circleci/deploy.sh 
and .circleci/process_submitted_data.py
"""


import os
import sys
import yaml
import pandas
import pathlib
import logging

from ..defaults import DATA_TYPES
from ..defaults import CONFIG_FILE_EXTENSION
from ..defaults import BULK_EXTENSION


def match_magic():
    """https://www.garykessler.net/library/file_sigs.html"""
    pass


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


def predefine_yaml():
    """A git-status-independent function used for testing"""
    return pathlib.Path('upload_data_example.yml')


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
            raise Exception(f' {yaml_file} is not an acceptible yaml file')

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
    import urllib

    weburl = urllib.request.urlopen(path)

    if weburl.getcode() == 200:
        return True
    else:
        return False


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
        if member.name.startswith(DATA_TYPES):
            return True
        else:
            logging.error(f' data_type not indicated in dir tree')
            sys.exit(1)


class Yamler:
    """A collection of utilities to parse the yaml file, check metadata
    and trigger cimr processing of the contributed file
    """

    def __init__(self, yaml_data):
        self.yaml_data = yaml_data
        self.data_type = None
        self.keys = None
        self.hash = None
        self.outdir = None
        self.downloaded_file = None


    def pick_keys(self):
        """List keys for the dictionarized yaml data and store in self.
        The following keys are expected:
        ['defined_as', 'data_file', 'contributor', 'data_info', 'method']
        """
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
    

    def check_hash(self):
        """Compare md5 of the downloaded file to the provided value"""
        if validate_hash(self.downloaded_file, self.hash):
            logging.info(f' data is ready for cimr processing.')
        else:
            raise ValueError(' provided md5 hash didn\'t match.')


    def download(self):
        """Check if provided weblink to the file exists. 
        Download if verified.
        """
        path = self.yaml_data['data_file']['location']['url']
        self.downloaded_file = path.split('/')[-1]

        self.outdir_root = 'submitted_data/'
        pathlib.Path(self.outdir_root).mkdir(exist_ok=True)
        self.outdir = self.outdir_root + str(self.data_type) + '/'
        pathlib.Path(self.outdir).mkdir(exist_ok=True)

        if verify_weblink(path):
            logging.info(f' starting download')
            download_file(path, self.outdir)
            self.hash = self.yaml_data['data_file']['location']['md5']
            self.downloaded_file = self.outdir + self.downloaded_file
            self.check_hash()
        else:
            logging.error(f' file unavailable')
            sys.exit(1)


    def bulk_download(self):
        """Check if file is a recognized tarfile, compare against the 
        provided md5 hash value and download
        """
        import tarfile

        if tarfile.is_tarfile(self.downloaded_file):
            tarred_data = tarfile.open(
                self.downloaded_file, 
                mode='r:*'
            )
        else:
            raise Exception(' invalid archive file for upload_bulk')

        if self.data_type == 'multiple':
            verify_dir(tarred_data)
            tarred_data.extractall(path=self.outdir_root)
        else:
            for member in tarred_data.getmembers():
                if member.isreg():
                    member.name = os.path.basename(member.name)
                    tarred_data.extract(member, path=self.outdir)
                    
        tarfile_dir = 'submitted_data/downloaded_archive/'
        pathlib.Path(tarfile_dir).mkdir(exist_ok=True)
        new_path = tarfile_dir + self.downloaded_file.split('/')[-1]
        os.rename(
            self.downloaded_file,
            tarfile_dir + self.downloaded_file.split('/')[-1]
        )
        self.downloaded_file = new_path


    def check_defined(self):
        """Check whether the submitted data is a single file"""
        if self.yaml_data['defined_as'] == 'upload':
            self.download()
        elif self.yaml_data['defined_as'] == 'upload_bulk':
            self.bulk_download()
        else:
            logging.error(f' not an acceptible \'defined_as\' variable')
            sys.exit(1)


    def check_data_file(self):
        """Standard set of Yamler functions to check information on the
        contributed data file for ci cimr processing.
        """
        self.set_data_type()
        self.check_defined()


if __name__ == '__main__':

    logging.basicConfig(level='DEBUG')

    if len(sys.argv) == 1:
        yaml_files = find_yaml_in_dir()
    else:
        yaml_files = pathlib.Path(sys.argv[1])

    try:
        for yaml_file in yaml_files:
            yaml_file = pathlib.Path(yaml_file)
            yaml_file_path = yaml_file.resolve(strict=True)
            logging.info(f' processing metadata {yaml_file_path}')
            yaml_data = load_yaml(yaml_file)
            y = Yamler(yaml_data)
            y.check_data_file()
    except FileNotFoundError:
        logging.info(f' no new yaml file found to process')
        sys.exit(0)


