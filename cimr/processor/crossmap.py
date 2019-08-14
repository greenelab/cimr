
import os
import sys
import subprocess


infile = sys.argv[1]

clone_cmd = 'git clone git@github.com:liguowang/CrossMap.git'
build_cmd = 'python3 setup.py build'
install_cmd = 'python3 setup.py install'
xmap_cmd = 'CrossMap.py bed /work/drug/refsnp/ensembl_chain/GRCh38_to_GRCh37.chain.gz /work/drug/refsnp/test_chr22.bed /work/drug/refsnp/test_chr22_hg19_inpy.bed'
clean_cmd = 'rm -rf CrossMap'


def run_cmd(cmd):
    running = subprocess.Popen(cmd.split(), stdout = subprocess.PIPE)
    output, error = running.communicate()
    print('output is : ', output)
    print('error is : ', error)

run_cmd(clone_cmd)

os.chdir('CrossMap')
run_cmd(build_cmd)
run_cmd(install_cmd)

run_cmd(xmap_cmd)

os.chdir('../')
run_cmd(clean_cmd)
