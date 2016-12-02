import sys
import os
import subprocess
import time
from datetime import datetime

RUN=True

def shell_run(command, hide=False):
    if not RUN:
        time.sleep(3.5)
        print(command)
    else:
        print(command)
        if hide:  # hide output
            FNULL = open(os.devnull, 'w')
            subprocess.call(command, shell=True, stdout=FNULL, stderr=subprocess.STDOUT)
            # subprocess.call(command, shell=True, stdout=FNULL)
            FNULL.close()
        else:
            subprocess.call(command, shell=True)

jellyfish_count = 'jellyfish count -m 63 -L 3 -s 3G -t 10 -C '
jellyfish_dump = 'jellyfish dump '
fasta_dir = '/storage/home/cxs1031/sbt_nobackup/test_data/'
output_dir = '/storage/home/cxs1031/sbt_nobackup/bft/'

kmerlist_filename = output_dir + 'kmerlist.txt'

kmerlist_file = open(kmerlist_filename, 'w')

with open(sys.argv[1]) as f:
    for line in f.readlines():
        output_name = line.strip()
        jf_name = output_name + '.jf'
        kmer_name = output_name + '.kmer'
        jf_file = output_dir + jf_name
        kmer_file = output_dir + kmer_name
        fasta_file = fasta_dir + output_name + '.fa'
        count_command = jellyfish_count + ' -o ' + jf_file + ' ' + fasta_file
        dump_command = jellyfish_dump + ' ' + jf_file + " | grep -v '^>' > " + kmer_file
        shell_run(count_command)
        shell_run(dump_command)
        kmerlist_file.write(kmer_file)

kmerlist_file.close()