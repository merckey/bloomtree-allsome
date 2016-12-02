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


def check_command(command):
    """
    check if corresponding command available
    """
    if os.path.isfile(command):
        return True

    for cmdpath in os.environ['PATH'].split(':'):
        if os.path.isdir(cmdpath) and command in os.listdir(cmdpath):
            return True
    return False

btquery = '/gpfs/cyberstar/pzm11/backup/cxs1031/code/bloomtree-dev/src/bt query-original'

tree_file_list = ['/gpfs/cyberstar/pzm11/backup/sbt/compressedSBT/SBT_list.txt', '/gpfs/cyberstar/pzm11/backup/sbt/clusteringSBT/clustering.rrr.sbt']

#query_file_name_list = []
query_file_name_list = ['thousand', 'hundred.0', 'hundred.1', 'hundred.2']

for i in range(10):
    query_file_name_list.append('ten.'+str(i))

query_file_name_list = [x+'.seq' for x in query_file_name_list]

print query_file_name_list

query_directory = '/gpfs/cyberstar/pzm11/nobackup/sbt/test_data/gencode_transcript/'

output_directory = '/gpfs/cyberstar/pzm11/nobackup/cxs1031/sbt_private/query_stat/'

tree_type_list = ['solomon', 'clustered']

for i in range(2):
	tree_type = tree_type_list[i]
	tree_file = tree_file_list[i]
	for query_file_name in query_file_name_list:
		output_file_name = tree_type + '.' + query_file_name + '.query'
		log_file_name = output_file_name + '.log'
		output_file_name = output_file_name + '.out'
		query_file = query_directory + query_file_name
		output_file = output_directory + output_file_name
		log_file = output_directory + log_file_name
		query_command = 'time ' + btquery + ' ' + tree_file + ' ' + query_file + ' ' + output_file + ' > ' + log_file
		print '----before run :'
		print datetime.now()
		shell_run(query_command)
		print datetime.now()
		print '----run finished.'

