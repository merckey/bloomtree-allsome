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

with open(sys.argv[1]) as list_f:
    for line in list_f:
        line = line.strip()
        filename = line.split()[-1]
        delete_command = 'rm ' + filename
        shell_run(delete_command)