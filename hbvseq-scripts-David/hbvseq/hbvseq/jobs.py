"""Run jobs either locally or on a HPC cluster using slurm"""

import os
import tempfile
from subprocess import call

from .settings import hbvseqSettings

def executeCommand(commStr, runLocally, cpu=8):
    """Run command either locally or on a cluster"""
    if(runLocally):
        call(commStr, shell=True)
    else:
        submit = ('sbatch --wrap="{comm}"'
                  '--cpus-per-task={cpu} --mem-per-cpu=4000').format(comm=commStr,
                                            cpu=cpu)
        call(submit, shell=True)

def executeCommandsOnCluster(commStrs, cpu=8):
    """Run a list of commands on a cluster"""
    submitComms = []
    for str in commStrs:
        submit = ('sbatch --wrap="{comm}" '
                         '--cpus-per-task={cpu}').format(comm=str,
                                                         cpu=cpu)
        submitComms.append(submit)
        call(submit, shell=True)
    return submitComms

