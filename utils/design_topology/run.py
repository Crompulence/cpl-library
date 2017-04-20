#!/usr/bin/env python

import subprocess as sp
import os
import atexit
import signal

#To ensure processes killed if this script is killed
def kill_sp(pids):
    for pid in pids:
        if pid is None:
            pass
        else:
            try:
                os.kill(pid, signal.SIGTERM)
            except OSError:
                pass


class cd:
    """Context manager for changing the current working directory"""
    def __init__(self, newPath):
        self.newPath = os.path.expanduser(newPath)

    def __enter__(self):
        self.savedPath = os.getcwd()
        os.chdir(self.newPath)

    def __exit__(self, etype, value, traceback):
        os.chdir(self.savedPath)

def setup():

    pass

def run():


    md = sp.Popen("mpiexec -n 4 ./fortran/md", shell=True)
    cfd = sp.Popen("mpiexec -n 1 python ./CFD.py", shell=True)

    atexit.register(kill_sp, [md.pid, cfd.pid])

    cfd.wait()
    md.wait()

#setup()
run()

