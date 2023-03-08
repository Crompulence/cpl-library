import subprocess as sp
import os


class cd:
    """Context manager for changing the current working directory"""
    def __init__(self, newPath):
        self.newPath = os.path.expanduser(newPath)

    def __enter__(self):
        self.savedPath = os.getcwd()
        os.chdir(self.newPath)

    def __exit__(self, etype, value, traceback):
        os.chdir(self.savedPath)


def get_subprocess_error(e):
    print("subprocess ERROR")
    import json
    error = json.loads(e[7:])
    print((error['code'], error['message']))


def buildstring(code, exename):

    bldstr = ("mpif90 " + code 
             + "-I" + os.environ["CPL_PATH"] 
             + "/include  -L" + os.environ["CPL_PATH"] + "/lib  " 
             + "-Wl,-rpath=" + os.environ["CPL_PATH"]  + "/lib/ -lcpl"
             + " -o " + exename)

    return bldstr

def runcmd(cmd, raiseerror=False):
    try:
        p = sp.Popen(cmd, stdout=sp.PIPE, stderr=sp.PIPE, shell=True)
        pcommout = p.communicate()
        output = pcommout[0].decode("utf-8")
        error = pcommout[1].decode("utf-8")
        if p.returncode != 0: 
            print("returncode", p.returncode)
            print("Stdout = ", output)
            print("Stderror = ", error)
    except sp.CalledProcessError as e:
        print("Stdout = ", e.stdout)
        print("Stderror = ", e.stderr)
        if e.output.startswith(b'error: {'):
            get_subprocess_error(e.output)
        if raiseerror:
            raise sp.CalledProcessError

    return output+error

