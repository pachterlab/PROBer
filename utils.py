# utils module

import sys
import subprocess
import argparse

class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def _split_lines(self, text, width):
        if text.startswith("L:"):
            return text[2:].splitlines()
        return argparse.ArgumentDefaultsHelpFormatter._split_lines(self, text, width)

class MyParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write("{}: error: {}\n".format(sys.argv[0], message))
        self.print_help()
        sys.exit(-1)

demo = False
        
def runProg(command):
    """ Run command using a subprocess """

    if demo:     
        print(command)
        return None
    
    try:
        subprocess.check_call(command, shell = True)
    except subprocess.CalledProcessError as error:
        print("Command \"{}\" failed!".format(error.cmd))
        sys.exit(-1)
