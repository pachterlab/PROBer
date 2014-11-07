# utils module

import os
import sys
import subprocess
import argparse
import textwrap

class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def _split_lines(self, text, width):
        if text.startswith("L:"):
            return text[2:].splitlines()
        return argparse.ArgumentDefaultsHelpFormatter._split_lines(self, text, width)
    
    def _fill_text(self, text, width, indent):
        if text.startswith("L:"):
            return ''.join([indent + line for line in text[2:].splitlines(True)])
        else:
            text = self._whitespace_matcher.sub(' ', text).strip()
            return textwrap.fill(text, width, initial_indent=indent,
                                           subsequent_indent=indent)

class MyParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write("{}: error: {}\n".format(os.path.basename(sys.argv[0]), message))
        self.print_help()
        sys.exit(-1)

def nargs_range(n_min, n_max):
    """ Require number of arguments between n_min and n_max """
    class _StoreConstraintAction(argparse.Action):
        def __call__(self, parser, namespace, values, option_string = None):
            if len(values) > n_max or len(values) < n_min:
                raise argparse.ArgumentTypeError("{} needs {} ~ {} arguments".format(self.dest, n_min, n_max))
            setattr(namespace, self.dest, values)
    return _StoreConstraintAction
        
demo = False
        
def runProg(command):
    """ Run command using a subprocess """

    print(command)

    if demo:     
        return None
    
    try:
        subprocess.check_call(command, shell = True)
    except subprocess.CalledProcessError as error:
        print("Command \"{}\" failed!".format(error.cmd))
        sys.exit(-1)
