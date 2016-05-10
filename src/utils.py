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

def nargs_range(list_of_range):
    """ Require number of arguments between n_min and n_max """
    class _StoreConstraintAction(argparse.Action):
        def __call__(self, parser, namespace, values, option_string = None):
            if not (len(values) in list_of_range):
                raise argparse.ArgumentTypeError("{} needs {} arguments".format(self.dest, list_of_range))
            setattr(namespace, self.dest, values)
    return _StoreConstraintAction
        
def expand(input):
    """ Expand input string to remove ~ and resovle symbolic link"""
    return os.path.realpath(os.path.expanduser(input))

def expandAll(input):
    """ input is a list separated by comma """
    inputs = input.split(',')
    res = []
    for afile in inputs:
        res.append(os.path.realpath(os.path.expanduser(afile)))
    return ",".join(res)

demo = False
        
def runProg(command, command2 = None, catch_stderr = None):
    """ Run command using a subprocess, if command2 != None, use Pipe """

    commandStr = " ".join(command) + (" 2> {}".format(catch_stderr) if catch_stderr != None else "") + (" | " + " ".join(command2) if command2 != None else "")
    print(commandStr)

    if demo:     
        return None
    
    try:
        if command2 == None:
            subprocess.check_call(command)
        else:
            fd = open(catch_stderr, "w") if catch_stderr != None else None
            p1 = subprocess.Popen(command, stdout = subprocess.PIPE, stderr = fd)
            subprocess.check_call(command2, stdin = p1.stdout)
            p1.stdout.close()
            if fd != None:
                fd.close()

    except subprocess.CalledProcessError as error:
        print("Command \"{}\" failed!".format(commandStr))
        sys.exit(-1)

