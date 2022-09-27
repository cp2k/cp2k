#!/usr/bin/env python3


import sys
import argparse


class Parser:
    "Parser for files with IF_XYZ(A|B) constructs"

    def __init__(self, switches):
        self.mSwitches = switches

    def Switches(self):
        "outputs the switches used by parser"
        return self.mSwitches

    def SetSwitch(self, key, val):
        "Set or add a (key,val) switch"
        self.mSwitches[key] = val

    def ParseSingleIf(self, string, switch):
        """
        switch should be a tuple (key, val)
        ParseSingleIf(string, switch) will replace in string the
        first occurrence of IF_key(A|B) with A if val = True;
        B if val = False
        """
        init = string.find("IF_" + switch[0])
        start = init
        end = len(string)
        mark = end
        counter = 0
        ind = start
        # determine the correct location for '|'
        for cc in string[init:]:
            if cc == "(":
                if counter == 0:
                    start = ind
                counter += 1
            elif cc == ")":
                counter -= 1
                if counter == 0:
                    end = ind
                    break
            elif cc == "|" and counter == 1:
                mark = ind
            ind += 1
        # resolve the option
        if switch[1]:
            result = (
                string[0:init]
                + string[start + 1 : mark]
                + string[end + 1 : len(string)]
            )
        else:
            result = (
                string[0:init] + string[mark + 1 : end] + string[end + 1 : len(string)]
            )
        return result

    def ParseIf(self, string, switch):
        """
        ParseIf(string, switch) will replace in string recursively
        occurance of IF_key(A|B) statements with A if val = True;
        B if val = False
        """
        result = string
        while result.find("IF_" + switch[0]) > -1:
            result = self.ParseSingleIf(result, switch)
        return result

    def ParseString(self, string):
        """
        ParseString(string) will parse in string recursively
        all of IF_key(A|B) statements for all (key, val) pairs
        in dictionary self.mSwitches
        """
        result = string
        for switch in self.mSwitches.items():
            result = self.ParseIf(result, switch)
        return result

    def ParseDocument(self, input_file, output_file):
        """
        ParseDocument(input_file, output_fiel) will replace recursively
        all of IF_key(A|B) statements for all (key, val) pairs
        in dictionary self.mSwitches in the input_file stream and output
        to the output_file stream.
        """
        output = []
        for line in input_file:
            output.append(self.ParseString(line))

        if input_file == output_file:
            output_file.seek(0)
            output_file.truncate()

        # write to the same file
        for line in output:
            output_file.write(line)


# ------------------------------------------------------------------------
# main program
# ------------------------------------------------------------------------
if __name__ == "__main__":
    argparser = argparse.ArgumentParser(
        description="Resolve IF_*() macros based on given tags"
    )
    argparser.add_argument(
        "-f",
        "--file",
        metavar="FILENAME",
        type=str,
        help="read from  given file instead of stdin",
    )
    argparser.add_argument(
        "-i",
        "--inplace",
        action="store_true",
        help="do in-place replacement for the given file instead of echo to stdout",
    )
    argparser.add_argument("--selftest", action="store_true", help="run self test")
    argparser.add_argument(
        "flags",
        metavar="FLAG",
        type=str,
        nargs="*",
        help="specified flags are set to true",
    )

    args = argparser.parse_args()

    # default list of switches used by the parser
    switches = {
        "MPI": False,
        "CUDA": False,
        "HIP": False,
        "OPENCL": False,
        "WARNALL": False,
        "DEBUG": False,
        "ASAN": False,
        "STATIC": False,
        "VALGRIND": False,
        "COVERAGE": False,
    }

    parser = Parser(switches)

    # set the list of switches given on the command line to True
    for flag in args.flags:
        parser.SetSwitch(flag, True)

    if args.selftest:
        sys.exit(0)  # TODO implement selftest

    # do parsing

    if not args.file:
        parser.ParseDocument(sys.stdin, sys.stdout)
    else:
        if args.inplace:
            with open(args.file, mode="r+") as fhandle:
                parser.ParseDocument(fhandle, fhandle)
        else:
            with open(args.file, mode="r") as fhandle:
                parser.ParseDocument(fhandle, sys.stdout)
