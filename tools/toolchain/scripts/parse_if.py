import sys

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
first occurance of IF_key(A|B) with A if val = True;
B if val = False
        """
        init = string.find('IF_' + switch[0])
        start = init
        end = len(string)
        mark = end
        counter = 0
        ind = start
        # determine the correct location for '|'
        for cc in string[init:]:
            if cc == '(':
                if counter == 0: start = ind
                counter += 1
            elif cc == ')':
                counter -= 1
                if counter == 0: end = ind; break
            elif cc == "|" and counter == 1:
                mark = ind
            ind += 1
        # resolve the option
        if switch[1]:
            result = string[0         : init] + \
                     string[start + 1 : mark] + \
                     string[end + 1   : len(string)]
        else:
            result = string[0         : init] + \
                     string[mark + 1  : end]  + \
                     string[end + 1   : len(string)]
        return result
    def ParseIf(self, string, switch):
        """
ParseIf(string, switch) will replace in string recursively
occurance of IF_key(A|B) statements with A if val = True;
B if val = False
        """
        result = string
        while result.find('IF_' + switch[0]) > -1:
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
    def ParseDocument(self, filename):
        """
ParseDocument(filename) will replace in a file recursively
all of IF_key(A|B) statements for all (key, val) pairs
in dictionary self.mSwitches
        """
        input_file = open(filename, "r")
        output=[]
        for line in input_file:
            output.append(self.ParseString(line))
        input_file.close()
        # write to the same file
        output_file = open(filename, "w")
        for line in output:
            output_file.write(line)
        output_file.close()

# ------------------------------------------------------------------------
# main program
# ------------------------------------------------------------------------

if (len(sys.argv) < 2) or \
   ("--help" in sys.argv) or \
   ("-h" in sys.argv):
    sys.stderr.write("usage: parse_if <filename> MPI OMP CUDA ...\n")
    sys.exit()

# default list of switches used by the parser
switches = {"MPI"      : False, \
            "OMP"      : False, \
            "CUDA"     : False, \
            "WARNALL"  : False, \
            "DEBUG"    : False, \
            "VALGRIND" : False, \
            "COVERAGE" : False}

parser=Parser(switches)

# set the list of switches given on the command line to True
for ii in sys.argv[2:]:
    parser.SetSwitch(ii, True)

# do parsing
parser.ParseDocument(sys.argv[1])
