#!/apps/python/bin/python

import sys
from sys import argv

def instantiateTemplate(infile,outfile,subs):
  import re
  # this ignores a unbalanced '[' or ']' and ignores a comment
  # like [ ! ], but well,...
  directive=re.compile("(\[[ \t]*([^\]]*)[ \t]*\])|(^[ \t]*![ \t]*\[[ \t]*template[ \t]*\(([^\)]*)\)\][ \t]*$)|(!)")
  print "instantiating '",infile.name,"' to '",outfile.name,"'"
  lineNr=0
  while 1:
    line=infile.readline()
    lineNr=lineNr+1
    if not line: break
    inComment=0
    try:
      while 1:
	match=directive.search(line)
	if not match:
	  outfile.write(line)
	  break
	if match.groups()[4]: # it is a comment
	  outfile.write(line[:match.end()])
	  line=line[match.end():]
	  inComment=1
	elif match.groups()[2]: # the template description
	  outfile.write(line)
	  outfile.write("! ARGS:\n")
	  for arg in match.groups()[3].split(","):
	    arg=arg.strip()
	    if not subs.has_key(arg):
	      print "ERROR: missing required argument:",arg
	      outfile.write("! ERROR argument '"+arg+"' missing\n") 
          for arg in subs.keys():
	    sost=subs[arg].split("\n")
            if (len(sost)>1):
              outfile.write('!  '+arg+' = ')
              outfile.write('!    "'+sost[0])
              for sostLine in sost[1:]:
                outfile.write('\n!     '+sostLine)
              outfile.write('!     "\n')
            else:
              outfile.write('!  '+arg+' = "'+sost[0]+'"\n')
          outfile.write("\n")
	  break
	elif not match.groups()[1]:
	  if not inComment:
	    print "WARNING, ingnoring empty group at line ",lineNr
	  outfile.write(line[:match.end()])
	  line=line[match.end():]
	else:
	  if subs.has_key(match.groups()[1]):
	    outfile.write(line[:match.start()]+subs[match.groups()[1]])
	    line=line[match.end():]
	  else:
	    if not inComment:
	      print "WARNING ignoring unknown token '",match.groups()[1],"' at line ",lineNr
	    outfile.write(line[:match.end()])
	    line=line[match.end():]
    except:
      import sys
      print "error in ",infile.name," at line ",lineNr
      print sys.exc_type, sys.exc_value
  outfile.close()
  infile.close()


if len(sys.argv)<2:
    print "usage:", sys.argv[0]," template1.instantiation [template2.instantiation ...]"
else:
  for name in sys.argv[1:]:
    try:
      input = open(name,'r')
      subst = eval(input.read())
      for substitution in subst:
        extension=".instantiation"
        if not name[-len(extension):]==extension :
          print "ERROR input",name,"is not a ",extension," file!!"
          break
        inName=input.name.replace(extension, ".template")
        try: infile=open(inName,"r")
        except:
          print("ERROR opening template '",inName,"'")
          print sys.exc_type, sys.exc_value
        outName=name[:-len(extension)]
        for token in substitution.keys():
          outName=outName.replace("_"+token+"_",
                                  substitution[token])
        outName=outName+".F"
        try: outfile=open(outName,'w')
        except:
          print "ERROR opening template '",outName,"'"
          print sys.exc_type, sys.exc_value
        instantiateTemplate(infile,outfile,substitution)
    except:
      import sys
      print "error evaluating substitutions from file ",name
      print sys.exc_type, sys.exc_value

# Local Variables:
# py-indent-offset: 2
# End:
