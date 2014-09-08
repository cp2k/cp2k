#! /usr/bin/env python

import sys, re
from sys import argv
import prettify
from formatting import normalizeFortranFile
from formatting import replacer
from formatting import addSynopsis

def instantiateTemplate(infile,outfile,subs,logFile=sys.stdout):
  import re
  # this ignores a unbalanced '[' or ']' and ignores a comment
  # like [ ! ], but well,...
  directive=re.compile("(\[[ \t]*([^\]]*)[ \t]*\])|(^[ \t]*![ \t]*\[[ \t]*template[ \t]*\(([^\)]*)\)\][ \t]*$)|(!)")
  logFile.write("instantiating '"+infile.name+"' to '"+outfile.name+"'\n")
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
	      logFile.write("ERROR: missing required argument:"+arg+"\n")
	      outfile.write("! ERROR argument '"+arg+"' missing\n") 
          kList=subs.keys()
          kList.sort()
          for arg in kList:
            sost=subs[arg].split("\n")
            if (len(sost)>1):
              outfile.write('!  '+arg+' = \n')
              outfile.write('!    "'+sost[0])
              for sostLine in sost[1:]:
                outfile.write('\n!     '+sostLine)
              outfile.write('"\n')
            else:
              outfile.write('!  '+arg+' = "'+sost[0]+'"\n')
          outfile.write("\n")
	  break
	elif not match.groups()[1]:
	  if not inComment:
	    logFile.write("WARNING, ingnoring empty group at line %d\n"%lineNr)
	  outfile.write(line[:match.end()])
	  line=line[match.end():]
	else:
	  if subs.has_key(match.groups()[1]):
	    outfile.write(line[:match.start()]+subs[match.groups()[1]])
	    line=line[match.end():]
	  else:
	    if not inComment:
	      logFile.write("WARNING ignoring unknown token '%s' at line %d\n"%
                            (match.groups()[1],lineNr))
	    outfile.write(line[:match.end()])
	    line=line[match.end():]
    except:
      logFile.write("error in '%s' at line %d\n"%(infile.name,lineNr))
      outfile.close()
      infile.close()
      os.rename(outfile.name,outfile.name+".err")
      raise
  outfile.close()
  infile.close()

def evaluateInstantiationFile(instantiationFile,logFile=sys.stdout,outDir=None):
    import os
    generatedFiles=[]
    try:
      extension=".instantiation"
      if not instantiationFile.endswith(extension):
          logFile.write("ERROR input '"+ instantiationFile+"' is not a "+extension+" file!!\n")
          raise
      input = open(instantiationFile,'r')
      subst = eval(input.read())
      errors=0
      for substitution in subst:
        inName = instantiationFile.replace(extension, ".template")
        if substitution.has_key("template"):
            inName = substitution["template"]
        try: infile=open(inName,"r")
        except:
          logFile.write("ERROR opening template '"+inName+"'\n")
          raise
        outName=instantiationFile[:-len(extension)]
        id_ext=0
        for token in substitution.keys():
          if token=="ext" :
            id_ext=1
          outName=re.sub("(?<![a-zA-Z0-9])_"+token+"_(?![a-zA-Z0-9])",
                         substitution[token],outName)
        tmpName=outName+".F"
        if id_ext==1 :
          tmpName=outName+substitution["ext"]
        outName=tmpName
        if outDir:
          outName=os.path.join(outDir,os.path.basename(outName))
        try: outfile=open(outName,'w')
        except:
          logFile.write("ERROR opening template '"+outName+"'\n")
          raise
        instantiateTemplate(infile,outfile,substitution,logFile)
        prettify.prettfyInplace(outName,logFile=logFile)
        generatedFiles.append(outName)
    except:
        import sys
        logFile.write("error evaluating substitutions from file '"+
                      instantiationFile+"'\n")
        import sys, traceback
        logFile.write('-'*60+"\n")
        traceback.print_exc(file=logFile)
        logFile.write('-'*60+"\n")
    return generatedFiles

if __name__ == '__main__':
    if len(sys.argv)<2:
        print "usage:", sys.argv[0]," template1.instantiation [template2.instantiation ...]"
    else:
        for name in sys.argv[1:]:
            evaluateInstantiationFile(name,sys.stdout)
# Local Variables:
# py-indent-offset: 2
# End:
