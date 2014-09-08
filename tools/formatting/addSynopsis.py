#! /usr/bin/env python

import sys
from sys import argv

def parseInterface(iFile,logFile=sys.stdout):
    """Returns a dictionary with the interface contained in the file iFile.
    
    iFile should be a file object (stream).
    The result is a dictionary of the following kind: 
    {'interface':{
      'interface1':{
        'name':'interface1','kind':'interface','expansion':
            'Interface Interface1\n  Module procedure func1\nEnd Interface'
        },'interface2':{...}
      },
      'type':{...},'function':,{...},'subroutine':,{...},'module':,{...}
    }
    """
    import re
    import exceptions
    # removed the |(!?[ \t]*contains)
    nameRe=re.compile("[ \t]*(((module)[ \t]*:?|subroutine|function|interface|(type)[ \t]*:*)[ \t]*([a-zA-Z0-9_]+)|(interface)|(end)[ \t]*([a-zA-Z_0-9]*)[ \t]*([a-zA-Z0-9_]*))",flags=re.IGNORECASE)
    (kindPos, modulePos, typePos, namePos, emptyInterfacePos,endPos, endKind) = (1,2,3,4,5,6,7)
    expansion=""
    lineNr=0
    result={"interface":{},"type":{},"function":{},"subroutine":{},"module":{}}
    stack=[]
    while 1:
	line=iFile.readline()
	lineNr=lineNr+1
	if not line: break
	try:
	    match=nameRe.match(line)
	    if not match:
		if stack:
		    stack[0]["expansion"]=stack[0]["expansion"]+line
	    elif match.groups()[namePos]:
		if match.groups()[modulePos]: # module def
		    if match.groups()[namePos] and match.groups()[namePos].lower()=="procedure": #module procedure
			if not stack: raise SyntaxError("empty stack")
			stack[0]["expansion"]=stack[0]["expansion"]+line
		    else:
			if stack: raise SyntaxError("stack non empty and module start")
			stack=[{"name":match.groups()[namePos].lower(), "kind":"module",
			"expansion":line}]
		elif match.groups()[typePos]: # type def
		    stack=[{"name":match.groups()[namePos].lower(), "kind":"type",
			"expansion":line}]+stack
		else: # other def start
		    if not match.groups()[kindPos]: raise SyntaxError("undefined kind")
		    stack=[{"name":match.groups()[namePos].lower(), "kind":match.groups()[kindPos].lower(),
			"expansion":line}]+stack
	    elif match.groups()[emptyInterfacePos]: # empty interfce begin
		stack=[{"name":"", "kind":"interface","expansion":line}]+stack
	    elif match.groups()[endPos]: # end of a def
		if not stack: raise SyntaxError("empty stack")
		if not match.groups()[endKind]:
		    raise SyntaxError("end without kind") # ignore this error??
		elif not match.groups()[endKind].lower()==stack[0]["kind"].lower():
		    raise SyntaxError("mismatched end") # ignore this error??
		stack[0]["expansion"]=stack[0]["expansion"]+line
		if stack[0]["name"]:
		    if not result.has_key(stack[0]["kind"]): result[stack[0]["kind"]]={}
		    kind_table=result[stack[0]["kind"]]
		    if kind_table.has_key(stack[0]["name"]):
			logFile.write("\nWARNING double definition of '"+
                                      stack[0]["name"]+"' at line "+
                                      str(lineNr)+", ignoring\n")
                    else:
                        kind_table[stack[0]["name"]]=stack[0] # put only the expansion??
		else:
		    logFile.write("ignoring group "+str(stack[0])+"\n")
		del stack[0]
	    else:
		raise SyntaxError("inconsistent parse '"+str(match.groups())+"'")
	except SyntaxError, e:
	    e.lineno=lineNr
	    e.text=line
	    logFile.write("stack="+str(stack)+"\n")
	    raise
	except Warning, w:
	    logFile.write("ignoring warning at line %d of file %s\n" %
                          (lineNr,iFile.name))
	    logFile.write(str(w))
    return result

def writeSynopsis(objDef,outfile,logFile=sys.stdout):
    """Writes out the synopsis (i.e. objDef['expansion']) removing empty lines.
    """
    if objDef.has_key("expansion") and objDef["expansion"]:
	outfile.write("!!   SYNOPSIS\n")
	wasEmptyLine=1
	for line in objDef["expansion"].splitlines(1):
	    if not (wasEmptyLine and line.isspace()):
                wasEmptyLine=line.isspace()
                while len(line)>=74:
                    import string
                    posToBreak=string.rfind(line," ",30,74)
                    if posToBreak<0: posToBreak=74
                    outfile.write("!! "+line[:posToBreak]+"&\n")
                    line="       "+line[posToBreak:]
                outfile.write("!! "+line)
            else:
                wasEmptyLine=line.isspace()
	if not wasEmptyLine:
	    outfile.write("!!\n")
    else:
	logFile.write("no expansion for "+str(objDef)+"\n")

def insertSynopsis(defs,infile,outfile,logFile=sys.stdout):
    """Writes the file infile to outfile inserting the synopsis defined in defs.
    
    infile and outfile should be file objects (streams).
    defs is a dictionary of the following form:
    {
        'funct1':{'expansion':
            'Function funct1()\n  integer :: funct1\nEnd function\n'}
        'funct2':{...}
        ...
    }
    
    If there is a new synopsis, the old one is removed from the robDoc comment.
    """
    import re
    import exceptions
    roboCommentRe=re.compile("[ \t]*!!(\*\*\*(\*?)|[ \t]*(NAME|COPYRIGHT|SYNOPSIS|USAGE|SOURCE|FUNCTION|DESCRIPTION|PURPOSE|AUTHOR|CREATION[ \t]+DATE|HISTORY|MODIFICATION[ \t]+HISTORY|INPUTS|ARGUMENTS|PARAMETERS|OUTPUT|SIDE[ \t]+EFFECTS|SWITCHES|RESULT|RETURN[ \t]+VALUE|EXAMPLE|OPTIONS|NOTES|DIAGNOSTICS|WARNINGS|ERRORS|BUGS|TODO|IDEAS|PORTABILITY|SEE[ \t]+ALSO|METHODS|NEW[ \t]+METHODS|ATTRIBUTES|NEW[ \t]+ATTRIBUTES|TAGS|COMMANDS|DERIVED[ \t]+FROM|DERIVED[ \t]+BY|USES|CHILDREN|USED[ \t]+BY|PARENTS|SOURCE)[ \t]*$|[ \t]*([a-zA-Z0-9_]+)|(.?))")
    (startStopPos,directivePos,namePos,genericCommentPos)=(1,2,3,4)
    # 0: not in comment, 1: in roboDoc comment, synopsis not yet inserted 2: reading name,
    # 3: in roboDoc comment, synopsis inserted, 4: skipping section (synopsis inserted)
    status=0
    lineNr=0
    name=""
    while 1:
	line= infile.readline()
	lineNr=lineNr+1
	if not line: break
	try:
	    match= roboCommentRe.match(line)
	    if not match:
		outfile.write(line)
	    elif match.groups()[startStopPos]!=None:
		if status==1: # already in comment: end
		    if name:
			if not defs.has_key(name):
			    e=SyntaxWarning("ignoring unknown object with name '"+name+"'")
			    outfile.write(line.lstrip())
			    status=0
			    raise e
			#elif defs[name]:
			#    writeSynopsis(defs[name],outfile,logFile)
		    status=0
		elif status==2: # searching for name
		    status=1
		    outfile.write(line.lstrip())
		    raise SyntaxWarning("name not found")
		elif status==3 or status==4: # synopsis inserted
		    status=0
		elif status==0: # not in comment
		    status=1
		    name=""
		    if not match.groups()[startStopPos]: # only 3 *, start should have 4
			outfile.write(line.lstrip())
			raise SyntaxWarning("incorrect start of comment")
		else:
		    status=0
		    outfile.write(line.lstrip())
		    raise SyntaxWarning("internal error, status="+str(status))
		outfile.write(line.lstrip())
	    elif match.groups()[directivePos]:
		if status==4: status=3
		if status==1 and name and defs.has_key(name) and defs[name]:
		    # writeSynopsis(defs[name],outfile)
		    status=3
		if match.groups()[directivePos]=="NAME":
		    status=2
		    outfile.write(line.lstrip())
		elif match.groups()[directivePos]=="SYNOPSIS":
		    # print "found synopsis at line",lineNr
		    if status==3 or status==4:
			status=4
			# print "skipping synopsis"
		    elif status==2:
			status=1
			outfile.write(line.lstrip())
			raise SyntaxWarning("did not find the name")
		    elif status==1:
			outfile.write(line.lstrip())
		    else:
			raise SyntaxWarning("unknown status "+str(status))
		else:
		    outfile.write(line.lstrip())
	    elif match.groups()[namePos]:
		if status==2:
		    name=match.groups()[namePos]
		    status=1
		if status!=4:
		    outfile.write(line.lstrip())
	    elif match.groups()[genericCommentPos]!=None:
		if status!=4:
		    outfile.write(line.lstrip())
	    else:
		raise SyntaxError("unknown match: "+str(match))
	except SyntaxError, e:
	    e.lineno=lineNr
	    e.text=line
	    raise
	except Warning, w:
	    logFile.write("ignoring warning at line %d of file '%s'\n"%
                          (lineNr,infile.name))
	    logFile.write(str(w)+"\n")

def addSynopsisToFile(ifile,sfile,outfile,logFile=sys.stdout):
    """Adds the synopsis to the file sfile, reading the interface from
    ifile, and writing the result to outfile"""
    rawDefs=parseInterface(ifile,logFile)
    defs={}
    #print "rawDefs=",rawDefs
    for kind in ["function","subroutine"]:
        if rawDefs.has_key(kind):
            for key in rawDefs[kind].keys():
                if (defs.has_key(key)):
                    raise Exception("duplicate name")
                else:
                    defs[key]=rawDefs[kind][key]
    #print "def=",defs
    insertSynopsis(defs,sfile,outfile,logFile)
    

def addSynopsisInDir(interfaceDir, outDir,filePaths,logFile=sys.stdout):
    """Add the synopsis to the files in filePaths, reading the interfaces from interfaceDir, and writing them to outDir.
    
    interfaceDir should be the path of the directory where the .int interface files are.
    outDir should be the path of an existing directory.
    filepath should be a list the the path of the source files to be processed (with extension)
    """
    fileMapping={}
    for sourceFilePath in filePaths:
        import string, os.path
        try:
            fileName=os.path.basename(sourceFilePath)
            fileName=fileName[:string.rfind(fileName,".")]
            interfaceFilePath=interfaceDir+"/"+fileName+".int"
            outFilePath=os.path.join(outDir,os.path.basename(sourceFilePath))
            logFile.write("=== began brocessing of"+sourceFilePath+"\n")
            ifile = open(interfaceFilePath,'r')
            sfile = open(sourceFilePath,'r')
            outfile= open(outFilePath,'w')
            addSynopsisToFile(ifile,sfile,outfile,logFile=logFile)
            fileMapping[sfile]=outfile
        except:
            import sys, traceback
            logFile.write('-'*60+"\n")
            traceback.print_exc(file=logFile)
            logFile.write('-'*60+"\n")
            if os.access(outFilePath,os.F_OK):
                try:
                    os.rename(outFilePath,outFilePath+".err")
                except:
                    logFile.write("+++ error renaming"+outFilePath+"\n")
    return fileMapping

if __name__ == '__main__':
    if len(sys.argv)<4:
        print "usage:", sys.argv[0]," interface_dir out_dir sourcefile1.F [sourcefile2.F ...]"
    else:
        interfaceDir=sys.argv[1]
        outDir=sys.argv[2]
        addSynopsisInDir(interfaceDir, outDir,sys.argv[3:],sys.stdout)

