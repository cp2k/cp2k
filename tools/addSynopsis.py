#! /usr/bin/env python

import sys
from sys import argv

def parseInterface(iFile):
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
			e=SyntaxWarning("double definition of '"+name+"' ignored")
			del stack[0]
			raise e
		    kind_table[stack[0]["name"]]=stack[0] # put only the expansion??
		else:
		    print "ignoring group ",stack[0]
		del stack[0]
	    else:
		raise SyntaxError("inconsistent parse '"+str(match.groups())+"'")
	except SyntaxError, e:
	    e.lineno=lineNr
	    e.text=line
	    print "stack=",stack
	    raise
	except Warning, w:
	    print "ignoring warning at line "+lineNr+" of file "+iFile.name
	    print w
    return result

def writeSynopsis(objDef,outfile):
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
	print "no expansion for ",objDef

def insertSynopsis(defs,infile,outfile):
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
			elif defs[name]:
			    writeSynopsis(defs[name],outfile)
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
		    writeSynopsis(defs[name],outfile)
		    status=3
		if match.groups()[directivePos]=="NAME":
		    status=2
		    outfile.write(line.lstrip())
		elif match.groups()[directivePos]=="SYNOPSIS":
		    print "found synopsis at line",lineNr
		    if status==3 or status==4:
			status=4
			print "skipping synopsis"
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
	    print "ignoring warning at line "+str(lineNr)+" of file "+ infile.name
	    print w


if len(sys.argv)<4:
    print "usage:", sys.argv[0]," interface_dir out_dir sourcefile1.F [sourcefile2.F ...]"
else:
    interfaceDir=sys.argv[1]
    outDir=sys.argv[2]
    for sourceFilePath in sys.argv[3:]:
        import string, os.path
        try:
            fileName=os.path.basename(sourceFilePath)
            fileName=fileName[:string.rfind(fileName,".")]
            interfaceFilePath=interfaceDir+"/"+fileName+".int"
            outFilePath=outDir+"/"+os.path.basename(sourceFilePath)
            print "=== began brocessing of",sourceFilePath
            ifile = open(interfaceFilePath,'r')
            sfile = open(sourceFilePath,'r')
            outfile= open(outFilePath,'w')
            rawDefs=parseInterface(ifile)
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
            insertSynopsis(defs,sfile,outfile)
        except:
            import sys, traceback
            print '-'*60
            traceback.print_exc(file=sys.stdout)
            print '-'*60
            if os.access(outFilePath,os.F_OK):
                try:
                    os.rename(outFilePath,outFilePath+".err")
                except:
                    print "+++ error renaming",outFilePath
            #for obj in sys.exc_info():
            #    print str(obj)
            #tb=sys.exc_info()[2]
            #while 1:
            #    if tb==None: break
            #    print "at line "+str(tb.tb_lineno)+" in"
            #    tb=tb.tb_next
    
