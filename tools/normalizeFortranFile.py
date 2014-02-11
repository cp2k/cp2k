#! /usr/bin/env python

import sys
import re
import string
from sys import argv
from cStringIO import StringIO

rUse=0
rVar=0
varRe=re.compile(r" *(?P<var>[a-zA-Z_0-9]+) *(?P<rest>(?:\((?P<param>(?:[^()]+|\((?:[^()]+|\([^()]*\))*\))*)\))? *(?:= *(?P<value>(:?[^\"',()]+|\((?:[^()\"']+|\([^()\"']*\)|\"[^\"]*\"|'[^']*')*\)|\"[^\"]*\"|'[^']*')+))?)? *(?:(?P<continue>,)|\n?) *",re.IGNORECASE)
useParseRe=re.compile(
    r" *use +(?P<module>[a-zA-Z_][a-zA-Z_0-9]*)(?P<only> *, *only *:)? *(?P<imports>.*)$",
    flags=re.IGNORECASE)
commonUsesRe=re.compile("^#include *\"cp_common_uses.h\"")
localNameRe=re.compile(" *(?P<localName>[a-zA-Z_0-9]+)(?: *= *> *[a-zA-Z_0-9]+)? *$")

def readFortranLine(infile):
    """Reads a group of connected lines (connected with &)
    returns a touple with the joined line, and a list with the original lines.
    Doesn't support multiline character constants!"""
    lineRe=re.compile(
        r"(?:(?P<preprocessor>#.*\n?)| *(&)?(?P<core>(?:!\$|[^&!\"']+|\"[^\"]*\"|'[^']*')*)(?P<continue>&)? *(?P<comment>!.*)?\n?)",#$
        re.IGNORECASE)
    joinedLine=""
    comments=None
    lines=[]
    continuation=0
    while 1:
        line=infile.readline().replace("\t",8*" ")
        if not line: break
        lines.append(line)
        m=lineRe.match(line)
        if not m or m.span()[1]!=len(line):
            raise SyntaxError("unexpected line format:"+repr(line))
        if m.group("preprocessor"):
            if len(lines)>1:
                raise SyntaxError("continuation to a preprocessor line not supported "+repr(line))
            comments=line
            break
        coreAtt=m.group("core")
        joinedLine = joinedLine.rstrip("\n") + coreAtt
        if coreAtt and not coreAtt.isspace(): continuation=0
        if m.group("continue"): continuation=1
        if m.group("comment"):
            if comments:
                comments+="\n"+m.group("comment")
            else:
                comments=m.group("comment")
        if not continuation: break
    return (joinedLine,comments,lines)

def parseRoutine(inFile):
    """Parses a routine"""
    startRe=re.compile(r" *(?:recursive +|pure +|elemental +)*(?:subroutine|function)",re.IGNORECASE)
    endRe=re.compile(r" *end (?:subroutine|function)",re.IGNORECASE)
    startRoutineRe=re.compile(r" *(?:recursive +|pure +|elemental +)*(?P<kind>subroutine|function) +(?P<name>[a-zA-Z_][a-zA-Z_0-9]*) *(?:\((?P<arguments>[^()]*)\))? *(?:result *\( *(?P<result>[a-zA-Z_][a-zA-Z_0-9]*) *\))? *(?:bind *\([^()]+\))? *\n?",re.IGNORECASE)#$
    typeBeginRe=re.compile(r" *(?P<type>integer(?: *\* *[0-9]+)?|logical|character(?: *\* *[0-9]+)?|real(?: *\* *[0-9]+)?|complex(?: *\* *[0-9]+)?|type)",
                           re.IGNORECASE)
    typeRe=re.compile(r" *(?P<type>integer(?: *\* *[0-9]+)?|logical|character(?: *\* *[0-9]+)?|real(?: *\* *[0-9]+)?|complex(?: *\* *[0-9]+)?|type) *(?P<parameters>\((?:[^()]+|\((?:[^()]+|\([^()]*\))*\))*\))? *(?P<attributes>(?: *, *[a-zA-Z_0-9]+(?: *\((?:[^()]+|\((?:[^()]+|\([^()]*\))*\))*\))?)+)? *(?P<dpnt>::)?(?P<vars>[^\n]+)\n?",re.IGNORECASE)#$
    attributeRe=re.compile(r" *, *(?P<attribute>[a-zA-Z_0-9]+) *(?:\( *(?P<param>(?:[^()]+|\((?:[^()]+|\([^()]*\))*\))*)\))? *",re.IGNORECASE)
    ignoreRe=re.compile(r" *(?:|implicit +none *)$",re.IGNORECASE)
    interfaceStartRe=re.compile(r" *interface *$",re.IGNORECASE)
    interfaceEndRe=re.compile(r" *end +interface *$",re.IGNORECASE)
    routine={'preRoutine':[],
             'core':[],
             'strippedCore':[],
             'begin':[],
             'end':[],
             'preDeclComments':[],
             'declarations':[],
             'declComments':[],
             'parsedDeclarations':[],
             'postRoutine':[],
             'kind':None,'name':None,'arguments':None,'result':None,
             'interfaceCount':0,
             'use':[]
             }
    includeRe=re.compile(r"#? *include +[\"'](?P<file>.+)[\"'] *$",re.IGNORECASE)
    while 1:
        (jline,comments,lines)=readFortranLine(inFile)
        if len(lines)==0: break
        if startRe.match(jline):break
        routine['preRoutine'].extend(lines)
        m=includeRe.match(lines[0])
        if m:
            try:
                subF=file(m.group('file'))
                while 1:
                    (subjline,subcomments,sublines)=readFortranLine(subF)
                    if not sublines:
                        break
                    routine['strippedCore'].append(subjline)
                subF.close()
            except:
                import traceback
                print "error trying to follow include ",m.group('file')
                print "warning this might lead to the removal of used variables"
                traceback.print_exc()
    if jline:
        routine['begin']=lines
        m=startRoutineRe.match(jline)
        if not m or m.span()[1]!=len(jline):
            raise SyntaxError("unexpected subroutine start format:"+repr(lines))
        routine['name']=m.group('name')
        routine['kind']=m.group('kind')
        if (m.group('arguments') and m.group('arguments').strip()):
            routine['arguments']=map(lambda x: x.strip(),
                                     m.group('arguments').split(","))
        if (m.group('result')):
            routine['result']=m.group('result')
        if (not routine['result'])and(routine['kind'].lower()=="function"):
            routine['result']=routine['name']
    while 1:
        (jline,comments,lines)=readFortranLine(inFile)
        if len(lines)==0: break
        if not ignoreRe.match(jline):
            if typeBeginRe.match(jline):
                m=typeRe.match(jline)
                if (m.group('type').lower()=='type' and
                    not m.group('parameters')):
                    break
                if not m or m.span()[1]!=len(jline):
                    raise SyntaxError("unexpected type format:"+repr(jline))
                decl={'type':m.group("type"),
                      'parameters':None,
                      'attributes':[],
                      'vars':[]}
                if m.group('parameters'):
                    decl['parameters']=(m.group("parameters").replace(" ","").
                                        replace(",",", "))
                str=m.group("attributes")
                while(str):
                    m2=attributeRe.match(str)
                    if not m2:
                        raise SyntaxError("unexpected attribute format "+
                                          repr(str)+" in "+repr(lines))
                    decl['attributes'].append(m2.group().replace(" ","").
                                              replace(",",", ")[2:])
                    str=str[m2.span()[1]:]
                str=m.group("vars")
                while 1:
                    m2=varRe.match(str)
                    if not m2:
                        raise SyntaxError("unexpected var format "+
                                          repr(str)+" in "+repr(lines))
                    var=m2.group("var")
                    if m2.group("param"):var+="("+m2.group("param")+")"
                    if m2.group("value"):
                        var+=" = "
                        var+=m2.group("value")
                    decl['vars'].append(var)
                    str=str[m2.span()[1]:]
                    if not m2.group("continue"):
                        if str:
                            raise SyntaxError("error parsing vars (leftover="+
                                              repr(str)+") in "+repr(lines))
                        break
                routine['parsedDeclarations'].append(decl)
            elif interfaceStartRe.match(jline):
                istart=lines
                interfaceDeclFile=StringIO()
                while 1:
                    (jline,comments,lines)=readFortranLine(inFile)
                    if interfaceEndRe.match(jline):
                        iend=lines
                        break
                    interfaceDeclFile.writelines(lines)
                interfaceDeclFile=StringIO(interfaceDeclFile.getvalue())
                iroutines=[]
                while 1:
                    iroutine=parseRoutine(interfaceDeclFile)
                    if not iroutine['kind']:
                        if len(iroutines)==0:
                            interfaceDeclFile.seek(0)
                            raise SyntaxError("error parsing interface:"+
                                              repr(interfaceDeclFile.read()))
                        iroutines[-1]['postRoutine'].extend(iroutine['preRoutine'])
                        break
                    iroutines.append(iroutine)
                for iroutine in iroutines:
                    routine['interfaceCount']+=1
                    decl={'type':'z_interface%02d'%(routine['interfaceCount']),
                          'parameters':None,
                          'attributes':[],
                          'vars':[iroutine['name']],
                          'iroutine':iroutine,
                          'istart':istart,
                          'iend':iend
                          }
                    routine['parsedDeclarations'].append(decl)
            elif useParseRe.match(jline):
                routine['use'].append("".join(lines))
            else:
                break
        routine['declarations'].append("".join(lines))
        if (len(routine['parsedDeclarations'])==0 and len(routine['use'])==0 and
            not re.match(" *implicit +none *$",jline,re.IGNORECASE)):
            routine['preDeclComments'].append("".join(lines))
        elif comments:
            routine['declComments'].append(comments)
    containsRe=re.compile(r" *contains *$",re.IGNORECASE)
    while len(lines)>0:
        if endRe.match(jline):
            routine['end']=lines
            break
        routine['strippedCore'].append(jline)
        routine['core'].append("".join(lines))
        if containsRe.match(lines[0]):
            break
        m=includeRe.match(lines[0])
        if m:
            try:
                subF=file(m.group('file'))
                while 1:
                    (subjline,subcomments,sublines)=readFortranLine(subF)
                    if not sublines:
                        break
                    routine['strippedCore'].append(subjline)
                subF.close()
            except:
                import traceback
                print "error trying to follow include ",m.group('file')
                print "warning this might lead to the removal of used variables"
                traceback.print_exc()
        (jline,comments,lines)=readFortranLine(inFile)
    return routine

def findWord(word,text,options=re.IGNORECASE):
    """Returns the position of word in text or -1 if not found.
    A match is valid only if it is a whole word (i.e. findWord('try','retry')
    returns false)"""
    wordRe=re.compile("(?<![a-zA-Z_0-9%])"+word+"(?![a-zA-Z_0-9])|(?<=[0-9.]_)"+word+"(?![a-zA-Z_0-9])",options)
    m=wordRe.search(text)
    if m:
        pos=m.span()[0]
    else:
        pos=-1
    return pos

def enforceDeclDependecies(declarations):
    """enforces the dependencies between the vars
    and compacts the declarations, returns the variables needed by other variables"""
    idecl=0
    ii=0
    while idecl<len(declarations):
        typeParam="".join(declarations[idecl]['attributes'])
        if declarations[idecl]['parameters']:
            typeParam+=" "+declarations[idecl]['parameters']
        typeParam=typeParam.lower()

        ivar=0
        while ivar<len(declarations[idecl]['vars']):
            moved=0
            m=varRe.match(declarations[idecl]['vars'][ivar])
            if not m:
                raise SyntaxError('could not match var '+repr(declarations[idecl]['vars'][ivar]))
            rest=m.group("rest")
            rest=rest.lower()
            if rest:
                for ivar2 in range(ivar+1,len(declarations[idecl]['vars'])):
                    m=varRe.match(declarations[idecl]['vars'][ivar2])
                    if findWord(m.group('var').lower(),rest)!=-1:
                        moved=ivar2+1
            if moved:
                declarations[idecl]['vars'][moved:moved]=[
                    declarations[idecl]['vars'][ivar]]
                del declarations[idecl]['vars'][ivar]
            else:
                for idecl2 in range(idecl+1,len(declarations)):
                    for ivar2 in range(len(declarations[idecl2]['vars'])):
                        ii+=1
                        if ii>100000:
                            raise StandardError,"could not enforce all constraints"
                        m=varRe.match(declarations[idecl2]['vars'][ivar2])
                        if (ivar==0 and
                            findWord(m.group('var').lower(),typeParam)!=-1):
                            declarations.insert(idecl2+1,declarations[idecl])
                            del declarations[idecl]
                            ivar=0
                            moved=1
                            break
                        if rest and findWord(m.group('var').lower(),rest)!=-1:
                            if len(declarations[idecl]['vars'])>1:
                                newDecl={}
                                newDecl.update(declarations[idecl])
                                newDecl['vars']=[declarations[idecl]['vars'][ivar]]
                                declarations.insert(idecl2+1,newDecl)
                                del declarations[idecl]['vars'][ivar]
                            else:
                                declarations.insert(idecl2+1,
                                                    declarations[idecl])
                                del declarations[idecl]
                                ivar=0
                            moved=1
                            break
                    if moved:
                        break
            if not moved: ivar+=1
        idecl+=1

    for i in range(len(declarations)-1,0,-1):
        if (declarations[i]['normalizedType'].lower()==
            declarations[i-1]['normalizedType'].lower()):
            declarations[i-1]['vars'].extend(declarations[i]['vars'])
            del declarations[i]

def sortDeclarations(declarations):
    """sorts, compacts declarations and respects dependencies
    normalizedType has to be defined for the declarations"""

    declarations.sort(lambda x,y:cmp(x['normalizedType'].lower(),
                                  y['normalizedType'].lower()))

    for i in range(len(declarations)-1,0,-1):
        if (declarations[i]['normalizedType'].lower()==
            declarations[i-1]['normalizedType'].lower()):
            declarations[i-1]['vars'].extend(declarations[i]['vars'])
            del declarations[i]

    for decl in declarations:
        decl['vars'].sort(lambda x,y:cmp(x.lower(),y.lower()))
    enforceDeclDependecies(declarations)

def writeRoutine(routine, outFile):
    """writes the given routine to outFile"""
    outFile.writelines(routine["preRoutine"])
    outFile.writelines(routine["begin"])
    outFile.writelines(routine["declarations"])
    outFile.writelines(routine["core"])
    outFile.writelines(routine["end"])
    outFile.writelines(routine["postRoutine"])

def writeInCols(dLine,indentCol,maxCol,indentAtt,file):
    """writes out the strings (trying not to cut them) in dLine up to maxCol
    indenting each newline with indentCol.
    The '&' of the continuation line is at maxCol.
    indentAtt is the actual intent, and the new indent is returned"""
    strRe=re.compile(r"('[^'\n]*'|\"[^\"\n]*\")")
    nonWordRe=re.compile(r"(\(/|/\)|[^-+a-zA-Z0-9_.])")
    maxSize=maxCol-indentCol-1
    tol=min(maxSize/6,6)+indentCol
    for fragment in dLine:
        if indentAtt+len(fragment)<maxCol:
            file.write(fragment)
            indentAtt+=len(fragment)
        elif len(fragment.lstrip())<=maxSize:
            file.write("&\n"+(" "*indentCol))
            file.write(fragment.lstrip())
            indentAtt=indentCol+len(fragment.lstrip())
        else:
            sPieces=strRe.split(fragment)
            for sPiece in sPieces:
                if sPiece and (not (sPiece[0]=='"' or sPiece[0]=="'")):
                    subPieces=nonWordRe.split(sPiece)
                else:
                    subPieces=[sPiece]
                for subPiece in subPieces:
                    if indentAtt==indentCol:
                        file.write(subPiece.lstrip())
                        indentAtt+=len(subPiece.lstrip())
                    elif indentAtt<tol or indentAtt+len(subPiece)<maxCol:
                        file.write(subPiece)
                        indentAtt+=len(subPiece)
                    else:
                        file.write("&\n"+(" "*indentCol))
                        file.write(subPiece.lstrip())
                        indentAtt=indentCol+len(subPiece.lstrip())
    return indentAtt

def writeCompactDeclaration(declaration,file):
    """Writes a declaration in a compact way"""
    d=declaration
    if d.has_key('iroutine'):
        file.writelines(d['istart'])
        writeRoutine(d['iroutine'],file)
        file.writelines(d['iend'])
    else:
        dLine=[]
        if len(d['vars'])>0:
            dLine.append("    "+d['type'])
            if d['parameters']: # do not drop empty parameter lists?
                dLine.append(d['parameters'])
            if d['attributes']:
                for a in d['attributes']:
                    dLine[-1:]=[dLine[-1]+", "]
                    dLine.append(a)
            dLine[-1:]=[dLine[-1]+" :: "]
            decl=[]
            decl.extend(dLine)
            lVars=0
            for var in d['vars'][:-1]:
                if lVars>600:
                    dLine[-1]=dLine[-1][:-2]
                    writeInCols(dLine,6,79,0,file)
                    file.write("\n")
                    lVars=0
                    dLine=[]
                    dLine.extend(decl)
                dLine.append(var+", ")
                lVars+=len(var)+2
            dLine.append(d['vars'][-1])

        writeInCols(dLine,6,79,0,file)
        file.write("\n")

def writeExtendedDeclaration(declaration,file):
    """Writes a declaration in a nicer way (using more space)"""
    d=declaration
    if len(d['vars'])==0: return
    if d.has_key('iroutine'):
        file.writelines(d['istart'])
        writeRoutine(d['iroutine'],file)
        file.writelines(d['iend'])
    else:
        dLine=[]
        dLine.append("    "+d['type'])
        if d['parameters']: # do not drop empty parameter lists?
            dLine.append(d['parameters'])
        if d['attributes']:
            for a in d['attributes']:
                dLine[-1:]=[dLine[-1]+", "]
                dLine.append(a)

        indentAtt=writeInCols(dLine,6,45,0,file)
        file.write(" "*(44-indentAtt))
        file.write(" :: ")
        indentAtt=48

        dLine=[]
        for var in d['vars'][:-1]:
            dLine.append(var+", ")
        dLine.append(d['vars'][-1])

        writeInCols(dLine,48,79,indentAtt,file)
        file.write("\n")

def writeDeclarations(parsedDeclarations,file):
    """Writes the declarations to the given file"""
    for d in parsedDeclarations:
        maxLenVar=0
        totalLen=0
        for v in d['vars']:
            maxLenVar=max(maxLenVar,len(v))
            totalLen+=len(v)
        if maxLenVar>30 or totalLen>75:
            writeCompactDeclaration(d,file)
        else:
            writeExtendedDeclaration(d,file)

def cleanDeclarations(routine,logFile=sys.stdout):
    """cleans up the declaration part of the given parsed routine
    removes unused variables"""
    global rVar
    containsRe=re.compile(r" *contains *$",re.IGNORECASE)
    if routine['core']:
        if containsRe.match(routine['core'][-1]):
            logFile.write("*** routine %s contains other routines ***\n*** declarations not cleaned ***\n"%
                (routine['name']))
            return
    commentToRemoveRe=re.compile(r" *! *(?:interface|arguments|parameters|locals?|\** *local +variables *\**|\** *local +parameters *\**) *$",re.IGNORECASE)
    nullifyRe=re.compile(r" *nullify *\(([^()]+)\) *\n?",re.IGNORECASE|re.MULTILINE)

    if not routine['kind']: return
    if (routine['core']):
        if re.match(" *type *[a-zA-Z_]+ *$",routine['core'][0],re.IGNORECASE):
            logFile.write("*** routine %s contains local types, not fully cleaned ***\n"%
                      (routine['name']))
        if re.match(" *import+ *$",routine['core'][0],re.IGNORECASE):   
            logFile.write("*** routine %s contains import, not fully cleaned ***\n"%
                      (routine['name']))
    if re.search("^#","".join(routine['declarations']),re.MULTILINE):
        logFile.write("*** routine %s declarations contain preprocessor directives ***\n*** declarations not cleaned ***\n"%(
            routine['name']))
        return
    try:
        rest="".join(routine['strippedCore']).lower()
        nullifys=",".join(nullifyRe.findall(rest))
        rest=nullifyRe.sub("",rest)
        paramDecl=[]
        decls=[]
        for d in routine['parsedDeclarations']:
            d['normalizedType']=d['type']
            if d['parameters']:
                d['normalizedType']+=d['parameters']
            if (d["attributes"]):
                d['attributes'].sort(lambda x,y:cmp(x.lower(),y.lower()))
                d['normalizedType']+=', '
                d['normalizedType']+=', '.join(d['attributes'])
            if "parameter" in map(str.lower,d['attributes']):
                paramDecl.append(d)
            else:
                decls.append(d)

        sortDeclarations(paramDecl)
        sortDeclarations(decls)
        has_routinen=0
        pos_routinep=-1
        for d in paramDecl:
            for i in xrange(len(d['vars'])):
                v=d['vars'][i]
                m=varRe.match(v)
                lowerV=m.group("var").lower()
                if lowerV=="routinen":
                    has_routinen=1
                    d['vars'][i]="routineN = '"+routine['name']+"'"
                elif lowerV=="routinep":
                    pos_routinep=i
                    d['vars'][i]="routineP = moduleN//':'//routineN"
            if not has_routinen and pos_routinep>=0:
                d['vars'].insert(pos_routinep,"routineN = '"+routine['name']+"'")


        if routine['arguments']:
            routine['lowercaseArguments']=map(lambda x:x.lower(),routine['arguments'])
        else:
            routine['lowercaseArguments']=[]
        if routine['result']: routine['lowercaseArguments'].append(routine['result'].lower())
        argDeclDict={}
        localDecl=[]
        for d in decls:
            localD={}
            localD.update(d)
            localD['vars']=[]
            argD=None
            for v in d['vars']:
                m=varRe.match(v)
                lowerV=m.group("var").lower()
                if lowerV in routine['lowercaseArguments']:
                    argD={}
                    argD.update(d)
                    argD['vars']=[v]
                    if argDeclDict.has_key(lowerV):
                        raise SyntaxError(
                            "multiple declarations not supported. var="+v+
                            " declaration="+str(d)+"routine="+routine['name'])
                    argDeclDict[lowerV]=argD
                else:
                    pos=findWord(lowerV,rest)
                    if (pos!=-1):
                        localD['vars'].append(v)
                    else:
                        if findWord(lowerV,nullifys)!=-1:
                            if not rmNullify(lowerV,routine['core']):
                                raise SyntaxError(
                                    "could not remove nullify of "+lowerV+
                                    " as expected, routine="+routine['name'])
                        logFile.write("removed var %s in routine %s\n" %
                                      (lowerV,routine['name']))
                        rVar+=1
            if (len(localD['vars'])):
                localDecl.append(localD)
        argDecl=[]
        for arg in routine['lowercaseArguments']:
            if argDeclDict.has_key(arg):
                argDecl.append(argDeclDict[arg])
            else:
                print "warning, implicitly typed argument '",arg,"' in routine",routine['name']
        if routine['kind'].lower()=='function':
            aDecl=argDecl[:-1]
        else:
            aDecl=argDecl

        # try to have arg/param/local, but checks for dependencies arg/param and param/local
        argDecl.extend(paramDecl)
        enforceDeclDependecies(argDecl)
        splitPos=0
        for i in xrange(len(argDecl)-1,-1,-1):
            if not 'parameter' in map(str.lower,argDecl[i]['attributes']):
                splitPos=i+1
                break
        paramDecl=argDecl[splitPos:]
        argDecl=argDecl[:splitPos]
        paramDecl.extend(localDecl)
        enforceDeclDependecies(paramDecl)
        splitPos=0
        for i in xrange(len(paramDecl)-1,-1,-1):
            if 'parameter' in map(str.lower,paramDecl[i]['attributes']):
                splitPos=i+1
                break
        localDecl=paramDecl[splitPos:]
        paramDecl=paramDecl[:splitPos]

        newDecl=StringIO()
        for comment in routine['preDeclComments']:
            if not commentToRemoveRe.match(comment):
                newDecl.write(comment)
        newDecl.writelines(routine['use'])
        writeDeclarations(argDecl,newDecl)
        if argDecl and paramDecl:
            newDecl.write("\n")
        writeDeclarations(paramDecl,newDecl)
        if (argDecl or paramDecl) and localDecl:
            newDecl.write("\n")
        writeDeclarations(localDecl,newDecl)
        if argDecl or paramDecl or localDecl:
            newDecl.write("\n")
        wrote=0
        for comment in routine['declComments']:
            if not commentToRemoveRe.match(comment):
                newDecl.write(comment)
                newDecl.write("\n")
                wrote=1
        if wrote:
            newDecl.write("\n")
        routine['declarations']=[newDecl.getvalue()]
    except:
        if routine.has_key('name'):
            logFile.write("**** exception cleaning routine "+routine['name']+" ****")
        logFile.write("parsedDeclartions="+str(routine['parsedDeclarations']))
        raise

def rmNullify(var,strings):
    removed=0
    var=var.lower()
    nullifyRe=re.compile(r" *nullify *\(", re.IGNORECASE)
    nullify2Re=re.compile(r"(?P<nullif> *nullify *\()(?P<vars>[^()!&]+)\)",re.IGNORECASE)

    for i in xrange(len(strings)-1,-1,-1):
        line=strings[i]
        comments=[]
        if nullifyRe.match(line) and findWord(var,line)!=-1:
            core=""
            comments=[]
            for l in line.splitlines():
                pos=l.find("&")
                pos2=l.find("!")
                if pos==-1:
                    if pos2==-1:
                        core+=l
                    else:
                        core+=l[:pos2]
                        comments.append(l[pos2:]+"\n")
                else:
                    core+=l[:pos]
                    if pos2!=-1:
                        comments.append(l[pos2:]+"\n")
            m=nullify2Re.match(core)
            if not m:
                raise SyntaxError("could not match nullify to "+repr(core)+
                                  "in"+repr(line))
            allVars=[]
            vars=m.group("vars")
            v= map(string.strip,vars.split(","))
            removedNow=0
            for j in xrange(len(v)-1,-1,-1):
                if findWord(var,v[j].lower())!=-1:
                    del v[j]
                    removedNow=1
            if removedNow:
                if len(v)==0:
                    if not comments:
                        del strings[i]
                    else:
                        strings[i]="".join(comments)
                else:
                    for j in xrange(len(v)-1):
                        v[j]+=", "
                    v[-1]+=")"
                    newS=StringIO()
                    v.insert(0,m.group("nullif"))
                    writeInCols(v,len(v[0])-len(v[0].lstrip())+5,77,0,newS)
                    newS.write("\n")
                    if comments:
                        for c in comments:
                            newS.write(c)
                    strings[i]=newS.getvalue()
                removed+=1
    return removed

def parseUse(inFile):
    """Parses the use statements in inFile
    The parsing stops at the first non use statement.
    Returns something like:
    ([{'module':'module1','only':['el1','el2=>el3']},...],
     '! comment1\\n!comment2...\\n',
     'last line (the line that stopped the parsing)')
    """
    lineNr=0
    preComments=[]
    modules=[]
    origLines=[]
    commonUses=""
    while 1:
        (jline,comments,lines)=readFortranLine(inFile)
        lineNr=lineNr+len(lines)
        if not lines: break
        origLines.append("".join(lines))
        # parse use
        m=useParseRe.match(jline)
        if m:
            if comments:
                print "jline",jline,"lines",lines
            useAtt={'module':m.group('module'),'comments':[]}

            if m.group('only'):
                useAtt['only']=map(string.strip,
                                   string.split(m.group('imports'),','))
            else:
                useAtt['renames']=map(string.strip,
                                      string.split(m.group('imports'),','))
                if useAtt['renames']==[""]: del useAtt['renames']
            if comments : useAtt['comments'].append(comments)
            # add use to modules
            modules.append(useAtt)
        elif jline and not jline.isspace():
            break
        else:
            if comments and commonUsesRe.match(comments):
                commonUses="".join(lines)
            elif len(modules)==0:
                preComments.append(("".join(lines)))
            elif comments:
                modules[-1]['comments'].append(comments)

    return {'modules':modules,'preComments':preComments,'commonUses':commonUses,
            'postLine':"".join(lines),'origLines':origLines[:-1]}

def normalizeModules(modules):
    """Sorts the modules and their export and removes duplicates.
    renames aren't sorted correctly"""
    # orders modules
    modules.sort(lambda x,y:cmp(x['module'],y['module']) )
    for i in range(len(modules)-1,0,-1):
        if modules[i]['module']==modules[i-1]['module']:
            if not (modules[i-1].has_key('only') and
                    modules[i].has_key('only')):
                raise SyntaxError('rejoining of module '+
                                  str(modules[i]['module'])+
                                  ' failed as at least one of the use is not a use ...,only:')
            modules[i-1]['only'].extend(modules[i]['only'])
            del modules[i]
    # orders imports
    for m in modules:
        if m.has_key('only'):
            m['only'].sort()
            for i in range(len(m['only'])-1,0,-1):
                if m['only'][i-1]==m['only'][i]: del m['only'][i]

def writeUses(modules,outFile):
    """Writes the use declaration using a long or short form depending on how
    many only statements there are"""
    for m in modules:
        if m.has_key('only') and len(m['only'])>8:
            writeUseShort(m,outFile)
        else:
            writeUseLong(m,outFile)

def writeUseLong(m,outFile):
    """Writes a use declaration in a nicer, but longer way"""
    if m.has_key('only'):
        outFile.write("  USE "+m['module']+","+
                      string.rjust('ONLY: ',38-len(m['module'])))
        if m['only']: outFile.write(m['only'][0])
        for i in range(1,len(m['only'])):
            outFile.write(",&\n"+string.ljust("",45)+m['only'][i])
    else:
        outFile.write("  USE "+m['module'])
        if m.has_key('renames') and m['renames']:
            outFile.write(","+string.ljust("",38)+
                          m['renames'][0])
            for i in range(1,len(m['renames'])):
                outFile.write(",&\n"+string.ljust("",45)+m['renames'][i])
    if m['comments']:
        outFile.write("\n")
        outFile.write('\n'.join(m['comments']))
    outFile.write("\n")

def writeUseShort(m,file):
    """Writes a use declaration in a compact way"""
    uLine=[]
    if m.has_key('only'):
        file.write("  USE "+m['module']+","+
                   string.rjust('ONLY: &\n',40-len(m['module'])))
        for k in m['only'][:-1]:
            uLine.append(k+", ")
        uLine.append(m['only'][-1])
        uLine[0]=" "*7+uLine[0]
    elif m.has_key('renames') and m['renames']:
        uLine.append("  USE "+m['module']+", ")
        for k in m['renames'][:-1]:
            uLine.append(k+", ")
        uLine.append(m['renames'][-1])
    else:
        uLine.append("  USE "+m['module'])
    writeInCols(uLine,7,79,0,file)
    if m['comments']:
        file.write("\n")
        file.write('\n'.join(m['comments']))
    file.write("\n")

def prepareImplicitUses(modules):
    """Transforms a modulesDict into an implictUses (dictionary of module names
    each containing a dictionary with the only, and the special key '_WHOLE_'
    wich is true if the whole mosule is implicitly present"""
    mods={}
    for m in modules:
        m_name=m['module'].lower()
        if (not mods.has_key(m_name)):
            mods[m['module']]={'_WHOLE_':0}
        m_att=mods[m_name]
        if m.has_key('only'):
            for k in m['only']:
                m=localNameRe.match(k)
                if not m:
                    raise SyntaxError('could not parse use only:'+repr(k))
                impAtt=m.group('localName').lower()
                m_att[impAtt]=1
        else:
            m_att['_WHOLE_']=1
    return mods

def cleanUse(modulesDict,rest,implicitUses=None,logFile=sys.stdout):
    """Removes the unneded modules (the ones that are not used in rest)"""
    global rUse
    exceptions={}
    modules=modulesDict['modules']
    rest=rest.lower()
    for i in range(len(modules)-1,-1,-1):
        m_att={}
        m_name=modules[i]['module'].lower()
        if implicitUses and implicitUses.has_key(m_name):
            m_att=implicitUses[m_name]
        if m_att.has_key('_WHOLE_') and m_att['_WHOLE_']:
            rUse+=1
            logFile.write("removed USE of module "+m_name+"\n")
            del modules[i]
        elif modules[i].has_key("only"):
            els=modules[i]['only']
            for j in range(len(els)-1,-1,-1):
                m=localNameRe.match(els[j])
                if not m:
                    raise SyntaxError('could not parse use only:'+repr(els[j]))
                impAtt=m.group('localName').lower()
                if m_att.has_key(impAtt):
                    rUse+=1
                    logFile.write("removed USE "+m_name+", only: "+repr(els[j])+"\n")
                    del els[j]
                elif not exceptions.has_key(impAtt):
                    if findWord(impAtt,rest)==-1:
                        rUse+=1
                        logFile.write("removed USE "+m_name+", only: "+repr(els[j])+"\n")
                        del els[j]
            if len(modules[i]['only'])==0:
                if modules[i]['comments']:
                    modulesDict['preComments'].extend(
                        map(lambda x:x+"\n",modules[i]['comments']))
                del modules[i]

def resetModuleN(moduleName,lines):
    "resets the moduleN variable to the module name in the lines lines"
    moduleNRe=re.compile(r".*:: *moduleN *= *(['\"])[a-zA-Z_0-9]+\1",
                         flags=re.IGNORECASE)
    for i in xrange(len(lines)):
        lines[i]=moduleNRe.sub(
            "  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = '"+moduleName+"'",
            lines[i])

def rewriteFortranFile(inFile,outFile,logFile=sys.stdout,orig_filename=None):
    """rewrites the use statements and declarations of inFile to outFile.
    It sorts them and removes the repetitions."""
    import os.path
    moduleRe=re.compile(r" *(?:module|program) +(?P<moduleName>[a-zA-Z_][a-zA-Z_0-9]*) *(?:!.*)?$",
                        flags=re.IGNORECASE)
    commonUsesIncludeFilepath=os.path.join(
        os.path.split(os.path.abspath(inFile.name))[0],"cp_common_uses.h")
    coreLines=[]
    while 1:
        line=inFile.readline()
        if not line: break
        if line[0]=='#':
            coreLines.append(line)
        outFile.write(line)
        m=moduleRe.match(line)
        if m:
            if not orig_filename: orig_filename=inFile.name
            fn = os.path.basename(orig_filename).rsplit(".",1)[0]
            if (m.group('moduleName')!=fn) :
                raise SyntaxError("Module name is different from filename ("+
                                  m.group('moduleName')+"!="+fn+")")
            break
    try:
        modulesDict=parseUse(inFile)
        routines=[]
        coreLines.append(modulesDict['postLine'])
        routine=parseRoutine(inFile)
        coreLines.extend(routine['preRoutine'])
        if m:
            resetModuleN(m.group('moduleName'),routine['preRoutine'])
        routines.append(routine)
        while routine['kind']:
            routine=parseRoutine(inFile)
            routines.append(routine)
        map(lambda x:cleanDeclarations(x,logFile),routines)
        for routine in routines:
            coreLines.extend(routine['declarations'])
            coreLines.extend(routine['strippedCore'])
        rest="".join(coreLines)
        nonStPrep=0
        for line in modulesDict['origLines']:
            if (re.search('^#',line) and not commonUsesRe.match(line)):
                print 'noMatch',repr(line)
                nonStPrep=1
        if nonStPrep:
            logFile.write("*** use statements contains preprocessor directives, not cleaning ***")
            outFile.writelines(modulesDict['origLines'])
        else:
            implicitUses=None
            if modulesDict['commonUses']:
                try:
                    f=file(commonUsesIncludeFilepath)
                    implicitUsesRaw=parseUse(f)
                    f.close()
                    implicitUses=prepareImplicitUses(implicitUsesRaw['modules'])
                except:
                    print ("ERROR trying to parse use statements contained in common",
                           "uses precompiler file ", commonUsesIncludeFilepath)
                    raise
            cleanUse(modulesDict,rest,implicitUses=implicitUses,logFile=logFile)
            normalizeModules(modulesDict['modules'])
            outFile.writelines(modulesDict['preComments'])
            writeUses(modulesDict['modules'],outFile)
            outFile.write(modulesDict['commonUses'])
            if modulesDict['modules']:
                outFile.write('\n')
        outFile.write(modulesDict['postLine'])
        for routine in routines:
            writeRoutine(routine,outFile)
    except:
        import traceback
        logFile.write('-'*60+"\n")
        traceback.print_exc(file=logFile)
        logFile.write('-'*60+"\n")

        logFile.write("Processing file '"+inFile.name+"'\n")
        raise

if __name__ == '__main__':
    import os.path
    if len(sys.argv)<2:
        print "usage:", sys.argv[0]," out_dir file1 [file2 ...]"
    else:
        outDir=sys.argv[1]
        if not os.path.isdir(outDir):
            print "out_dir must be a directory"
            print "usage:", sys.argv[0]," out_dir file1 [file2 ...]"
        else:
            for fileName in sys.argv[2:]:
                try:
                    print "normalizing",fileName
                    infile=open(fileName,'r')
                    outfile=open(os.path.join(outDir,
                                              os.path.basename(fileName)),'w')
                    rewriteFortranFile(infile,outfile)
                    infile.close()
                    outfile.close()
                except:
                    print "error for file", fileName
            print "*** "*6
            print "removedUse=",rUse
            print "removedVar=",rVar
                # print "done"
