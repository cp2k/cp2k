import sys,re,os.path,StringIO
import normalizeFortranFile
from normalizeFortranFile import readFortranLine
"""script to convert robodoc headers to doxygen headers

2008 fawzi
"""

def parseRoboDoc(lines,inFile):
    info={
        'special':{},'order':[],'trash':[],'name':[],'brief':[],
        'date':[],'author':[],'note':[],'todo':[],'params':[],'bug':[],
        'literature':[],'history':[],'see':[],'pname':[],'version':[],
        'warning':[],'origLines':[]}
    sAtt=info['brief']
    transl={
        'NAME':'name','FUNCTION':'brief','AUTHOR':'author',
        'CREATION DATE':'date','DESCRIPTION':'brief','PURPOSE':'brief',
        'INPUTS':'params','ARGUMENTS':'params','PARAMETERS':'params',
        'ATTRIBUTES':'params','BUGS':'bug','TODO':'todo',
        'NOTE':'note','NOTES':'note','LITERATURE':'literature',
        'REFERENCES':'literature','OUTPUTS':'params','AUTHORS':'author',
        'AURHOR':'author','WARNINGS':'warning','WARNING':'warning',
        'MODIFICATION HISTORY':'history','HISTORY':'history',
        'SEE ALSO':'see','SOURCE':'trash','SYNOPSIS':'trash'}
    roboDocRe=re.compile(" *!!?(.*)")
    headerRe=re.compile(r" *!!\*+[cdfhmstuv]?\** *(?P<name>[a-zA-Z_0-9/]+) *(?:\[(?P<version>[0-9.a-zA-Z_ ]+)\])?")
    labelRe=re.compile(" +(?P<label>NAME|COPYRIGHT|USAGE|FUNCTION|DESCRIPTION|PURPOSE|AUTHOR|CREATION DATE|MODIFICATION HISTORY|HISTORY|INPUTS|ARGUMENTS|OPTIONS|PARAMETERS|SWITCHES|OUTPUT|SYNOPSIS|SIDE EFFECTS|RESULT|RETURN VALUE|EXAMPLE|NOTES?|WARNINGS?|ERROR|DIAGNOSTICS|BUGS|TODO|IDEAS|PORTABILITY|SEE ALSO|METHODS|ATTRIBUTES|SOURCE|LITERATURE|TAGS|USED BY|[A-Z]+ ?[A-Z]* ?[A-Z]*) *$")
    fluffRe=re.compile(r" *([-+*#!= ]{3,})? *$")
    preprocessorRe=re.compile(r" *#")
    m=headerRe.match(lines[0])
    if m:
        if m.group('version'):
            info['version'].append(m.group('version'))
        if m.group('name'):
            info['pname'].append(m.group('name'))
        info['origLines'].append(lines[0])
        lines=lines[1:]
    while 1:
        for line in lines:
            info['origLines'].extend(lines)
            m=roboDocRe.match(line)
            if not m:
                raise SyntaxError('unexpectly out of robodoc')
            l=m.groups()[0]+'\n'
            m2=labelRe.match(l)
            if m2:
                label=m2.group('label')
                if transl.has_key(label):
                    sAtt=info[transl[label]]
                    info['order'].append(transl[label])
                else:
                    if info['order'] and info['order'][-1]=='author':
                        # print ('treating "%s"as author'%(label))
                        sAtt.append(l)
                    else:
                        label=label.title()
                        info['order'].append(label)
                        print ('WARNING, found odd label "%s", added to special'%(label))
                        if info['special'].has_key(label):
                            sAtt=info['special'][label]
                        else:
                            sAtt=[]
                            info['special'][label]=sAtt
            else:
                m=fluffRe.match(l)
                if m:
                    info['trash'].append(l)
                    if len(sAtt) and not sAtt[-1].isspace():
                        sAtt.append('\n')
                else:
                    sAtt.append(l.replace('@author',''))
        while 1:
            (jline,comments,lines)=readFortranLine(inFile)
            if not lines:
                break
            if jline!=None and jline!="" and not jline.isspace() or preprocessorRe.match(lines[0]):
                break
            if roboDocRe.match(lines[0]):
                break
            else:
                for l in lines:
                    if not l.isspace():
                        print "rmv",repr(l)
            info['origLines'].extend(lines)
        if not lines:
            break
        if jline!=None and jline!="" and not jline.isspace() or preprocessorRe.match(lines[0]):
            break
    if not info['name'] and info['pname']:
        info['name'].append(os.path.basename(info['pname'][0]))
    return (jline,comments,lines,info)

def non_empty(linee):
    for l in linee:
        if not l.isspace():
            return 1
    return 0
def drop_trailing(linee):
    lines=list(linee)
    for i in xrange(len(lines)-1,-1,-1):
        if lines[i].isspace():
            del lines[i]
        else:
            break
    while len(lines)>0 and lines[0].isspace():
        del lines[0]
    return lines

def printRoboDoc(info,outF=sys.stdout,isspace=0):
    # print info['name'],info['pname']
    toRm=['trash','pname','name']
    pre=['brief','special','params']
    fixedOrder=['note','todo','warning','bug',
        'literature','history','author','date','see']
    normal=[]
    normal.extend(toRm)
    normal.extend(pre)
    normal.extend(fixedOrder)
    rest=[]
    for o in info['order']:
        if not o in normal or o=='special':
            rest.append(o)
    res=[]
    pre="!> "
    if not info['name']:
        print 'ignoring robodoc like comment'
        for l in info['origLines']:
            print repr(l)
            outF.write(l)
        return 0
    #if info['name'][0].strip()[-5:]=='_type':
    #    outF.write("!> \struct %s\n"%(info['name'][0].strip()))
    lines=drop_trailing(info['brief'])
    if len(lines)>0:
        outF.write(pre)
        outF.write('\\brief ')
        i0=len(lines)
        for i,line in enumerate(lines):
            if not line.isspace():
                outF.write(line.lstrip())
                i0=i+1
                break
        for line in lines[i0:]:
            outF.write(pre)
            outF.write(line)
    for label in rest:
        lines=drop_trailing(info['special'][label])
        if lines:
            outF.write(pre)
            outF.write('\\par ')
            outF.write(label)
            outF.write('\n')
            for line in lines:
                outF.write(pre)
                if line.isspace():
                    outF.write('\\par\n')
                else:
                    outF.write(line)
        else:
            print "removing empty special label",label
    paramRe=re.compile(r" *- *(?P<param>[a-zA-Z_]+) *:? *(?P<rest>.+)")
    param2Re=re.compile(r" *(?P<param>[a-zA-Z_]+) *: *(?P<rest>.+)")
    lines=drop_trailing(info['params'])
    if len(lines)>0:
        outF.write(pre)
        outF.write("\param ")
        i0=len(lines)
        for i,line in enumerate(lines):
            if not line.isspace():
                m=paramRe.match(line)
                if not m: m=param2Re.match(line)
                if m:
                    outF.write(m.group('param'))
                    outF.write(' ')
                    outF.write(m.group('rest'))
                    outF.write('\n')
                else:
                    print "failed matching of param in ",repr(line)
                    outF.write(line.lstrip())
                i0=i+1
                break
        for line in lines[i0:]:
            outF.write(pre)
            m=paramRe.match(line)
            if not m : m=param2Re.match(line)
            if m:
                outF.write("\param %s %s\n"%(m.group('param'),m.group('rest')))
            else:
                outF.write(line)
    for label in fixedOrder:
        lines=info[label]
        lines=drop_trailing(lines)
        if len(lines)==0: continue
        if label in ['literature','history']:
            outF.write(pre)
            outF.write('\\par ')
            outF.write(label.title())
            outF.write('\n')
            for line in lines:
                outF.write(pre)
                if line.isspace():
                    outF.write('\\par\n')
                else:
                    outF.write(line)
        else:
            outF.write(pre)
            outF.write('\\')
            outF.write(label)
            if label in ['note','todo']:
                outF.write('\n')
                i0=0
            else:
                outF.write(' ')
                outF.write(lines[0].lstrip())
                i0=1
            for line in lines[i0:]:
                outF.write(pre)
                outF.write(line)
    basicFluffRe=re.compile(r" *(\** *\**|-*!?) *$")
    for l in info['trash']:
        if l and not basicFluffRe.match(l):
            print 'rmvd ',repr(l)
    return 1
            
if __name__=='__main__':
    inFile=file(sys.argv[1])
    outF=file(sys.argv[1][:-1]+'f90','w')
    roboDocRe=re.compile(" *!!")
    # fluffRe=re.compile(r" *![+=* ]+$")
    basicFluffRe=re.compile(r" *! *\** *\*{20,} *$")
    basicFluff2Re=re.compile(r" *!! *\** *\** *$")
    separatePreRe=re.compile(r" *(?:(?:end +)?(?:(?P<module>module)|program) +[a-zA-Z]|type +[a-zA-Z]|(?:recursive +|pure +|elemental +)*(?:subroutine|function) +[a-zA-Z])",re.IGNORECASE)
    #separatePostRe=re.compile(r" *end +(?:subroutine|function|type|module)",re.IGNORECASE)
    bodyRe=re.compile(r" *(?P<end>end +)(?:type *(?![ ()])|subroutine|function).")
    isspace=1
    didDoc=0
    inBody=0
    didModule=0
    while 1:
        (jline,comments,lines)=readFortranLine(inFile)
        if not lines: break
        m=bodyRe.match(lines[0])
        if m:
            if m.group('end'):
                inBody=None
            else:
                inBody=bodyRe.group()
        if roboDocRe.match(lines[0]):
            if basicFluff2Re.match(lines[0]):
                continue
            (jline,comments,lines,info)=parseRoboDoc(lines,inFile)
            output = StringIO.StringIO()
            didDocAtt=printRoboDoc(info,output,isspace)
            doc=output.getvalue()
            if not doc:
                didDocAtt=0
            if didDocAtt:
                if not isspace:
                    outF.write('\n')
                outF.write('! '+(77*'*')+'\n')
            outF.write(doc)
            if didDocAtt:
                didDoc=1
                if inBody:
                    print '* documentation in body',inBody
            isspace=1
        if separatePreRe.match(lines[0]):
            if separatePreRe.match(lines[0]).group('module'):
                if not didModule:
                    didModule=1
                    outF.write('! '+(77*'*')+'\n')
                elif not didDoc :
                    print "* module without documentation, documentation at wrong place?"
            else:
                outF.write('! '+(77*'*')+'\n')
            for l in lines:
                outF.write(l)
            isspace=0
            continue
        #if separatePostRe.match(lines[0]):
        #    for l in lines:
        #        outF.write(l)
        #    outF.write('! '+(77*'*')+'\n')
        #    isspace=1
        #    continue
        for l in lines:
            if not basicFluffRe.match(l):
                if not isspace or not l.isspace():
                    outF.write(l)
                isspace=l.isspace()
    outF.close()