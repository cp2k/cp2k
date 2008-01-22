import sys,re
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
        'warning':[]}
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
    roboDocRe=re.compile(" *!!(.*)")
    headerRe=re.compile(r" *!!\*+[cdfhmstuv]?\** *(?P<name>[a-zA-Z_0-9/]+) *(?:\[(?P<version>[0-9.a-zA-Z_ ]+)\])?")
    labelRe=re.compile(" +(?P<label>NAME|COPYRIGHT|USAGE|FUNCTION|DESCRIPTION|PURPOSE|AUTHOR|CREATION DATE|MODIFICATION HISTORY|HISTORY|INPUTS|ARGUMENTS|OPTIONS|PARAMETERS|SWITCHES|OUTPUT|SYNOPSIS|SIDE EFFECTS|RESULT|RETURN VALUE|EXAMPLE|NOTES?|WARNINGS?|ERROR|DIAGNOSTICS|BUGS|TODO|IDEAS|PORTABILITY|SEE ALSO|METHODS|ATTRIBUTES|SOURCE|LITERATURE|TAGS|USED BY|[A-Z]+ ?[A-Z]* ?[A-Z]*) *$")
    fluffRe=re.compile(r" *([-+*#!= ]{3,})? *$")
    m=headerRe.match(lines[0])
    if m:
        if m.group('version'):
            info['version'].append(m.group('version'))
        if m.group('name'):
            info['pname'].append(m.group('name'))
        lines=lines[1:]
    while 1:
        for line in lines:
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
        (jline,comments,lines)=readFortranLine(inFile)
        if not lines:
            break
        if jline!=None and jline!="" and not jline.isspace():
            break
        if not roboDocRe.match(lines[0]):
            break
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

def printRoboDoc(info,outF=sys.stdout):
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
    if info['brief']==None:
        print 'brief0',info
    lines=drop_trailing(info['brief'])
    if lines==None:
        print 'lines0',info
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
            
if __name__=='__main__':
    inFile=file(sys.argv[1])
    outF=file(sys.argv[1]+'.doxy','w')
    roboDocRe=re.compile(" *!!")
    while 1:
        (jline,comments,lines)=readFortranLine(inFile)
        if not lines: break
        if roboDocRe.match(lines[0]):
            (jline,comments,lines,info)=parseRoboDoc(lines,inFile)
            printRoboDoc(info,outF)
        for l in lines:
            outF.write(l)
    outF.close()