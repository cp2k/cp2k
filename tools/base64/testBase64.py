#!/usr/bin/env python

import sys,commands,os,os.path,re

def compareFiles(path1,path2,logFile,small=0.5**1021):
    """compares two files and prints out a summery on how different thy are"""
    if (not os.path.exists(path1) or not os.path.exists(path2)):
        logFile.write("*** ERROR comparing "+repr(path1)+" with "+repr(path2)+
                      " ***\n")
        if (not os.path.exists(path1)):
            logFile.write(repr(path1)+" does not exist.\n")
        if (not os.path.exists(path2)):
            logFile.write(repr(path2)+" does not exist.\n")
        return float('inf')
#    try:
    f1=file(path1)
    f2=file(path2)
    diff=comparenum(f1,f2,small)
    if diff!=None:
        if diff['diff']>1.e-6:
            logFile.write('O ')
        elif diff['diff']>1.e-15:
            logFile.write('* ')
        else:
            logFile.write('- ')
        logFile.write("err:%8g %20s:%-5i vs %20s:%-5i chunk:%5i %s!=%s\n"%
                      (diff['diff'],repr(path1),diff['line1'],repr(path2),
                       diff['line2'],diff['chunk'],str(diff['e1']),
                       str(diff['e2'])))
        res=diff['diff']
    else:
        logFile.write("- err:%8g %25s vs %25s -\n" %
                      (0.,repr(path1),repr(path2)))
        res=0.
#    except:
#        import traceback
#        print "** INTERNAL ERROR trying to compare "+repr(path1)+" with "+repr(path2)+" **"
#        traceback.print_stack(file=logFile)
#        res=float('inf')
    return res

def comparenum(file1,file2,small=0.5**(1021)):
    """ compares two files of real numbers """
    longNRe=re.compile(r"([+-]?[0-9]*\.?[0-9]*)[dD]?([+-][0-9]+)$")
    l1=[]
    l2=[]
    diffR=0.
    diffRMax=None
    i=0
    il1=0
    il2=0
    while 1:
        if not l1:
            while 1:
                linea=file1.readline()
                if not linea:break
                il1+=1
                l1+=linea.split()
                if l1: break
        if not l2:
            while 1:
                linea=file2.readline()
                if not linea:break
                il2+=1
                l2+=linea.split()
                if l2: break
        if not l1:
            if l2:
                diffRMax={'line1':il1,'line2':il2,'chunk':i,'e1':None,
                          'e2':l2[0],'diff':float('inf')}
            break
        if not l2:
            diffRMax={'line1':il1,'line2':il2,'chunk':i,'e1':l1[0],
                      'e2':None,'diff':float('inf')}
            break
        i+=1
        e1=l1[0]
        l1=l1[1:]
        e2=l2[0]
        l2=l2[1:]
        m=longNRe.match(e1)
        if m:
            e1=m.groups()[0]+"e"+m.groups()[1]
        m=longNRe.match(e2)
        if m:
            e2=m.groups()[0]+"e"+m.groups()[1]
        comp=(e1==e2)
        if not comp:
            try:
                r1=float(e1)
                r2=float(e2)
                if r1==0. and r2==0.:
                    myDiff=0.
                elif r1!=r1 or r2!=r2:
                    if r1!=r1 and r2!=r2:
                        myDiff=0.
                    else:
                        myDiff=float('inf')
                else:
                    myDiff=abs(r1-r2)/max(small,(abs(r1)+abs(r2)))
                if diffR<myDiff:
                    diffRMax={'line1':il1,'line2':il2,'chunk':i,'e1':r1,
                              'e2':r2,'diff':myDiff}
                    diffR=myDiff
                comp=1
            except ValueError:
                None
        if not comp:
            diffRMax={'line1':il1,'line2':il2,'chunk':i,'e1':e1,'e2':e2,
                      'diff':float('inf')}
            break
    return diffRMax

def testB64(filename,logFile,small=0.5**1021):
    if not os.path.exists(filename):
        logFile.write("ERROR testfile "+repr(filename)+" does not exist\n")
        return float('inf')
    diff=0
    for suffix in ['N_1','N_2','_1','_2']:
        if os.path.exists(filename+suffix):
            os.remove(filename+suffix)
    r=commands.getoutput(conv_exe+' --from-base64 '+filename+' '+filename+'N_1')
    if not os.path.exists(filename+"N"):
        logFile.write("WARNING reference result "+repr(filename+"N")+" does not exist\n")
    else:
        diff=max(diff,compareFiles(filename+"N",filename+"N_1",logFile=logFile,small=small))
    r2=commands.getoutput(conv_exe+' --to-base64 '+filename+'N_1 '+filename+'_1')
    r3=commands.getoutput(conv_exe+' --from-base64 '+filename+'_1 '+filename+'N_2')
    diff2=compareFiles(filename+"N_1",filename+"N_2",logFile=logFile,small=small)
    diff=max(diff,diff2)
    if (diff2!=0.):
        logFile.write("WARNING some number are unexpectedly different\n")
        logFile.write(" the compiler might have problems generating exactly the same numbers\n")
        logFile.write(" when reading back what it wrote... ignoring result of next test\n")
    r4=commands.getoutput(conv_exe+' --to-base64 '+filename+'N_2 '+filename+'_2')
    diff3=compareFiles(filename+"_1",filename+"_2",logFile=logFile,small=small)
    if diff2==0.:
        diff=max(diff,diff3)
    return diff

if __name__=="__main__":
    doc=os.path.basename(sys.argv[0])+""" base64Converter.x
    This script tests the conversion of reals to base64.
    You have to give the executable obtained by compiling convertBase64.F90
    as argument to this script"""

    numDataRe=re.compile(" *digits *(?P<digits>[0-9]+) *maxexp *(?P<maxexp>[0-9]+) *minexp *(?P<minexp>-[0-9]+) *radix *(?P<radix>[0-9]+)")

    if len(sys.argv)!=2:
        print doc
        sys.exit(1)
    conv_exe=sys.argv[1]
    logFile=sys.stdout
    print 20*"="+ " converter self tests results "+20*"="
    r=commands.getoutput(conv_exe)
    print r
    d=filter(lambda x:numDataRe.match(x),r.splitlines())
    if d:
        m=numDataRe.match(d[0])
        # failures in the denormalized range deemed not important
        small=float(m.group("radix"))**(int(m.group("minexp"))+int(m.group("digits"))-1)
    else:
        small=0.5**(-1021)
    print 60*"="

    for fname in ['ref_32','ref_64','normalNrs']:
        print "** convert",fname,"**"
        diff=testB64(fname,logFile,small)

    print "\n======== ERROR:",diff,"=========\n"
    


