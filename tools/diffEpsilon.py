#! /usr/bin/env python

import sys

def compareNr(n1,n2):
  return abs(n1-n2)/max(abs(n1),abs(n2))

def diffEpsilon(str1, str2,incomparable_val=1):
    """retuns the difference between two strings, parsing numbers and confronting them."""
    import re
    nrRe=re.compile("[-+]?[0-9]*\\.?[0-9]+([eEdD][-+]?[0-9]+)?$")
    tokens1=str1.split()
    tokens2=str2.split()
    distance=0.0
    if len(tokens1)!=len(tokens2):
        return incomparable_val
    i=0
    for t1 in tokens1:
        t2=tokens2[i]
        i=i+1
        if (t1!=t2):
            if nrRe.match(t1) and nrRe.match(t2):
                (f1,f2)=(float(t1),float(t2))
                distance=max(distance, compareNr(f1,f2))
            else:
                return incomparable_val
    return distance

def getCoreLine(file,lineNr):
  import re
  bannerRe=re.compile("([ *]+PROGRAM | CP2K\\| Input file name| CP2K\\| Program compiled)")
  timingRe=re.compile(" *- +T I M I N G +- *$")
  scfLineRe=re.compile(" *([0-9]+) +([a-zA-Z]+) +([-+]?[0-9]*\\.?[0-9]+[EedD]?[-+]?[0-9]*) +([-+]?[0-9]*\\.?[0-9]+[EedD]?[-+]?[0-9]*) +([-+]?[0-9]*\\.?[0-9]+[EedD]?[-+]?[0-9]*) +([-+]?[0-9]*\\.?[0-9]+[EedD]?[+-]?[0-9]*) *$")
  line=file.readline()
  while bannerRe.match(line):
    lineNr=lineNr+1
    line=file.readline()
  if scfLineRe.match(line):
    # removing timing and convergence info (too system dependent)
    # print "found scfLine=",`line`
    groups=scfLineRe.match(line).groups()
    line=groups[0]+" "+groups[1]+" "+groups[2]+" x x "+groups[5]
  if timingRe.match(line):
    while not line=="":
      lineNr=lineNr+1
      line=file.readline()
  if line!="":
    lineNr=lineNr+1
  return (line,lineNr)

def compareCp2kOutput(file1,file2, epsilon=0.01,logFile=sys.stdout):
    lineNr1=lineNr2=0
    distance=0.0
    try:
      while 1:
        (line1,lineNr1)= getCoreLine(file1,lineNr1)
        (line2,lineNr2)= getCoreLine(file2,lineNr2)
#        print "line1=",`line1`,"line2=",`line2`
        if line1=="" and line2=="": break
        lineDiff=diffEpsilon(line1,line2)
        if lineDiff>epsilon:
          logFile.write("diff line "+`lineNr1`+" vs line "+`lineNr2`+" = "+`lineDiff`
                        +":\n"+`line1`+"\n"+`line2`+"\n")
        distance=max(distance, lineDiff)
    except:
      import sys, traceback
      logFile.write('-'*60+"\n")
      traceback.print_exc(file=logFile)
      logFile.write(file1.name+":"+`lineNr1`+", "+file2.name+":"+`lineNr2`+"\n")
      logFile.write('-'*60+"\n")
    return distance

if __name__ == '__main__':
    if len(sys.argv)<3 or len(sys.argv)>4:
      print "usage:", sys.argv[0]," file2 out_dir file2 [epsilon]"
    else:
        file1=open(sys.argv[1],"r")
        file2=open(sys.argv[2],"r")
        epsilon=0.01
        if len(sys.argv)==4: epsilon=float(sys.argv[3])
        print "distance=",compareCp2kOutput(file1,file2,epsilon,sys.stdout)
        
