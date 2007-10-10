#!/usr/bin/env python
# copies inFile to outFile converting high precision numbers to low precision ones
# useful with .restart files because some compilers cannot read back with read * what
# they wrote with very high precision (16 decimals)
# cp2k 2007 (fawzi)
import re,sys,os.path

defaultHighPrecision=re.compile('[-+]?[0-9]*\.[0-9]{10,20}[EeDd][+-][0-9]+')

def transcodeToLowerPecision(f_in,f_out,highPrecision=defaultHighPrecision,outPrec='%10.10g'):
    'writes file f to f_out converting high precision numbers to low precision ones'
    while 1:
        line=f_in.readline()
        if not line: break
        ipos=0
        for m in highPrecision.finditer(line):
            f_out.write(line[ipos:m.start()])
            f_out.write(outPrec % (float(m.group())))
            ipos=m.end()
        f_out.write(line[ipos:])


if __name__=='__main__':
    if (len(sys.argv)!=3):
        usage='usage: '+os.path.basename(sys.argv[0])+""" inFile outFile
        copies inFile to outFile converting high precision numbers to low precision ones\
        useful with .restart files of compilers that cannot read back what they wrote..."""
        print usage
        sys.exit(1)
    
    f_in=file(sys.argv[1])
    f_out=file(sys.argv[2],'w')
    transcodeToLowerPecision(f_in,f_out)
    f_in.close()
    f_out.close()
    
