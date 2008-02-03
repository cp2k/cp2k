#!/usr/bin/env python
import glob,commands,os,shutil,re,os.path

print """prepares for doxygen using the .F and .f90  from ../src and
transferring them to ../docF.

should be started from the tools or src directory (so that the relative
paths are correct).

It converts the extensions to .f90 and slightly changes them,
so that doxygen can be easily applied to them."""

doxyRe=re.compile(' *!>')
sepRe=re.compile(r" *! *\*{40,} *$")
includeRe=re.compile(r"^( *#?include *['\"].+\.f90)(['\"])",re.IGNORECASE)

# remove old files
shutil.rmtree('../docF')

os.mkdir('../docF')
os.mkdir('../docF/lib')
# copy new files
fromTo={'../src':'../docF/','../src/lib':'../docF/lib'}
for fromDir,toDir in fromTo.iteritems():
    for fromExt,toExt in {'.F':'.f90','.f90':'.f90.inc'}.iteritems():
        for f in glob.glob(os.path.join(fromDir,"*"+fromExt)):
            inF=file(f)
            outF=file(os.path.join(toDir,os.path.basename(f)[:-len(fromExt)]+toExt),'w')
            outF.write('# 1 "'+f+'"\n')
            wasDoxy=0
            while 1:
                line=inF.readline()
                if not line: break
                if not (wasDoxy and sepRe.match(line)):
                    outF.write(includeRe.sub(r"\1.inc\2",line))
                wasDoxy=doxyRe.match(line)
            inF.close()
            outF.close()

print 'done'