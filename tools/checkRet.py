import sys

if (len(sys.argv)!=3):
    print """ usage: checkRet.py object trace_file
    search in trace file for mismatched retain/release of object
    you should add to the trace instruction like these:
      PRINT *, "XXX created object ",object%id_nr
      PRINT *, "XXX retained object ",object%id_nr,object%ref_count
      PRINT *, "XXX released object ",object%id_nr,object%ref_count
    """
      
f=file(sys.argv[2])

sects={}

l1=" XXX created "+sys.argv[1]
l2=" XXX retained "+sys.argv[1]
l3=" XXX released "+sys.argv[1]
while 1:
    line=f.readline()
    if not line: break
    if line[:len(l1)]==l1:
        ll=line.split()
        sects[int(ll[3])]=1
    if line[:len(l2)]==l2:
        ll=line.split()
        sects[int(ll[3])]+=1
        if(sects[int(ll[3])]!=int(ll[4])):
            print "error 1 ",ll[3]
    if line[:len(l3)]==l3:
        ll=line.split()
        sects[int(ll[3])]-=1
        if(sects[int(ll[3])]!=int(ll[4])):
            print "error 2 ",ll[3]
        if sects[int(ll[3])]==0:
            del sects[int(ll[3])]
print "REST:"
for n in sects.keys():
    print  sys.argv[1],n," ref_count ",sects[n]

        
