#!/usr/bin/python
# -*- coding: utf-8 -*-
#
# author: Ole Schuett

import sys
from os import path

#===============================================================================
def main():
    if(len(sys.argv) != 3):
        print("Usage: bcast_tmpl_coverage.py <lcov_in.info> <lcov_out.info>")
        print("Broadcasts coverage stats across all template instances.")
        sys.exit(1)

    fn_in, fn_out = sys.argv[1:]

    print("Reading "+fn_in)
    content = open(fn_in).read()

    # pass 1: find all records
    records = dict()
    curr_record = []
    curr_fn = None
    for line in content.split("\n"):
        curr_record.append(line)
        if(line.startswith("SF:")):
            curr_fn = line[3:]
            #print curr_fn
        elif(line.startswith("end_of_record")):
            if(not records.has_key(curr_fn)):
                records[curr_fn] = list()
            records[curr_fn] += curr_record
            curr_record = []
            curr_fn = None

    # pass 3: sort records
    content_sorted = []
    for k in sorted(records.keys()):
        content_sorted += records[k]

    # pass 4: join records into broadcast groups
    groups = []
    prev_fn = prev_line_nums = None
    curr_fn = curr_line_nums = None
    curr_group = []
    for line in content_sorted:
        curr_group.append(line)

        if(line.startswith("SF:")):
            prev_fn = curr_fn
            prev_line_nums = curr_line_nums
            curr_fn = line[3:]
            curr_line_nums = []
            #print curr_fn
        elif(line[:3] in ("FN:", "DA:")):
            lineno = int(line[3:].split(",")[0])
            curr_line_nums.append(lineno)
        elif(line.startswith("end_of_record")):
            if(similar(curr_fn, prev_fn) and curr_line_nums==prev_line_nums):
                groups[-1] += curr_group
                print "Broadcasting: "+path.basename(prev_fn), "<->",path.basename(curr_fn)
            else:
                groups.append(curr_group)
            curr_group = [] # start a new group

    output = []
    for group in groups:
        FNs = dict()
        FNDAs = dict()
        DAs = dict()
        # pass 5: broadcast coverage stats within groups
        for line in group:
            if(line.startswith("FN:")):
                lno, name = line[3:].split(",")
                FNs[name] = lno
            elif(line.startswith("FNDA:")):
                count, name = line[5:].split(",")
                uid = FNs[name] if(name in FNs) else name
                FNDAs[uid] = max(int(count), FNDAs.get(uid, 0))
            elif(line.startswith("DA:")):
                lno, count = line[3:].split(",")
                DAs[lno] = max(int(count), DAs.get(lno, 0))

        # pass 6: write new records
        for line in group:
            if(line.startswith("FNDA:")):
                count, name = line[5:].split(",")
                uid = FNs[name] if(name in FNs) else name
                output.append( "FNDA:%d,%s"%(FNDAs[uid], name) )
            elif(line.startswith("DA:")):
                lno, count = line[3:].split(",")
                output.append( "DA:%s,%d"%(lno,DAs[lno]) )
            elif(line.startswith("LH:")):
                count = len([v for v in DAs.values() if v>0])
                output.append( "LH:%d"%count )
            elif(line.startswith("FNH:")):
                count = len([v for v in FNDAs.values() if v>0])
                output.append( "FNH:%d"%count )
            else:
                output.append(line)


    print("Writting "+fn_out)
    f = open(fn_out, "w")
    f.write("\n".join(output))
    f.close()

#===============================================================================
def similar(a, b):
    """ Checks if two filenames are very similar """
    if(a==None or b==None): return(False)
    if(path.dirname(a) != path.dirname(b)): return(False)
    a_base, a_ext = path.basename(a).rsplit(".", 1)
    b_base, b_ext = path.basename(b).rsplit(".", 1)
    if(a_ext != b_ext): return(False)
    a_parts = a_base.split("_")
    b_parts = b_base.split("_")
    if(len(a_parts) != len(b_parts)): return(False)
    diffs = len(a_parts) - sum([i==j for i,j in zip(a_parts, b_parts)])
    return(diffs ==  1) # at most one field may differ


main()
#EOF
