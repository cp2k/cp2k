#!/bin/bash -e

arr=("before" "after")


for i in 0 1
do
for ncpus in 1 2
do
for t in ${arr[@]}
do
echo "edge index ${t}:           ${i}"
grep "edge index ${t}:           ${i}" out${ncpus}.out | awk '{print $(NF-1), $NF}' | sort -k1 -n  >  edge_index_${t}_${i}_${ncpus}cpu.dat
grep "edge index ${t}:" out${ncpus}.out | awk '{print $(NF-1), $NF}' | sort -k1 -n >   edge_index_full_${t}_${ncpus}cpu.dat
done
done
done
