set terminal png
set output "plot.png"
set grid xtics lc "grey"
set grid ytics lc "grey"
set xlabel "batch size"
set ylabel "GFLOP/s"
set autoscale fix
set mytics 2

plot "plot.txt" using 2:3 notitle smooth sbezier lc "grey", \
     "" using 2:3 notitle with points pt 7
