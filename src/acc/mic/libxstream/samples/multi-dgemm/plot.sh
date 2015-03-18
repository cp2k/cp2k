#!/bin/bash

FILE="benchmark.txt"

grep -A1 "Running " ${FILE} | tr "\n" " " | sed \
  -e "s/Running //g" \
  -e "s/ batche*s* of//g" \
  -e "s/items*... Performance: //g" \
  -e "s/ GFLOPS\/s//g" \
  -e "s/ -- /\n/g" \
> plot.txt

if [ "${OS}" != "Windows_NT" ] ; then
  gnuplot plot.plt
else
  export GDFONTPATH=/cygdrive/c/Windows/Fonts

  if [[ -f /cygdrive/c/Program\ Files/gnuplot/bin/wgnuplot ]] ; then
    /cygdrive/c/Program\ Files/gnuplot/bin/wgnuplot plot.plt
  else
    /cygdrive/c/Program\ Files\ \(x86\)/gnuplot/bin/wgnuplot plot.plt
  fi
fi
