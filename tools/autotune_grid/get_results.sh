#!/usr/local/bin/bash
#
#   get_results
#   ==============================================
#       Extract the best combinations from all the
#          directories.
#
#    @author: Ruyman Reyes (rreyesc@epcc.ed.ac.uk)
#
source config.in

# Number of combinations with the existing options
Ncomb=$((2 ** $Nopt))

# Error tolerance
err_tol="1.0E-10"

# For both collocate (0) and integrate (1)
for icollo in 0 1; do
  echo "$icollo"
  echo "l - Time - Best Opt"
  echo "l Time Opt" > best_timings
  for l in $(seq 0 $lmax); do
    bsf=999999.0
    # First optimisation starts at ONE (all_options[0] is not initialised and generate does not produces valid code)
    for iopt in $(seq 1 $Ncomb); do
      # Extract the error
      tmp=$(cat out_${l}_${iopt}/out_test_${l}_${iopt}_F_${icollo} | tr -s ' ' | grep "largest error" | cut -f 5 -d ' ')
      test "$tmp" || continue
      error=$(python -c "print $err_tol>$tmp")
      # Skip if error is too big
      if [[ $error == "False" ]]; then
        #   echo "($error) Error $tmp was too big for $iopt"
        continue
      fi
      # echo "($error) Error was $tmp , accepted $iopt"
      # Extract the time of the combination
      tmp=$(cat out_${l}_${iopt}/out_test_${l}_${iopt}_T_${icollo} | tr -s ' ' | grep "best time" | cut -f 5 -d ' ')

      test "$tmp" || continue
      dif=$(python -c "print $bsf>$tmp")

      if [[ $dif == "True" ]]; then
        # echo " New bsf is $tmp from $bsf "
        bsf=$tmp
        bsf_iopt=$iopt
        #else
        # echo " $bsf is lower than $tmp ($dif, print $bsf>$tmp)"
      fi
    done # for iopt

    bsf_global[$l]=$bsf
    bsf_iopt_global[$l]=$bsf_iopt

    echo $l -- ${bsf_global[$l]} -- ${bsf_iopt_global[$l]}
    echo $l ${bsf_global[$l]} ${bsf_iopt_global[$l]} >> best_timings

  done # for l

  if [[ $icollo -eq 0 ]]; then
    bsf_global_0=("${bsf_global[@]}")
    bsf_iopt_global_0=("${bsf_iopt_global[@]}")
  else
    bsf_global_1=("${bsf_global[@]}")
    bsf_iopt_global_1=("${bsf_iopt_global[@]}")
  fi

done # for icollo

# Generate the input file for generate.x with the optimal combinations for each problem size
echo " Generating optimal combination "
(
  echo $lmax
  for l in $(seq 0 $lmax); do
    echo ${bsf_iopt_global_0[$l]}
  done
  echo $lmax
  for l in $(seq 0 $lmax); do
    echo ${bsf_iopt_global_1[$l]}
  done
) > generate_best
