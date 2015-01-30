#!/bin/sh
# ----------------------------------------------------
# CP2K (remote) regtester 
# Developed by Matthias Krack, Marko Misic
# 
# Update paths to correctly point to www directory, cp2k directory, and regtest directory 
#
# ----------------------------------------------------

# Owner mail address
owner_mail="marko.misic@etf.bg.ac.rs"

# Directory structure
# - regtestdir  - where all the scripts are
# - cp2kdir     - where cp2k and do_regtest script reside
# - wwwtestdir  - all output files are generated here and then copied to wwwdir

cp2kdir=/disk/cp2k-www0/home/cp2ktest/cp2k/
regtestdir=/disk/cp2k-www0/home/cp2ktest/regtest/revisited/
wwwtestdir=/disk/cp2k-www0/home/cp2ktest/regtest/revisited/wwwtest

# Command line configuration for do_regtest script
# .conf file defines environment variables for do_regtest script
conf_file="${cp2kdir}/conf/regtest.popt.conf"
cmd_line="-config ${conf_file} -nosvn -noemptycheck"

# Must be specified to produce memory leaks report
# Only if g95 is the compiler!
memcheck=false

# Execution type: local or remote
exec_type="local"

# This function should be rewritten to support different batchin systems

function do_regression_test_local() {

 ${cp2kdir}/do_regtest ${cmd_line} >>${reglog} 

# Testing line
# Enable this one, and disable the one above if you want to test the regtester
# echo "GREPME 99 0 2335 0 2434 0" >>${reglog}

}

function do_regression_test_remote() {

# Submit the job through the batching system
  job_id=`bsub < regtest.sopt.sh | sed 's/[^0-9]*<\([0-9]*\)>[^0-9]*/\1/'`

# Wait some time for the results
  while true
  do
    status=`bjobs ${job_id} | grep -E 'DONE|EXIT'`
    if [[ ${status} ]]
    then
      break
    else
      sleep 600
    fi
  done
 
# Do not forget to copy output of do_regtest from batch job 
  cat "=========== Standard output ===========" >>${reglog}
  cat regtest_sopt_bsub.out | head -n -7 >>${reglog}
  cat "=========== Standard error ============" >>${reglog}
  cat regtest_sopt_bsub.err >>${reglog}
  cat "=========== End of output =============" >>${reglog}

}

function do_regression_test {
  if [[ ${exec_type} == "remote" ]] 
  then
    do_regression_test_remote
  else
    do_regression_test_local  
  fi
}

# ----------------------------------------------------
# From here, no need to change anything
# ----------------------------------------------------

# The script takes two or three arguments (the third is config file) 
# Name of the tarball to pack the results
# To force or not regtest

# Check the number of command line arguments
if [[ $# -ne 2 && $# -ne 3 ]]
then
  echo "Command line: cp2k_regtester [config_file] tarball_name (YES|NO)"
  exit 500
fi

# If configuration file is provided, load it to override default settings
if [[ $# == 3 ]]
then
  source $1
  regtest_tarball=$2
  forceregtest=$3
else
  regtest_tarball=$1
  forceregtest=$2
fi

# Script that should handle any PATH or LD_LIBRARY_PATH issues
source ${regtestdir}/setup_environment.sh

datetime=$(date --utc +%F_%H-%M-%S)
date=${datetime%_*}
time=${datetime#*_}

# Try to find out the configuration
dir_triplet=`cat ${conf_file} | grep "dir_triplet" | sed 's/.*=\(.*\)/\1/'`
cp2k_version=`cat ${conf_file} | grep "cp2k_version" | sed 's/.*=\(.*\)/\1/'`
dir_last=${cp2kdir}/LAST-${dir_triplet}-${cp2k_version}

# Memlogs are not used unless g95 is the compiler
memlog=${wwwtestdir}/${datetime}_memory_leaks_log.txt
reglog=${wwwtestdir}/${datetime}_regtester_log.txt

# Create those files
touch ${memlog}
touch ${reglog}

# Fill in memlog with some text, if leaks checking is not supported
echo "================== CP2K automated regression tester ==================" >>${memlog}
if [[ ${memcheck} == false ]] 
then
  echo "Memory leaks checking is not supported for this compiler/configuration." >> ${memlog}
  echo "If you need memory checking compiler, try g95." >> ${memlog}
fi

echo "================== CP2K automated regression tester ==================" >>${reglog}
echo "Date: ${date}" >>${reglog}
echo "Time: ${time//-/:}" >>${reglog}
echo "Configuration: ${dir_triplet}-${cp2k_version}" >>${reglog}

# Do SVN update
${regtestdir}/cp2k_svn_update.sh ${conf_file} ${cp2kdir} ${wwwtestdir} >>${reglog}
exitcode=$?
if [[ ${exitcode} -eq 1 ]]
then
  echo "ERROR: Problem with SVN update!" >>${reglog}
  svn_update=false
  mail -s "ERROR: Problem with SVN update!" ${owner_mail} <<EOF
Check the status of the automatic regression tester!
The regression tester tries to continue, since this might only be a temporary problem.
EOF
elif [[ ${exitcode} -eq 2 ]]
then
  svn_update=false
else
  svn_update=true
fi    

if [[ ${forceregtest} == "NO" && ${svn_update} == false ]]
then
  rm -f ${reglog}
  rm -f ${memlog}
  exit 1
fi

# Regression testing
cd ${cp2kdir}

echo "Starting regression testing at ${datetime}!"

do_regression_test

echo "Ended regression testing from ${datetime}!"

# Copy ChangeLog, arch file, update archive with reference results and do some cleanup
  tar -cjf ${wwwtestdir}/trunk_LAST-${dir_triplet}-${cp2k_version}.tar.bz2 ${dir_last}
  
  cp ${dir_last}/ChangeLog ${wwwtestdir}/ChangeLog

  cp ${cp2kdir}/cp2k/arch/${dir_triplet}.${cp2k_version} ${wwwtestdir}

  if [[ ${memcheck} == true ]] 
  then
    mv memory_leaks.txt ${memlog}
  fi
 
  # Cleanup - only last 10 TEST directories are kept
  removedir=`ls | grep "TEST-${dir_triplet}-${cp2k_version}" | head -n -10`
  if [[ -n ${removedir} ]]
  then
    for r in ${removedir} 
    do
      rm -rf $r
    done
  fi

# Produce big tar file with the results
  cd ${wwwtestdir}
  tar -czf ${regtest_tarball} ${datetime}_memory_leaks_log.txt ${datetime}_regtester_log.txt ChangeLog ChangeLog.diff ${dir_triplet}.${cp2k_version} trunk_LAST-${dir_triplet}-${cp2k_version}.tar.bz2

exit 0
