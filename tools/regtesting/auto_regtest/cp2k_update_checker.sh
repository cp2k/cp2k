#!/bin/sh
# ----------------------------------------------------
# CP2K automated regtester 
# Developed by Matthias Krack (2012), Marko Misic (2013)
# 
# Update paths to correctly point to www directory, cp2k directory, and regtest directory 
# Update e-mail account of the owner of automated regtester
#
# ----------------------------------------------------

# Owner mail address
owner_mail="marko.misic@etf.bg.ac.rs"
# If this switch is enabled, emails are sent to developers
# Owner of the automated regtester is always mailed if there are any problems
email_developers=false

# Directory structure
# - regtestdir  - where all the scripts are
# - wwwdir      - actual www dir 
# - wwwtestdir  - all output files are generated here and then copied to wwwdir

regtestdir=/disk/cp2k-www0/home/cp2ktest/regtest/remote-testing/indy-sopt/
wwwdir=/disk/cp2k-www0/websites/cp2k-www.epcc.ed.ac.uk/indy/sopt/
wwwtestdir=/disk/cp2k-www0/home/cp2ktest/regtest/remote-testing/indy-sopt/www-indy-sopt/

# Regtest filenames 
# - html - file name of the main page
# - plot - file name of the plot
# - par  - grace input file with plot configuration 
# - dat  - data for the plot
regtesthtml="regtest.html"
regtestplot="regplot.png"
regtestpar="regplot.par"
regtestdat="regplot.dat"

# Must be specified to produce memory leaks report
# Only if g95 is the compiler!
memcheck=false

# Maximum number of errors to be displayed in the plot
(( maxvalue = 40 ))

# Sleep time
(( sleep_time = 600 ))

# Test configuration
regtest_tarball="regtest-indy-sopt.tar.gz"

# Test command should call cp2k_regtester script either localy or remotely
#test_command="${regtestdir}/cp2k_regtester.sh ${regtest_tarball}"
test_command="ssh cp2ktest@indy0 \"/home/w02/cp2ktest/regtest/cp2k_regtester.sh ${regtest_tarball}\""

# Copy command should copy the tarball produced by cp2k_regtester script
# For local execution, keep it empty, and make sure wwwtestdir is the same for both scripts
#copy_command=""
copy_command="scp cp2ktest@indy0:/home/w02/cp2ktest/regtest/www-indy-sopt/${regtest_tarball} ${wwwtestdir}"

# Tests are run depending on several factors:
# - if a change is commited to SVN
# - FORCE_AUTO_REGTEST file forces regtest to proceed
# - SUSPEND_AUTO_REGTEST file suspends regtesting

# ----------------------------------------------------
# From here, no need to change anything
# ----------------------------------------------------

# If configuration file is provided, load it to override default settings
if [[ $# == 1 ]]
then
  source $1
fi

# Spin and test
while true; do

  send_email=false

# Enter wwwtestdir  to update links and generate new plot     
    cd ${wwwtestdir}

# Check if test is suspended or forced
  forceregtest=false
  if [[ -f SUSPEND_AUTO_REGTEST ]]; then
    echo "Testing has been manually suspended"
    sleep ${sleep_time}
    continue
  elif [[ -f FORCE_AUTO_REGTEST ]]; then
    forceregtest=true    
    rm FORCE_AUTO_REGTEST
  fi
   
# Execute remote commands
  if [[ ${forceregtest} == false ]]
  then
    eval "${test_command} NO"
    exitcode=$?
  else
    echo "Testing has been forced, even if SVN has not changed"
    eval "${test_command} YES"
    exitcode=$?
    forceregtest=false
  fi

  if [[ ${exitcode} -eq 1 ]]
  then
    echo "CP2K has not changed - nothing to test!"
    sleep ${sleep_time}
    continue
  elif [[ ${exitcode} -eq 255 ]]
  then
    echo "Execution of the (remote) command failed!"
  fi

# Copy the results   
  eval ${copy_command}
  exitcode=$?

  if [[ $? -eq 255 ]]
  then
    echo "Execution of the (remote) command failed!"
  fi
  
  tar -xzmf ${wwwtestdir}/${regtest_tarball}
  
  rm -f ${wwwtestdir}/${regtest_tarball}

  reglog=`ls ${wwwtestdir} | grep regtester_log.txt | tail -n 1`
  memlog=`ls ${wwwtestdir} | grep memory_leaks_log.txt | tail -n 1`

  date=`cat ${wwwtestdir}/${reglog} | grep "Date: " | sed 's/Date: \([0-9\-]*\)/\1/'`
  time=`cat ${wwwtestdir}/${reglog} | grep "Time: " | sed 's/Time: \([0-9\:]*\)/\1/'`
  config=`cat ${wwwtestdir}/${reglog} | grep "Configuration: " | sed 's/Configuration: \(.*\)/\1/'`

  dir_triplet=`cat ${wwwtestdir}/${reglog} | grep "Configuration: " | sed 's/Configuration: \(.*\)-.*/\1/'`
  cp2k_version=`cat ${wwwtestdir}/${reglog} | grep "Configuration: " | sed 's/Configuration: .*-\(.*\)/\1/'`
  
  new_results=false
 
  if [[ $(grep "error happened : make realclean VERSION=" ${reglog} ) ]]
  then 
    new_results=true
    send_email=true
  elif [[ $(grep "error happened : make VERSION=" ${reglog} ) ]]
  then 
    new_results=true
    send_email=true
  else
    new_results=true
  fi
  
  if [[ ${new_results} == true ]]
  then 

# First we update links of the last 10 regtests
    for ((i=0;i<10;i++)); do
      ln -sf $(ls -1 *_regtester_log.txt    | tail -n 10 | head -n $((i+1)) | tail -1) regtest-$((10-i-1))
      ln -sf $(ls -1 *_memory_leaks_log.txt | tail -n 10 | head -n $((i+1)) | tail -1) memory_leaks-$((10-i-1))
    done

# Update arch link
    ln -sf ${dir_triplet}.${cp2k_version} regtest-arch

# Update regtest-results.html
    cat ${regtestdir}/regtest-template.html | sed "s/{{TAG_CP2K_CONFIG}}/${config}/" > ${wwwtestdir}/regtest.html
    sed -i "s/{{TAG_CP2K_REF_OUT}}/trunk_LAST-${dir_triplet}-${cp2k_version}.tar.bz2/" ${wwwtestdir}/regtest.html
 
# Check revision number
    if [[ $(grep -m 1 "At revision" ${reglog}) ]]; then
      revision=$(grep -m 1 "At revision" ${reglog} | cut -c13-17)
      revision_update=true
    elif [[ $(grep -m 1 "Updated to revision" ${reglog}) ]]; then
      revision=$(grep -m 1 "Updated to revision" ${reglog} | cut -c21-25)
      revision_update=true
    else
      revision_update=false
    fi
  
# Update plot  
    if [[ $(grep "GREPME" ${reglog}) ]]; then
      set -- $(grep "GREPME" ${reglog})
      runtime_errors=$2
      wrong_results=$3
      correct_tests=$4
      new_inputs=$5
      number_of_test_inputs=$6
      memory_leaks=$7
      if [[ ${memory_leaks} == "X" ]]
      then
        (( memory_leaks = 0 ))
      fi
      sed -i "s/\(Results of the last run from\) [0-9-]* [0-9:]*/\1 ${date} ${time//-/:}/" ${regtestpar}
      sed -i "s/\(Number of test inputs:\) [0-9]*/\1 ${number_of_test_inputs}/" ${regtestpar}
      sed -i "s/\(Runtime errors:\) [0-9]*/\1 ${runtime_errors}/" ${regtestpar}
      sed -i "s/\(Wrong results:\) [0-9]*/\1 ${wrong_results}/" ${regtestpar}
      sed -i "s/\(New test inputs:\) [0-9]*/\1 ${new_inputs}/" ${regtestpar}
      sed -i "s/\(Memory leaks:\) [0-9]*/\1 ${memory_leaks}/" ${regtestpar}
      (( runtime_errors > maxvalue )) && (( runtime_errors = maxvalue ))
      (( wrong_results > maxvalue )) && (( wrong_results = maxvalue ))
      (( new_inputs > maxvalue )) && (( new_inputs = maxvalue ))
      (( memory_leaks > maxvalue )) && (( memory_leaks = maxvalue ))
      if [[ -z $(grep "${revision} ${runtime_errors} ${wrong_results} ${new_inputs} ${memory_leaks}" ${regtestdat}) ]]; then
        echo "${revision} ${runtime_errors} ${wrong_results} ${new_inputs} ${memory_leaks}" >>${regtestdat}
      fi
      set -- $(head -2 ${regtestdat} | tail -1)
      (( xmin = ($1 - 1)/10*10 ))
      set -- $(tail -1 ${regtestdat})
      (( xmax = ($1 + 10)/10*10 ))
      (( xmax - 50 > xmin )) && (( xmin = xmax - 50 ))
      (( ymax = maxvalue + 20 ))
      sed -i "s/\(  world\) [0-9]*, [0-9]*, [0-9]*, [0-9]*/\1 ${xmin}, 0, ${xmax}, ${ymax}/" ${regtestpar}
      gracebat -autoscale none -param ${regtestpar} -hdevice PNG -printfile "${regtestplot}" -nxy ${regtestdat}

      (( runtime_errors > 0 )) && send_email=true
      (( wrong_results > 0 )) && send_email=true
      (( memory_leaks > 0 )) && send_email=true
    
    elif [[ ${revision_update} == false ]] 
    then
        send_email=true
    fi

# Copyback everything to www folder if it is accessible    
    if [[ -r ${wwwtestdir}/${regtesthtml} ]]; then
      cp ${wwwtestdir}/* ${wwwdir}
    else
      echo "ERROR: www directory is not accessible!" &>>${reglog}
    fi

  fi
  
  if [[ ${send_email} == true ]] 
  then
    echo "Errors found: sending emails..."
  fi

# Send emails to the villains and owner of the regtester 
  if ${send_email}; then
    tolist=""
    if [[ ${email_developers} == true ]]  
    then
      for developer in $(egrep "[0-9]{4}-[0-9]{2}-[0-9]{2} [0-9]{2}:[0-9]{2}" ${wwwtestdir}/ChangeLog.diff | awk '{print $NL}' | uniq); do
        tolist="${tolist} $(grep ${developer} ${regtestdir}/addressbook.dat | awk '{print $2}')"
      done
      tolist=$(echo ${tolist} | tr -s [:space:] \\n | sort | uniq)
    fi
      mail -s "${config} ${date} ${time}" -a ${reglog} -a ${memlog} ${tolist} ${owner_mail} <<EOF

Dear CP2K developer,

$(echo ${tolist})

the automatic gfortran regtester has detected at least one issue due to a recent commit:

$(cat ${wwwtestdir}/ChangeLog.diff)

Please, check the attached regtester outputs for detailed information or visit

 http://cp2k-www.epcc.ed.ac.uk/

For complains please mailto: $(echo ${owner_mail})
EOF
  fi
  sleep ${sleep_time}
done
