#/bin/bash
# cp2g_page_update
# Updates main page of the CP2K regtester

# Those lines should point to the correct locations
regtestdir=/disk/cp2k-www0/home/cp2ktest/regtest/remote-testing/
wwwdir=/disk/cp2k-www0/websites/cp2k-www.epcc.ed.ac.uk/
wwwtestdir=/disk/cp2k-www0/home/cp2ktest/regtest/remote-testing/www/

# Page refresh time
sleep_time="600"

# Spin endlessly
while true
do

# Timestamp for page update
  datetime=$(date "+%F %H:%M:%S")

# Copy the template
  cp ${regtestdir}/regtest-main-template.html ${wwwtestdir}/index.html

# Read the list of testers one by one
  while read -r t
  do
    
# If configuration line begins with a comment, skip it
    if [[ ${t:0:1} ==  '#' ]]
    then
      continue
    fi

# Parse configuration line
    name=`echo $t | cut -d":" -f1`
    localdir=`echo $t | cut -d":" -f2`
    remotedir=`echo $t | cut -d":" -f3`

    localdir="${regtestdir}/${localdir}"

# reglog and memlog files
    reglog=`ls ${localdir} | grep regtester_log.txt | tail -n 1`
    memlog=`ls ${localdir} | grep memory_leaks_log.txt | tail -n 1`

# Extracting info from regtester log
    date=`cat ${localdir}/${reglog} | grep "Date: " | sed 's/Date: \([0-9\-]*\)/\1/'`
    time=`cat ${localdir}/${reglog} | grep "Time: " | sed 's/Time: \([0-9\:]*\)/\1/'`
    config=`cat ${localdir}/${reglog} | grep "Configuration: " | sed 's/Configuration: \(.*\)/\1/'`

    if [[ $(grep -m 1 "At revision" ${localdir}/${reglog}) ]]; then
      revision=$(grep -m 1 "At revision" ${localdir}/${reglog} | cut -c13-17)
    elif [[ $(grep -m 1 "Updated to revision" ${localdir}/${reglog}) ]]; then
      revision=$(grep -m 1 "Updated to revision" ${localdir}/${reglog} | cut -c21-25)
    else
      revision="Unknown"
    fi

    (( runtime_errors = 0 ))
    (( wrong_results = 0 ))
    (( correct_tests = 0 ))
    (( new_inputs = 0 ))
    (( number_of_test_inputs = 0 ))
    (( memory_leaks = 0 ))

    if [[ $(grep "GREPME" ${localdir}/${reglog}) ]]; then
      set -- $(grep "GREPME" ${localdir}/${reglog})
      runtime_errors=$2
      wrong_results=$3
      correct_tests=$4
      new_inputs=$5
      number_of_test_inputs=$6
      memory_leaks=$7
      regtester_problem=false
    else
      regtester_problem=true
    fi 

    if [[ ${memory_leaks} == "X" ]]
    then
      memory_leaks="Not tested"
    fi

# Link to regtest and arch file pages
    regtest_link="<a href=\"${remotedir}\/regtest.html\">${name}<\/a>"

    config_link="<a href=\"${remotedir}\/regtest-arch\">${config}<\/a>"

    if [[ ${regtester_problem} == true ]]
    then
      regtest_link="${regtest_link}*"
    else
      regtest_link="${regtest_link}"
    fi

    table_row="<tr align=\"center\"><td align="left">${regtest_link}<\/td><td align="left">${config_link}<\/td><td>${number_of_test_inputs}<\/td><td>${correct_tests}<\/td><td>${wrong_results}<\/td><td>${new_inputs}<\/td><td>${runtime_errors}<\/td><td>${memory_leaks}<\/td><td>${revision}<\/td><td>${date} ${time}<\/td><\/tr>"
    
    sed -i "/{{TAG_TESTERS_ROW}}/i ${table_row}" ${wwwtestdir}/index.html

  done < cp2k_regtesters_list.conf

# Get results from http://people.web.psi.ch/krack/cp2k/regtest/regtest.html
  
  wget -q -N -P ${regtestdir}/mkrack-popt/www/ http://people.web.psi.ch/krack/cp2k/regtest/regtest-0
  if [[ -f ${regtestdir}/mkrack-popt/www/regtest-0 ]]
  then
    date=`cat ${regtestdir}/mkrack-popt/www/regtest-0 | grep "Date: " | sed 's/Date: \([0-9\-]*\)/\1/'`
    time=`cat ${regtestdir}/mkrack-popt/www/regtest-0 | grep "Time: " | sed 's/Time: \([0-9\:]*\)/\1/'`
    
    (( runtime_errors = 0 ))
    (( wrong_results = 0 ))
    (( correct_tests = 0 ))
    (( new_inputs = 0 ))
    (( number_of_test_inputs = 0 ))
    (( memory_leaks = 0 ))

    if [[ $(grep "GREPME" ${regtestdir}/mkrack-popt/www/regtest-0) ]]; then
      set -- $(grep "GREPME" ${regtestdir}/mkrack-popt/www/regtest-0)
      runtime_errors=$2
      wrong_results=$3
      correct_tests=$4
      new_inputs=$5
      number_of_test_inputs=$6
      memory_leaks=$7
      regtester_problem=false
    else
      regtester_problem=true
    fi 
  
    regtest_link="<a href=\"http://people.web.psi.ch/krack/cp2k/regtest/regtest.html\">PSI, Switzerland<\/a>"

    config_link="Linux-x86-64-g95-popt"

    table_row="<tr align=\"center\"><td align="left">${regtest_link}<\/td><td align="left">${config_link}<\/td><td>${number_of_test_inputs}<\/td><td>${correct_tests}<\/td><td>${wrong_results}<\/td><td>${new_inputs}<\/td><td>${runtime_errors}<\/td><td>${memory_leaks}<\/td><td>${revision}<\/td><td>${date} ${time}<\/td><\/tr>"
    
    sed -i "/{{TAG_TESTERS_ROW}}/i ${table_row}" ${wwwtestdir}/index.html
  fi

# Finally, strip {{TAG_TESTERS_ROW}}
  
  sed -i "/{{TAG_TESTERS_ROW}}/d" ${wwwtestdir}/index.html

# And update page date and time
  sed -i "s/{{TAG_DATE_TIME}}/${datetime}/" ${wwwtestdir}/index.html

# Copy back index.html to www dir
  cp ${wwwtestdir}/index.html ${wwwdir}

  sleep ${sleep_time}
done
