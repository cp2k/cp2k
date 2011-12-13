#!/bin/sh
# clean_cwd: Cleans the current working directory and all its subdirectories from scratch files,
#            i.e. files which are not committed to SVN
echo "The directory $PWD is cleaned from scratch files"
if [[ $1 == "-f" ]]; then
  echo "NOTE: You requested to remove all files without being prompted"
fi
for file in $(svn status | grep '? ' | awk '{print $2}'); do
  if [[ -d ${file} ]]; then
    if [[ $1 == "-f" ]]; then
      rm -rf ${file}
      echo Directory ${file} removed
    elif [[ $1 == "-i" ]]; then
      rm -ir ${file}
    else
      rm -r ${file}
      echo Directory ${file} removed
    fi
  elif [[ -f ${file} ]]; then
    if [[ $1 == "-f" ]]; then
      rm -f ${file}
      echo File ${file} removed
    elif [[ $1 == "-i" ]]; then
      rm -i ${file}
    else
      if [[ ${file##*.} == "inp" ]]; then
        echo "This files could be a new input file. Do you really want to delete it?"
        rm -i ${file}
      elif [ ${file##*.} == "sopt" ] ||
           [ ${file##*.} == "popt" ] ||
           [ ${file##*.} == "sdbg" ] || 
           [ ${file##*.} == "pdbg" ] ||
           [ ${file##*.} == "ssmp" ] ||
           [ ${file##*.} == "psmp" ]; then
        echo "This files could be a new arch file. Do you really want to delete it?"
        rm -i ${file}
      else
        rm -f ${file}
        echo File ${file} removed
      fi
    fi
  else
    echo Unknown object ${file} found
    exit
  fi
done
