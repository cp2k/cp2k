#!/bin/bash -e 

# Test if all templates are instantiated and prettify was run
# author: Ole Schuett

echo -n "Date: "
date --utc --rfc-3339=seconds

svn up
svn info

find ./src/ -type f                        -not -path "*/.svn/*" -print0 | xargs -0 md5sum > checksums.md5
find ./src/ -type f -name "*instantiation" -not -path "*/.svn/*" -print0 | xargs -0 ./tools/instantiateTemplates.py

cd makefiles
make --jobs=20 pretty
cd ..

nfiles=`wc -l checksums.md5 | cut -f 1 -d " "`
summary="Checked $nfiles files."
status="OK"

if grep -r "UNMATCHED_PROCEDURE_ARGUMENT" ./src/* ; then
  summary="Found UNMATCHED_PROCEDURE_ARGUMENT"
  status="FAILED"
fi

if ! md5sum --quiet --check checksums.md5 ; then
  summary="Code not invariant under instantiateTemplates.py or prettify.py"
  status="FAILED"
fi

rm checksums.md5

echo "Summary:" $summary
echo "Status:" $status

#EOF
