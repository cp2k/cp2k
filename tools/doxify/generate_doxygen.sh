#!/bin/bash -e

# author: Ole Schuett

DOXYGENDIR="./doxygen"

while IFS= read -r -d '' SRC_FN; do
  echo "Preprocessing ${SRC_FN}"
  DEST_FN="${DOXYGENDIR}/${SRC_FN}"
  mkdir -p "$(dirname "${DEST_FN}")"
  if [[ "${SRC_FN}" == *.F ]]; then
    ./tools/build_utils/fypp "${SRC_FN}" "${DEST_FN}" &> /dev/null
  else
    cp "${SRC_FN}" "${DEST_FN}"
  fi
done < <(find ./src/ -type f -print0)

cp tools/logo/cp2k_logo_100.png ${DOXYGENDIR}/cp2k_logo.png

REVISION=$(./tools/build_utils/get_revision_number ./src)
sed -e "s/#revision#/${REVISION}/" ./tools/doxify/Doxyfile.template > ${DOXYGENDIR}/Doxyfile

echo "Running Doxygen..."
cd ${DOXYGENDIR}
doxygen ./Doxyfile

# Remove CP2K from header bars because we already have the logo.
cd html
sed -i 's/"projectname">CP2K/"projectname">/' ./*.html ./*/*.html ./*/*/*.html

#EOF
