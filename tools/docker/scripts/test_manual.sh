#!/bin/bash -e

# author: Ole Schuett

# shellcheck disable=SC1091
source /opt/cp2k-toolchain/install/setup

echo -e "\n========== Compiling CP2K =========="
cd /opt/cp2k
echo -n "Compiling cp2k... "
if make -j VERSION=psmp &> make.out; then
  echo "done."
else
  echo -e "failed.\n\n"
  tail -n 100 make.out
  mkdir -p /workspace/artifacts/
  cp make.out /workspace/artifacts/
  echo -e "\nSummary: Compilation failed."
  echo -e "Status: FAILED\n"
  exit 0
fi

echo -e "\n========== Installing Dependencies =========="
apt-get update -qq
apt-get install -qq --no-install-recommends \
  default-jre-headless \
  libsaxonhe-java \
  python3 \
  python3-pip \
  python3-venv
rm -rf /var/lib/apt/lists/*

# Create and activate a virtual environment for Python packages.
python3 -m venv /opt/venv
export PATH="/opt/venv/bin:$PATH"

# install python packages
pip3 install -r /opt/cp2k/docs/requirements.txt

echo -e "\n========== Generating Manual =========="

mkdir -p /workspace/artifacts/manual
cd /workspace/artifacts/manual

/opt/cp2k/exe/local/cp2k.psmp --version
/opt/cp2k/exe/local/cp2k.psmp --xml

set +e # disable error trapping for remainder of script
(
  set -e # abort if error is encountered
  /opt/cp2k/docs/generate_input_reference.py ./cp2k_input.xml
  echo ""
  sphinx-build /opt/cp2k/docs/ /workspace/artifacts/manual -W -n --keep-going --jobs 16
  /opt/cp2k/docs/fix_github_links.py /workspace/artifacts/manual
  rm -rf /workspace/artifacts/manual/.doctrees
)
EXIT_CODE=$?
if ((EXIT_CODE)); then
  echo -e "\nSummary: Sphinx build failed"
  echo -e "Status: FAILED\n"
  exit
fi

echo -e "\n========== Generating VIM Plugin =========="
(
  set -e                            # abort if error is encountered
  sed -i 's/\x0/?/g' cp2k_input.xml # replace null bytes which would crash saxon
  SAXON="java -jar /usr/share/java/Saxon-HE.jar"
  $SAXON -o:cp2k.vim ./cp2k_input.xml /opt/cp2k/tools/input_editing/vim/vim.xsl
)
EXIT_CODE=$?
if ((EXIT_CODE)); then
  echo -e "\nSummary: Saxon run failed."
  echo -e "Status: FAILED\n"
  exit
fi

echo -e "\nSummary: Manual generation works fine."
echo -e "Status: OK\n"

#EOF
