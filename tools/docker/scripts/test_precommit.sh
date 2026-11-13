#!/bin/bash

# author: Ole Schuett

echo -e "\n========== Starting Precommit Server =========="
cd ./tools/precommit/ || exit 1
export REVISION="unknown revision"
./precommit_server.py --port=8080 &> /var/tmp/precommit_server.logs &
sleep 3
cat /var/tmp/precommit_server.logs

echo -e "\n========== Running Precommit Checks =========="
export CP2K_PRECOMMIT_SERVER="http://127.0.0.1:8080"
./precommit.py --no-cache --plot-file-count --progressbar-wait=10 -j 32
EXIT_CODE=$?

if ((EXIT_CODE)); then
  echo ""
  echo "********************************************"
  echo "*                                          *"
  echo "*  O/     Try running ./make_pretty.sh     *"
  echo "* /|                                       *"
  echo "* / \                                      *"
  echo "********************************************"
fi

echo -e "\n"
exit 0 # Prevent CI from overwriting precommit.py's summary message.

#EOF
