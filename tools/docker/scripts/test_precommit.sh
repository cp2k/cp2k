#!/bin/bash

# author: Ole Schuett

echo -e "\n========== Starting Precommit Server =========="
cd ./tools/precommit/ || exit 1
export REVISION="unknown revision"
gunicorn --bind=:8080 --workers=1 --threads=8 --timeout=0 precommit_server:app &> /var/tmp/precommit_server.logs &
sleep 3
cat /var/tmp/precommit_server.logs

echo -e "\n========== Running Precommit Checks =========="
export CP2K_PRECOMMIT_SERVER="http://127.0.0.1:8080"
./precommit.py --no-cache --progressbar-wait=10
echo -e "\n"

exit 0 # Prevent CI from overwriting precommit.py's summary message.

#EOF
