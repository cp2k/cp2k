#!/bin/bash -e

# author: Ole Schuett

# Install clang-format, shellcheck, and other Ubuntu packages.
# https://github.com/koalaman/shellcheck
# https://clang.llvm.org/docs/ClangFormat.html
export DEBIAN_FRONTEND=noninteractive DEBCONF_NONINTERACTIVE_SEEN=true
apt-get update -qq
apt-get install -qq --no-install-recommends \
  ca-certificates \
  clang-format \
  git \
  less \
  nano \
  python3 \
  python3-venv \
  python3-pip \
  python3-wheel \
  python3-setuptools \
  shellcheck \
  shfmt
rm -rf /var/lib/apt/lists/*

# Create and activate a virtual environment for Python packages.
python3 -m venv /opt/venv
export PATH="/opt/venv/bin:$PATH"

# TODO Add a pylock.toml file (https://peps.python.org/pep-0751/)
# Install Python packages. Upgrade via:
#   pip3 install black flask gunicorn cmakelang fortitude-lint \
#      mdformat mdformat-gfm mdformat-myst
#   pip3 freeze > requirements.txt
pip3 install --quiet -r requirements.txt

#EOF
