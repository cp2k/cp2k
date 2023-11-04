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
  wget
rm -rf /var/lib/apt/lists/*

# Create and activate a virtual environment for Python packages.
python3 -m venv /opt/venv
export PATH="/opt/venv/bin:$PATH"

# Install Python packages. Upgrade via:
#   pip3 install black flask gunicorn mdformat-gfm cmake-format
#   pip3 freeze > requirements.txt
pip3 install --quiet -r requirements.txt

# Install shfmt.
# https://github.com/mvdan/sh
wget -q https://github.com/mvdan/sh/releases/download/v3.2.2/shfmt_v3.2.2_linux_amd64
echo '3a32a69286a19491a81fcd854154f0d886c379ff28d99e32d5594490b8bbef4b shfmt_v3.2.2_linux_amd64' | sha256sum --check
chmod +x shfmt_v3.2.2_linux_amd64
ln -s /opt/cp2k-precommit/shfmt_v3.2.2_linux_amd64 /usr/bin/shfmt

#EOF
