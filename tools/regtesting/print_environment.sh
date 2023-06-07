#!/bin/bash -e

if (($# != 2)); then
  echo "Usage: print_environment.sh <ARCH> <VERSION>"
  echo "Always run from the cp2k repository root directory."
  exit 1
fi

ARCH=$1
VERSION=$2

echo "--------------------------- GIT ------------------------------------------"

git_sha="<N/A>"
submodule_sha=()

if [[ -d ".git" ]]; then
  head_ref=$(< ".git/HEAD")

  if [[ ${head_ref} =~ ^ref:\ (.*) ]]; then
    ref_file=".git/${BASH_REMATCH[1]}"

    if [[ -f "${ref_file}" ]]; then
      git_sha=$(< "${ref_file}")
    fi
  else
    # this is a detached HEAD, no further deref needed
    git_sha="${head_ref}"
  fi

  modules_git_dir=".git/modules/"

  # find submodule indexes, if any
  while IFS= read -r -d $'\0'; do
    submodule_basedir="${REPLY%/*}"
    submodule_hash=$(< "${submodule_basedir}/HEAD")
    submodule_dir="${submodule_basedir#"${modules_git_dir}"}"
    submodule_sha+=("${submodule_dir}:${submodule_hash}")
  done < <(find "${modules_git_dir}" -name index -print0)
fi

echo "CommitSHA: ${git_sha}"

for submodule in "${submodule_sha[@]}"; do
  echo "Commmit SHA of submodule in ${submodule%%:*}: ${submodule##*:}"
done

# add some more information about that last commit if `git` is available
if [[ "${git_sha}" != "<N/A>" ]] && command -v git > /dev/null 2>&1; then
  GIT_WORK_TREE="." GIT_DIR=".git" \
    git --no-pager \
    log -n 1 \
    --format="CommitTime: %ci%nCommitAuthor: %an%nCommitAuthorEmail: %ae%nCommitSubject: %s%n" \
    "${git_sha}"
fi

echo "------------------------- System information -----------------------------"
uname -mrs
case "$(uname -s)" in
  Linux)
    ldd --version | head -n 1
    ;;
  Darwin)
    sw_vers
    ;;
esac
echo ""

echo "--------------------------- Resource limits ------------------------------"
case "$(uname -s)" in
  Linux)
    prlimit
    ;;
  Darwin)
    launchctl limit
    ;;
esac
echo ""

echo "--------------------------- SELinux --------------------------------------"
if [[ -f /usr/sbin/getenforce ]]; then
  echo "SELinux is installed and is $(/usr/sbin/getenforce)"
else
  echo "No SELinux installation found"
fi
echo ""

echo "--------------------------- ARCH-file ------------------------------------"
cat "./arch/${ARCH}.${VERSION}"

if [ -x "$(command -v nvidia-smi)" ]; then
  echo "--------------------------- NVIDIA-SMI -----------------------------------"
  nvidia-smi || true
  echo ""
fi

if [ -x "$(command -v rocm-smi)" ]; then
  echo "--------------------------- ROCm-SMI -------------------------------------"
  rocm-smi
  echo ""
fi

echo "-------------------------- Build-Tools -----------------------------------"
make toolversions ARCH="${ARCH}" VERSION="${VERSION}"

echo "----------------------- External Modules ---------------------------------"
make extversions ARCH="${ARCH}" VERSION="${VERSION}"

echo "---------------------------- Modules -------------------------------------"
if [ "$(type -t module)" = 'function' ]; then
  module list 2>&1
else
  echo "Module system not installed."
fi

#EOF
