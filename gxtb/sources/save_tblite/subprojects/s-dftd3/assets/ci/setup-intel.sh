set -ex
if [ $(uname) = Darwin ]; then
  OUT=webimage-base.dmg
  URL=https://registrationcenter-download.intel.com/akdlm/irc_nas/17969/m_BaseKit_p_2021.3.0.3043.dmg
  COMPONENTS=intel.oneapi.mac.mkl.devel
  curl --output $OUT --url "$URL" --retry 5 --retry-delay 5
  hdiutil attach $OUT
  if [ -z "$COMPONENTS" ]; then
    sudo /Volumes/"$(basename "$URL" .dmg)"/bootstrapper.app/Contents/MacOS/bootstrapper -s --action install --eula=accept --continue-with-optional-error=yes --log-dir=.
    installer_exit_code=$?
  else
    sudo /Volumes/"$(basename "$URL" .dmg)"/bootstrapper.app/Contents/MacOS/bootstrapper -s --action install --components="$COMPONENTS" --eula=accept --continue-with-optional-error=yes --log-dir=.
    installer_exit_code=$?
  fi
  hdiutil detach /Volumes/"$(basename "$URL" .dmg)" -quiet

  URL=https://registrationcenter-download.intel.com/akdlm/irc_nas/17890/m_HPCKit_p_2021.3.0.3226_offline.dmg
  COMPONENTS=all
  curl --output $OUT --url "$URL" --retry 5 --retry-delay 5
  hdiutil attach $OUT
  if [ -z "$COMPONENTS" ]; then
    sudo /Volumes/"$(basename "$URL" .dmg)"/bootstrapper.app/Contents/MacOS/bootstrapper -s --action install --eula=accept --continue-with-optional-error=yes --log-dir=.
    installer_exit_code=$?
  else
    sudo /Volumes/"$(basename "$URL" .dmg)"/bootstrapper.app/Contents/MacOS/bootstrapper -s --action install --components="$COMPONENTS" --eula=accept --continue-with-optional-error=yes --log-dir=.
    installer_exit_code=$?
  fi
  hdiutil detach /Volumes/"$(basename "$URL" .dmg)" -quiet
  exit $installer_exit_code
else
  KEY=GPG-PUB-KEY-INTEL-SW-PRODUCTS-2023.PUB
  wget https://apt.repos.intel.com/intel-gpg-keys/$KEY
  sudo apt-key add $KEY
  rm $KEY
  echo "deb https://apt.repos.intel.com/oneapi all main" | sudo tee /etc/apt/sources.list.d/oneAPI.list
  sudo apt-get update
  sudo apt-get install \
    intel-oneapi-compiler-fortran-2021.2.0 \
    intel-oneapi-compiler-dpcpp-cpp-and-cpp-classic-2021.2.0
fi
