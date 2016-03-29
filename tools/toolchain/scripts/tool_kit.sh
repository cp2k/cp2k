# A set of tools used in the toolchain installer, intended to be used
# by sourcing this file inside other scipts.

SYS_INCLUDE_PATH=${SYS_INCLUDE_PATH:-"/usr/local/include:/usr/include"}
SYS_LIB_PATH=${SYS_LIB_PATH:-"/user/local/lib64:/usr/local/lib:/usr/lib64:/usr/lib:/lib64:/lib"}
INCLUDE_PATHS=${INCLUDE_PATHS:-"CPATH SYS_INCLUDE_PATH"}
LIB_PATHS=${LIB_PATHS:-"LIBRARY_PATH LD_LIBRARY_PATH LD_RUN_PATH SYS_LIB_PATH"}

# report a warning message with script name and line number
report_warning() {
    if [ $# -gt 1 ] ; then
        local __lineno=", line $1"
        local __message="$2"
    else
        local __lineno=''
        local __message="$1"
    fi
    echo "WARNING: (${SCRIPT_NAME}${__lineno}) $__message" >&2
}

# report an error message with script name and line number
report_error() {
    if [ $# -gt 1 ] ; then
        local __lineno=", line $1"
        local __message="$2"
    else
        local __lineno=''
        local __message="$1"
    fi
    echo "ERROR: (${SCRIPT_NAME}${__lineno}) $__message" >&2
}

# error handler for line trap from set -e
error_handler() {
    local __lineno="$1"
    report_error $1 "Non-zero exit code detected."
    exit 1
}

# source a file if it exists, otherwise do nothing
load() {
    if [ -f "$1" ] ; then
        source "$1"
    fi
}

# A more portable command that will give the full path, removing
# symlinks, of a given path. This is more portable than readlink -f
# which does not work on Mac OS X
realpath() {
    local __path="$1"
    if [ "x$__path" = x ] ; then
       return 0
    fi
    local __basename=$(basename "$__path")
    if [ -e "$__path" ] ; then
        echo $(cd "$(dirname "$__path")" ; pwd -P)/"$__basename"
        return 0
    else
        return 1
    fi
}

# given a list, outputs a list with duplicated items filtered out
unique() (
    # given a list, outputs a list with duplicated items filtered out.
    # If -d <delimiter> option exists, then output the list delimited
    # by <delimiter>; note that this option does not effect the input.
    local __result=''
    local __delimiter=' '
    local __item=''
    if [ "$1" = "-d" ] ; then
        shift
        __delimiter="$1"
        shift
    fi
    # It is essential that we quote $@, which makes it equivalent to
    # "$1" "$2" ...  So this works if any of the arguments contains
    # space.  And we use \n to separate the fields in the
    # __result for now, so that fields that contain spaces are
    # correctly grepped.
    for __item in "$@" ; do
        if [ x"$__result" = x ] ; then
            __result="${__item}"
        # Note that quoting $__result after echo is essential to
        # retain the \n in the variable from the output of echo.  Also
        # remember grep only works on a line by line basis, so if
        # items are delimited by newlines, then for grep search it
        # should be delimited by ^ and $ (beginning and end of line)
        elif ! (echo "$__result" | \
                grep -s -q -e "^$__item\$") ; then
            __result="${__result}
${__item}"
        fi
    done
    __result="$(echo "$__result" | paste -s -d "$__delimiter" -)"
    # quoting $__result below is again essential for correct
    # behaviour if IFS is set to be the same $__delimiter in the
    # parent shell calling this macro
    echo "$__result"
)

# reverse a list
reverse() (
    # given a list, output a list with reversed order. If -d
    # <delimiter> option exists, then output the list delimited by
    # <delimiter>; note that this option does not effect the input.
    local __result=''
    local __delimiter=' '
    local __item=''
    if [ "$1" = "-d" ] ; then
        shift
        __delimiter="$1"
        shift
    fi
    for __item in "$@" ; do
        if [ x"$__result" = x ] ; then
            __result="$__item"
        else
            __result="${__item}${__delimiter}${__result}"
        fi
    done
    echo "$__result"
)

# get the number of processes avaliable for compilation
get_nprocs() {
    if $(command -v nproc >&- 2>&-) ; then
        echo $(nproc --all)
    else
        echo 1
    fi
}

# convert a list of paths to -L<dir> ... used by ld
paths_to_ld() {
    # need to define the POSIX default IFS values here, cannot just do
    # __ifs=$IFS first, because IFS can be unset, and if so __ifs will
    # becomes an empty string (null) and NOT unset, so later when IFS
    # is set to __ifs it becomes null rather than unset, and thus
    # causing wrong behaviour.  So if IFS is unset, __ifs should be
    # the POSIX default value.  Further more, due to shell
    # automatically remove the tailing "\n" in a string during
    # variable assignment, we need to add x after \n and then remove
    # it.
    local __paths=$@
    local __name=''
    local __raw_path=''
    local __dir=''
    local __lib_dirs=''
    # set default IFS first
    local __ifs=$(printf " \t\nx"); __ifs="${__ifs%x}"
    [ "$IFS" ] && __ifs="$IFS"
    for __name in $__paths ; do
        eval __raw_path=\$"$__name"
        # change internal field separator to :
        IFS=':'
        # loop over all dirs in path, and filter out duplications
        for __dir in $__raw_path ; do
            if ! [ x"$__dir" = x ] ; then
                if ! [[ "$__lib_dirs" =~ (^|[[:space:]])"-L'$__dir'"($|[[:space:]]) ]] ; then
                    __lib_dirs="$__lib_dirs -L'$__dir'"
                fi
            fi
        done
        IFS="$__ifs"
    done
    echo $__lib_dirs
}

# Find a file from directories given in a list of paths, each has the
# same format as env variable PATH. If the file is found, then echos
# the full path of the file. If the file is not found, then echos
# __FALSE__. The file name can also contain wildcards that are
# acceptable for bash, and in that case the full path of the first
# matching file will be echoed.
find_in_paths() {
    local __target=$1
    shift
    local __paths=$@
    local __name=''
    local __raw_path=''
    local __dir=''
    local __file=''
    local __files=''
    # use the IFS variable to take care of possible spaces in file/dir names
    local __ifs="$(printf " \t\nx")"; __ifs="${__ifs%x}"
    [ "$IFS" ] && __ifs="$IFS"
    for __name in $__paths ; do
        eval __raw_path=\$"$__name"
        # fields in paths are separated by :
        IFS=':'
        for __dir in $__raw_path ; do
            # files in possible glob expansion are to be delimited by "\n\b"
            IFS="$(printf "\nx")"; IFS="${IFS%x}"
            for __file in $__dir/$__target ; do
                if [ -e "$__file" ] ; then
                    echo $(realpath "$__file")
                    # must remember to change IFS back when exiting
                    IFS="$__ifs"
                    return 0
                fi
            done
            IFS=':'
        done
        IFS=$__ifs
    done
    echo "__FALSE__"
}

# search through a list of given paths, try to find the required file
# or directory, and if found then add full path of dirname file, or
# directory, to the -I include list for CFLAGS and append to a user
# specified variable (__cflags_name). If not found, then nothing is
# done. If the option -p is present, then if the search target is a
# directory, then the parent directory of the direcotry is used for -I
# instead.  The search target accepts bash wildcards, and in this case
# the first match will be used.
add_include_from_paths() {
    local __parent_dir_only=false
    if [ $1 = "-p" ] ; then
        __parent_dir_only=true
        shift
    fi
    local __cflags_name=$1
    shift
    local __search_target=$1
    shift
    local __paths=$@
    local __found_target=""
    local __cflags=""
    __found_target="$(find_in_paths "$__search_target" \
                                    $__paths)"
    if [ "$__found_target" != "__FALSE__" ] ; then
        if [ -f "$__found_target" ] || $__parent_dir_only ; then
           __found_target="$(dirname "$__found_target")"
        fi
        echo "Found include directory $__found_target"
        eval __cflags=\$"${__cflags_name}"
        __cflags="${__cflags} -I'${__found_target}'"
        # remove possible duplicates
        __cflags="$(unique $__cflags)"
        # must escape all quotes again before the last eval, as
        # otherwise all quotes gets interpreted by the shell when
        # assiging to variable because eval will reduce one escape
        # level
        __cflags="${__cflags//'/\\'}"
        eval $__cflags_name=\"$__cflags\"
    fi
}

# search through a list of given paths, try to find the required file
# or directory, and if found then add full path of dirname file, or
# directory, to the -L library list (including -Wl,-rpath) for LDFLAGS
# and append to a user specified variable (__ldflags_name). If not
# found, then nothing is done. If the option -p is present, then if
# the search target is a directory, then the parent directory of the
# direcotry is used for -L instead.  The search target accepts bash
# wildcards, and in this case the first match will be used.
add_lib_from_paths() {
    local __parent_dir_only=false
    if [ $1 = "-p" ] ; then
        __parent_dir_only=true
        shift
    fi
    local __ldflags_name=$1
    shift
    local __search_target=$1
    shift
    local __paths=$@
    local __found_target=""
    local __ldflags=""
    __found_target="$(find_in_paths "$__search_target" \
                                    $__paths)"
    if [ "$__found_target" != "__FALSE__" ] ; then
        if [ -f "$__found_target" ] || $__parent_dir_only ; then
           __found_target="$(dirname "$__found_target")"
        fi
        echo "Found lib directory $__found_target"
        eval __ldflags=\$"${__ldflags_name}"
        __ldflags="${__ldflags} -L'${__found_target}' -Wl,-rpath='${__found_target}'"
        # remove possible duplicates
        __ldflags="$(unique $__ldflags)"
        # must escape all quotes again before the last eval, as
        # otherwise all quotes gets interpreted by the shell when
        # assiging to variable because eval will reduce one escape
        # level
        __ldflags="${__ldflags//'/\\'}"
        eval $__ldflags_name=\"$__ldflags\"
    fi
}

# check if environment variable is assigned and non-empty
require_env() {
    local __env_var_name=$1
    local __env_var="$(eval echo \"\$$__env_var_name\")"
    if [ -z "$__env_var" ] ; then
        report_error "requires environment variable $__env_var_name to work"
        return 1
    fi
}

# check if a command is available
check_command() {
    local __command=$1
    if [ $# -eq 1 ] ; then
        local __package=$1
    elif [ $# -gt 1 ] ; then
        local __package=$2
    fi
    if $(command -v $__command >&- 2>&-) ; then
        echo "path to $__command is " $(command -v $__command)
    else
        report_error "cannot find $__command, please check if $__package is installed or in system search path"
        return 1
    fi
}

# check if directory exists
check_dir() {
    local __dir=$1
    if [ -d "$__dir" ] ; then
        echo "Found directory $__dir"
    else
        report_error "Cannot find $__dir"
        return 1
    fi
}

# check if a library can be found by ld, library names should in the
# format -lname, which would then referred to libname.a or libname.so
# by ld
check_lib() {
    local __libname="${1#-l}"
    if [ $# -eq 1 ] ; then
        local __package=lib"$__libname"
    elif [ $# -gt 1 ] ; then
        local __package=$2
    fi
    # Note that LD_LIBRARY_PATH is NOT used by ld during linking
    # stage, and is only used for searching to the shared libraries
    # requred by the executable AFTER it has already been compiled, to
    # override its internal search paths built into the binary when it
    # was compiled. Here, we explicitly include the commonly defined
    # library search paths---including LD_LIBRARY_PATH---in the -L
    # search paths of ld.  This is the only way ld can include
    # non-standard directories in its search path.  If we use gcc
    # instead of ld for linker then we can use LIBRARY_PATH, which IS
    # used during link stage. However, I think using ld is more
    # general, as in most systems LIBRARY_PATH is rarely defined, and
    # we would have to reply on gcc.
    local __search_engine="ld -o /dev/null"
    local __search_paths="$LIB_PATHS"
    # convert a list of paths to -L<dir> list used by ld
    __search_engine="$__search_engine $(paths_to_ld $__search_paths)"
    # needed the eval to interpret the quoted directories correctly (somehow)
    if (eval $__search_engine -l$__libname 2>&1 | grep -q -s "\-l$__libname") ; then
        # if library not found, ld will return error message
        # containing the library name
        report_error \
        "ld cannot find -l$__libname, please check if $__package is installed or in system search path"
        return 1
    else
        # if library is found, then ld will return error message about
        # not able to find _start or _main symbol
        echo "lib$__libname is found in ld search path"
    fi
}

# check if a module is available for the current version of gfortran,
# returns 0 if available and 1 if not
check_gfortran_module() {
    local __module_name=$1
    local __FC=${FC:-gfortran}
    cat <<EOF | $__FC -c -o /dev/null -xf95 -ffree-form - >&- 2>&-
PROGRAM check_gfortran_module
USE ${__module_name}
IMPLICIT NONE
PRINT *, "PASS"
END PROGRAM check_gfortran_module
EOF
}

# check if a flag is allowed for the current version of
# gfortran. returns 0 if allowed and 1 if not
check_gfortran_flag() {
    local __flag=$1
    local __FC=${FC:-gfortran}
    # no need to do a full compilation, just -E -cpp would do for
    # checking flags
    cat <<EOF | $__FC -E -cpp $__flag -xf95 -ffree-form - >&- 2>&-
PROGRAM test_code
  IMPLICIT NONE
  PRINT *, "PASS"
END PROGRAM test_code
EOF
}

# check if a flag is allowed for the current version of
# gcc. returns 0 if allowed and 1 if not
check_gcc_flag() {
    local __flag=$1
    local __CC=${CC:-gcc}
    # no need to do a full compilation, just -E -cpp would do for
    # checking flags
    cat <<EOF | $__CC -E -cpp $__flag -xc - >&- 2>&-
#include <stdio.h>
int main() {
  printf("PASS\n");
}
EOF
}

# check if a flag is allowed for the current version of
# g++. returns 0 if allowed and 1 if not
check_gxx_flag() {
    local __flag=$1
    local __CXX=${CXX:-g++}
    # no need to do a full compilation, just -E -cpp would do for
    # checking flags
    cat <<EOF | $__CC -E -cpp $__flag -xc - >&- 2>&-
#include <stdio.h>
int main() {
  printf("PASS\n");
}
EOF
}

# given a list of flags, only print out what is allowed by the current
# version of gfortran
allowed_gfortran_flags() {
    local __flags=$@
    local __flag=''
    local __result=''
    for __flag in $__flags ; do
        if (check_gfortran_flag $__flag) ; then
            [ -z "$__result" ] && __result="$__flag" || __result="$__result $__flag"
        fi
    done
    echo $__result
}

# given a list of flags, only print out what is allowed by the current
# version of gcc
allowed_gcc_flags() {
    local __flags=$@
    local __flag=''
    local __result=''
    for __flag in $__flags ; do
        if (check_gcc_flag $__flag) ; then
            [ -z "$__result" ] && __result="$__flag" || __result="$__result $__flag"
        fi
    done
    echo $__result
}

# given a list of flags, only print out what is allowed by the current
# version of g++
allowed_gxx_flags() {
    local __flags=$@
    local __flag=''
    local __result=''
    for __flag in $__flags ; do
        if (check_gxx_flag $__flag) ; then
            [ -z "$__result" ] && __result="$__flag" || __result="$__result $__flag"
        fi
    done
    echo $__result
}

# prepend a directory to a given path
prepend_path() {
    # prepend directory to $path_name and then export path_name. If
    # the directory already exists in path, bring the directory to the
    # front of the list.
    local __path_name=$1
    local __directory=$2
    if eval [ x\"\$$__path_name\" = x ] ; then
        eval $__path_name=\"$__directory\"
    else
        eval $__path_name=\"$__directory:\$$__path_name\"
    fi
    eval $__path_name=\"$(IFS=:; eval unique -d ':' \$$__path_name)\"
    export $__path_name
}

# append a directory to a given path
append_path() {
    # append directory to $path_name and then export path_name. If
    # the directory already exists in path, bring the directory to the
    # back of the list.
    #
    # $1 is path_name
    # $2 is directory
    local __path_name=$1
    local __directory=$2
    if eval [ x\"\$$__path_name\" = x ] ; then
        eval $__path_name=\"$__directory\"
    else
        eval $__path_name=\"\$$__path_name:$__directory\"
    fi
    # here, we reverse the list, apply unique, and then reverse back,
    # so that the last instance of the repeated directory is kept
    eval $__path_name=\"$(IFS=':'; \
                reverse -d ':' \
                        $(unique -d ':' \
                                 $(eval reverse -d ':' \$$__path_name)))\"
    export $__path_name
}

# helper routine for reading --enable=* input options
read_enable() {
    local __input_var="${1#*=}"
    case $__input_var in
        $1)
            # if there is no "=" then treat as "yes"
            echo "__TRUE__"
            ;;
        yes)
            echo "__TRUE__"
            ;;
        no)
            echo "__FALSE__"
            ;;
        *)
            echo "__INVALID__"
            ;;
    esac
}

# helper routine for reading --with=* input options
read_with() {
    local __input_var="${1#--with*=}"
    case $__input_var in
        $1)
            # if there is no "=" then treat as "install"
            echo '__INSTALL__'
            ;;
        install)
            echo '__INSTALL__'
            ;;
        system)
            echo '__SYSTEM__'
            ;;
        no)
            echo '__DONTUSE__'
            ;;
        *)
            echo "${__input_var//\~/$HOME}"
            ;;
    esac
}

# helper routine to check integrity of downloaded files
checksum() {
   local __filename=$1
   local __checksums=$2
   local __shasum_command='sha256sum'
   # check if we have sha256sum command, Mac OS X does not have
   # sha256sum, but has an equivalent with shasum -a 256
   command -v "$__shasum_command" >&- 2>&- || \
       __shasum_command="shasum -a 256"
   if eval "grep $__filename $__checksums | ${__shasum_command} --check" ; then
       echo "Checksum of $__filename Ok"
   else
      rm -v ${__filename}
      report_error "Checksum of $__filename could not be verified, abort."
      return 1
   fi
}

# downloader for the package tars, excludes checksum
download_pkg_no_checksum() {
    # usage: download_pkg_no_checksum [-n] [-o output_filename] url
    local __wget_flags=''
    local __url=''
    while [ $# -ge 1 ] ; do
        case "$1" in
            -n)
                local __wget_flags="$__wget_flags --no-check-certificate"
                ;;
            -o)
                shift
                __wget_flags="$__wget_flags -O $1"
               ;;
            *)
                __url="$1"
                ;;
        esac
        shift
    done
    # download
    if ! wget $__wget_flags $__url ; then
        report_error "failed to download $__url"
        return 1
    fi
}

# downloader for the package tars, includes checksum
download_pkg() {
    # usage: download_pkg [-n] [-o output_filename] url
    local __wget_flags=''
    local __filename=''
    local __url=''
    while [ $# -ge 1 ] ; do
        case "$1" in
            -n)
                __wget_flags="$__wget_flags --no-check-certificate"
                ;;
            -o)
                shift
                __wget_flags="$__wget_flags -O $1"
                __filename="$1"
               ;;
            *)
                __url="$1"
                ;;
        esac
        shift
    done
    if [ "$__filename" = "" ] ; then
        __filename="$(basename $__url)"
    fi
    # env variable for checksum file must be provided
    require_env SHA256_CHECKSUMS
    # download
    if ! wget $__wget_flags $__url ; then
        report_error "failed to download $__url"
        return 1
    fi
    # checksum
    checksum "$__filename" "$SHA256_CHECKSUMS"
}
