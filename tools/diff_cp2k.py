#!/usr/bin/env python3

# Compare CP2K outputs
# Author: Alfio Lazzaro
# Email: alfio.lazzaro@mat.ethz.ch
# Year: 2016

# Example 1: show timings for a CP2K output
#       > diff_cp2k.py <name_file>
#    It shows the values of the timings for the MAXIMUM SELF
#    timings as extracted from the final table of timings of
#    CP2K output. The values are sorted (only values >0).
#    You can use the option
#       -f <1 || 2 || 3 || 4>
#    to change between AVERAGE SELF (1), MAX SELF (2), AVERAGE TOTAL (3) or MAX TOTAL (4).
#    The last line CP2K_Total refers always to the MAXIMUM TOTAL TIME.
#    There is also the possibility to filter between the SUBROUTINE names
#    by using the options:
#       -g <name> : partial comparison (similar to linux command grep)
#       -e <name> : exact comparison
#    (regexp are not implemented)
#
# Example 2: compare two (or more) CP2K outputs
#       > diff_cp2k.py <list of files>
#    You can use wild cards (for example *.out).
#    It shows the timings from all outputs, sorted by the values
#    of the first file, which is the reference for the comparison.
#    It also shows the relative difference (in percentage) with respect to
#    the reference values. Colors/bold are used to easy spot the larger discrepancies:
#       blue: smaller than reference
#       blue bold: smaller than reference > 100%
#       green: bigger than reference
#       green bold: bigger than reference > 100%
#    A set of dashes "-------" are reported for SUBROUTINES that
#    are only in the reference file, while the SUBROUTINES
#    that are only in the other files are reported for each file at the end.
#    You can use the option
#       -b <#>
#    to change the file used as reference (default is 1).
#    The other options mentioned in Example 1 are still valid here.
#    It is possible to replace the SUBROUTINE names. This feature allows, for example,
#    to compare SUBROUTINEs with different names belonging to different files.
#    Create a file, called diff_cp2k_keys.py, where you declare
#    the SUBROUTINE names and their replacements, e.g.
#        special_keys={'cannon_multiply_low_rma_metroc':'cannon_multiply_low_metrocomm1' ,
#                      'cannon_multiply_low_rma':'cannon_multiply_low'}
#    In this case the SUBROUTINE with name cannon_multiply_low_rma_metroc will be
#    replaced by the name cannon_multiply_low_metrocomm1.
#    The file is automatically loaded from the local directory where
#    you run the script or from the home directory. Alternatively it is possible
#    to use the option
#       -k <file keys>
#    to specify a different file.
#
# Example 3: grep for some other values
#    As described in Example 2, create a file, called diff_cp2k_keys.py,
#    where you declare the keywords that you want to grep from the output, e.g.
#       stats_keys={'flops total':[0],'average stack size':[1,2]}
#    The script splits the line by the keyword in two parts and reports
#    the field at the given position of the right part.
#    The file is automatically loaded from the local directory where
#    you run the script or from the home directory. Alternatively it is possible
#    to use the option
#       -k <file keys>
#    to specify a different file.
#    The values will appear under "Stats report".
#

import sys
import argparse
import operator
import os
import imp


def read_file(filename, field, special_keys, stats_keys):
    try:
        nline = 0
        nstats = 0
        dict_values = {}
        dict_stats = {}
        nameout = ["", ""]
        with open(filename, "r") as f:
            for line in f:
                # start reading
                if "NAMEOUT=" in line:
                    nameout[0] = line.split("=", 2)[1].strip()
                    continue
                if "ENERGY| Total FORCE_EVAL ( QS ) energy (a.u.):" in line:
                    nameout[1] = line.split(":", 2)[1].strip()
                    continue
                if "DBCSR STATISTICS" not in line and nstats == 0:
                    continue
                nstats = 1
                for stats_key in stats_keys:
                    if stats_key in line:
                        for index in stats_keys[stats_key]:
                            index_key = stats_key.strip() + " [" + str(index) + "]"
                            dict_stats[index_key] = line.split(stats_key, 2)[1].split()[
                                index
                            ]
                        break
                if "T I M I N G" not in line and nline == 0:
                    continue
                nline += 1
                if nline < 6:
                    continue
                # end reading
                if "-----" in line:
                    nline = 0
                    continue
                values = line.split()
                # filter
                if float(values[3 + field]) <= 0.001 and values[0] != "CP2K":
                    continue
                if values[0] in special_keys:
                    values[0] = special_keys[values[0]]
                # take only he first timing of duplicate special_keys
                if values[0] in dict_values:
                    continue
                if values[0] == "CP2K":
                    dict_values[values[0] + "_Total"] = float(values[6])
                else:
                    dict_values[values[0]] = float(values[3 + field])

        f.closed
        return dict_values, dict_stats, nameout
    except IOError:
        print("Cannot open " + filename)
        print("Exit")
        sys.exit(-1)


def print_value(ref, value):
    if ref > 0:
        comp = (value - ref) / ref * 100
    else:
        comp = float("Inf")
    color = "\033[0m"
    endc = "\033[0m"
    if comp > 0:
        color = "\033[92m"
    elif comp < 0:
        color = "\033[94m"
    if abs(comp) > 100:
        color += "\033[1m"
    sys.stdout.write(color + "%10.3f" % value + "%5.0f" % comp + endc)


#################
# Main function #
#################


def main():
    parser = argparse.ArgumentParser(description="Comparison of CP2K output timings.")
    parser.add_argument("file_lists", nargs="+", help="list of files")
    parser.add_argument(
        "-f",
        metavar="field",
        type=int,
        dest="field",
        choices=range(1, 5),
        default=2,
        help="which field to show (default is 2)",
    )
    parser.add_argument(
        "-b",
        metavar="base",
        type=int,
        dest="base",
        default=1,
        help="which file to use as base for the comparison (default is 1)",
    )
    parser.add_argument(
        "-g",
        metavar="grep",
        nargs="+",
        dest="grep",
        default="",
        help="Fields to grep (check the inclusion correspondance of the words)",
    )
    parser.add_argument(
        "-e",
        metavar="filter",
        nargs="+",
        dest="filter",
        default="",
        help="Fields to grep (check the exact correspondance of the words)",
    )
    parser.add_argument(
        "-k", metavar="file_keys", dest="file_keys", default="", help="File of keys"
    )
    args = parser.parse_args()

    # Empty keys by default
    special_keys = {}
    stats_keys = {}

    # Check for keys file
    file_keys = []
    if len(args.file_keys) > 0:
        file_keys.append(os.path.abspath(args.file_keys))
    else:
        # if not file_keys is provided, then look for it in the local directory and home
        file_keys.append(os.getcwd() + "/diff_cp2k_keys.py")
        file_keys.append(os.path.expanduser("~") + "/diff_cp2k_keys.py")

    for filename in file_keys:
        try:
            module = imp.load_source("*", filename)
            special_keys = module.special_keys
            stats_keys = module.stats_keys
        except IOError:
            if len(args.file_keys) > 0:
                print("Cannont open file keys " + filename + "!")
                print("Exit")
                sys.exit(-1)

    if args.base < 1 or args.base > len(args.file_lists):
        print(
            "Value for -b option out-of-bounds! Allowed values are between 1 and "
            + str(len(args.file_lists))
        )
        print("Exit")
        sys.exit(-1)

    dict_values = {}
    dict_stats = {}
    files = {}
    for filename in args.file_lists:
        dict_values[filename], dict_stats[filename], files[filename] = read_file(
            filename, args.field - 1, special_keys, stats_keys
        )

    print("===== Timings report =====")

    # sorted by first file timings
    sorted_values = sorted(
        dict_values[args.file_lists[args.base - 1]].items(), key=operator.itemgetter(1)
    )
    for key in sorted_values:
        # Apply filtering
        if key[0] != "CP2K_Total" and (
            (len(args.grep) > 0 and any(s not in key[0] for s in args.grep))
            or (len(args.filter) > 0 and key[0] not in args.filter)
        ):
            continue
        sys.stdout.write(key[0].ljust(30) + "%10.3f" % key[1])
        for filename in args.file_lists:
            if filename == args.file_lists[args.base - 1]:
                continue
            if key[0] not in dict_values[filename]:
                sys.stdout.write(("-" * 10).rjust(15))
                continue
            print_value(key[1], dict_values[filename][key[0]])
            del dict_values[filename][key[0]]
        print("")

    print("")

    ref = 0
    if len(files[args.file_lists[args.base - 1]][1]) > 0:
        ref = float(files[args.file_lists[args.base - 1]][1])
    color = "\033[0m"
    endc = "\033[0m"
    for filename in args.file_lists:
        if len(files[filename][1]) > 0:
            comp = (float(files[filename][1]) - ref) / ref
            if abs(comp) > 1e-14:
                color = "\033[91m"
            else:
                color = "\033[0m"
            print(
                ("{0} ==> {1} : {2} : " + color + "{3}" + endc).format(
                    files[filename][0],
                    filename,
                    files[filename][1],
                    (float(files[filename][1]) - ref) / ref,
                )
            )
        else:
            print(("{0} ==> {1} : ").format(files[filename][0], filename)),
            sys.stdout.write(("-" * 20).rjust(20))
            print("")

    print("")

    for filename in args.file_lists:
        if filename == args.file_lists[args.base - 1]:
            continue
        print("Remaining entries in " + files[filename][0] + " ==> " + filename)
        sorted_values = sorted(
            dict_values[filename].items(), key=operator.itemgetter(1)
        )
        count = 0
        for key in sorted_values:
            # Apply filtering
            if (len(args.grep) > 0 and any(s not in key[0] for s in args.grep)) or (
                len(args.filter) > 0 and key[0] not in args.filter
            ):
                continue
            print(key[0].ljust(30) + "%10.3f" % key[1])
            count += 1
        if count == 0:
            print("<None>")
        print("")

    print("===== Stats report =====")

    if len(stats_keys) > 0:
        for stats_key in stats_keys:
            for index in stats_keys[stats_key]:
                index_key = stats_key.strip() + " [" + str(index) + "]"
                sys.stdout.write(index_key.ljust(35))
                for filename in args.file_lists:
                    if index_key not in dict_stats[filename]:
                        sys.stdout.write(("-" * 18).ljust(20))
                        continue
                    sys.stdout.write(dict_stats[filename][index_key].ljust(20))
                print("")
    else:
        print("<None>")

    print("")


# ===============================================================================
main()
