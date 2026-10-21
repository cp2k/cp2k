#!/usr/bin/env python3

# Compare CP2K outputs
# Author: Alfio Lazzaro
# Email: alfio.lazzaro@hpe.com
# Year: 2016-2022

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
#       green: faster than reference (beyond threshold)
#       red: slower than reference (beyond threshold)
#       plain: within +/- threshold of the reference (neutral)
#       bold: at least 2x faster or slower than reference
#    You can use the option
#       -t <percent>
#    to change the threshold (default is 5%).
#    Colors are disabled if the output is not a terminal or NO_COLOR is set.
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
#    The values will appear under "Stats report" (omitted if stats_keys is empty).
#    Only lines after the "DBCSR STATISTICS" banner are considered.
#

import sys
import argparse
import operator
import os
import importlib.util
from importlib.machinery import SourceFileLoader


def read_file(filename, field, special_keys, stats_keys):
    try:
        nline = 0
        nstats = 0
        dict_values = {}
        dict_stats = {}
        nameout = ["", ""]
        with open(filename, "r", encoding="utf8") as f:
            for line in f:
                # start reading
                if "NAMEOUT=" in line:
                    nameout[0] = line.split("=", 2)[1].strip()
                    continue
                if "ENERGY| Total FORCE_EVAL ( QS ) energy" in line:
                    split = line.split()
                    nameout[1] = split[-1]  # last item
                    continue
                if "DBCSR STATISTICS" in line:
                    nstats = 1
                if nstats != 0:
                    for stats_key in stats_keys:
                        if stats_key in line:
                            fields = line.split(stats_key, 1)[1].split()
                            for index in stats_keys[stats_key]:
                                index_key = stats_key.strip() + " [" + str(index) + "]"
                                if index < len(fields):
                                    dict_stats[index_key] = fields[index]
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

        return dict_values, dict_stats, nameout
    except IOError:
        print("Cannot open " + filename)
        print("Exit")
        sys.exit(-1)


def ansi(code, colors):
    return "\033[" + code + "m" if colors else ""


def print_value(ref, value, show_comp, threshold, colors):
    if ref > 0:
        comp = (value - ref) / ref * 100
    else:
        comp = float("Inf")
    color = ansi("0", colors)  # reset
    nocol = ansi("0", colors)  # reset
    # timings: lower is better
    if comp < -threshold:
        color = ansi("92", colors)  # green
    elif comp > threshold:
        color = ansi("91", colors)  # red
    # at least 2x faster or slower
    if abs(comp) > threshold and (comp <= -50 or comp >= 100):
        color += ansi("1", colors)  # bold
    if show_comp:
        sys.stdout.write(color + "%10.3f" % value + "%5.0f" % comp + nocol)
    else:
        sys.stdout.write(color + "%10.3f" % value + nocol)


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
        "-t",
        metavar="threshold",
        type=float,
        dest="threshold",
        default=5.0,
        help="relative difference in percent considered neutral (default is 5)",
    )
    parser.add_argument(
        "-k", metavar="file_keys", dest="file_keys", default="", help="File of keys"
    )
    parser.add_argument(
        "--show_comp",
        action="store_true",
        default=False,
        help="show comparison values",
    )
    args = parser.parse_args()
    colors = sys.stdout.isatty() and not os.environ.get("NO_COLOR")

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

    # first file found takes precedence (local directory before home)
    for filename in file_keys:
        try:
            loader = SourceFileLoader("diff_cp2k_keys", filename)
            spec = importlib.util.spec_from_file_location(
                "diff_cp2k_keys", filename, loader=loader
            )
            module = importlib.util.module_from_spec(spec)
            loader.exec_module(module)
            special_keys = getattr(module, "special_keys", {})
            stats_keys = getattr(module, "stats_keys", {})
            break
        except IOError:
            if len(args.file_keys) > 0:
                print("Cannot open file keys " + filename + "!")
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
    ndash = 10 if args.show_comp else 5
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
                sys.stdout.write(("-" * ndash).rjust(5 + ndash))
                continue
            print_value(
                key[1],
                dict_values[filename][key[0]],
                args.show_comp,
                args.threshold,
                colors,
            )
            del dict_values[filename][key[0]]
        print("")

    print("")

    ref = 0
    if len(files[args.file_lists[args.base - 1]][1]) > 0:
        ref = float(files[args.file_lists[args.base - 1]][1])
    nocol = ansi("0", colors)  # reset
    for filename in args.file_lists:
        if len(files[filename][1]) > 0 and ref != 0:
            comp = (float(files[filename][1]) - ref) / ref
            if abs(comp) > 1e-14:
                color = ansi("91", colors)  # red
            else:
                color = ansi("0", colors)  # reset
            print(
                ("{0} ==> {1} : {2} : " + color + "{3}" + nocol).format(
                    files[filename][0],
                    filename,
                    files[filename][1],
                    (float(files[filename][1]) - ref) / ref,
                )
            )
        else:
            print(("{0} ==> {1} : ").format(files[filename][0], filename) + "-" * 20)

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

    if len(stats_keys) > 0:
        print("===== Stats report =====")
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
        print("")


# ===============================================================================
main()
