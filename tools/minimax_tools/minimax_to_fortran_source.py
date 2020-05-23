#!/usr/bin/env python3

import re
import sys
import os.path

# 1) parse data #
#################

reldir = r"www.mis.mpg.de/scicomp/EXP_SUM/1_x/"

txt = open(reldir + "tabelle", "r")

re_header = re.compile("k =  \| *(?P<kval>([\d]+ *)+)")
re_body = re.compile("(?P<Rc>\dE\d\d) \|(?P<err>[\w\.\- ]+)")
re_sep = re.compile("[\-]{50,100}")

read_header = True
read_sep1 = False
read_body = False
read_sep2 = False

k = []
Rc = []
err = []
coeff_file = []
missing = []
a = []
w = []

start_line = 46
end_line = 288

for l, line in enumerate(txt.readlines()):
    if l < start_line - 1:
        continue
    if l > end_line - 1:
        break

    if read_header:
        if re_header.match(line):
            kvals = [int(_) for _ in re_header.match(line).group("kval").split()]
            read_sep1 = True
            read_header = False
            continue

    if read_sep1:
        if re_sep.match(line):
            read_body = True
            read_sep1 = False
            continue

    if read_body:
        if re_body.match(line):
            Rcval = float(re_body.match(line).group("Rc"))
            Rcstr = re_body.match(line).group("Rc")
            # slight change of notation e.g. 6E09 --> 6E9
            Rcstr = re.sub("E0(?=\d)", "E", Rcstr)
            errvals = [
                float(_) if _ != "--" else 0.0
                for _ in re_body.match(line).group("err").split()
            ]

            for i, errval in enumerate(errvals):
                if errval > 0.0:
                    # hack: slight change of notation e.g. 1 --> 01
                    kstr = str(kvals[i])
                    if len(kstr) == 1:
                        kstr = "0" + kstr
                    filename = reldir + "1_xk" + kstr + "_" + Rcstr
                    if not os.path.isfile(filename):
                        missing.append(filename)
                        continue

                    # read parameters
                    txt2 = open(filename, "r")
                    ww = []
                    aa = []
                    for l2, line2 in enumerate(txt2.readlines()):
                        if l2 < kvals[i]:
                            ww.append(float(line2.split()[0]))
                        else:
                            aa.append(float(line2.split()[0]))
                    a.append(aa)
                    w.append(ww)
                    coeff_file.append(filename)
                    Rc.append(Rcval)
                    err.append(errval)
                    k.append(kvals[i])

            continue
        else:
            read_body = False
            read_sep2 = True

    if read_sep2:
        if re_sep.match(line):
            read_body = False
            read_header = True
            read_sep2 = False
            continue

# 2) sort all data w.r.t. 1) Rc, 2) k #
#######################################

data_2_sort = zip(k, Rc, err, a, w)
data_sorted = sorted(data_2_sort)

k = [_[0] for _ in data_sorted]
Rc = [_[1] for _ in data_sorted]
err = [_[2] for _ in data_sorted]
a = [_[3] for _ in data_sorted]
w = [_[4] for _ in data_sorted]

print("missing files")
print(missing)

txt.close()

# pointers to k
k_p = []
my_k = 0
for i, kk in enumerate(k):
    if my_k < kk:
        my_k = kk
        k_p.append(i + 1)

k_p.append(len(k) + 1)

# 3) generate fortran file #
############################

out = open("../../src/minimax/minimax_exp_k53.F", "w")
out.write(
    "!--------------------------------------------------------------------------------------------------!\n\
!   CP2K: A general program to perform molecular dynamics simulations                              !\n\
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !\n\
!--------------------------------------------------------------------------------------------------!\n\n"
)

out.write(
    "! **************************************************************************************************\n\
!> \\brief Routines to calculate the minimax coefficients in order to\n\
!>        approximate 1/x as a sum over exponential functions\n\
!>        1/x ~ SUM_{i}^{K} w_i EXP(-a_i * x) for x belonging to [1:Rc].\n\
!>        This module contains coefficients for minimax approximations with 1 <= k <= 53.\n\
!>        Generated from data from\n\
!>        http://www.mis.mpg.de/scicomp/EXP_SUM/1_x\n\
!>        This module should not be modified manually and should not be used anywhere\n\
!>        except in main minimax module.\n\
!>        This file was created using the scripts in cp2k/tools/minimax_tools.\n\
! **************************************************************************************************\n\n"
)

out.write("MODULE minimax_exp_k53\n")
out.write("USE kinds, ONLY: dp\n")

out.write("IMPLICIT NONE\n")
out.write("PRIVATE\n")
out.write(
    "PUBLIC :: R_max, R_mm, err_mm, get_minimax_coeff_low, k_max, k_mm, k_p, n_approx\n\n"
)

out.write("INTEGER, PARAMETER :: n_approx = {}\n".format(len(k)))
out.write("INTEGER, PARAMETER :: n_k = {}\n".format(len(k_p) - 1))
out.write("INTEGER, PARAMETER :: k_max = {}\n".format((max(k))))
out.write("REAL(KIND=dp), PARAMETER :: R_max = {:.1E}_dp\n\n".format((max(Rc))))

out.write("INTEGER, PARAMETER, DIMENSION(n_k+1) :: k_p = &\n[ ")
for i, kkp in enumerate(k_p):
    if i % 5 == 0 and i > 0:
        out.write("&\n  ")
    out.write("{:3}".format(kkp))
    if not i + 1 == len(k_p):
        out.write(",")
    out.write(" ")
out.write("]\n\n")

out.write("INTEGER, PARAMETER, DIMENSION(n_approx) :: k_mm = &\n[ ")
for i, kk in enumerate(k):
    if i % 5 == 0 and i > 0:
        out.write("&\n  ")
    out.write("{:2}".format(kk))
    if not i + 1 == len(k):
        out.write(",")
    out.write(" ")
out.write("]\n\n")

out.write("REAL(KIND=dp), PARAMETER, DIMENSION(n_approx) :: R_mm = &\n[ ")
for i, RR in enumerate(Rc):
    if i % 5 == 0 and i > 0:
        out.write("&\n  ")
    out.write("{:.1E}_dp".format(RR))
    if not i + 1 == len(Rc):
        out.write(",")
    out.write(" ")
out.write("]\n\n")

out.write("REAL(KIND=dp), PARAMETER, DIMENSION(n_approx) :: err_mm = &\n[ ")
for i, EE in enumerate(err):
    if i % 5 == 0 and i > 0:
        out.write("&\n  ")
    out.write("{:.3E}_dp".format(EE))
    if not i + 1 == len(err):
        out.write(",")
    out.write(" ")
out.write("]\n\n")

out.write("CONTAINS\n\n")
out.write("SUBROUTINE get_minimax_coeff_low(i, aw)\n")
out.write("INTEGER, INTENT(IN) :: i\n")
out.write("REAL(KIND=dp), DIMENSION(k_mm(i)*2), INTENT(OUT) :: aw\n\n")

out.write("SELECT CASE(i)\n\n")
for i, (kk, RR, EE, CC) in enumerate(zip(k, Rc, err, coeff_file)):

    out.write("CASE({})\n\n".format(i + 1))
    # out.write(''.format(i))
    out.write("  aw(:) = & ! a\n[ ")

    offset = len(a[i])
    for j, aa in enumerate(a[i]):
        if j % 3 == 0 and j > 0:
            out.write("&\n  ")
        out.write("{}_dp, ".format(repr(aa)))
    out.write("&\n            ! w\n  ")
    for j, ww in enumerate(w[i]):
        if j % 3 == 0 and j > 0:
            out.write("&\n  ")
        out.write("{}_dp".format(repr(ww)))
        if not j + 1 == offset:
            out.write(",")
        out.write(" ")
    out.write("]\n\n")

out.write("END SELECT\n\n")
out.write("END SUBROUTINE get_minimax_coeff_low\n")
out.write("END MODULE minimax_exp_k53\n")
out.close()
