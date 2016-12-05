#!/usr/bin/env python
# -*- coding: utf-8 -*-

#: maximum tensor rank
maxrank = 4

#: datatypes
dtype_float_prec = ['real_8', 'real_4', 'real_8', 'real_4']
dtype_float_type = ['REAL(kind=real_8)', 'REAL(kind=real_4)', 'COMPLEX(kind=real_8)', 'COMPLEX(kind=real_4)']
dtype_float_suffix = ['r_dp', 'r_sp', 'c_dp', 'c_sp']
dtype_float_param = ['dbcsr_type_real_8', 'dbcsr_type_real_4', 'dbcsr_type_complex_8', 'dbcsr_type_complex_4']

dtype_int_type = ['INTEGER']
dtype_int_suffix = ['i']
dtype_int_param = ['dbcsr_type_int_4']

dtype_all_type = dtype_float_type + dtype_int_type
dtype_all_suffix = dtype_float_suffix + dtype_int_suffix
dtype_all_param = dtype_float_param + dtype_int_param

dtype_float_list = list(zip(dtype_float_param, dtype_float_type, dtype_float_suffix))
dtype_int_list = list(zip(dtype_int_param, dtype_int_type, dtype_int_suffix))
dtype_all_list = list(zip(dtype_all_param, dtype_all_type, dtype_all_suffix))

def arrlist(name, n):
    """ expand array into list of elements "name(1), name(2), ..., name(n)"
    """
    return ", ".join([name + "(" + str(i) +")" for i in range(1, n+1)])

def varlist(name, n):
    """ create variable list "name_1, name_2, ..., name_n"
    """
    return ", ".join([name + "_" + str(i) for i in range(1, n+1)])

def shape_colon(n):
    """ repeated colon ':' for e.g. assumed shape array notation
    """
    return ','.join([':']*n)

def uselist(list_in):
    """ comma-separated list of unique entries of list_in.
    """
    return ", ".join(list(set(list_in)))
