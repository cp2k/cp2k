/*****************************************************************************
 *  CP2K: A general program to perform molecular dynamics simulations        *
 *  Copyright (C) 2000 - 2015 the CP2K developers group                      *
 *****************************************************************************/
/* Hans Pabst (Intel Corp.)
******************************************************************************/
#ifndef LIBMICSMM_HPP
#define LIBMICSMM_HPP

#if defined(__ACC) && defined(__ACC_MIC) && defined(__DBCSR_ACC)

#include "libmicsmm.h"


template<typename T, bool Complex> struct dbcsr_elem  { static const dbcsr_elem_type type = DBCSR_ELEM_UNKNOWN;
                                                        static const char* name() { return "unknown"; } };
template<> struct dbcsr_elem<float,false>             { static const dbcsr_elem_type type = DBCSR_ELEM_F32;
                                                        static const char* name() { return "f32"; } };
template<> struct dbcsr_elem<double,false>            { static const dbcsr_elem_type type = DBCSR_ELEM_F64;
                                                        static const char* name() { return "f64"; } };
template<> struct dbcsr_elem<float,true>              { static const dbcsr_elem_type type = DBCSR_ELEM_C32;
                                                        static const char* name() { return "c32"; } };
template<> struct dbcsr_elem<double,true>             { static const dbcsr_elem_type type = DBCSR_ELEM_C64;
                                                        static const char* name() { return "c64"; } };

#endif // defined(__ACC) && defined(__ACC_MIC) && defined(__DBCSR_ACC)
#endif // LIBMICSMM_HPP
