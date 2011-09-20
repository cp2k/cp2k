/*****************************************************************************
*   CP2K: A general program to perform molecular dynamics simulations
*   Copyright (C) 2000 * 2011 Christiane Ribeiro and the CP2K developers group
******************************************************************************/
#ifdef __LIBNUMA
#ifndef __MA_LINUX_H
#define __MA_LINUX_H

#define ALL 999
#define OS  0
#define LOCAL 1
#define INTERLEAVE 2


/*******************************************************************************
 LINUX based functions
 - to map process and threads to the machine processing units
 - to set the memory policy for memory allocation
*******************************************************************************/
void linux_set_proc_core(int core);
int linux_proc_core();

void linux_set_my_core(int core);
int linux_my_core();

int linux_get_nodeid();
void linux_set_mempol(int mempol);

#endif
#endif
