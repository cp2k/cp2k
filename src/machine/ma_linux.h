/*****************************************************************************
*   CP2K: A general program to perform molecular dynamics simulations
*   Copyright (C) 2000 * 2011 Christiane Ribeiro and the CP2K developers group
******************************************************************************/
#ifdef __LIBNUMA
#ifndef MA_LINUX_H
#define MA_LINUX_H

#define ALL 999
#define OS  0
#define LOCAL 1
#define INTERLEAVE 2
#define MANUAL -1
/*
 * * Describes the machine components
 * * Used as an interface between Fortran and C
 * *
 * * nnodes - number of numa nodes
 * * nsocktes - total number of sockets
 * * ncores - total number of cores
 * * ncaches - total number of caches per socket
 * * nshared_caches - total number of shared caches between cores
 * * nsiblings - number of direct sibling cores
 * * nmemcontroller - number of memory banks per socket
 * */
struct arch_topology
{
 int nnodes;
 int nsockets;
 int ncores;
 int npus;
 int ngpus;
 int ncaches;
 int nshared_caches;
 int nsiblings;
 int nmemcontroller;
 int nnetcards;
};

/*
 * * Describes the components of a node
 * * Internal structure of the machine library
 * *
 * * id - node id
 * * ncores - number of core per node
 * * mycores - physical ids of the node cores
 * * memory - the amount of memory per node
 * * nneighbors - number of direct neighbors of the node
 * * neighbors_id - the ids of the neighbors of the node
 * */
struct node{
 unsigned id;
 int ncores;
 int ngpus;
 size_t memory;
 unsigned *mycores;
 unsigned *mygpus;
 int nneighbors;
 int *neighbors_id; //TODO: get the node neighbors of each NUMA node
};


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
void linux_set_mempol(int mempol, int node);

#endif
#endif
