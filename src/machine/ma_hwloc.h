/*****************************************************************************
 *  CP2K: A general program to perform molecular dynamics simulations        *
 *  Copyright (C) 2000 - 2015  CP2K developers group                         *
 *****************************************************************************/

//Describes the machine topology
/*
* nnodes - number of numa nodes
* nsocktes - total number of sockets
* ncores - total number of cores
* ncaches - total number of caches per socket
* nshared_caches - total number of shared caches between cores
* nsiblings - number of direct sibling cores
*/
#ifdef __HWLOC

struct arch_topology
{
 int nnodes;
 int nsockets;
 int ncores;
 int npus;
 int ncaches;
 int nshared_caches;
 int nsiblings;
 int nmemcontroller;
 int len_machine;
};

struct core{
 unsigned id;
 size_t *caches;
 int *shared_caches;
 int nsiblings;
 unsigned *siblings_id;
 double freq;  // TODO: get frequence of each core
};

struct node{
 unsigned id;
 int ncores;
 unsigned *mycores;
 size_t memory;
 int nneighbors;
 int *neighbors_id; //TODO: get the node neighbors of each NUMA node
};

void set_phys_siblings(int index, unsigned myid, hwloc_obj_t obj, int ncores, int nsiblings, int type);
#endif
