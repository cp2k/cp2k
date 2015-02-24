/*****************************************************************************
 *  CP2K: A general program to perform molecular dynamics simulations        *
 *  Copyright (C) 2000 - 2015  CP2K developers group                         *
 *****************************************************************************/

#ifdef __HWLOC
#ifndef MA_COMPONENTS_H
#define MA_COMPONENTS_H

#define ALL 999
#define OS  0
#define LOCAL 1
#define INTERLEAVE 2
#define MANUAL -1

/*
* Describes the machine components
* Used as an interface between Fortran and C
*
* nnodes - number of numa nodes
* nsocktes - total number of sockets
* ncores - total number of cores
* ncaches - total number of caches per socket
* nshared_caches - total number of shared caches between cores
* nsiblings - number of direct sibling cores
* nmemcontroller - number of memory banks per socket
*/
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
* Console output
*/
struct machine_output{
 char console_output[8192];
 int len;
};

/*
* Describes the components of a core
* Internal structure of the machine library
*
* id - physical id of the core
* caches - size of the caches levels for each core
* shared_caches - cache memories levels that are shared between cores
* siblings_id - physical ids of the core neighbors \
* freq - core frequency
*/
struct core{
 unsigned id;
 double freq;  // TODO: get frequence of each core
 int numaNode;
 size_t *caches;
 int *shared_caches;
 int nsiblings;
 unsigned *siblings_id;
};

/*
* Describes the components of a node
* Internal structure of the machine library
*
* id - node id
* ncores - number of core per node
* mycores - physical ids of the node cores
* memory - the amount of memory per node
* nneighbors - number of direct neighbors of the node
* neighbors_id - the ids of the neighbors of the node
*/
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


/******************************************************************************
 HWLOC based functions 
 - to extract the machine topology
 - to map process and threads to the machine processig units
******************************************************************************/

//Get the machine hierarchy
int hw_topology_init (struct arch_topology *topo);
int hw_topology_destroy (struct arch_topology *topo);

void get_phys_id(hwloc_topology_t topology, int ncores,int cur_core);  
void set_phys_siblings(int index, unsigned myid, hwloc_obj_t obj, int ncores, 
                       int nsiblings,int type);
void set_node_cores(hwloc_topology_t topology, hwloc_obj_t obj, int numa_node, int *num_core);

void hw_phys_pu_topology(struct machine_output *ma_out);
void hw_mem_topology (struct machine_output *ma_out);
void hw_pu_topology (struct machine_output *ma_out);
void hw_machine_topology (struct machine_output *ma_out);
void hw_high_level_show(struct machine_output *ma_out);

//Create a string with the machine topology - by object type
void print_machine(hwloc_topology_t topology, hwloc_obj_t obj, int depth);
void print_children_mem(hwloc_topology_t topology, hwloc_obj_t obj, int depth);
void print_children_physical(hwloc_topology_t topology, hwloc_obj_t obj);
void print_children_pu(hwloc_topology_t topology, hwloc_obj_t obj, int depth);
void print_machine_branch(hwloc_topology_t topology, hwloc_obj_t obj, int depth, 
                          int obj_type);

//Process related functions - set/get core, numa node, mempol
int hw_get_proc_node();
int hw_get_proc_core();
void hw_set_proc_node(int node);
void hw_set_proc_core(int core);
int hw_set_mempol(int mempol, int node);
void hw_get_mempol(int *node, int *mempol);
int hw_get_firstcore();

//Thread related functions - set/get core and numa node
int hw_my_node();
int hw_my_core();
unsigned hw_get_myid();
void hw_set_my_core(int cpu);
void hw_set_thread_cores(int distance,int core);

/*******************************************************************************
 LINUX based functions
 - to map process and threads to the machine processing units
 - to set the memory policy for memory allocation
*******************************************************************************/
void linux_set_proc_core(int core);
int linux_proc_core();
void linux_set_my_core(int core);
int linux_my_core();
#endif
#endif
