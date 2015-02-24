/*****************************************************************************
 *  CP2K: A general program to perform molecular dynamics simulations        *
 *  Copyright (C) 2000 - 2015  CP2K developers group                         *
 *****************************************************************************/

#ifdef __LIBNUMA
#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#include <sys/syscall.h>
#include <sched.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <numa.h>
#include <numaif.h>
#include <dirent.h>

#include "ma_linux.h"
#include "ma_components.h"


static struct node *machine_nodes;
static struct arch_topology *local_topo;


static inline int mask_isset(unsigned long *a, unsigned long sa, int node) {
        if ((unsigned int)node >= sa*8) return 0;
        return bit_is_marked(a,node);
}

int linux_topology_init(struct arch_topology *topo)
{
  int count, i, j, error,k,tmpNode;
  
#ifdef  __DBCSR_ACC
  int nDev;
  ma_get_ndevices_cu(&nDev);
  topo->ngpus = nDev;
#endif

  local_topo = malloc(sizeof(struct arch_topology));

  topo->nnodes = linux_get_nnodes();
  local_topo->nnodes = topo->nnodes;
  topo->ncores = linux_get_ncores();
  local_topo->ncores = topo->ncores;
  topo->npus = topo->ncores;
  local_topo->npus = topo->npus;

  //libnuma has no support for I/O devices
  topo->nnetcards = 0;
  local_topo->nnetcards = 0;
  topo->nsockets = linux_get_nsockets();
  local_topo->nsockets = topo->nsockets;
   //Compute number of memory controlers per socket
   //basically the number of NUMA nodes per socket
  if (topo->nnodes > topo->nsockets)
    topo->nmemcontroller = topo->nnodes/topo->nsockets;
  else
    topo->nmemcontroller = 1;
                  
  topo->ncaches = linux_get_ncaches(); 

   local_topo->nmemcontroller = topo->nmemcontroller;
   local_topo->ncaches = topo->ncaches;

  topo->nshared_caches = linux_get_nshared_caches();
  topo->nsiblings = linux_get_nsiblings();  

  local_topo->nshared_caches =  topo->nshared_caches;
  local_topo->nsiblings = topo->nsiblings;

  //Machine node and core representation
  machine_nodes = (struct node*) malloc (topo->nnodes*sizeof(struct node));

  int ncore_node = topo->ncores/topo->nnodes;

 for (i = 0; i < topo->nnodes ; i++)
   {
        machine_nodes[i].id = i;
        machine_nodes[i].memory = 0;
        machine_nodes[i].ncores = ncore_node;
#ifdef  __DBCSR_ACC
       ma_get_nDevcu(i,&nDev);
       machine_nodes[i].mygpus = malloc (nDev*sizeof(int));
       ma_get_cu(i,machine_nodes[i].mygpus);
#endif
   }

   if (topo->nnodes == -1 || topo->ncores == -1 || topo->npus == -1 ||
       topo->nsockets == -1)
        return -1;
   else
        return 0;

}

int linux_get_nshared_caches()
{
   int i, j, nshared=0;
   char cpu_list[256], temp[12], filename[512];
   FILE *file;

   cpu_list[0] = '\0'; 

   

   for (i=2; i <= local_topo->ncaches ; i++)
   {
      strcpy(filename,"/sys/devices/system/cpu/cpu0/cache/index");
      temp[0]='\0';
      sprintf(temp,"%i", i);
      strcat(filename,temp);
      strcat(filename,"/shared_cpu_list");
      file = fopen(filename,"r");
      if(file != NULL){
        fscanf(file,"%s",cpu_list);
        strcat(cpu_list,"\0");
        for (j = 0; j < strlen(cpu_list); j++)
          if((cpu_list[j] == '-') || (cpu_list[j] == ','))
            nshared++;
      } 
   }
   
  return nshared;
}

/*
TODO: get thread id - portable
*/
unsigned linux_get_myid()
{
  return syscall(SYS_gettid);
}

/*
* Set the core where the current process will run
* return the core
*/
void linux_set_proc_core(int core)
{
  cpu_set_t set;

  CPU_ZERO(&set);
  CPU_SET(core, &set);
  if(sched_setaffinity(0,sizeof(set),&set)!=0)
     printf("\n Error binding threads");
}

/*
* Get the core where the current process is running
* return the core
*/
int linux_proc_core()
{
  int core;
  core = sched_getcpu();
  return core;
}

/*
* Set the core where the current thread will run
* return the core
*/
void linux_set_my_core(int core)
{
  cpu_set_t set;

  CPU_ZERO(&set);
  CPU_SET(core, &set);

  if(sched_setaffinity(syscall(SYS_gettid),sizeof(set),&set)!=0)
     printf("\n Error binding threads");
}

/*
* Set the core where the current thread will run
* return the core
*/
void linux_set_proc_cores(int distance, int core)
{
  cpu_set_t set;
  int i;

  CPU_ZERO(&set);
  for(i=core;i<core+distance;i++)
       CPU_SET(i,&set);
  
  if(sched_setaffinity(syscall(SYS_gettid),sizeof(set),&set)!=0)
     printf("\n Error binding threads");
}

/*
* Get the core where the current thread is running
* return the core
*/
int linux_my_core()
{
  int core;
  core = sched_getcpu();
  return core;
}

int linux_is_numa()
{
 return numa_available();
}

int linux_get_nsiblings()
{
        int             i,socket_id=0,nsiblings=0;
        char            dirnamep[256];
        char            temp[33];
        struct dirent   *dirent;
        FILE             *file;
        int maxcores = 0;

        temp[0]=dirnamep[0]='\0';

        do{
        strcpy(dirnamep,"/sys/devices/system/cpu/cpu");
        temp[0]='\0';
        sprintf(temp, "%i", maxcores);
        strcat(dirnamep,temp);
        strcat(dirnamep,"/topology/physical_package_id");

        file = fopen(dirnamep,"r");
           if(file != NULL){
             fscanf(file,"%d",&socket_id);
             nsiblings++;
          }
          maxcores++;
        } while((socket_id == 0) && (file != NULL));

        fclose(file);

 return nsiblings-1;

}

int linux_get_nsockets()
{
        int             i,socket_id=0, cur_socket_id=0,nsockets=0;
        char            dirnamep[256];
        char            temp[33];
        struct dirent   *dirent;
        FILE             *file;
        int maxcores = local_topo->ncores,cores=0;
       

        temp[0]=dirnamep[0]='\0';

        do{
        strcpy(dirnamep,"/sys/devices/system/cpu/cpu");
        temp[0]='\0';
        sprintf(temp, "%i", cores);
        strcat(dirnamep,temp);
        strcat(dirnamep,"/topology/physical_package_id");
        file = fopen(dirnamep,"r");
           if(file != NULL){
             fscanf(file,"%d",&cur_socket_id);
            if(cur_socket_id != socket_id)  
             {   nsockets++; socket_id = cur_socket_id;}
          }
          cores++;
        } while(cores <  maxcores);

        fclose(file);

 return nsockets+1;

}

int linux_get_ncaches()
{
        int i;
        int             node=-1;
        char            dirnamep[256];
        char            temp[33];
        struct dirent   *dirent;
        DIR             *dir;
        int maxcaches = 0;

        temp[0]=dirnamep[0]='\0';

        do{
        strcpy(dirnamep,"/sys/devices/system/cpu/cpu0/cache/index");
        temp[0]='\0';
        sprintf(temp, "%i", maxcaches);
        strcat(dirnamep,temp);
        dir = opendir(dirnamep);
           if(dir != NULL)
             maxcaches++;
        } while(dir != NULL);

        closedir(dir);

 return maxcaches-1;

}


/*get inumber of nodes
return 0 if the machine is not NUMA
*/
int linux_get_nnodes()
{
        int i;
        int             node=-1;
        char            dirnamep[256];
        char            temp[33];
        struct dirent   *dirent;
        DIR             *dir;
        int maxnodes = 0;

        temp[0]=dirnamep[0]='\0';

        do{
        strcpy(dirnamep,"/sys/devices/system/node/node");
        temp[0]='\0';
        sprintf(temp, "%i", maxnodes);
        strcat(dirnamep,temp);
        dir = opendir(dirnamep);
           if(dir != NULL)
             maxnodes++;
        } while(dir != NULL);

        closedir(dir);

 return maxnodes;

}

/*get inumber of cores
*/
int linux_get_ncores()
{
        int i;
        int             node=-1;
        char            dirnamep[256];
        char            temp[33];
        struct dirent   *dirent;
        DIR             *dir;
        int maxcores = 0;

        temp[0]=dirnamep[0]='\0';

        do{
        strcpy(dirnamep,"/sys/devices/system/cpu/cpu");
        temp[0]='\0';
        sprintf(temp, "%i", maxcores);
        strcat(dirnamep,temp);
        dir = opendir(dirnamep);
           if(dir != NULL)
             maxcores++;
        } while(dir != NULL);

        closedir(dir);

  return maxcores;
}


/*get node id of a cpu/core
return: 0 if UMA 
*/
int linux_get_nodeid()
{
   	int             cpu,i;
        int             node=-1;
        char            dirnamep[256];
        char            cpuid[7],temp[33];
        struct dirent   *dirent;
        DIR             *dir;
        int maxnodes;
        maxnodes  =  local_topo->nnodes;
        cpu = sched_getcpu();
      
      if( maxnodes > 1){  
	  temp[0]=cpuid[0]=dirnamep[0]='\0';
          strcpy(cpuid,"cpu");
          sprintf(temp, "%i", cpu);
          strcat(cpuid,temp);
          strcat(cpuid,"\0");


        for(i=0;i<maxnodes && node == -1;i++){
        strcpy(dirnamep,"/sys/devices/system/node/node");
        temp[0]='\0';
        sprintf(temp, "%i", i);
        strcat(dirnamep,temp);
        dir = opendir(dirnamep);
        if (dir == NULL) {
                return 0;
        }
        while ((dirent = readdir(dir)) != 0) {
                if (strcmp(cpuid, dirent->d_name)==0) {
                        node = i;
                        break;
                }
        }
        closedir(dir);
     }
    }
    else
      node = 0;

  return node;
}

int linux_get_nodeid_cpu(int cpu)
{
        int             i;
        int             node=-1;
        char            dirnamep[256];
        char            cpuid[8],temp[33];
        struct dirent   *dirent=NULL;
        DIR             *dir;
        int maxnodes;
        maxnodes  =  local_topo->nnodes;

      if( maxnodes > 1){
          temp[0]=cpuid[0]=dirnamep[0]='\0';
          strcpy(cpuid,"cpu");
          sprintf(temp, "%i", cpu);
          strcat(cpuid,temp);
          strcat(cpuid,"\0");


        for(i=0;i<maxnodes && node == -1;i++){
        strcpy(dirnamep,"/sys/devices/system/node/node");
        temp[0]='\0';
        sprintf(temp, "%i", i);
        strcat(temp,"\0");
        strcat(dirnamep,temp);
        dir = opendir(dirnamep);
        if (dir == NULL) {
               return 0;
        }
        dirent = readdir(dir);
        while (dirent != NULL) {
          if (strcmp(cpuid, dirent->d_name)==0) {
                        node = i;                   
           }
         dirent = readdir(dir);
        }
     }
       closedir(dir);
    }
    else
      node = 0;

  return node;
}

void linux_set_proc_node(int node)
{
   cpu_set_t set;
   int i;

   CPU_ZERO(&set);

   for(i=0;i<local_topo->ncores;i++)
    if(linux_get_nodeid_cpu(i) == node)
        CPU_SET(i,&set);

   if(sched_setaffinity(syscall(SYS_gettid),sizeof(set),&set)!=0)
         printf("\n Error binding threads");

}

/*
* Set the memory policy for data allocation for a process
*/
void linux_set_mempol(int mempol, int realnode)
{

  int node = local_topo->nnodes;
  int core, i, error;
  unsigned long mask;

  switch(mempol)
  {
   case OS:
    if(node > 0)
      error = set_mempolicy(MPOL_DEFAULT,NULL,0); 
   break;
   case LOCAL:   
    if(node >0)
    {  
     core = sched_getcpu();
     realnode = linux_get_nodeid(core);
     mask = 1 << realnode;
     error = set_mempolicy(MPOL_BIND,&mask,sizeof(unsigned long)+1);
    }
   break;
   case INTERLEAVE:
    if(node >0){
      for(i=0;i<node+1;i++)
        mask = mask + (1 << i);
        error = set_mempolicy(MPOL_INTERLEAVE,&mask,sizeof(mask)+1);      
    }
   break;
   case MANUAL:   
    if(node >0)
    {  
     mask = 1 << realnode;
     error = set_mempolicy(MPOL_BIND,&mask,sizeof(unsigned long)+1);
    }
   break;
   default:
    if(node > 0)
      error = set_mempolicy(MPOL_DEFAULT,NULL,0); 
   break;
  }

  if (error < 0)
    printf("\nWARNING: Memory binding not supported or can not be enforced\n");

}

/*
* Get the memory policy of a process
*/
void linux_get_mempol(int *node, int *mem_pol)
{
  unsigned long nset;
  unsigned long addr;
  int mempol=-1,nnodes,i, error;
    
  nnodes = local_topo->nnodes;

 if (nnodes != 0 ){
  (*node) = 0;

  error = get_mempolicy (&mempol,NULL,0,0,0); 

  switch(mempol)
  {
   case MPOL_DEFAULT:
     (*mem_pol) = OS;
   break;
   case MPOL_BIND:
     (*mem_pol) = LOCAL;
     (*node)   =  linux_get_nodeid();
   break;
   case MPOL_INTERLEAVE:
     (*mem_pol) = INTERLEAVE;
     (*node)   = -1;
   break;
   default:
     (*mem_pol) = -1;
     (*node) = -1;
   break;
  }
 }
 else
  (*mem_pol) = -1;

 if (error < 0)
   printf("\nWARNING: Can not retrieve memory binding.\n");

}



#ifdef  __DBCSR_ACC

/*
 *Build the GPUs list for a core
 *param the core Id
 *return the list of gpus
 */
void linux_my_gpuList(int nodeId,int *devList)
{

  if( local_topo->nnodes == 0)
    ma_get_uma_closerDev_cu(devList);
  else
    ma_get_closerDev_cu(nodeId,devList);
  
}

/*
 *Get the GPU ID for a MPI process
 *param core id where the process is running
 *      rank MPI rank number
 *return the GPU ID to connect to           
 * */
int linux_my_gpu(int coreId, int myRank, int nMPIs)
{
  int *devList,i,nDev,nodeId;

  int ngpus = local_topo->ngpus;
  devList = malloc(ngpus*sizeof(int));

  nodeId = linux_get_nodeid_cpu(coreId);


  linux_my_gpuList(nodeId,devList);

/*
 *The algorithm: 
 * if UMA machine - similar costs to access devices
 *   just make a round-robin distribution of the GPUs
 * if NUMA machine - not similar costs to access devices
 *   if number of MPI process equal to number of GPUs
 *     just give one GPU for each MPI process
 *   if the MPI is in a NUMA node that has no GPU
 *     just make a round-robin distribution of the GPUs
 *   if number of MPIs is larger than number of GPUs 
 *      try to balance the GPU usage    	  
 * */

  if( local_topo->nnodes == 0 )
   return devList[myRank%ngpus];
  else {
   ma_get_nDevcu(nodeId,&nDev);
   if (nDev == 0)
    return devList[myRank%ngpus];
   else  
    return devList[myRank%nDev];  
 }

}

int linux_get_gpu_node(gpu)
{
  int node=0;
  ma_get_NUMAnode_cu(gpu, &node);
  return node;
}
#endif


#endif
