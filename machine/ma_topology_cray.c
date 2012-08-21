/*****************************************************************************
*   CP2K: A general program to perform molecular dynamics simulations         
*   Copyright (C) 2000 * 2011 Christiane Ribeiro and the CP2K developers group
******************************************************************************/

#include "ma_topology_cray.h"  

//Supported CRAY architectures: XT4,XT5,XT6, XE6 and XK6
#if defined __GEMINI || __SEASTAR

//Extract the machine architecture from the CRAY database
//Shoudl be called by only one MPI process
//Return a topology file
int extract_topology()
{
  //get information form database - not the best way to do it 
  // but it is the simplest one 
  FILE *fpexec;

  fpexec = fopen("topo_cray.sh", "w");
  fprintf (fpexec, "#!/bin/bash\n");
  fprintf (fpexec, "xtdb2proc\n");
  fprintf (fpexec, "tr ',' ' ' < processor > processor_tmp\n");
  fprintf (fpexec,"tr '=' ' ' < processor_tmp > processor\n");
  fprintf (fpexec,"sed -e '/#/d' processor > processor_tmp\n");
  fprintf (fpexec,"cat processor_tmp | awk '{print $30\"\t\"$14\"\t\"$16\"\t\"$18}' > processor\n");
  fprintf (fpexec,"rm processor_tmp\n");
  fprintf (fpexec,"N=`cat processor |wc -l`\n");
  fprintf (fpexec,"echo $N > topology\n");
  fprintf (fpexec,"cat processor >> topology\n");
  fprintf (fpexec,"rm processor\n");
  fclose(fpexec);

  char command[]="chmod +x topo_cray.sh; ./topo_cray.sh; rm topo_cray.sh";
  int ret = system(command);
  return ret;
}

int remove_topology()
{
  char command[]="if [ -f topology ]; then rm topology; fi";
  int ret = system(command);
  return ret;
}

#endif


