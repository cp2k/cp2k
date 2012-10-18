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
  int pid;
  char name[24], name_dir[256];
  char pidstr[12];

  name[0] = pidstr[0] = '\0';

  pid = getpid();
  sprintf(pidstr,"%d",pid);  

  strcpy(name_dir,"topo_");
  strcat(name_dir,pidstr);
  strcat(name_dir,"\0");

  char mk_dir[256];
  strcpy(mk_dir,"mkdir topo_");
  strcat(mk_dir,pidstr);
  strcat(mk_dir,"\0");
  int ret = system(mk_dir);

  strcpy(name,"topo_");
  strcat(name,pidstr);
  strcat(name,"/topo_cray.sh\0");

  fpexec = fopen(name, "w");
  fprintf (fpexec, "#!/bin/bash\n");
  fprintf (fpexec, "xtdb2proc\n");
  fprintf (fpexec, "tr ',' ' ' < processor > processor_tmp\n");
  fprintf (fpexec,"tr '=' ' ' < processor_tmp > processor\n");
  fprintf (fpexec,"sed -e '/#/d' processor > processor_tmp\n");
  fprintf (fpexec,"cat processor_tmp | awk '{print $30\"\t\"$14\"\t\"$16\"\t\"$18}' > processor\n");
  fprintf (fpexec,"rm processor_tmp\n");
  fprintf (fpexec,"N=`cat processor | awk 'BEGIN{max=0}{if($1>max)max=$1}END{print max}'`\n");
  fprintf (fpexec,"echo $N > topology\n");
  fprintf (fpexec,"cat processor >> topology\n");
  fprintf (fpexec,"rm processor\n");
  fclose(fpexec);

  char command[256];
   
  strcpy(command,"chmod +x ./");
  strcat(command,name);
  strcat(command,";cd ");
  strcat(command,name_dir);
  strcat(command,"; ./topo_cray.sh");
  strcat(command,"; rm topo_cray.sh");
  strcat(command,";cd ../");


  ret = system(command);
  return ret;
}

int remove_topology()
{
  char pidstr[12];
  char name_dir[256];
  int pid;

  pid = getpid();
  sprintf(pidstr,"%d",pid);  

  strcpy(name_dir,"topo_");
  strcat(name_dir,pidstr);

  char command[256];
  
  strcpy(command,"rm -rf ");
  strcat(command,name_dir);
  strcat(command,";");

  int ret = system(command);
  return ret;
}

int get_job_alloc_strategy()
{
  int strategy = -1;
  char nid_order;
  FILE *fpexec;
 
  char command[]="less /etc/sysconfig/alps | grep ALPS_NIDORDER | cut -c 18";
  fpexec = popen(command,"r");
  if (fpexec != NULL){
     fscanf(fpexec,"%c",&nid_order);
     if (nid_order == 'n' || nid_order == 'x' || nid_order == 'y' ||
         nid_order == 'r' || nid_order == 'r')
          strategy = 1;                    
  }  

  return strategy;
}

#endif


