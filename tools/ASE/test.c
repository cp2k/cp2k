#include <stdio.h>
#include "cp2k_c_bridge.h"

int main(int argc,char *argv[])
{
  int env_id,i,j;
  static const int nel=9;
  f_real pos[nel];
  f_real force[nel];
  f_real e_pot;
  if (argc<3) {
    printf("you need to specify input and output\n");
    return 1;
  }

  for (j=0;j<nel/3;++j)
    for (i=0;i<3;++i) pos[i]=(f_real)j;
  for (i=0;i<nel;++i) force[i]=0.;
  printf("1\n");
  if (ccp_init_cp2k(1)) {
    printf("init error\n");
    return 1;
  }
  printf("2\n");
  if (ccp_create_fenv(&env_id,argv[1],argv[2])) {
    printf("create_fenv error\n");
    return 1;
  }
  printf("3\n");
  if (ccp_set_pos(env_id,pos,nel)) {
    printf("set pos error\n");
    return 1;
  }
  printf("4\n");
  if (ccp_calc_energy_force(env_id,1)) {
    printf("calc e_f error\n");
    return 1;
  }
  printf("5\n");
  if (ccp_get_energy(env_id,&e_pot)) {
    printf("get_energy error\n");
    return 1;
  }
  printf("6\n");
  printf("energy=%f\n",e_pot);
  if (ccp_get_force(env_id,force,nel)) {
    printf("get_force error\n");
    return 1;
  }
  printf("7\n");
  printf("force=");
  for (i=0;i<nel;i++) printf(" %f",force[i]);
  printf("\n");
  if (ccp_finalize_cp2k(1)) {
    printf("finalize error\n");
    return 1;
  }
  printf("*** finished! ***\n");
  return 0;
}
