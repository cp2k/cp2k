#include <math.h>
#include "cp2k_shell.h"
#include <stdio.h>
#define myErrCheck(cond,msg) if (cond) { printf(msg);printf("\n%s:%d\n",__FILE__,__LINE__);exit(2);}

/**
 * \brief very simple testing of cp2k_shell interface
 *
 * \author fawzi
 *
 * Copyright (C) 2008, 2009  CP2K developers group
 */

int main (int argc, const char * argv[]) {
	SubT st;
	int env_id,i,natoms,ninfo;
	double *pos,*f,e1,e2,max_err;
	const char *cp2k_shellEXE="/Volumes/Users/fawzi/cp2k/exe/Darwin-IntelMacintosh-g95/cp2k_shell.sopt",
		*cp2k_test_dir="/Volumes/Users/fawzi/cp2k/tests/QS/regtest-gpw-1",
		*cp2k_input_file="Ar-2.inp";
	char *info;
	if (argc>1) cp2k_shellEXE=argv[1];
	if (argc>2) cp2k_test_dir=argv[2];
	if (argc>3) cp2k_input_file=argv[3];
    fm_spawnSubtask(cp2k_shellEXE,
		cp2k_test_dir,&st);
	printf("t1\n"); 
	myErrCheck(fm_getready(&st),"task did not get ready, start error?");
	printf("t2\n"); 
	myErrCheck(fm_load(&st,cp2k_input_file,&env_id),"error loading input");
	printf("loaded %i\n",env_id);
	myErrCheck(fm_natom(&st,env_id,&natoms),"error getting natoms");
	f=malloc(sizeof(double)*3*natoms);
	pos=malloc(sizeof(double)*3*natoms);	
	myErrCheck(fm_getpos(&st,env_id,3*natoms,pos),"error getting pos");
	printf("getPos\n");
	myErrCheck(fm_calc_e(&st,env_id,&e1),"error calculating e");
	printf("calcE %g\n",e1);
	pos[1]=0.2;
	pos[2]=0.4321;
	myErrCheck(fm_setpos(&st,env_id,3*natoms,pos,&max_err),"error setting pos");
	printf("setPos %g\n",max_err);
	myErrCheck(fm_calc_ef(&st,env_id,&e2,3*natoms,f),"error calculating ef");
	printf("calcEF %g\n",e2);
	printf("Translation non invariance: %g, mean energy: %g,%g,%g\n",fabs(e1-e2),0.5*(e1+e2),e1,e2);
	printf("force:");
	for (i=0;i<3*natoms;i++) {
		printf(" %24.13g",f[i]);
		if (i%3==2) printf("\n");
	}
// async test
	pos[1]=0.0;
	pos[2]=0.0;
	myErrCheck(fm_setpos(&st,env_id,3*natoms,pos,&max_err),"error setting pos");
	printf("setPos %g\n",max_err);
	myErrCheck(fm_eval_e(&st,env_id),"error evaluating e");
	printf("evalE\n");
	myErrCheck(fm_get_e(&st,env_id,&e1),"error getting e");
	printf("getE %g\n",e1);
	pos[1]=0.2;
	pos[2]=0.4321;
	myErrCheck(fm_setpos(&st,env_id,3*natoms,pos,&max_err),"error setting pos");
	printf("setPos %g\n",max_err);
	myErrCheck(fm_eval_ef(&st,env_id),"error evaluating ef");
	printf("evalEF\n");
	myErrCheck(fm_get_e(&st,env_id,&e2),"error getting e");
	printf("getE %g\n",e2);
	myErrCheck(fm_get_f(&st,env_id,3*natoms,f),"error getting f");
	printf("getF\n");
	printf("Translation non invariance: %g, mean energy: %g,%g,%g\n",fabs(e1-e2),0.5*(e1+e2),e1,e2);
	printf("force:");
	for (i=0;i<3*natoms;i++) {
		printf(" %24.13g",f[i]);
		if (i%3==2) printf("\n");
	}
	// bg load
	myErrCheck(fm_bg_load(&st,cp2k_input_file),"error loading input");
	printf("bg_load\n");
	myErrCheck(fm_last_env_id(&st,&env_id),"error getting env id");
	printf("last_env_id %d\n",env_id);
	myErrCheck(fm_calc_e(&st,env_id,&e1),"error calculating e");
	printf("calcE %g\n",e1);
	// info
	info=NULL;
	myErrCheck(fm_get_info(&st,0,info,&ninfo),"error getting info length");
	info=malloc(ninfo*sizeof(char));
	myErrCheck(fm_get_info(&st,(ninfo<100?ninfo:100),info,&ninfo),"error getting short info");
	printf("got info start %s\n",info);
	myErrCheck(fm_get_info(&st,ninfo,info,&ninfo),"error getting long info");
	printf("got info \n%s\n",info);
	free(info);
	printf("\n*** SUCCESS ***\n");
	
	free(f);
	free(pos);
	fclose(st.outF);
	fclose(st.inF);
    return 0;
}
