#include <math.h>
#include "cp2k_shell.h"
#include <stdio.h>

/* very simple testing of cp2k_shell interface
 *
 * fawzi
 */

int main (int argc, const char * argv[]) {
	SubT st;
	int env_id,i;
	double pos[3],f[3],e1,e2;
    fm_spawnSubtask("/Volumes/Users/fawzi/cp2k/exe/Darwin-IntelMacintosh-g95/cp2k_shell.sdbg",
		"/Volumes/Users/fawzi/cp2k/tests/QS/regtest-gpw-1",&st);
	printf("t1\n"); 
	if (fm_getready(st)) {
		printf("task did not get ready, start error?");
		exit(1);
	}
	printf("t2\n"); 
	if (fm_load(st,"Ar-2.inp",&env_id)){
		printf("error loading input");
		exit(2);
	}
	printf("loaded %i\n",env_id);
	if (fm_getpos(st,env_id,3,pos)){
			printf("error getting pos");
		exit(2);
	}
	printf("getPos\n");
	if (fm_calc_e(st,env_id,&e1)){
			printf("error calculating e");
		exit(2);
	}
	printf("calcE %g\n",e1);
	pos[1]=0.2;
	pos[2]=0.4321;
	if (fm_setpos(st,env_id,3,pos)){
			printf("error setting pos");
		exit(2);
	}
	printf("setPos\n");
	if (fm_calc_ef(st,env_id,&e2,3,f)){
			printf("error calculating ef");
		exit(2);
	}
	printf("calcEF %g\n",e2);
	fclose(st.outF);
	fclose(st.inF);
	printf("Translation non invariance: %g, mean energy: %g,%g,%g\n",fabs(e1-e2),0.5*(e1+e2),e1,e2);
	printf("force:");
	for (i=0;i<3;i++) printf(" %24.13g",f[i]);
	printf("\n");
    return 0;
}
