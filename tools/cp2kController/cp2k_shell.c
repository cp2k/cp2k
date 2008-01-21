/*
 *  cp2k_shell.c
 *  cp2kController
 *
 *  Created by Fawzi Mohamed on 1/20/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

#include "cp2k_shell.h"
#include <unistd.h>
#include <stdlib.h>

/* executes (execvp) task with arguments taskArgs in the directory taskDir as subtask.
 * returns file descriptors taskIn, taskOut to std in and std out of the task
 */
void fm_spawnSubtask(const char *task,const char *taskDir, SubT *st)
{
        int     fdIn[2],fdOut[2];
        pid_t   childpid;

        pipe(fdIn);
        pipe(fdOut);
        
        if((childpid = fork()) == -1)
        {
                perror("fork");
                exit(1);
        }

        if(childpid == 0)
        {
                /* Child process */
				close(fdOut[0]);
                close(fdIn[1]);
				chdir(taskDir);
                dup2(fdIn[0],STDIN_FILENO);
                dup2(fdOut[1],STDOUT_FILENO);
                execlp("sh","sh","-c",task,NULL);
        }
        else
        {
                /* Parent process */
                close(fdIn[0]);
                close(fdOut[1]);
				st->inF=fdopen(fdIn[1],"w");
				st->outF=fdopen(fdOut[0],"r");
        }
}

int fm_getready(const SubT st){
	int c,err;
	while ((c=getc(st.outF))!=EOF){
		if (c=='*') break;
	}
	if (c==EOF) return -1;
	c=-1;
	err=fscanf(st.outF," READY%n",&c)<0;
	return (err || c<0);
}

int fm_load(const SubT st,const char *inputF,int *env_id){
	int err;
	err=fprintf(st.inF,"load %s\n",inputF)<1;
	if (!err) err=fflush(st.inF);
	if (!err) err=fscanf(st.outF,"%i",env_id)!=1;
	if (!err) err=*env_id<1;
	if (!err) err=fm_getready(st);
	return err;
}

int fm_natom(const SubT st,int env_id,int *natom){
	int err;
	if (env_id>0) {
		err=fprintf(st.inF,"natom %d\n",env_id)<1;
	} else {
		err=fprintf(st.inF,"natom\n")<1;
	}
	if (!err) err=fflush(st.inF);
	if (!err) err=fscanf(st.outF,"%i",env_id)!=1;
	if (!err) err=fm_getready(st);
	return err;
}
	
int fm_setpos(const SubT st,int env_id,int npos,const double pos[]){
	int err,i;
	if (env_id>0) {
		err=fprintf(st.inF,"setpos %d\n",env_id)<1;
	} else {
		err=fprintf(st.inF,"setpos\n")<1;
	}
	if (!err) err=fprintf(st.inF,"%d\n",npos)<1;
	for (i=0;i<npos;i++){
		if (!err) err=fprintf(st.inF," %23.14g",npos)<1;
	}
	if (!err) err=fprintf(st.inF,"\n* END\n")<1;
	if (!err) err=fflush(st.inF);
	if (!err) err=fm_getready(st);
	return err;
}

int fm_getpos(const SubT st,int env_id,int npos,double pos[]){
	int err, npos2,i,c;
	if (env_id>0) {
		err=fprintf(st.inF,"getpos %d\n",env_id)<1;
	} else {
		err=fprintf(st.inF,"getpos\n")<1;
	}
	if (!err) err=fflush(st.inF);
	if (!err) err=fscanf(st.outF,"%d",&npos2)!=1;
	if (!err) err=npos!=npos2;
	for (i=0;i<npos;i++){
		if (!err) err=fscanf(st.outF,"%lg",&pos[i])!=1;
	}
	while ((c=getc(st.outF))!=EOF){
		if (c=='*') break;
		putchar(c);
	}
	c=-1;
	if (!err) err=fscanf(st.outF," END%n",&c);
	if (!err) err=c<0;
	if (!err) err=fm_getready(st);
	return err;
}

int fm_calc_ef(const SubT st,int env_id,double *energy,int nforce,double force[]){
	int err, nforce2, i,c;
	if (env_id>0) {
		err=fprintf(st.inF,"calc_ef %d\n",env_id)<1;
	} else {
		err=fprintf(st.inF,"calc_ef\n")<1;
	}
	if (!err) err=fflush(st.inF);
	if (!err) err=fscanf(st.outF,"%lg",energy)!=1;
	if (!err) err=fscanf(st.outF,"%d",&nforce2)!=1;
	printf("nforce %d, %d\n",nforce,nforce2);
	if (!err) err=nforce!=nforce2;
	for (i=0;i<nforce;i++){
		if (!err) err=fscanf(st.outF,"%lg",&force[i])!=1;
	}
	while ((c=getc(st.outF))!=EOF){
		if (c=='*') break;
		putchar(c);
	}
	c=-1;
	if (!err) err=fscanf(st.outF," END%n",&c);
	if (!err) err=c<1;
	if (!err) err=fm_getready(st);
	return err;
}

int fm_calc_e(const SubT st,int env_id,double *energy){
	int err;
	if (env_id>0) {
		err=fprintf(st.inF,"calc_e %d\n",env_id)<1;
	} else {
		err=fprintf(st.inF,"calc_e\n")<1;
	}
	if (!err) err=fflush(st.inF);
	if (!err) err=fscanf(st.outF,"%lg",energy)!=1;
	if (!err) err=fm_getready(st);
	return err;
}
