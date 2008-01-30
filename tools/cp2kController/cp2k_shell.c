/**
 *  \file cp2k_shell.c
 *
 *  Functions to control an external cp2k_shell task through std_in/std_out.
 *
 *  The functions returns an int as error value.
 *  
 * \note this file is available also with a BSD style license
 *
 * \copyright (C) 2008  CP2K developers group
 * \author Fawzi Mohamed
 */

#include "cp2k_shell.h"
#include <unistd.h>
#include <stdlib.h>

/**
 * \brief executes (execvp) task with arguments taskArgs in the directory taskDir as subtask.
 *
 * returns file handles std in and std out of the task in \param st
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
				st->isReady=0;
        }
}

/**
 * \brief reads the readiness signal of the task in not already done
 * (asynchronous tasks do not read readiness signal before returning)
 */
int fm_getready_if_needed(SubT *st){
	if (! st->isReady) {
		int err=0;
		err=fm_getready(st);
		if (err) return 200;
	}
	return 0;
}

/**
 * \brief reads the readiness signal of the task
 */
int fm_getready(SubT *st){
	int c,err,writeout=0;
	err=0;
	while ((c=getc(st->outF))!=EOF){
		if (c=='*') break;
		if (c!=' ' && c!='\n' && c!='\t') writeout=1;
		if (writeout) putchar(c); // write out skipped chars
	}
	if (c==EOF) err=-1;
	c=-1;
	if (!err) err=fscanf(st->outF," READY%n",&c)<0;
	if (!err) err=c<0;
	st->isReady= !err;
	return err;
}

/**
 * \brief loads an input file and return the env_id of it
 */
int fm_load(SubT *st,const char *inputF,int *env_id){
	int err;
	err=fm_getready_if_needed(st);
	if (!err) err=fprintf(st->inF,"load %s\n",inputF)<1;
	if (!err) err=fflush(st->inF);
	if (!err) err=fscanf(st->outF,"%d",env_id)!=1;
	if (!err) err=*env_id<1;
	if (!err) err=fm_getready(st);
	return err;
}

/**
 * \brief Begins loading an input file.
 *
 * Asynchronous, does wait or for completion.
 */
int fm_bg_load(SubT *st,const char *inputF){
	int err;
	err=fm_getready_if_needed(st);
	if (!err) err=fprintf(st->inF,"bg_load %s\n",inputF)<1;
	if (!err) err=fflush(st->inF);
	st->isReady=0;
	return err;
}

/**
 * \brief returns the env_id of the last file loaded
 */
int fm_last_env_id(SubT *st,int *env_id){
	int err;
	err=fm_getready_if_needed(st);
	if (!err) err=fprintf(st->inF,"last_env_id\n")<1;
	if (!err) err=fflush(st->inF);
	if (!err) err=fscanf(st->outF,"%d",env_id)!=1;
	if (!err) err=*env_id<1;
	if (!err) err=fm_getready(st);
	return err;
}

/**
 * \brief returns number of atoms in the given env_id
 */
int fm_natom(SubT *st,int env_id,int *natom){
	int err;
	err=fm_getready_if_needed(st);
	if (!err) {
		if (env_id>0) {
			err=fprintf(st->inF,"natom %d\n",env_id)<1;
		} else {
			err=fprintf(st->inF,"natom\n")<1;
		}
	}
	if (!err) err=fflush(st->inF);
	if (!err) err=fscanf(st->outF,"%d",natom)!=1;
	if (!err) err=fm_getready(st);
	return err;
}
	
/**
 * \brief changes the position of the atoms of the given env_id
 */
int fm_setpos(SubT *st,int env_id,int npos,const double pos[]){
	int err,i;
	err=fm_getready_if_needed(st);
	if (!err) {
		if (env_id>0) {
			err=fprintf(st->inF,"setpos %d\n",env_id)<1;
		} else {
			err=fprintf(st->inF,"setpos\n")<1;
		}
	}
	if (!err) err=fprintf(st->inF,"%d\n",npos)<1;
	for (i=0;i<npos;i++){
		if (!err) err=fprintf(st->inF," %22.13g",pos[i])<1;
	}
	if (!err) err=fprintf(st->inF,"\n* END\n")<1;
	if (!err) err=fflush(st->inF);
	if (!err) err=fm_getready(st);
	return err;
}

/**
 * \brief returns the position of the atoms of the given env_id
 */
int fm_getpos(SubT *st,int env_id,int npos,double pos[]){
	int err, npos2,i,c;
	err=fm_getready_if_needed(st);
	if (!err){
		if (env_id>0) {
			err=fprintf(st->inF,"getpos %d\n",env_id)<1;
		} else {
			err=fprintf(st->inF,"getpos\n")<1;
		}
	}
	if (!err) err=fflush(st->inF);
	if (!err) err=fscanf(st->outF,"%d",&npos2)!=1;
	if (!err) err=npos!=npos2;
	for (i=0;i<npos;i++){
		if (!err) err=fscanf(st->outF,"%lg",&pos[i])!=1;
	}
	while ((c=getc(st->outF))!=EOF){
		if (c=='*') break;
	}
	c=-1;
	if (!err) err=fscanf(st->outF," END%n",&c);
	if (!err) err=c<0;
	if (!err) err=fm_getready(st);
	return err;
}

/**
 * \brief begins the computation of e and f.
 * 
 * Asynchronous, does not wait for completion
 */
int fm_eval_ef(SubT *st,int env_id){
	int err;
	err=fm_getready_if_needed(st);
	if (!err) {
		if (env_id>0) {
			err=fprintf(st->inF,"eval_ef %d\n",env_id)<1;
		} else {
			err=fprintf(st->inF,"eval_ef\n")<1;
		}
	}
	if (!err) err=fflush(st->inF);
	st->isReady=0;
	return err;
}

/**
 * \brief gets the force of the last force computation.
 * \warning if you calculated just the energy, the content of the array are undefined.
 */
int fm_get_f(SubT *st,int env_id,int nforce,double force[]){
	int err, nforce2, i, c;
	err=fm_getready_if_needed(st);
	if (!err) {
		if (env_id>0) {
			err=fprintf(st->inF,"get_f %d\n",env_id)<1;
		} else {
			err=fprintf(st->inF,"get_f\n")<1;
		}
	}
	if (!err) err=fflush(st->inF);
	if (!err) err=fscanf(st->outF,"%d",&nforce2)!=1;
	if (!err) err=nforce!=nforce2;
	for (i=0;i<nforce;i++){
		if (!err) err=fscanf(st->outF,"%lg",&force[i])!=1;
	}
	while ((c=getc(st->outF))!=EOF){
		if (c=='*') break;
	}
	c=-1;
	if (!err) err=fscanf(st->outF," END%n",&c);
	if (!err) err=c<1;
	if (!err) err=fm_getready(st);
	return err;
}


/**
 * \brief Calculates the energy and forces and returns them
 */
int fm_calc_ef(SubT *st,int env_id,double *energy,int nforce,double force[]){
	int err,nforce2,i,c;
	err=fm_getready_if_needed(st);
	if (!err) {
		if (env_id>0) {
			err=fprintf(st->inF,"calc_ef %d\n",env_id)<1;
		} else {
			err=fprintf(st->inF,"calc_ef\n")<1;
		}
	}
	if (!err) err=fflush(st->inF);
	if (!err) err=fscanf(st->outF,"%lg",energy)!=1;
	if (!err) err=fscanf(st->outF,"%d",&nforce2)!=1;
	if (!err) err=nforce!=nforce2;
	for (i=0;i<nforce;i++){
		if (!err) err=fscanf(st->outF,"%lg",&force[i])!=1;
	}
	while ((c=getc(st->outF))!=EOF){
		if (c=='*') break;
	}
	c=-1;
	if (!err) err=fscanf(st->outF," END%n",&c);
	if (!err) err=c<1;
	if (!err) err=fm_getready(st);
	return err;
}

/**
 * \brief begins the computation of e.
 * 
 * Asynchronous, does not wait for completion
 * \note the force of the environment is undefined after this call
 */
int fm_eval_e(SubT *st,int env_id){
	int err;
	err=fm_getready_if_needed(st);
	if (!err) {
		if (env_id>0) {
			err=fprintf(st->inF,"eval_e %d\n",env_id)<1;
		} else {
			err=fprintf(st->inF,"eval_e\n")<1;
		}
	}
	if (!err) err=fflush(st->inF);
	st->isReady=0;
	return err;
}

/**
 * \brief Calculates the energy and returns it
 *
 * the force of the environment is undefined after this call
 */
int fm_calc_e(SubT *st,int env_id,double *energy){
	int err;
	if (env_id>0) {
		err=fprintf(st->inF,"calc_e %d\n",env_id)<1;
	} else {
		err=fprintf(st->inF,"calc_e\n")<1;
	}
	if (!err) err=fflush(st->inF);
	if (!err) err=fscanf(st->outF,"%lg",energy)!=1;
	if (!err) err=fm_getready(st);
	return err;
}

/**
 * \brief returns the energy of the last computation
 */
int fm_get_e(SubT *st,int env_id,double *energy){
	int err;
	err=fm_getready_if_needed(st);
	if (!err) {
		if (env_id>0) {
			err=fprintf(st->inF,"get_e %d\n",env_id)<1;
		} else {
			err=fprintf(st->inF,"get_e\n")<1;
		}
	}
	if (!err) err=fflush(st->inF);
	if (!err) err=fscanf(st->outF,"%lg",energy)!=1;
	if (!err) err=fm_getready(st);
	return err;
}

/**
 * \brief returns some informations about the underlying cp2k task
 *
 * \note assumes that no line begins with '*'
 */
int fm_get_info(SubT *st,int sizeBuf,char *buf,int *totalSize){
	int err,nread=0,atStart=1,c,totc=1;
	err=fm_getready_if_needed(st);
	if (!err) err=fprintf(st->inF,"info\n")<1;
	if (!err) err=fflush(st->inF);
	while ((c=getc(st->outF))!=EOF){
		if (atStart && c=='*') break;
		if (c!=' ' && c!='\t') atStart=0;
		if (c=='\n' || c=='\r') atStart=1;
		if (nread<sizeBuf) {
			buf[nread]=(char)c;
			nread++;
		}
		totc++;
	}
	if (nread<sizeBuf) nread++;
	if (nread<sizeBuf) buf[nread]='\0';
	if (c==EOF) err=-1;
	c=-1;
	if (!err) err=fscanf(st->outF," READY%n",&c)<0;
	if (!err) err=c<0;
	st->isReady= !err;
	*totalSize=totc;
	return err;
}
