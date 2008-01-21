/*
 *  cp2k_shell.h
 *  cp2kController
 *
 *  Created by Fawzi Mohamed on 1/20/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */
#ifndef __CP2K_SHELL_H
#define __CP2K_SHELL_H
#include <stdio.h>
#include <sys/types.h>
#include <stdlib.h>

struct Subtask{
	FILE *inF;
	FILE *outF;
};
typedef struct Subtask SubT;

void fm_spawnSubtask(const char *task,const char *taskDir, SubT *st);
int fm_getready(const SubT st);
int fm_load(const SubT st,const char *inputF,int *env_id);
int fm_natom(const SubT st,int env_id,int *natom);
int fm_setpos(const SubT st,int env_id,int npos,const double pos[]);
int fm_getpos(const SubT st,int env_id,int npos,double pos[]);
int fm_calc_ef(const SubT st,int env_id,double *energy,int nforce,double force[]);
int fm_calc_e(const SubT st,int env_id,double *energy);
#endif