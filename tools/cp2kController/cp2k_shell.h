/**
 *  \file cp2k_shell.h
 *
 *  Functions to control an external cp2k_shell task through std_in/std_out.
 *  
 * \note this file is available also with a BSD style license
 *
 * \copyright (C) 2008  CP2K developers group
 * \author Fawzi Mohamed
 */
#ifndef __CP2K_SHELL_H
#define __CP2K_SHELL_H
#include <stdio.h>
#include <sys/types.h>
#include <stdlib.h>

struct Subtask{
	FILE *inF;
	FILE *outF;
	int isReady;
};
typedef struct Subtask SubT;

void fm_spawnSubtask(const char *task,const char *taskDir, SubT *st);
int fm_getready_if_needed(SubT *st);
int fm_getready(SubT *st);
int fm_load(SubT *st,const char *inputF,int *env_id);
int fm_bg_load(SubT *st,const char *inputF);
int fm_last_env_id(SubT *st,int *env_id);
int fm_natom(SubT *st,int env_id,int *natom);
int fm_setpos(SubT *st,int env_id,int npos,const double pos[]);
int fm_getpos(SubT *st,int env_id,int npos,double pos[]);
int fm_eval_ef(SubT *st,int env_id);
int fm_eval_e(SubT *st,int env_id);
int fm_calc_ef(SubT *st,int env_id,double *energy,int nforce,double force[]);
int fm_calc_e(SubT *st,int env_id,double *energy);
int fm_get_f(SubT *st,int env_id,int nforce,double force[]);
int fm_get_e(SubT *st,int env_id,double *energy);
int fm_get_info(SubT *st,int sizeBuf,char *buf,int *fullSize);
#endif
