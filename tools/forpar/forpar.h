#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>   /* needed for variable num. of parameters to w_fehler */
#include "slist.h"    /* sorted lists of char* */
#if defined ( __DECC )
#include <alloca.h>
#endif

#define YYSTYPE char*

#define NI wws='n';
#define YK wws='b';

int wws=0, fixed=0, chkint=0, ok=0, lineno=0, filedep=1;
int submod=1, cont=0, ext=1, internal=0, interf=0, wstat='\0';

int private=0, public=1;
int on=1;

#define KB *1024
#define MB *1048576

int indent=0, dind=2, nargbuf=-1, nfunname=-1;
slist *pub_l=NULL, *pri_l=NULL, *pl;
slist *use_l=NULL;
char stack[16 KB], *sp=stack, *spmax=stack;
char heap[64 KB], *hp=heap, *hpmax=heap;
char line[4 KB], inpline[16 KB], argbuf[3][4 KB], *args;
char modname[64], typename[64], intername[64]; /*module, type, interface name*/
char proname[64], buffunname[3][64], *funname; /* program, funct-subroutine */
char dimension[1 KB]="";	/* dimension of array */
char semp[1]="";
char *command;			/* name of command (*argv[0]) */
FILE *flout;
char fnin[256]="", fnout[256]="", fnold[256]="", path[256]="";
char *pc;                       /* temporary pointer to char */

void parerr(char *meldung, ...);
void parwarn(char *meldung, ...);
char* find(char* pat,char* in);
char* present(char* tok, char* list, char* delim);
char* spr(int m, char* format,...);
void printw(int m, char* format,...);
int ispub(char* name);
int check_ext(int* c);
FILE *newfile(char *fnin, char *fnout, int open, int addpath);
