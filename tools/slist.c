#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "slist.h"


/* finds element "newele" in list "list" (NULL: not in list) */
slist* find_lele(slist **list, char* newele){
  slist *p;
  int f;

  for( p=*list; p && (f=strcmp(newele,p->ele))>0; p=p->nxt );
  if( p && f!=0 ) p=NULL;

  return p;
}


/* delete element newele from list */
void del_lele(slist **list, char* newele){
  slist *p,*pp;
  int f;

  for( p=*list,pp=NULL; p && (f=strcmp(newele,p->ele))>0; pp=p,p=p->nxt );
  if( !p || f!=0 ) return;	/* element not in list */
  if( pp )
    pp->nxt=p->nxt;		/* remove element p from list */
  else
    *list=p->nxt;
  free( p->ele );
  free( p );

  return;
}


/* delete and dellocate all elements from list */
void zap_lele(slist **list){
  slist *p=*list,*pp=NULL;

  while( p ){ pp=p; p=p->nxt; free(pp); }
  *list=NULL;

  return;
}


/*   insert element "newele" into sorted list "list".
 *   "list" is entry-point of sorted list.
 *   pp and p point to the element before and after newele
 */
void ins_lele(slist **list, char* newele){
  slist *p,*pn,*pp;
  int f;
  char *s;

  for( p=*list,pp=NULL; p && (f=strcmp(newele,p->ele))>0; pp=p,p=p->nxt );
  if( p && f==0 ) return; /* newele already exists */
  if( !(pn=malloc(sizeof(slist))) )
    { fprintf(stderr,"malloc error in ins_lele (pn)\n"); exit(1); }
  if( !pp )			/* put at head of list */
    { pn->nxt=*list; *list=pn; }
  else				/* insert pn between pp and p */
    { pn->nxt=p; pp->nxt=pn; }
  if( !(s=malloc(strlen(newele)+1)) )
    { fprintf(stderr,"malloc error in ins_lele (s)\n"); exit(1); }
  strcpy(s,newele);
  pn->ele=s;
}


