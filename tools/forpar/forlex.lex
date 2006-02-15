/* the input routine wget() takes care of most skipping to be done
   it garantees that each command will be on a single line ending
   with '\n', without comments and repeated blanks, and that quotes
   come in pairs
*/

%option stack

%{

#define DEBUG
#undef  DEBUG

#define YYSTYPE char*
#define YY_INPUT(buf,result,max_size) { result=wget(buf,max_size); }

extern int wws, submod, fixed, ext, lineno;
extern char line[], inpline[];

#define MAX_INCLUDE_DEPTH 10
YY_BUFFER_STATE include_stack[MAX_INCLUDE_DEPTH];
int include_stack_ptr = 0, nstrt=0;
void lexerr(char *meldung, ...);
void wdebug(char *meldung, ...);
void strt_ni();
void mm();
void un_strt();
void zr_strt();
void loc(char* s);

#include "forpar.tab.h"
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>   /* needed for variable num. of parameters to w_fehler */

extern char stack[], *sp, *spmax;
char* push(char* s);
int getstr(char* str, char* lim);
int wget(char buf[],unsigned max_size);
char get1c(char* s);
void getid(char* s, char* id);
char c, wts[256], *pc;
int modsub=INITIAL;

%}




/* Start conditions:
   INITIAL: reckognize keywords
   NURID:   do not reckognize keywords
   M:       in module after CONTAINS statement, beginning of line
   MM:      in module after CONTAINS statement, not at beginning of line
*/

%x NURID
%x M
%x MM
%x SKIP

yy_push_state(
/* ws: at least one white-space, maybe more; ow: optional ws */
ws      " "+
ow      " "*
/* id: identifier */
id      [A-Za-z][A-Za-z_0-9]*
/* cid: composite id = identifier of element of defined type */
cid     {id}({ow}"%"{ow}{id})+
/* strc: string constant */
strc   \'([^'\n]|\'\')*\'|\"([^"\n]|\"\")*\"
/* logc: logical constant */
logc   ".true."|".false."
/* realc: real constant */
realc  (([0-9]+(\.[0-9]*)?)|(\.[0-9]+))([ed][+-]?[0-9]+)?(_(([0-9]+)|{id}))?
/* s: separate keywords */
s      [^A-Za-z_0-9]
/* ope: any operator */
ope: "."[a-zA-Z]+"."|[+\-/\*><]|"**"|"//"|"=="|"/="|"<="|">="



%%

       yylval=NULL;
       switch(wws){
       case 'n': strt_ni(); wws='p'; break;
       case 'b': un_strt(); wws=0;   break;
       case 'p': /* do nothing */;   break;
       }

<*>^\n { /* empty: skip empty lines */ }
<*>\n { zr_strt(); if( wws=='p' ); wws=0;
        sp=stack; /* at end-of-command reset stack-pointer sp */
        return EOC;
       }

<INITIAL>^"program"/{s}              { strt_ni(); return PROGRAM; }
<INITIAL>^"module"/{s}               { strt_ni(); return MODULE; }
<INITIAL,M>^"end module"/{s}         { zr_strt(); BEGIN(INITIAL);
                                       strt_ni(); submod=0;
                                       return ENDMOD; }
<INITIAL,M>^"end"{ow}\n              { return END; }
<INITIAL>^"contains"/{s}             { zr_strt(); BEGIN(M); return CONTAINS; }
<INITIAL,M>^"use"/{s}                { mm(); strt_ni(); return USE; }
<INITIAL,M>^"integer"/{s}            { mm(); return INTEGER; }
<INITIAL,M>^"real"/{s}               { mm(); return REAL; }
<INITIAL,M>^"double precision"/{s}   { mm(); return DPREC; }
<INITIAL,M>^"complex"/{s}            { mm(); return COMPLEX; }
<INITIAL,M>^"character"/{s}          { mm(); return CHARACTER; }
<INITIAL,M>^"logical"/{s}            { mm(); return LOGICAL; }
<INITIAL>^"type"/{s}                 { strt_ni(); return TYPE; }
<M>^"type"{ow}/"("                   { mm(); strt_ni(); return TYPE; }
<INITIAL>^"end"{ow}"type"/{s}        { strt_ni(); return ENDTYPE; }
<INITIAL,M>^"implicit none"/{s}      { mm(); return IMPLNONE; }
<INITIAL>^"include"/{s}              { strt_ni(); return INCLUDE; }
<INITIAL>^"parameter"/{s}            { strt_ni(); return PARAMETER; }
<INITIAL>^"public"/{s}               { strt_ni(); return PUBLIC; }
<INITIAL>^"private"/{s}              { strt_ni(); return PRIVATE; }
<INITIAL,M>^"interface"/{s}          { mm(); strt_ni(); return INTERFACE; }
<INITIAL>^"module procedure "/{id}   { strt_ni(); return MODPROC; }
<INITIAL,M>^"end"{ow}"interface"/{s} { mm(); strt_ni(); return ENDINTER; }
<INITIAL,M>^"recursive"/{s}          { mm(); return RECURSIVE; }
<INITIAL,M>^"pure"/{s}               { mm(); return PURE; }
<INITIAL,M>^"elemental"/{s}          { mm(); return ELEMENTAL; }

<M>^"function"/{s}                   { mm(); strt_ni(); return FUNCTION; }
<M>^"subroutine"/{s}                 { mm(); strt_ni(); return SUBROUTINE; }
<M>^"end function"/{s}               { strt_ni(); return ENDFUN; }
<M>^"end subroutine"/{s}             { strt_ni(); return ENDSUB; }
<M>^" "*\n                           { }                 /* skip empty line */
<M>^[0-9]{0,5}.                      { BEGIN(SKIP); }    /* skip line */
<SKIP>[^\n]*\n                       { BEGIN(M); }       /* skip line */


<INITIAL,MM>"function"/{s}           { mm(); strt_ni(); return FUNCTION; }
<INITIAL,MM>"result"/{s}             { strt_ni(); return RESULT; }
<INITIAL,MM>"subroutine"/{s}         { mm(); strt_ni(); return SUBROUTINE; }
<INITIAL,MM>"end function"/{s}       { strt_ni(); return ENDFUN; }
<INITIAL,MM>"end subroutine"/{s}     { strt_ni(); return ENDSUB; }
<INITIAL,MM>"include"/{s}            { return INCLUDE; }
<INITIAL,MM>"parameter"/{s}          { return PARAMETER; }
<INITIAL,MM>"public"/{s}             { return PUBLIC; }
<INITIAL,MM>"private"/{s}            { return PRIVATE; }
<INITIAL,MM>"sequence"/{s}           { return SEQUENCE; }
<INITIAL,MM>"allocatable"/{s}        { return ALLOCATABLE; }
<INITIAL,MM>"dimension"/{s}          { return DIMENSION; }
<INITIAL,MM>"external"/{s}           { return EXTERNAL; }
<INITIAL,MM>"intent"/{s}             { return INTENT; }
<INITIAL,MM>"intrinsic"/{s}          { return INTRINSIC; }
<INITIAL,MM>"optional"/{s}           { return OPTIONAL; }
<INITIAL,MM>"pointer"/{s}            { return POINTER; }
<INITIAL,MM>"save"/{s}               { return SAVE; }
<INITIAL,MM>"target"/{s}             { return TARGET; }
<INITIAL,MM>"in"/{s}                 { return IN; }
<INITIAL,MM>"out"/{s}                { return OUT; }
<INITIAL,MM>"inout"/{s}              { return INOUT; }





<NURID>{id}     { loc(yytext); yylval=push(yytext); return ID; }
<NURID>{cid}    { zapw(yytext); loc(yytext); yylval=push(yytext); return CID; }
<NURID>{ow}"only"{ow}":"{ow}                        { return ONLY; }
<NURID>{ow}"operator"{ow}/"("{ow}.[a-zA-Z]+{ow}.{ow}")" { return OPERATOR; }
<NURID>{ow}"operator"{ow}/"("{ow}[+\-/\*><]{ow}")"      { return OPERATOR; }
<NURID>{ow}"operator"{ow}/"("{ow}"**"{ow}")"            { return OPERATOR; }
<NURID>{ow}"operator"{ow}/"("{ow}"//"{ow}")"            { return OPERATOR; }
<NURID>{ow}"operator"{ow}/"("{ow}"=="{ow}")"            { return OPERATOR; }
<NURID>{ow}"operator"{ow}/"("{ow}"/="{ow}")"            { return OPERATOR; }
<NURID>{ow}"operator"{ow}/"("{ow}"<="{ow}")"            { return OPERATOR; }
<NURID>{ow}"operator"{ow}/"("{ow}">="{ow}")"            { return OPERATOR; }
<NURID>{ow}"assignment"{ow}/"("{ow}"="{ow}")"           { return ASSIGNMENT; }

<*>{realc} { yylval=push(yytext); return REALC; /* real const */ }
<*>{strc}  { yylval=push(yytext); return STRC;  /* string const */ }
<*>{logc}  { yylval=push(yytext); return LOGC;  /* logical const */ }
<*>"."[a-zA-Z]+"." { yylval=push(yytext); return DEFOPER ; /* operator */ }
<*>[+\-/><] { yylval=push(yytext); return OPER ; /* operator */ }
<*>"*"      { return '*'; }
<*>"**"|"//"|"=="|"/="|"<="|">=" { yylval=push(yytext); return OPER ; }
<*>" " { /* skip blanks */ }
<*>{ow}"::"{ow}        { return DCOLON; }
<*>{ow}"=>"{ow}        { return POINT; }
<*>{ow}"("{ow}"/"{ow}  { return ARRIN; }
<*>{ow}"/"{ow}")"{ow}  { return ARROUT; }
<*>{ow}[:(),=]{ow}     { for( pc=yytext; *pc==' '; ++pc ); return *pc; }
<*>{ow}"*"{ow}"("{ow}"*"{ow}")"{ow}  { return ARBCHARLEN; }

<INITIAL>{ow}"#include"{ow}\"/[^"\n]+\"\n  {
  if( getstr(wts,"\"\n")==EOF || input()!='"' )
    lexerr("illegal file-name \"%s\" in include",wts);
  if ( include_stack_ptr >= MAX_INCLUDE_DEPTH )
     lexerr("Includes nested too deeply");
  getstr(NULL,"\n"); input();	/* skip to end of include command */
  include_stack[include_stack_ptr]=YY_CURRENT_BUFFER;
  yyin=fopen( wts, "r" );
  include_stack_ptr++;
  if( ! yyin ) lexerr("could not open file \"%s\"\n",wts);
  yy_switch_to_buffer(yy_create_buffer(yyin, YY_BUF_SIZE));
}
<*><<EOF>> {
  if ( --include_stack_ptr < 0 )  yyterminate();
  else {
    yy_delete_buffer(YY_CURRENT_BUFFER);
    yy_switch_to_buffer(include_stack[include_stack_ptr]);
  }
}



<*>. { yylval=push(yytext);
       wdebug("UNKNOWN is >>>%s<<<, in line:\n \"%s\"",yylval,line);
       return UNKNOWN;
     }

%%


int yywrap(void){ return 1; }

void lexerr(char *meldung, ...){
  va_list args;

  va_start(args,meldung);	/* make args point to first argument of list */
  fprintf(stderr,"lexerr: error in line %d:\n",lineno);
  vfprintf(stderr,meldung,args);
  fprintf(stderr,"\n");
  va_end(args);			/* clear-up arg-list */

  exit(1);
}


void wdebug(char *meldung, ...){
  va_list args;
#ifdef DEBUG
  va_start(args,meldung);	/* make args point to first argument of list */
  vfprintf(stderr,meldung,args);
  fprintf(stderr,"\n");
  va_end(args);			/* clear-up arg-list */
#endif
  return;
}


int getc_safe(FILE *fp){
  static int eof=0;
  int c;
  if( eof ) return EOF;
  c=getc(fp);
  if( c==EOF ) eof=1;
  return c;
}


int getcc(void){}
int getcs(void){}


/* get first non-blank char in s */
char get1c(char* s){
  char* cp;
  for( cp=yytext; *cp!=' '&&*cp!='\0'; ++cp );
  return *cp;
}


/* get identifier from s and put it in id */
void getid(char* s, char* id){
  char* cp=s;
  while(*cp && *cp!=' ') ++cp;
  if( !isalpha(*cp) ) lexerr("illegal identifier in \"%s\"",s);
  *s=*cp; 
  while( *cp && (isalnum(c)||c=='_') ) *++s=*++cp;
  return;
}


/* reads string str from input until one of the chars in lim or EOF
   is encountered
   returns the last character encountered
   if str==NULL skip characters instead of putting them in str;
*/
int getstr(char* str, char* lim){
  int l,n,j,t,e;
  char c;

  l=strlen(lim); n=0;
  do {
    c=t=input();
    for( j=0,e=0; j<n; ++j ) e |= (c==lim[j]);
    e |= (t==EOF);
    if( str!=NULL ) str[n++]=c; else ++n;
  } while( !e );
  if( t!=EOF) unput(c);
  if( str!=NULL ) str[--n]='\0'; else --n;

  return t;
}

/* skip whitespace: return first non-ws char */
char* skpws(char *s){
  while( *s==' ' || *s=='\t' ) ++s;
  return s;
}

/* read a line and store it in zeile; if line ends with EOF, put '\n';
   returns number of chars read
*/
int readzeile(char *zeile){
  int n=0, c;
  while( (c=getc(yyin))!='\n' && c!=EOF )
       { *zeile++=c; ++n; }
  ++lineno;
  if( n==0 && c==EOF ) return 0;
  *zeile='\n'; ++n; *++zeile='\0';
  return n;
}



/* read a free format command, skipping double blanks and comments and
  joining lines */
int getfree(char buf[],unsigned max_size){
  int pc, c, cc, n, j, l;
  int istr=0;
  char *p, *pp;

  do {     /* skip empty lines and comments */
    if( (l=readzeile(inpline))==0 ) return YY_NULL;
  } while( (c=*skpws(inpline))=='\n' || c=='!' );

  pc='\0'; n=0; p=inpline;
  do {
    c=*p;
    if( c=='\0' ) lexerr("illegal char (null-char) encountered");
    if( istr ){
      if( c=='&' && *skpws(p+1)=='\n' ){    /* continue line in string */
        do {
          if( (l=readzeile(inpline))==0 ) lexerr("EOF in string");
        } while( (cc=*skpws(inpline))=='\n' || cc=='!' );
        p=inpline; pp=skpws(inpline); c='\0';
        if( *pp=='&' ) p=pp+1;
      }
      else if( c==istr ){ buf[n++]=c; istr=0; ++p; } /* end of string */
      else if( c=='\n' ) lexerr("end of line in string");
      else { buf[n++]=c; ++p; }                      /* char in string */
    }
    else
    {
      if( c=='\'' || c=='\"' ){ buf[n++]=c; istr=c; ++p; }
      else if( c==' ' || c=='\t' ){ /* skip leading and repeated blanks */
        p=skpws(p);
        if( pc!=' ' && pc!='\n' && n>0 ) c=buf[n++]=' ';   /* no tabs remain */
      }
      else if( c=='!' || c=='\n' || c==';' ){ /* end command */
        if( pc==' ' || pc=='\n' ) --n; /* kill trailing whitespace */
	buf[n++]='\n';
        if( c!=';' ) c='\n';    /* ';' ends command, not line */
        ++p;
      }
      else if( c=='&' ){        /* continue line */
        c=*(skpws(p+1));
        if( c!='!' && c!='\n' && c!='\0' ) lexerr(
            "only comments or whitespace allowed after continuation mark (&)");
        do {                    /* skip inline-comments */
          if( (l=readzeile(inpline))==0 )
            lexerr("unexpected end of file after continuation mark (&)");
        } while( (cc=*skpws(inpline))=='!' || cc=='\n' || cc=='\0' );
        p=inpline; pp=skpws(inpline); c='\0';
        if( *pp=='&' ) p=pp+1;
      }
      else { buf[n++]=c; ++p; }
    }
    /* at this point p points to the next char to be read */
    if( c!='\0' ) pc=c;         /* c=='\0' means no char added to buf */
    if( n>=max_size ) lexerr("OOOOPS... line to long (>=%d)",max_size);
  } while( c!='\n' );
  
  return n;
}


/* read a line of at most 72 chars, skip rest of line, return last char read */
int rzeile(int r, char zeile[]){
  static int c='\0';
  int j;

  if( !r && *zeile ) return c;
  for( j=0; j<72 && (c=getc(yyin))!=EOF && c!='\n'; ) zeile[j++]=c;
  while( c!='\n' && c!=EOF ) c=getc(yyin); /* skip to end of line */
  zeile[j]='\0'; ++lineno;

  return c;
}


/* check if zeile is continuation line */
int contline(char zeile[]){
  int j;
  if( tolower(zeile[0])=='c' ) return 0;
  if( strlen(zeile)<6 ) return 0;
  for( j=0; j<6; ++j ) if( zeile[j]=='\t' ) return 0;
  if( zeile[5]==' ' ) return 0;
  return 1;
}


/* get a fixed-format command, skipping double blanks and comments
   and joining lines */
int getfixed(char buf[],unsigned max_size){
  static char zeile[73]="";
  int c,cc,n=0,pc=0,j;
  int istr=0;
  static int fine=0;

  cc=rzeile(0,zeile);		/* read only if not already present */
  while(    tolower(*zeile)=='c'   || *zeile=='*'       /* skip comments and */
         || (c=*skpws(zeile))=='!' || c=='\0'     ){    /* empty lines */
    cc=rzeile(1,zeile);
  }
  if( *zeile=='\0' && cc==EOF ) return YY_NULL;

  j=0;
  do{
    while( j<strlen(zeile) ){
      c=zeile[j];
      if( c=='\0' ) lexerr("illegal char (null-char) encountered");
      if( istr ){
        if( c==istr ){ buf[n++]=c; istr=0; } /* end of string */
        else buf[n++]=c;
      } else {
        if( c=='\'' || c=='\"' ){ buf[n++]=c; istr=c; }
        else if( c==' ' || c=='\t' ) /* skip leading and repeated blanks */
	  { if( n>0 && pc!=' ' ) buf[n++]=' '; } /* no tabs remain */
        else buf[n++]=c;
      }
      if( n>0 ) pc=buf[n-1];
      ++j;
    }
    cc=rzeile(1,zeile); j=6;
  } while( contline(zeile) );
  buf[n++]='\n'; buf[n]='\0';
  if( n>max_size ) lexerr("OOOOPS... line to long (>%d)",max_size);
  if( istr ) lexerr("line ended while reading string const");
  for( j=n-1; j>0 && buf[j]==' '; ){ buf[j]=buf[j+1]; buf[j+1]='\0'; --j,--n; }
  return n;
}


/* wget: input routine:
  get a command, skipping double blanks and comments and joining lines */
int wget(char buf[],unsigned max_size){
  int j,n;

  switch( fixed ){
  case 0: n=getfree(buf,max_size); break;
  case 1: n=getfixed(buf,max_size); break;
  default: parerr("PANICO: wrong value of fixed in wget()");
  }
  for( j=0; j<n; ++j ) line[j]=buf[j];
  if( n>0 && line[n-1]=='\n' ) line[n-1]='\0'; else line[n]='\0';
  wdebug("Read line: %s",buf);

  return n;
}


char* push(char* s){
  char* p=sp;
  strcpy(sp,s); sp+=strlen(sp)+1;
  if( sp>spmax ) spmax=sp;
  return p;
}


void strt_ni(){
  yy_push_state(NURID); ++nstrt; 
  wdebug("strt_ni: state now is %d (n=%d)",YY_START,nstrt); 
}


void mm(){
  if( !submod ) return;
  if( ext ){ ext=0; BEGIN(M); }
  yy_push_state(MM); ++nstrt;
  wdebug("strt_ni(mm): state now is %d (n=%d)",YY_START,nstrt); 
}


void un_strt(){
  while( nstrt>0 && yy_top_state()==YY_START){ yy_pop_state(); --nstrt;}
  if( nstrt>0 ){ yy_pop_state(); --nstrt;}
  wdebug("un_strt: state now is %d (n=%d)",YY_START,nstrt); 
}

void zr_strt(){
  wdebug("reset start-state stack");
  while(nstrt>0){ yy_pop_state(); --nstrt; } 
}

void loc(char* s){ while( *s && s ){ *s=tolower(*s); ++s; } }

/* zap whitespace */
int zapw(char* s){
  int j;
  char* p;
  for( p=s; *p!='\0'; ++p ){}
}

/*
int yywrap(void){ return 1; }


int main(void){
  while( yylex()!=EOF );
}
*/
