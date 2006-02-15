/*
 *      wns!:  non-standard feature
 *      +++:   to finish
 *      Test on QuickStep 2.2
 *      - fix string continuation ("text"\n  "text2" and "text"\n  &"text2"
 *                                 seem to be ok for IBMs only!!!)
 *      - errors with interfaces in "interfaces.f" and "pao.f"
 *      - character*(*)
 *      - 
 *      - 
 *
 *      TODO:
 *      - complete save command
 *      - external functions cannot be in the same file as modules or program
 *        (what if a module has the same name as the ext-fun file?)
 *      - write *.int file in same directory as *.f90 file, also if it
 *        is not in the current directory
 *      - parse executable part?
 *      - end statt end function
 *      - DIMENSION always after ID, not in attributes
 *      - specification statements
 *      - Befehle alphabetisch ordnen
 *      - Include?
 *      
 *      DONE:
 *      - when possible zero stack and heap
 *      - no nested functions/subroutines
 *      - parse program
 *      - parse more than one module in same file
 *      
 */


%{
#define YYDEBUG 1
#undef  YYDEBUG
#include "forpar.h"
extern FILE* yyin;
%}

%token MODULE ENDMOD CONTAINS PROGRAM
%token INTEGER REAL DPREC COMPLEX CHARACTER LOGICAL TYPE ENDTYPE
%token USE ONLY OPERATOR ASSIGNMENT
%token ID CID
%token POINT "=>"
%token ARBCHARLEN "*(*)"
%token DCOLON "::"
%token ARRIN "(/"
%token ARROUT "/)"
%token IMPLICIT IMPLNONE PARAMETER
%token OPER DEFOPER REALC STRC LOGC
%token INTERFACE ENDINTER MODPROC
%token RECURSIVE PURE ELEMENTAL FUNCTION RESULT ENDFUN SUBROUTINE ENDSUB

%token END EOC
%token UNKNOWN

%token PUBLIC PRIVATE SEQUENCE ALLOCATABLE DIMENSION EXTERNAL
%token INTENT IN OUT INOUT INTRINSIC OPTIONAL POINTER SAVE TARGET INCLUDE

%token FORMAT ENTRY MODULSUB

 

%%

/************ Grammar ************/

unit: program
    | extfunsub
    | module_l
    ;

program: program_s { parwarn("main program");
                     ok=1; exit(0);
                   }
       ;

program_s: PROGRAM ID EOC { strcpy(proname,$2); wstat='p';
                            if( !filedep ) flout=newfile(proname,fnout,1,1);
                            printw(-1,"PROGRAM %s\n\n",proname);
                            indent+=dind;
                            submod=public=private=ext=0;
                            if( !filedep ){ fclose(flout); flout=NULL; }
                          }
         ;
extfunsub: modulsub_l {
             if( !filedep ) { 
               fclose(flout); flout=NULL;
               if( chkint ) do_chkint(fnout,fnold);
             }
           }
         ;
module_l: module
        | module_l module
        ;
module: module_s specpart modsub_o endmod_s { public=private=0; }
      ;

module_s: MODULE ID EOC { strcpy(modname,$2); wstat='m';
                          if( !filedep ) flout=newfile(modname,fnout,1,1);
                          printw(-1,"MODULE: %s\n\n",modname);
                          indent+=dind;
                          submod=public=private=ext=0;
                        }
        ;
specpart: use_slp implpart_l declconstr_l
                  { printw(0,"\n");
		    if( private ) printw(-1,"Private\n");
		    if( public )  printw(-1,"Public\n");
		    if( pub_l ) printw(-1,"Public:: ");
		    for( pl=pub_l; pl; pl=pl->nxt )
		      printw(0,"%s %s",((pl==pub_l)? "" : ","), pl->ele);
		    if( pub_l ) printw(0,"\n");
		    if( pri_l ) printw(-1,"Private:: ");
		    for( pl=pri_l; pl; pl=pl->nxt )
		      printw(0,"%s %s",((pl==pri_l)? "" : ","), pl->ele);
		    if( pri_l ) printw(0,"\n");
		    printw(0,"\n");
                  }
        ;
use_s_l: /* empty */
       | use_s_l use_s                  { ins_lele(&use_l,$2); hp=heap; }
       ;
use_slp: use_s_l { for( pl=use_l; pl; pl=pl->nxt ) printw(0,pl->ele);
                   zap_lele(&use_l); }
       ;
implpart_l: /* empty */
          | implpart_l implpart         { hp=heap; }
          ;
declconstr_l: /* empty */
            | declconstr_l declconstr    { hp=heap; }
            ;
/******** use_s ********/
use_s: USE ID usearg EOC { $$=spr(-1,"Use %s%s\n",$2,$3); }
     ;
use_spl: /* empty */ { $$=semp; }
       | use_spl use_sp
       ;
use_sp: use_s { printw(0,$1); }
      ;
usearg: crename_l
      | ',' ONLY only_lo   { $$=spr(0,", Only: %s",$3); }
      ;
only_lo: /* empty */ { $$=semp; }
       | only_l
       ;
only_l: only              { $$=spr(0,"%s",$1); }
      | only_l ',' only   { $$=spr(0,"%s, %s",$1,$3); }
      ;
crename_l: /* empty */    { $$=semp; }
         | crename_l ',' rename  { $$=spr(0,", %s",$3); }
         ;
rename: ID "=>" ID       { $$=spr(0,"%s=>%s",$1,$3); }
      | opera "=>" opera { $$=spr(0,"%s=>%s",$1,$3); }
      ;

only: rename
    | ID
    | opera
/*    | genspec wns! */
    ;
opera: OPERATOR '(' oper ')' { $$=spr(0,"Operator(%s)",$3); }
     ;

/******** implpart ********/
implpart: implicit_s
        | include_s
        | parameter_s
        | format_s
        | entry_s
        ;
implicit_s: IMPLNONE EOC { printw(-1,"Implicit None\n"); }
            /* wns! implicit implicit-spec-list (R542) missing */
          ;
implicit_sl: /* empty */ { $$=semp; }
           | implicit_sl implicit_s
           ;
parameter_s: PARAMETER '(' nconstdef_l ')' EOC { printw(-1,"parameter\n"); }
           ;
nconstdef_l: nconstdef
           | nconstdef_l ',' nconstdef
           ;
nconstdef: ID '=' expr { $$=spr(0,"%s = %s",$1,$3); printw(-1,"%s\n",$$); }
         ;
const: srealc | STRC | LOGC | cmplxconst
     ;
srealc: REALC | '+' REALC | '-' REALC
      ;
expr: const
    | gid                 { $$=$1; }
/*    | gid '(' expr_l ')' { $$=spr(0,"%s(%s)",$1,$3); } */
    | gid '(' arfun_l ')' { $$=spr(0,"%s(%s)",$1,$3); }
    | oper expr           { $$=spr(0,"%s%s",$1,$2); }
    | expr oper expr      { $$=spr(0,"%s%s%s",$1,$2,$3); }
    | '(' expr ')'        { $$=spr(0,"(%s)",$2); }
    ;
gid: ID | CID;

arfun: expr { $$=$1; }
     | ':'  { $$=":"; } /* allow for unspecified index in arrays */
       ;
arfun_l: arfun
       | arfun_l ',' arfun { $$=spr(0,"%s,%s",$1,$3); }
       ;
expr_l: expr            { $$=$1; }
      | expr_l ',' expr { $$=spr(0,"%s,%s",$1,$3); }
      ;
cmplxconst: '(' srealc ',' srealc ')' { $$=spr(0,"(%s,%s)",$2,$4); }
          ;
oper: OPER
    | DEFOPER
    | '*'     { $$="*"; } 
    ;

/******** declconstr ********/
declconstr: dertypedef
          | interfaceblk
          | typedecl_s
          | specific_s
          | parameter_s
          | include_s
          | format_s       /* +++ Format statement */
          | entry_s        /* +++ Entry statement */
/*          | sfunc_s      /* +++ statement-function statement */
          ;
specific_s: access_s
          | save_s
/* +++         | allocatable_s    */
/* +++         | common_s     */
/* +++         | data_s       */
/* +++         | dimension_s  */
/* +++         | equivalence_s    */
/* +++         | external_s   */
/* +++         | intent_s     */
/* +++         | intrinsic_s  */
/* +++         | namelist_s   */
/* +++         | optional_s   */
/* +++         | pointer_s    */
/* +++         | target_s     */
          ;

dertypedef: dertype_s privseq_sl componentdef_l endtype_s
          ;

dertype_s: TYPE dcol_o ID EOC
                { strcpy(typename,$3); printw(-1,"Type:: %s\n",$3);
	                      indent+=dind; }
         | TYPE ',' {YK} accessspec "::" {NI} ID EOC { strcpy(typename,$7);
                         printw(-1,"Type, %s:: %s\n",$4,$7); indent+=dind; }
         ;
access_s: private_s
        | public_s
        ;

save_s: SAVE EOC { printw(-1,"Save\n"); } /* +++ should accept an arg-list */
      ;

include_s: INCLUDE '(' ID '.' ID ')' EOC { printw(-1,"Include '%s.%s'\n",$3,$5); }
         ;

privseq_sl: /* empty */ { $$=semp; }
          | privseq_sl privseq_s
          ;
privseq_s: PRIVATE  EOC { $$="Private"; }
         | SEQUENCE EOC { $$="Sequence"; }
         ;
componentdef_l: componentdef
              | componentdef_l componentdef
              ;
componentdef: typespec {NI} dcol_o entitydec_l EOC
                       { printw(-1,"%s:: %s\n",$1,$4); }
            | typespec {NI} ',' {YK} coattrspec_l "::" {NI} entitydec_l EOC {
                              printw(-1,"%s, %s:: %s\n",$1,$5,$8); }
            ;
coattrspec_l: coattrspec
            | coattrspec_l ',' coattrspec { $$=spr(0,"%s, %s",$1,$3); }
            ;
coattrspec: DIMENSION '(' {NI} arrspec {YK} ')' { 
                                         $$=spr(0,"Dimension(%s)",$4); }
          | POINTER          { $$="Pointer"; }
          ;

endtype_s: ENDTYPE EOC { indent-=dind; printw(-1,"End Type %s\n",typename); }
         | ENDTYPE ID EOC { indent-=dind; printw(-1,"End Type %s\n",$2);
	                    if( strcmp(typename,$2)!=0 )
			      parerr("type-name mismatch");
	                  }
         ;
private_s: PRIVATE EOC                 { private=1; }
         | PRIVATE dcol_o modent_l EOC {  }
         ;
public_s: PUBLIC EOC              { public=1; }
        | PUBLIC dcol_o pub_l EOC {  }
        ;
pub_l: modent            { ins_lele(&pub_l,$1); }
     | pub_l ',' modent  { $$=spr(0,"%s, %s",$1,$3); ins_lele(&pub_l,$3); }
        ;
pri_l: modent            { ins_lele(&pri_l,$1); }
     | pri_l ',' modent  { $$=spr(0,"%s, %s",$1,$3); ins_lele(&pri_l,$3); }
        ;
modent_l: modent
        | modent_l ',' modent { $$=spr(0,"%s, %s",$1,$3); }
        ;
modent: ID
      | opera
      | ASSIGNMENT '(' '=' ')' { $$=spr(0,"Assignment(=)",$3); }
      ;
modent_o: /* empty */ { $$=semp; }
        | modent
        ;

typedecl_s: typespec {NI}      entitydec_l EOC
                     { printw(-1,"%s:: %s\n",$1,$3); }
          | typespec {NI} "::" entitydec_l EOC
                     { printw(-1,"%s:: %s\n",$1,$4); }
          | typespec {NI} ',' {YK} attrspec_l "::" {NI} entitydec_l EOC
                     { printw(-1,"%s, %s:: %s\n",$1,$5,$8); }
          ;
entitydec_l: entitydec
           | entitydec_l ',' entitydec { $$=spr(0,"%s, %s",$1,$3); }
           ;
entitydec: ID arrspec_o lensel_o init_o { $$=spr(0,"%s%s%s%s",$1,$2,$3,$4); }
         ;
init_o: /* empty */     { $$=semp; }
      | '=' expr        { $$=spr(0," = %s",$2); }
      | '=' arrinit     { $$=spr(0," = %s",$2); }
      | "=>" ID '(' ')' { if( strcmp($2,"null")!=0 ) yyerror("expecting NULL");
                          $$=spr(0,"=>NULL()");
                        }
      ;
arrinit: ARRIN expr_l ARROUT   { $$=spr(0,"(/ %s /)",$2); }
       ;
arrspec_o: /* empty */       { $$=semp; }
         | '(' arrspec ')'   { $$=spr(0,"(%s)",$2); }
         ;
charlen_o: /* empty */ { $$=semp; }
         | '*' expr    { $$=spr(0,"*%s",$2); }
         ;
lensel_o: /* empty */  { $$=semp; }
        | lensel       { $$=$1; }
        ;
typespec: INTEGER {NI} kindsel        { $$=spr(0,"Integer%s",$3); }
        | REAL {NI} kindsel           { $$=spr(0,"Real%s",$3); }
        | DPREC {NI}                  { $$="Double Precision"; }
        | COMPLEX {NI} kindsel        { $$=spr(0,"Complex%s",$3); }
        | CHARACTER {NI} lensel       { $$=spr(0,"Character%s",$3); }
        | LOGICAL {NI} kindsel        { $$=spr(0,"Logical%s",$3); }
        | TYPE '(' ID ')'             { $$=spr(0,"Type(%s)",$3); }
        ;    /* wns! should be CHARACTER charsel */


kindsel: /* empty */         { $$=semp; }
       | '(' expr ')'        { $$=spr(0,"(Kind=%s)",$2); }
       | '(' ID '=' expr ')' { if( strcmp($2,"kind")!=0 )
                                 yyerror("expecting \"Kind\"\n"); 
                               $$=spr(0,"(Kind=%s)",$4);
                             }
       | '*' REALC           { $$=spr(0,"(Kind=%s)",$2);
                               for( pc=$2; *pc!='\0'; ++pc ) if( !isdigit(*pc) )
               parerr("\"Character*%s\" length-selector must be integer",$2);
                             }
       ;
lensel: /* empty */            { $$=semp; }
      | '(' charlen ')'        { $$=spr(0,"(Len=%s)",$2); }
      | '(' ID '=' charlen ')' { if( strcmp($2,"len")!=0 )
                                   yyerror("expecting \"Len\"\n"); 
                                 $$=spr(0,"(Len=%s)",$4);
                               }
      | '*' expr               { $$=spr(0,"(Len=%s)",$2); }
      | ARBCHARLEN             { $$="*(len=*)"; }
      ;
charlen: '*'  { $$="*"; }
       | expr { $$=$1; }
       ;
attrspec_l: attrspec
          | attrspec_l ',' attrspec { $$=spr(0,"%s, %s",$1,$3); }
          ;
attrspec: PARAMETER        { $$="Parameter"; }
        | PUBLIC           { $$="Public"; }
        | PRIVATE          { $$="Private"; }
        | ALLOCATABLE      { $$="Allocatable"; }
        | DIMENSION '(' {NI} arrspec {YK} ')' { $$=spr(0,"Dimension(%s)",$4); }
        | EXTERNAL         { $$="External"; }
        | accessspec
        | INTRINSIC        { $$="Intrinsic"; }
        | OPTIONAL         { $$="Optional"; }
        | POINTER          { $$="Pointer"; }
        | SAVE             { $$="Save"; }
        | TARGET           { $$="Target"; }
        ;
accessspec: INTENT '(' intentspec ')' { $$=spr(0,"Intent (%s)",$3); }
          ;
intentspec: IN { $$="IN"; } | OUT { $$="OUT"; } | INOUT { $$="INOUT"; }
          ;
arrs: expr ':' expr   { $$=spr(0,"%s:%s",$1,$3); }
    | expr ':'        { $$=spr(0,"%s:",$1); }
    | ':' expr        { $$=spr(0,":%s",$2); }
    | expr ':' '*'    { $$=spr(0,"%s:*",$1); }
    | expr            { $$=spr(0,"%s",$1); }
    | '*'             { $$=spr(0,"*"); }
    | ':'             { $$=spr(0,":"); }
    ;
arrspec: arrs
       | arrspec ',' arrs   { $$=spr(0,"%s,%s",$1,$3); }
       ;


interfaceblk: interface_s interfacespec_l endinter_s
            ;
interface_s: INTERFACE modent_o EOC
                       { strcpy(intername,$2);
		         if( !submod ){ printw(-1,"Interface");
 		           if( *$2 ) printw(0," %s",intername); printw(0,"\n");
			   indent+=dind;
			 }
			 interf=1;
		       }
           ;
endinter_s:  ENDINTER modent_o EOC
                     { if( !submod ){ indent-=dind;
		         printw(-1,"End Interface %s\n",intername);
		       }
		       if( *$2 && strcmp(intername,$2) )
			 parerr("Interface name mismatch");
		       interf=0;
		     }
           ;
interfacespec_l: interfacespec
               | interfacespec_l interfacespec
               ;
interfacespec: interbody
             | modproc_s
             ;
interbody: interfunspec
         | intersubspec
         ;
modproc_s: MODPROC procname_l EOC { printw(-1,"Module procedure %s\n",$2); }
         ;
procname_l: ID
          | procname_l ',' ID { $$=spr(0,"%s, %s",$1,$3); }
          ;
interfunspec: function_s   use_spl implicit_sl interspec_l endfun_s
            ;
intersubspec: subroutine_s use_spl implicit_sl interspec_l endsub_s
            ;
interspec_l: /* empty */
           | interspec_l interspec
           ;
interspec: itspec
/*         | specific_s */
         ;
itspec: typespec {NI}      interent_l EOC
                 { if( *$3 ) printw(-1,"%s:: %s\n",$1,$3); }
      | typespec {NI} "::" interent_l EOC
                 { if( *$4 ) printw(-1,"%s:: %s\n",$1,$4); }
      | typespec {NI} ',' {YK} attrspec_l "::" {NI} interent_l EOC
                 { if( *$8 ) printw(-1,"%s, %s:: %s\n",$1,$5,$8); }
      ;
interent_l: interent
          | interent_l ',' interent 
                { if( *$3 ){
		    if( *$1 ) $$=spr(0,"%s, %s",$1,$3);
		    else $$=$3;
		  } else $$=$1;
		}
          ;
interent: ID arrspec_o lensel_o init_o {
               if( present($1,args,", ") ) $$=spr(0,"%s%s%s%s",$1,$2,$3,$4);
	       else $$=semp;
             }
         ;



function_s: prefix FUNCTION ID '(' dummyarg_lo ')' {YK} result_o EOC {
                   check_ext(&wstat);
                   on=ispub($3);
		   funname=buffunname[++nfunname]; args=argbuf[++nargbuf];
		   strcpy(funname,$3); strcpy(args,$5);
		   if( submod && interf )
		     { printw(-1,"Interface\n"); indent+=dind; }
		   if( *$1 ) printw(-1,"%s ",$1); else printw(-1,"");
		   printw(0,"Function %s(%s)",$3,$5); strcat(args,", ");
		   if( *$8 ){ printw(0," Result(%s)",$8); strcat(args,$8); }
		   else strcat(args,$3);
		   indent+=dind; printw(0,"\n");
		 }
          ;
subroutine_s: prefix SUBROUTINE ID subargs_o EOC {
                   check_ext(&wstat);
		   on=ispub($3);
		   funname=buffunname[++nfunname]; args=argbuf[++nargbuf];
		   strcpy(funname,$3); strcpy(args,$4);
		   if( submod && interf )
		     { printw(-1,"Interface\n"); indent+=dind; }
		   if( *$1 ) printw(-1,"%s ",$1); else printw(-1,"");
		   printw(0,"Subroutine %s(%s)",$3,$4);
		   indent+=dind; printw(0,"\n");
		 }
            ;
endfun_s: endoffun
                   { indent-=dind;
		     if( *$1 && strcmp($1,funname)!=0 )
		       parerr("mismatching function name (%s<->%s)",
			      funname,$1);
		     printw(-1,"End Function %s\n",funname);
		     if( submod && interf )
		       { indent-=dind; printw(-1,"End Interface\n"); }
		     funname=buffunname[--nfunname]; args=argbuf[--nargbuf];
		     on=1; }
        ;
endsub_s: endofsub
                   { indent-=dind;
		     if( *$1 && strcmp($1,funname)!=0 )
		       parerr("mismatching subroutine name (%s<->%s)",
			      funname,$1);
		     printw(-1,"End Subroutine %s\n",funname);
		     if( submod && interf )
		       { indent-=dind; printw(-1,"End Interface\n"); }
		     funname=buffunname[--nfunname]; args=argbuf[--nargbuf];
		     on=1; }
        ;
endoffun: END { $$=semp; }
        | ENDFUN id_o EOC { $$=$2; }
        ;
endofsub: END { $$=semp; }
        | ENDSUB id_o EOC { $$=$2; }
        ;
result_o: /* empty */       { $$=semp; }
        | RESULT '(' ID ')' { $$=spr(0,"%s",$3); }
        ;
dummyarg_lo: /* empty */   { $$=semp; } 
           | dummyarg_l
           ;
subargs_o: /* empty */          { $$=semp; } 
         | '(' dummyarg_lo ')'  { $$=$2; }
         ;
dummyarg_l: dummyarg
          | dummyarg_l ',' dummyarg { $$=spr(0,"%s, %s",$1,$3); }
          ;
dummyarg: ID
        | '*' { $$="*"; }
        ;
prefix: /* empty */ { $$=semp; }
      | typespec
      | RECURSIVE   { $$=spr(0,"Recursive"); }
      | PURE        { $$=spr(0,"Pure"); }
      | ELEMENTAL   { $$=spr(0,"Elemental"); }
      ;
modsub_o: /* empty */
        | modsub
        ;
modsub: contain_s modulsub_l
      ;
modulsub_l: /* empty */ { $$=semp; }
          | modulsub_l modulsub       { hp=heap; }
          | modulsub_l modulfun       { hp=heap; }
          ;
modulfun: function_s   use_spl implicit_sl modspec_l internal_o endfun_s
        ;
modulsub: subroutine_s use_spl implicit_sl modspec_l internal_o endsub_s
        ;
fun: function_s   use_spl implicit_sl modspec_l endfun_s
   ;
sub: subroutine_s use_spl implicit_sl modspec_l endsub_s
   ;
modspec_l: /* empty */
         | modspec_l modspec
         ;
modspec: interspec
       | specific_s
       | modinter
       | entry_s        /* +++ Entry statement */
       ;
modinter: interface_s interfacespec_l endinter_s

internal: contain_s funsub_l { internal=0; }
        ;
funsub: fun
      | sub
      ;
internal_o: /* empty */
          | internal
          ;
funsub_l: funsub
        | funsub_l funsub
        ;



id_o: /* empty */ { $$=semp; }
    | ID          { $$=$1; }
    ;

dcol_o: /* empty */ { $$=semp; }
      | "::"        { $$="::"; }
      ;

/******** Dummies ********/
format_s: FORMAT;
entry_s: ENTRY;


contain_s: CONTAINS EOC
                    { if( submod ) internal=1;
  		      else {
		        indent-=dind; printw(0,"\n");
			printw(-1,"!CONTAINS\n\n"); indent+=dind;
			printw(0,"\n"); printw(-1,"Interface\n\n");
			indent+=dind; submod=cont=1;
		      }
		    }
         | /* empty */
         ;

endmod_s: ENDMOD id_o EOC { 
  if( *$2 && strcmp($2,modname)!=0 )
    parerr("end module: non-matching name (\"%s <-> %s\")",modname,$2);
  if( cont ){ indent-=dind; printw(0,"\n"); printw(-1,"End Interface\n\n"); }
  indent-=dind; printw(0,"\n"); printw(-1,"END MODULE %s\n",modname);
  cont=ext=0;
  if( indent!=0 ) parerr("problems in indentation (indent=%d)",indent);
  if( !filedep ) {
    fclose(flout); flout=NULL;
    if( chkint ) do_chkint(fnout,fnold);
  }
} ;

/************ End of grammar ************/

%%

int yyerror(char* s){
  fflush(stdout);
  fprintf(stderr,"%s: %s\n",command,s);
  fprintf(stderr,"   encountered error while reading line %d:\n   \"%s\"\n",
	  lineno,line);
  fflush(stderr);
  exit(1);
}

void init(){
  modname[0]='\0';
}



void parerr(char *meldung, ...){
  va_list args;

  va_start(args,meldung);	/* make args point to first argument of list */
  fflush(stdout);
  fprintf(stderr,"%s: ",command);
  vfprintf(stderr,meldung,args);
  fprintf(stderr,"\n");
  fflush(stderr);
  va_end(args);			/* clear-up arg-list */

  exit(1);
}


void parwarn(char *meldung, ...){
  va_list args;

  va_start(args,meldung); fflush(stdout);
  fprintf(stderr,"%s: ",command);
  vfprintf(stderr,meldung,args); fprintf(stderr,"\n"); fflush(stderr);
  va_end(args);			/* clear-up arg-list */
}


/* finds the first occurence of string pat in string in */
char* find(char* pat,char* in){
  int lin=strlen(in), lp=strlen(pat),dl=lin-lp,j;
  char *p=pat, *s;
  for( s=in; s<=in+dl; ++s ){
    for( j=0; j<lp; ++j ){ if( pat[j]!=s[j] ) break; }
    if( j==lp ) return s;
  }
  return NULL;
}


/* return pointer to first occurrence of token "tok" in list of takens "list",
   separated by delimeters in "delim"
*/
char* present(char* tok, char* list, char* delim){
  char c, *p=list;

  while( (p=find(tok,p)) ){
    if( p==list || strchr(delim,*(p-1)) ) /* check *p is beginning of token */
      if( (c=*(p+strlen(tok)))=='\0' || strchr(delim,c) ) return p;
    p+=strlen(tok);
  }
  return NULL;
}


char* spr(int m, char* format,...){
  va_list args;
  char* p=hp;
  int j;

  if( m==-1 ) m=indent;
  for( j=1; j<=m; ++j ) *hp++=' ';
  va_start(args,format);
  vsprintf(hp,format,args); hp+=strlen(hp)+1+m;
  if( hp>hpmax ) hpmax=hp;
  va_end(args);		/* clear-up arg-list */
  return p;
}


void printw(int m, char* format,...){
  va_list args;
  int j;

  if( !on ) return;
  if( m==-1 ) m=indent;
  va_start(args,format);
  for( j=1; j<=m; ++j ) fprintf(flout," ");
  vfprintf(flout,format,args);

  va_end(args);		/* clear-up arg-list */
}


int ispub(char* name){
  int on;
  slist *pub=NULL,*pri=NULL;
  if( internal ) return 0;
  if( submod && interf )	/* function is passed as an argument */
    if( present(name,args,", ") ) return 1; else return 0;
  if( public ) on=1; if( private ) on=0;
  if( submod && private && (pub=find_lele(&pub_l,name)) ) on=1;
  if( submod && public  && (pri=find_lele(&pri_l,name)) ) on=0;
  if( pub && pri )
    parerr("name \"%s\" appears in both PUBLIC and PRIVATE list",name);
  return on;
}


/* truncate string s to remove dot-separated suffix (e.g. file.for -> file) */
int delsuff(char* s){
  char *pc;
  for( pc=s+strlen(s)-1; pc!=s && *pc!='.' && *pc!='/'; --pc );
  if( pc==s || *pc!='.' ) return 0;
  *pc='\0';
  return 1;
}


int fex(char *n){
  FILE* fp;
  if( !(fp=fopen(n,"r")) ) return 0;
  fclose(fp);
  return 1;
}


void clearup(void){
  if( chkint && !ok ){
    remove(fnout); remove(fnold);
    parwarn("parsing of file %s FAILED\n",fnin);
  }
}


/* remove path from filename (e.g.: /home/users/paperino/file.for -> file.for */
char* base(char* name){
  char* p=name+strlen(name)-1;
  if( *p=='/' && p>name ) --p;
  while( p>name && *p!='/' ) --p;
  if( *p=='/' ) ++p;
  return p;
}


/* if file is external function file, open flout */
int check_ext(int* c){
  if( *c=='\0' && !filedep ){ /* external funcs & subs */
    if( flout!=NULL ) fclose(flout);
    flout=newfile(fnin,fnout,1,0);   /* do not add path */
    *c='e';
  }
}


/* get the path of file, including trailing "/" (path is "" if not specified)
   return length of path-string (0 if not specified) */
int dir_name(char *file, char *path){
  char *p, *pp;
  for( p=file+strlen(file)-1; *p!='/' && p>=file; --p );
  for( pp=file; pp<=p; ++pp,++path ) *path=*pp; *path='\0';
  return pp-file;
}


/* choose new filename; if it exists, move the old int-file 
   if moddep (i.e. !filedep), put path in front of fnout
   NB: very dirty: fnin and fnout are passed, fnold and path are global
*/
FILE *newfile(char *fnin, char *fnout, int open, int addpath){
  FILE *fl;
  if( addpath ) strcpy(fnout,path); else *fnout='\0';
  strcat(fnout,fnin); delsuff(fnout); strcpy(fnold,fnout);
  strcat(fnout,".int"); strcat(fnold,".oldint");
  remove(fnold);
  if( fex(fnout) )
    rename(fnout,fnold);
  if( open ){
    if( !(fl=fopen(fnout,"w")) )
      parerr("could not open file \"%s\" for writing\n",fnout);
    return fl;
  }
  return NULL;
}


int do_chkint(char *fnout, char *fnold){
  char cmd[512];
  char *filemod[2]={"module" /* filedep==0 */ ,"file" /* filedep==0 */  };
  int cmp;

  sprintf(cmd,"cmp -s %s %s",fnout,fnold);
  if( !fex(fnold) )
    parwarn("BUILT INTERFACE for %s %s\n",filemod[filedep],fnout);
  else if( (cmp=system(cmd))==0 ){
    parwarn("Interface of %s \"%s\" has NOT CHANGED\n",filemod[filedep],fnout);
    rename(fnold,fnout);
  } else parwarn("Interface of %s \"%s\" DID CHANGE\n",filemod[filedep],fnout);
  remove(fnold);
}


void prhelp(){
  printf("%s: create an interface for a fortran90 module\n",command);
  printf("  usage is: %s [options] <free-format-file>\n",command);
  printf("        or: %s [options] -fix <fixed-format-file>\n",command);
  printf("        or: %s [options] -autofix <any-format-file>\n",command);
  printf("  options are:\n");
  printf("    -h|-?:     print this help\n");
  printf("    -fix:      specify that following files are in fixed format\n");
  printf("    -autofix:  try to detect which files are fixed and wich are\n");
  printf("               free format\n");
  printf("    -o <file>: specify file on which to write interface\n");
  printf("    -chkint:   create interface and check if it has changed\n");
  printf("    -filedep:  (default) filename of interface is chosen based \n");
  printf("               on name of file, not name of module\n");
  printf("    -moddep:   filename of interface is chosen based on name of\n");
  printf("               module, not name of file\n");
  printf("\n");
  exit(1);
  return;
}


int main(int argc, char* argv[]){
  slist *p;

#ifdef YYDEBUG
  int yydebug=1;
#endif

  flout=stdout;
  init();			/* set defaults */
  command=base(argv[0]);

  /* check options */
  while( argc>1 && **++argv=='-' ){
    --argc;
    if( strcmp(*argv,"-?")==0 | strcmp(*argv,"-h")==0 ){ /* -?|-h print help */
      prhelp(); ++argv; --argc;
    } else if( strcmp(*argv,"-o")==0 ){       /* -o <file>: output to <file> */
      strcpy(fnout,*++argv); --argc;
    } else if( strcmp(*argv,"-fix")==0 ){     /* -fix:
                                                 parse fixed format file */
      fixed=1;
    } else if( strcmp(*argv,"-chkint")==0 ){  /* -chkint:
               handle interface files for dependencies */
      chkint=1;
    } else if( strcmp(*argv,"-filedep")==0 ){ /* -filedep:
               dependencies based on files (one .int file for each *.f) */
      filedep=1;
    } else if( strcmp(*argv,"-moddep")==0 ){  /* -moddep: (default)
               dependencies based on modules (one .int file for each module,
               + one for the program)          */
      filedep=0;
    } else parerr("illegal option (\"%s\")\n",*argv);
  }
  if( argc>2 ) parerr("specify at most one argument\n");
  if( chkint && argc!=2 )
    parerr("with the \"-chkint\" option you MUST specify the input file\n");
  if( chkint && *fnout )
    parerr("with the \"-chkint\" option you CANNOT specify the output file\n");
  
  if( argc==2 ){ strcpy(fnin,*argv); dir_name(fnin,path); }

  if( chkint ){
    atexit(clearup);
    if( filedep ) flout=newfile(fnin,fnout,1,0);
  }

  if( filedep && *fnout && !(flout=fopen(fnout,"w")) )
    parerr("could not open file \"%s\" for writing\n",fnout);
  if( *fnin && !(yyin=fopen(fnin,"r")) )
      parerr("could not open file \"%s\" for reading\n",fnin);



  yyparse(); ok=1;



  if( *fnout && flout ) { fclose(flout); flout=NULL; }
  if( chkint && filedep ) do_chkint(fnout,fnold);


#ifdef YYDEBUG
  parwarn("MAX OCCUPANCY OF:\n");
  parwarn("   stack: %d bytes (=%.1f kb = %.1f Mb)\n",
	 (spmax-stack),(spmax-stack)/1024.,(spmax-stack)/1024./1024.);
  parwarn("   heap:  %d bytes (=%.1f kb = %.1f Mb)\n",
	 (hpmax-heap),(hpmax-heap)/1024.,(hpmax-heap)/1024./1024.);
#endif

/*  fprintf(stderr,"read %d lines\n",lineno); */

  return 0;			/* everything ok */
}
