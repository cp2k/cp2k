BISON=bison
CFLAGS=-c -g -O

forlex.c: 	forlex.lex forpar.tab.h
		flex -i -oforlex.c forlex.lex

forpar.tab.c:	forpar.y forpar.h
		$(BISON) -d -v forpar.y

forpar.tab.h:	forpar.y
		$(BISON) -d -v forpar.y

clean:
	rm -f *.o *.x *.tab.* forpar.output forlex.c
