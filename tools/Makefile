.SUFFIXES: .o .c

### Machine dependent ###

CC	= cc
LDR	= cc

CDEFS	= -Dxxx
OPT	= -O
BFLAGS	= 
CFLAGS	= $(BFLAGS) $(OPT) $(CDEFS)
LFLAGS	= $(BFLAGS) $(OPT)

LIBS	= 
OBJECTS	= forlex.o  forpar.tab.o  slist.o

#################################
vpath %.c $(SRCDIR)

PROG	= forpar.x


all: $(PROG)

$(PROG): $(OBJECTS)
	$(LDR) $(LFLAGS) -o $(PROG) $(OBJECTS) $(LIBS)

%.o: %.c
	$(CC) -c -I$(SRCDIR) $(CFLAGS) $<

clean:
	rm -f *.o

realclean:
	rm -f $(PROG) *.o *~
