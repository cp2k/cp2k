#
# should be called automatically by the Makefile in cp2k/makefiles
# but it had better work if called manually
#

.SUFFIXES: .o .c
OBJECTS	= forlex.o  forpar.tab.o  slist.o

ifeq "$(SRCDIR)" ""
SRCDIR = .
endif
vpath %.c $(SRCDIR)

PROG	= forpar.x


all: $(PROG)

$(PROG): $(OBJECTS)
	$(CC) -O -o $(PROG) $(OBJECTS) 

%.o: %.c
	$(CC) -O -Dxxx -c -I$(SRCDIR) $<

clean:
	rm -f *.o

realclean:
	rm -f $(PROG) *.o *~
