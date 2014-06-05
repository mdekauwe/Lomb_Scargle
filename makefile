##############################################################################
#HOME     = /users/eow/mgdk
HOME     = /Users/mdekauwe
CFLAGS   = -g -Wall
#CFLAGS   = -O3 -Wall
FFLAGS   = -g
#ARCH     =  sun4u
ARCH     =  i386
INCLS    = -I./include -I$(HOME)/include/
LIBS     = -L$(HOME)/lib/ -lm -lmy_c_lib -lnrc
#INCLS    = -I./include -I$(HOME)/include/$(ARCH)
#LIBS     = -L$(HOME)/lib/$(ARCH) -lm -lmy_c_lib -lnrc
CC       =  gcc
PROGRAM  =  lombscargle

SOURCES  =  \
$(PROGRAM).c fft.c estimate_tau.c read_files.c \
parsar_and_usage.c trend_and_taper.c signif.c  \
outputs.c

OBJECTS = $(SOURCES:.c=.o)

MAN_PAGE =  man/man1/lombscargle.1
RM       =  rm -f
##############################################################################

# top level create the program...
all: 		$(PROGRAM)

# Compile the src file...
$(OBJECTS):	$(SOURCES)
		$(CC) ${INCLS} $(CFLAGS) -c $(SOURCES)

# Linking the program...
$(PROGRAM):	$(OBJECTS)
		$(CC) $(OBJECTS) $(LIBS) ${INCLS} $(CFLAGS) -o $(PROGRAM) 

# clean up...		
clean:
		$(RM) $(OBJECTS) $(PROGRAM)

install:
		cp $(PROGRAM) $(HOME)/bin/$(ARCH)/.
		cp $(MAN_PAGE) $(HOME)/man/man1/.	
		$(RM) $(OBJECTS)
##############################################################################
