
CC  	:= gcc
CPP	:= g++

CFLAGS 	:= -arch x86_64 -O2 -Wall 

# —sysroot=/Developer/SDKs/MacOSX10.6.sdk/

# these are machine-specific… change these two and should run ok on any machine
USRNAME := morroww
MACHINE := wrmorrow-MBP

# MAKE COMMANDS ######################################################################

INCLUDE  = -I./inc/ -I/usr/local/lib/ -I/usr/local/include/

SRCDIR   = .
OBJDIR   = .
EXCDIR   = .

DATAOBJS = $(OBJDIR)/mnl_data.o
		
RANLIBS  = -lunuran -lrngstreams 
BLASLIBS = -framework Accelerate 

# SNOPT SHORTCUTS ######################################################################

SNOPT_LIBDIR	:= /usr/local/lib
SNOPT_LIBS  	:= -lsnopt_c -lsnprint_c -lblas_c
SNOPT_AR     	:= libsnopt_c libsnprint_c libblas_c libsnoextras
SNOPT_AR_LIBS 	:= $(SNOPT_AR:%=$(SNOPT_LIBDIR)/%.a)
SNOPT_INC  	:= -I/usr/local/include/snopt 
SNOPT_OBJS 	:= $(OBJDIR)/snfilewrapper.o $(OBJDIR)/mnl_snopt.o

F2C_LIBDIR	:= /usr/local/lib
F2C_INCLUDE	:= /usr/local/include/snopt 

# SPQR SHORTCUTS ######################################################################

SPQR_LIB_DIR = /usr/local/lib
SPQR_LIBS  := -lcholmod -lamd -lcolamd -lcamd -lccolamd -lmetis -lspqr -lma57 -lgfortran

# MAKE COMMANDS ######################################################################

default: data

data:

	$(CC) $(CFLAGS) -c $(SRCDIR)/mnl_data.c -o $(OBJDIR)/mnl_data.o $(INCLUDE)

mnl_snopt: data snopt

	$(CC) $(CFLAGS) -c $(SRCDIR)/mnl_snopt.c -o $(OBJDIR)/mnl_snopt.o \
		$(INCLUDE) $(SNOPT_INC) 

test: snopt mnl_snopt
	
	$(CPP) $(CFLAGS) $(SRCDIR)/mnl_test.c -o $(EXCDIR)/test \
		$(INCLUDE) $(SNOPT_INC) \
		$(DATAOBJS) $(SNOPT_OBJS) $(SPQR_LIBS) $(RANLIBS) $(BLASLIBS) \
		-L$(SNOPT_LIBDIR) $(SNOPT_AR_LIBS) $(F2C_LIBDIR)/libf2c.a -lm

# SNOPT ############################################################################

# Fake target to remind people to set the F2C environment variable with SNOPT

$(F2CINCLUDE)/f2c.h:
	@echo "Could not find the f2c distribution."
	@echo "Set the following environment variables:"
	@echo "  F2CINCLUDE should be the path to f2c.h"
	@false

snfilewrapper.o:

	$(CC) $(CFLAGS) -c $(SRCDIR)/snfilewrapper.c -o $(OBJDIR)/snfilewrapper.o \
		$(INCLUDE) $(SNOPT_INC)

snopt: $(F2CINCLUDE)/f2c.h snfilewrapper.o 

