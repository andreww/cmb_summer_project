#
# Makefile for programs
# SccsID	=	@(#)Makefile	11.9	1/14/91
#

# *** FINAL PRODUCT MACROS ***
# The library name should be in the form libxxx.a that it can be 
# referenced as -lxxx (example:  -lpost)

PROG1		= remodl
PROG2		= setbrn
PROG3		= ttim


# *** DESTINATION DIRECTORY MACROS ***

DESTDIR		= 
DESTMAN		= 

# *** COMPONENT MACROS ***
# SRCS, OBJS: only those sources/objects making up the library
# (not including main program)

CSRCS1		= remodl.c
CSRCS2		= setbrn.c
CSRCS3		= ttim.c

OBJS1		= ${CSRCS1:%.c=%.o}
OBJS2		= ${CSRCS2:%.c=%.o}
OBJS3		= ${CSRCS3:%.c=%.o}

MANS		= 

HDRS		=

# *** AUXILIARY LIBRARY MACROS ***
# Libraries used for linking

STDLIBS		= -lm 

LIBS		= -ltau 

O		= -g

INCLUDES	= -I../include

LIBSEARCH	= -L../lib

# Intermediate-stage compiler flags
# -u REQUIRES explicitly defined variables in Fortran programs

CFLAGS		= $O $(INCLUDES)
FFLAGS		= -u $O
PCC		= pcc
PCCINC1		= $(ORACLE_HOME)/c/lib
PCCINC2		= ../../include
PCCFLAGS = include=$(PCCINC1) include=$(PCCINC2) ireclen=132 ltype=none oreclen=132 rebind=yes select_error=no


# Linking compiler flags

LDFLAGS		= $O $(LIBSEARCH)

# Lint preprocessor flags

LFLAGS		= -bchu $(INCLUDES)

# *** COMPILER MACROS ***
# use LD = $(F77) if main program is Fortran
# use LD = $(CC) if main program is C

CC		= cc
F77		= f77
LD		= $(CC)
LINT		= lint

# *** PERMISSION MACROS ***
# Permission bits for installed software: 
# LMODE = libraries, manual pages (rw_rw_r__)
# PMODE = programs (rwxrwxr_x)

LMODE		= 664
PMODE		= 775

# *** SCCS MACROS ***

GET		= sccs get
GFLAGS		=
SCLEAN		= sccs clean

#----------------------------------------------------------------------
# *** BEGIN TARGETS ***
#----------------------------------------------------------------------

# Default target (target executed if just "make" is typed) is all

all:		$(PROG1) $(PROG2) $(PROG3)

$(PROG1):	$(OBJS1) $(HDRS)
		$(LD) $(LDFLAGS) -o $@ $(OBJS1) $(LIBRARY) $(LIBS) $(STDLIBS)
		chmod +w $@

$(PROG2):	$(OBJS2) $(HDRS)
		$(LD) $(LDFLAGS) -o $@ $(OBJS2) $(LIBRARY) $(LIBS) $(STDLIBS)
		chmod +w $@

$(PROG3):	$(OBJS3) $(HDRS)
		$(LD) $(LDFLAGS) -o $@ $(OBJS3) $(LIBRARY) $(LIBS) $(STDLIBS)
		chmod +w $@

# *** INSTALLATION ***
# (1) Update program and manuals
# (2) Install program in destination directory
# (3) Strip program
# (4) Remove program from source directory

install:	all $(MANS)
		install -m $(PMODE) $(PROG1) $(DESTDIR)
		-rm -f $(PROG1)

# *** DIRECTORY MAINTENANCE ***
# clean: all .o files
# xclean: all .o files, executables, and SCCS-controlled source code
# 	not currently being edited

clean:
		-rm -f $(OBJS1) $(OBJS2) $(OBJS3) core

xclean:		clean
		-rm -f $(PROG1) $(PROG2) $(PROG3)
		-$(SCLEAN)

# *** MANUAL MAINTENANCE ***
# (1) Compare manuals in this directory ($$(@F)) w/ those in DESTDIR $(MANS)
# (2) Retrieve from SCCS
# (3) Copy to DESTDIR
# (4) Change permission of DESTDIR to allow copying in the future

manual:		$(MANS)

$(MANS):	$$(@F)
		sccs get $(@F)
		-rm -f $@
		cp $(@F) $@
		chmod $(LMODE) $@

lint:		$(CSRCS) $(PCSRCS) $(HDRS)
		$(LINT) $(LFLAGS) $(CSRCS) $(PCSRCS)

sccs:
		$(GET) $(GFLAGS) $(CSRCS) $(PCSRCS) Makefile


%.o:    	%.c
		-rm -f $@
		$(CC) -o $@ $(CFLAGS) -c $*.c
		chmod 664 $@
%.o:		%.pc
		-rm -f $@
		$(PCC) $(PCCFLAGS) iname=$*.pc
		$(CC) -o $@ $(CFLAGS) -c $*.c
		chmod 664 $@
%.o:		%.f
		-rm -f $@
		$(F77) -o $@ $(FFLAGS) -c $*.f
		chmod 664 $@

%.uid:		%.uil
		uil -I${XHOME}/include/uil -I${XHOME}/include/Uil $*.uil -o $@
