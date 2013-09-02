# Makefile for XTC analysis
#
# possible DEFINES are: -DDEBUG (for extra checks and extra output)
#                       -DXTC   (if xdrf/xtc routines are available)
#                       -DNOT32 (if OS is not 32 bit)

SHELL    = /bin/sh
CC       = g++
RM       = rm 
INCLUDES = -I/usr/include/xdrfile -I/use/include -I/usr/lib64
# INCLUDES = -I/usr/local/include/xdrfile -I/usr/local/lib64
LIBS     = -lm -lfftw3 
XTCLIB   = -lxdrfile
OPTFLAGS = 
DEFINES  = -DCPLUSPLUS
CFLAGS   = $(INCLUDES) $(LIBS) $(XTCLIB) $(OPTFLAGS) $(DEFINES) 

mainfile = main.cpp
distobj = sim_io.cpp grid.cpp random.cpp sim_libs.cpp sim_run.cpp \
	sim_field.cpp  \
	sim_surface93.cpp \
	sim_harmonicpin.cpp \
	sim_samsurface.cpp

interface: $(distobj)
	$(CC) $(CFLAGS) $(mainfile) $(distobj) -o makewaves
clean: 		
	$(RM) makewaves

