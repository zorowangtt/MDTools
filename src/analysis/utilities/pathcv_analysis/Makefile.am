#--------1---------2---------3---------4---------5---------6---------7---------8
#
#  File     : Makefile.am
#  Authors  : NT
#
#  (c) Copyright 2014 RIKEN. All rights reserved.
#
#--------1---------2---------3---------4---------5---------6---------7---------8

SUBDIRS := ../../libana .

bin_PROGRAMS = pathcv_analysis

pathcv_analysis_SOURCES = \
        pn_option_str.fpp \
        pn_option.fpp \
        pn_control.fpp \
        pn_setup.fpp \
        pn_analyze.fpp \
        pn_main.fpp

EXTRA_DIST = \
	Makefile.depends 

pathcv_analysis_ADDLIB = ../../../lib/lib.a
pathcv_analysis_ADDLIBANA = ../../libana/libana.a

SUFFIXES = .fpp
DEPENDS = Makefile.depends
MAINTAINERCLEANFILES = Makefile.in

CLEANFILES = $(bin_PROGRAMS) $(pathcv_analysis_OBJECTS) *.f90 *~ *.mod

.fpp.o:
if USEKCOMP
	cp $*.fpp $*.f90
	$(FPP) $(PPFLAGS) $(DEFS) $*.f90
	mv $*.cpp.f90 $*.f90
	$(FC)  $(DEFAULT_INCLUDES) $(INCLUDES) \
		-I../../../lib -I../../libana $(FCFLAGS) -c $*.f90
else
	$(FPP) $(PPFLAGS) $(DEFS) $*.fpp $*.f90
	$(FC)  $(DEFAULT_INCLUDES) $(INCLUDES) \
		-I../../../lib -I../../libana $(FCFLAGS) -c $*.f90
endif

pathcv_analysis$(EXEEXT) : $(pathcv_analysis_OBJECTS) \
		$(pathcv_analysis_ADDLIBANA) $(pathcv_analysis_ADDLIB)
	$(FC) -o $(bin_PROGRAMS) $(pathcv_analysis_OBJECTS) \
		$(pathcv_analysis_ADDLIBANA) $(pathcv_analysis_ADDLIB) $(LDFLAGS)

depend: clean_depend
	python ../../../../fortdep.py *.fpp > $(DEPENDS)

clean_depend:
	rm -f $(DEPENDS)

-include $(DEPENDS)
