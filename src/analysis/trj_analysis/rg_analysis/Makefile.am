#--------1---------2---------3---------4---------5---------6---------7---------8
#
#  File     : Makefile.am
#  Authors  : NT
#
#  (c) Copyright 2014 RIKEN. All rights reserved.
#
#--------1---------2---------3---------4---------5---------6---------7---------8

SUBDIRS := ../../libana .

bin_PROGRAMS = rg_analysis

rg_analysis_SOURCES = \
	rg_option_str.fpp \
	rg_option.fpp \
	rg_control.fpp \
	rg_setup.fpp \
	rg_analyze.fpp \
	rg_main.fpp

EXTRA_DIST = \
	Makefile.depends 

rg_analysis_ADDLIB = ../../../lib/lib.a
rg_analysis_ADDLIBANA = ../../libana/libana.a

SUFFIXES = .fpp
DEPENDS = Makefile.depends
MAINTAINERCLEANFILES = Makefile.in

CLEANFILES = $(bin_PROGRAMS) $(rg_analysis_OBJECTS) *.f90 *~ *.mod

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

rg_analysis$(EXEEXT) : $(rg_analysis_OBJECTS) \
		$(rg_analysis_ADDLIBANA) $(rg_analysis_ADDLIB)
	$(FC) -o $(bin_PROGRAMS) $(rg_analysis_OBJECTS) \
		$(rg_analysis_ADDLIBANA) $(rg_analysis_ADDLIB) $(LDFLAGS)

depend: clean_depend
	python ../../../../fortdep.py *.fpp > $(DEPENDS)

clean_depend:
	rm -f $(DEPENDS)

-include $(DEPENDS)
