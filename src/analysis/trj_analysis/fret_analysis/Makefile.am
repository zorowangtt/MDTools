#--------1---------2---------3---------4---------5---------6---------7---------8
#
#  File     : Makefile.am
#  Authors  : NT
#
#  (c) Copyright 2014 RIKEN. All rights reserved.
#
#--------1---------2---------3---------4---------5---------6---------7---------8

SUBDIRS := ../../libana .

bin_PROGRAMS = fret_analysis

fret_analysis_SOURCES = \
	fa_option_str.fpp \
	fa_option.fpp \
	fa_control.fpp \
	fa_setup.fpp \
	fa_analyze.fpp \
	fa_main.fpp

EXTRA_DIST = \
	Makefile.depends 

fret_analysis_ADDLIB = ../../../lib/lib.a
fret_analysis_ADDLIBANA = ../../libana/libana.a

SUFFIXES = .fpp
DEPENDS = Makefile.depends
MAINTAINERCLEANFILES = Makefile.in

CLEANFILES = $(bin_PROGRAMS) $(fret_analysis_OBJECTS) *.f90 *~ *.mod

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

fret_analysis$(EXEEXT) : $(fret_analysis_OBJECTS) \
		$(fret_analysis_ADDLIBANA) $(fret_analysis_ADDLIB)
	$(FC) -o $(bin_PROGRAMS) $(fret_analysis_OBJECTS) \
		$(fret_analysis_ADDLIBANA) $(fret_analysis_ADDLIB) $(LDFLAGS)

depend: clean_depend
	python ../../../../fortdep.py *.fpp > $(DEPENDS)

clean_depend:
	rm -f $(DEPENDS)

-include $(DEPENDS)
