#--------1---------2---------3---------4---------5---------6---------7---------8
#
#  File     : Makefile.am
#  Authors  : NT
#
#  (c) Copyright 2014 RIKEN. All rights reserved.
#
#--------1---------2---------3---------4---------5---------6---------7---------8

SUBDIRS := ../../libana .

bin_PROGRAMS = wham_analysis

wham_analysis_SOURCES = \
        wa_option_str.fpp \
        wa_option.fpp \
        wa_control.fpp \
        wa_setup.fpp \
        wa_analyze.fpp \
        wa_main.fpp

EXTRA_DIST = \
	Makefile.depends 

wham_analysis_ADDLIB = ../../../lib/lib.a
wham_analysis_ADDLIBANA = ../../libana/libana.a

SUFFIXES = .fpp
DEPENDS = Makefile.depends
MAINTAINERCLEANFILES = Makefile.in

CLEANFILES = $(bin_PROGRAMS) $(wham_analysis_OBJECTS) *.f90 *~ *.mod

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

wham_analysis$(EXEEXT) : $(wham_analysis_ADDLIBANA) $(wham_analysis_ADDLIB) \
							$(wham_analysis_OBJECTS)
	$(FC) -o $(bin_PROGRAMS) $(wham_analysis_OBJECTS) \
		$(wham_analysis_ADDLIBANA) $(wham_analysis_ADDLIB) $(LDFLAGS)

depend: clean_depend
	python ../../../../fortdep.py *.fpp > $(DEPENDS)

clean_depend:
	rm -f $(DEPENDS)

-include $(DEPENDS)
