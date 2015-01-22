# Makefile for the ROOT test programs.
# This Makefile shows nicely how to compile and link applications
# using the ROOT libraries on all supported platforms.
#
# Copyright (c) 2000 Rene Brun and Fons Rademakers
#
# Author: Fons Rademakers, 29/2/2000

include $(ROOTSYS)/etc/Makefile.arch
-include $(ROOTSYS)/MyConfig.mk
HdrSuf       = h
SOFLAGS       = #-----------------------------------------------------------------------------
# special programs
#-----------------------------------------------------------------------------
# special root
AVLOCROOTO   =  src/avlocroot.$(ObjSuf)
AVLOCROOT    =  bin/avlocroot$(ExeSuf)
# make summary ntuple
MAKENTUPLEO  =  src/make_ntuple.$(ObjSuf)
MAKENTUPLE   =  bin/make_ntuple$(ExeSuf)
# make a bunch of plots
MAKEPLOTSO   =  src/make_plots.$(ObjSuf)
MAKEPLOTS    =  bin/make_plots$(ExeSuf)
# make a quick flatmap
#FLATMAPO     =  src/plot_flatmap.$(ObjSuf)
#FLATMAP      =  bin/plot_flatmap$(ExeSuf)
#-----------------------------------------------------------------------------
# my object files for the library
#-----------------------------------------------------------------------------
AVLOCOBJS    =  src/AVLocTools.$(ObjSuf) src/AVLocProc.$(ObjSuf) \
		src/AVLocPlot.$(ObjSuf)
AVLOCHDRS    =  include/AVLocTools.$(HdrSuf) include/AVLocBasicProc.$(HdrSuf) \
		src/AVLocPlot.$(HdrSuf) 
AVLOCLIB     =  lib/libAVLoc.$(DllSuf)

#-----------------------------------------------------------------------------
# libraries to be included
#-----------------------------------------------------------------------------
ALLLIBS      =  $(LIBS) -L$(ROOTSYS)/lib -lMinuit -L$(RATROOT)/lib  -lRATEvent_Darwin
NEXTLIBS = -Wl,-rpath,$(CURDIR)/lib/ -L$(CURDIR)/lib/ -lAVLoc
#-----------------------------------------------------------------------------
.SUFFIXES: .$(SrcSuf) .$(ObjSuf) .$(DllSuf)
.PHONY:     READ

all:	 $(AVLOCSO) $(AVLOCLIB) $(MAKENTUPLE) $(MAKEPLOTS) $(AVLOCROOT)

clean:
		@rm -f $(AVLOCOBJS) core*  $(AVLOCLIB) $(AVLOCO) \
		$(AVLOCROOTO) $(MAKENTUPLEO) $(MAKEPLOTSO) $(FLATMAPO) \
		$(AVLOCROOT) $(MAKENTUPLE) $(MAKEPLOTS) $(FLATMAP)

$(AVLOCROOT):	$(AVLOCROOTO) $(AVLOCLIB)
		$(LD) $(LDFLAGS) $(ALLLIBS) $(NEXTLIBS) $(AVLOCROOTO) $(AVLOCLIB) \
		$(OutPutOpt) $@
		$(MT_EXE)
		@echo "$@ done"

$(MAKENTUPLE):	$(AVLOCLIB) $(MAKENTUPLEO)
			$(LD) $(LDFLAGS) $(MAKENTUPLEO) $(AVLOCLIB) $(ALLLIBS)  \
			$(OutPutOpt)$@
			$(MT_EXE)
			@echo "$@ done"

$(MAKEPLOTS):	$(AVLOCLIB) $(MAKEPLOTSO)
			$(LD) $(LDFLAGS) $(MAKEPLOTSO) $(AVLOCLIB) $(ALLLIBS) $(NEXTLIBS)\
			$(OutPutOpt)$@
			$(MT_EXE)
			@echo "$@ done"

$(FLATMAP):	$(AVLOCLIB) $(FLATMAPO)
			$(LD) $(LDFLAGS) $(FLATMAPO) $(AVLOCLIB) $(ALLLIBS) $(NEXTLIBS)\
			$(OutPutOpt)$@
			$(MT_EXE)
			@echo "$@ done"

$(AVLOCLIB):	$(AVLOCOBJS)
		$(LD)  -dynamiclib -single_module -install_name $(CURDIR)/$@ $(ALLLIBS) $^ $(OutPutOpt) $@
		@echo "$@ done"
		@echo "libs $(SOFLAGS)" 

.$(SrcSuf).$(ObjSuf):
	$(CXX)  $(CXXFLAGS) -I. -I$(RATROOT)/include -c $< -o $@








