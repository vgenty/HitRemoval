#
# This is an example GNUmakefile for my packages
#

# specific names for this package
DICT  = HitRemoval_NuTrackRemovalCint
SHLIB = libHitRemoval_NuTrackRemoval.so
SOURCES = $(filter-out $(DICT).cxx, $(wildcard *.cxx))
FMWK_HEADERS = LinkDef.h $(DICT).h
HEADERS = $(filter-out $(FMWK_HEADERS), $(wildcard *.h))
OBJECTS = $(SOURCES:.cxx=.o)

# include options for this package
INCFLAGS  = -I.                       #Include itself
INCFLAGS += $(shell larlite-config --includes)
INCFLAGS += $(shell basictool-config --includes)
INCFLAGS += -I$(LARLITE_USERDEVDIR)/HitRemoval/Base

# platform-specific options
OSNAME          = $(shell uname -s)
HOST            = $(shell uname -n)
OSNAMEMODE      = $(OSNAME)

# call kernel specific compiler setup
include $(LARLITE_BASEDIR)/Makefile/Makefile.${OSNAME}

# call the common GNUmakefile
LDFLAGS += $(shell larlite-config --libs)
LDFLAGS += $(shell basictool-config --libs)
LDFLAGS += -L$(LARLITE_LIBDIR) -lHitRemoval_Base
include $(LARLITE_BASEDIR)/Makefile/GNUmakefile.CORE


