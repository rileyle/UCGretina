# $Id: GNUmakefile,v 1.1 1999/01/07 16:05:40 gunter Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

# Use -D to define LHTARGET macro for the C preprocesssor
ifdef LHTARGET
  name := UCGretina_LH
  CPPFLAGS := "-DLHTARGET "
else
  name := UCGretina
endif

############################################
## uncomment for geant v. 4.7.1 and above
###########################################
CPPFLAGS += -DG4V47

############################################
## uncomment for geant v. 4.8.0 and above
###########################################
CPPFLAGS += -DG4V48

############################################
## uncomment for geant v. 4.9.0
###########################################
CPPFLAGS += -DG4V49


G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

.PHONY: all
all: lib bin

include $(G4INSTALL)/config/binmake.gmk
