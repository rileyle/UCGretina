# $Id: GNUmakefile,v 1.1 1999/01/07 16:05:40 gunter Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

# Need this for writing output files larger than 2GB on 32-bit Linux
# (untested on other systems)
CPPFLAGS := -D_FILE_OFFSET_BITS=64

# Use -D to define LHTARGET macro for the C preprocesssor
ifdef LHTARGET
  name := UCGretina_LH
  CPPFLAGS += -DLHTARGET
else
ifdef SCANNING
  name := UCGretina_Scan
  CPPFLAGS  += -DSCANNING
  CPPFLAGS  += -I/usr/local/CADMesh/include
  EXTRALIBS += -L/usr/local/CADMesh/lib -lcadmesh -lassimp -ltet
else
  name := UCGretina
endif
endif

ifdef HIGHMULT
  CPPFLAGS  += -DHIGHMULT
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

############################################
## uncomment for geant v. 4.9.6
###########################################
CPPFLAGS += -DG4V496

############################################
## uncomment for geant v. 4.10
###########################################
#CPPFLAGS += -DG4V10

G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

.PHONY: all
all: lib bin

include $(G4INSTALL)/config/binmake.gmk
