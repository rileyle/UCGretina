# $Id: GNUmakefile,v 1.1 1999/01/07 16:05:40 gunter Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

# Need this for writing output files larger than 2GB on 32-bit Linux
# (untested on other systems)
CPPFLAGS := -D_FILE_OFFSET_BITS=64

# Use -D to define LHTARGET, AD, SCANNING, NEUTRONS macros
# for the C preprocesssor
ifdef NEUTRONS
  CPPFLAGS += -DNEUTRONS
endif
ifdef LHTARGET
ifdef POL
  name := UCGretina_LH_Pol
  CPPFLAGS += -DLHTARGET
  CPPFLAGS += -DPOL
else
  name := UCGretina_LH
  CPPFLAGS += -DLHTARGET
endif
else # end LHTARGET / begin SCANNING
ifdef SCANNING
ifdef POL
  name := UCGretina_Scan_Pol
  CPPFLAGS  += -DSCANNING
  CPPFLAGS += -DPOL
else
  name := UCGretina_Scan
  CPPFLAGS  += -DSCANNING
endif
else # end SCANNING / begin UCGretina
ifdef POL
  name := UCGretina_Pol
  CPPFLAGS += -DPOL
else
  name := UCGretina
endif
endif
endif

ifdef HIGHMULT
  CPPFLAGS  += -DHIGHMULT
endif

############################################
## uncomment for geant v. 4.7.1 and above
###########################################
#CPPFLAGS += -DG4V47

############################################
## uncomment for geant v. 4.8.0 and above
###########################################
#CPPFLAGS += -DG4V48

############################################
## uncomment for geant v. 4.9.0
###########################################
#CPPFLAGS += -DG4V49

############################################
## uncomment for geant v. 4.9.6
###########################################
#CPPFLAGS += -DG4V496

############################################
## uncomment for geant v. 4.10
###########################################
#CPPFLAGS += -DG4V10

############################################
## uncomment for geant v. 4.11
###########################################
CPPFLAGS += -DG4V11

############################################
## uncomment for debugging
###########################################
#CPPFLAGS += -ggdb  

G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

.PHONY: all
all: lib bin

include $(G4INSTALL)/config/binmake.gmk
