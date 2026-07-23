# $Id: GNUmakefile,v 1.1 1999/01/07 16:05:40 gunter Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

# Need this for writing output files larger than 2GB on 32-bit Linux
# (untested on other systems)
CPPFLAGS := -D_FILE_OFFSET_BITS=64 -g

# Use -D to define LHTARGET, AD, SCANNING, NEUTRONS macros
# for the C preprocesssor
ifdef NEUTRONS
  CPPFLAGS += -DNEUTRONS
endif
ifdef LHTARGET
ifdef POL
  name := UCGretina_LH_Pol
  CPPFLAGS += -DLHTARGET -DPOL
else
  name := UCGretina_LH
  CPPFLAGS += -DLHTARGET
endif
else # end LHTARGET / begin SCANNING
ifdef SCANNING
ifdef POL
  name := UCGretina_Scan_Pol
  CPPFLAGS  += -DSCANNING -DPOL
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

ifdef CACHETEXT
  CPPFLAGS += -DCACHETEXT
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
CPPFLAGS += -DG4V10

G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

# Collect the git branch and commit hash.
GIT_HASH := $(shell git describe --always --abbrev=6 --exclude '*')
GIT_BRANCH := $(shell git rev-parse --abbrev-ref HEAD)
ifneq ("$(wildcard git_hash)","")
PREV_GIT_HASH := $(shell cat git_hash)
else
PREV_GIT_HASH := ""
endif

.PHONY: all
all: include/Git_Hash.hh lib bin

include/Git_Hash.hh: FORCE
ifneq ($(PREV_GIT_HASH),$(GIT_HASH))
	echo $(GIT_HASH) > git_hash
	echo "#ifndef Git_Hash_h" > include/Git_Hash.hh
	echo "#define Git_Hash_h" >> include/Git_Hash.hh
	echo "#define GIT_HASH \"$(GIT_HASH)\"" >> include/Git_Hash.hh
	echo "#define GIT_BRANCH \"$(GIT_BRANCH)\"" >> include/Git_Hash.hh
	echo "#endif" >> include/Git_Hash.hh
endif

FORCE:

.PHONY: test-smoke test-sources test-inbeam test-scanning test-background test-benchmark test test-functional test-baselines

test-smoke:
	python3 tests/benchmark.py --mode smoke

test-sources:
	python3 tests/benchmark.py --mode sources

test-inbeam:
	python3 tests/benchmark.py --mode inbeam

test-scanning:
	python3 tests/benchmark.py --mode scanning

test-background:
	python3 tests/benchmark.py --mode background

test-functional: test-sources test-inbeam test-background test-scanning

test-baselines:
	python3 tests/benchmark.py --update-baselines --events $(or $(BENCH_EVENTS),1000000)

include $(G4INSTALL)/config/binmake.gmk
