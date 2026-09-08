
# Makefile for GeoClaw Elliott Bay debris scenario.
CLAWMAKE = $(CLAW)/clawutil/src/Makefile.common

CLAW_PKG = geoclaw
EXE = xgeoclaw
SETRUN_FILE = setrun.py
OUTDIR = _output
SETPLOT_FILE = setplot.py
PLOTDIR = _plots

FFLAGS ?=

GEOLIB = $(CLAW)/geoclaw/src/2d/shallow
include $(GEOLIB)/Makefile.geoclaw

EXCLUDE_MODULES = \

EXCLUDE_SOURCES = \

MODULES = \

SOURCES = \
  $(CLAW)/riemann/src/rpn2_geoclaw.f \
  $(CLAW)/riemann/src/rpt2_geoclaw.f \
  $(CLAW)/riemann/src/geoclaw_riemann_utils.f \

include $(CLAWMAKE)

.PHONY: topo all
topo:
	python fetch_topo.py

all:
	$(MAKE) topo
	$(MAKE) .output
