# sample taken from
# INTEL(R) ADVISOR XE 2013.
# Copyright (C) 2009-2011 Intel Corporation. All rights reserved
# all additions are 
# DANIEL(R) JIMENEZ GO 1989.
# Copyleft (C) 2014 CICESE Quantum Optics lab. Not for distribution.
NAME:=correlation_g2

Arch=intel64
libdir=lib64
FCFLAGS+= -m64

FC = ifort
ADV_DIR = $(ADVISOR_XE_2013_DIR)

ANNOTATE_FLAGS = -I $(ADV_DIR)/include/$(Arch) -L$(ADV_DIR)/$(libdir) -ladvisor -ldl
MKL_VSL_FLAGS = -I$(MKLROOT)/include -i8 -L$(MKLROOT)/lib/intel64 -lmkl_intel_ilp64 -lmkl_core -lmkl_sequential -lpthread -lm

COHERENCE_DEPENDENCIES = source/simulation_driver.f90 source/generate_random_photons.f90 source/rotate_to_aperture_angle.f90 source/evaluate_wavefunction.f90

all: build
release: build
debug: build_debug

build: coherence
build_debug: rangen_debug


coherence: $(COHERENCE_DEPENDENCIES)
	$(FC) $(FCFLAGS) $(COHERENCE_DEPENDENCIES) -o coherence -r8 -O2 -xHost $(MKL_VSL_FLAGS)

coherence.o: source/simulation_driver.f90
	$(FC) $(FCFLAGS) source/simulation_driver.f90 -c -r8 -O2 -xHost $(MKL_VSL_FLAGS)
#rangen: source/rangen.f90
#	$(FC) $(FCFLAGS) source/rangen.f90 -o rangen -r8 -O2 $(MKL_VSL_FLAGS) #-g #$(ANNOTATE_FLAGS)

#rangen_debug : source/rangen.f90
#	$(FC) $(FCFLAGS) source/rangen.f90 -o rangen_debug -O0 $(MKL_VSL_FLAGS) -g -D_DEBUG

#1_nqueens_serial: nqueens_serial.f90
#	$(FC) $(FCFLAGS) nqueens_serial.f90 -o 1_nqueens_serial -O2 -g #$(ANNOTATE_FLAGS)
#1_nqueens_serial_debug: nqueens_serial.f90
#	$(FC) $(FCFLAGS) nqueens_serial.f90 -o 1_nqueens_serial_debug -O0 -g -D_DEBUG #$(ANNOTATE_FLAGS)

EXECUTABLES =  coherence

clean::
	rm -f $(EXECUTABLES) *.o

