# Copyright (c) Geoffrey Lentner 2015. All Rights Reserved.
# See LICENSE file (GPLv3)
#
# Makefile for GAIA-2.1.1 project.

# variables
compile = g++
debug   = 
ccflags = -fopenmp -std=c++11
exe     = Bin/Gaia
bin     = Bin
inc     = Include
lib     = Library
obj     = Object

# change this line to install to a different location!
install = /usr/bin/

# Primary `link`ing target
link: Main libTools libFramework libProfiles
	${compile} ${debug} -o ${exe} Main.o -I${inc} -L. -lTools \
		-lProfiles -lFramework ${ccflags}
	mv *.a *.o ${obj}/

Tools      = KernelFit.o Interpolate.o Random.o
Tools:       KernelFit.o Interpolate.o Random.o

Framework  = Simulation.o Parser.o Monitor.o FileManager.o PopulationManager.o
Framework:   Simulation.o Parser.o Monitor.o FileManager.o PopulationManager.o

Profiles   = ProfileBase.o ProfileManager.o
Profiles:    ProfileBase.o ProfileManager.o

# library construction template
lib%: %
	ar -rv $@.a ${$*}
	ranlib $@.a

# main has no header file
Main: ${lib}/Main.cc
	${compile} ${debug} -c ${lib}/Main.cc -I${inc} ${ccflags}

# all objects compiled with this template
%.o: ${lib}/%.cc ${inc}/%.hh
	${compile} ${debug} -c ${lib}/$*.cc -I${inc} ${ccflags}

.PHONY: clean

# clear all objects and archives
clean:
	rm -rf ${obj}/{*.o,*.a}

install:
	cp ${exe} ${install}
