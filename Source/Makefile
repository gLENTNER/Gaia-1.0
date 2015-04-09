# Copyright (c) Geoffrey Lentner 2015. All Rights Reserved.
# See LICENSE file (GPLv3)
#
# Makefile for GAIA-2.1.1 project.

# variables
compile = g++
debug   = -g
ccflags = -fopenmp
exe     = bin/GAIA
bin     = bin
inc     = Include
lib     = Library
obj     = Object

# change this line to install to a different location!
install = /usr/bin/

# included archives (code structure)
#Archives  = libFramework libProfiles libTools
#Framework = Simulation.o PopulationManager.o FileManager.o Parser.o Monitor.o
#Profiles  = ProfileManager.o ProfileBase.o
#Tools     = CurveFit.o Interpolation.o Random.o mt19937.o Vector.o

# Primary `link`ing target
link: Main libFramework
	${compile} ${debug} -o ${exe} Main.o -I${inc} -L. -lFramework ${ccflags}
	mv *.a *.o ${obj}/

Framework  = Parser.o Simulation.o Monitor.o
Framework:   Parser.o Simulation.o Monitor.o
	
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

