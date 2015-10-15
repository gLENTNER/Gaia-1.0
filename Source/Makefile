# Copyright (c) Geoffrey Lentner 2015. All Rights Reserved.
# See LICENSE file (GPLv3)
#
# Makefile for GAIA-2.1.1 project.

INSTALL   = /usr/local/bin/

CC        = g++-5
CCFLAGS   = -fopenmp -std=c++11
EXE       = gaia
INC       = Include
LIB       = Library
OBJ       = Objects
MAIN      = Library/Main

Tools     = KernelFit Interpolate Random
Framework = Simulation Parser Monitor FileManager PopulationManager
Profiles  = ProfileBase ProfileManager

Sources   = $(addprefix $(OBJ)/, $(Framework) $(Tools) $(Profiles))
Objects   = $(addsuffix .o, $(Sources))

# Primary `link`ing target
$(EXE): lib$(EXE).a $(MAIN).o
	$(CC) -o $(EXE) $(MAIN).o -I$(INC) -L. -l$(EXE) $(CCFLAGS)

# archive construction template
lib$(EXE).a: $(Objects)
	ar -rv $@ $(Objects)
	ranlib $@

# `main`
$(MAIN).o: $(MAIN).cpp
	$(CC) -c $< -o $(MAIN).o -I$(INC) $(CCFLAGS)

# all objects compiled with this template
$(OBJ)/%.o: $(LIB)/%.cpp $(INC)/%.hpp
	$(CC) -c $< -o $@ -I$(INC) $(CCFLAGS)

# copy to install directory
install:
	cp $(EXE) $(INSTALL)/

# clear all objects and archives
clean:
	rm -f $(Objects) lib$(EXE).a

.PHONY: clean
