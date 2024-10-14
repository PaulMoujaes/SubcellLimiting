# Copyright (c) 2010, Lawrence Livermore National Security, LLC. Produced at the
# Lawrence Livermore National Laboratory. LLNL-CODE-443211. All Rights reserved.
# See file COPYRIGHT for details.
#
# This file is part of the MFEM library. For more information and source code
# availability see http://mfem.org.
#
# MFEM is free software; you can redistribute it and/or modify it under the
# terms of the GNU Lesser General Public License (as published by the Free
# Software Foundation) version 2.1 dated February 1999.

# Use the MFEM build directory.

MFEM_DIR ?= ../../builds/MFEM
CONFIG_MK = $(MFEM_DIR)/share/mfem/config.mk

MFEM_LIB_FILE = mfem_is_not_built
-include $(CONFIG_MK)

# Enumerate code directories in logical order.

#PRO_DIR = systems
BUI_DIR = build
MET_DIR = methods
AUX_DIR = auxiliaries
OUT_DIR = output

DIRS = $(MET_DIR) $(AUX_DIR) 

#AUX_FILES = $(BUI_DIR)/dofs.o # $(BUI_DIR)/tools.o
MET_FILES = $(BUI_DIR)/feevol.met $(BUI_DIR)/loworder.met $(BUI_DIR)/clipandscale.met $(BUI_DIR)/convex_clipandscale.met $(BUI_DIR)/subcell_feevol.met $(BUI_DIR)/subcell_loworder.met $(BUI_DIR)/subcell_clipandscale.met
#PRO_FILES = $(BUI_DIR)/system.pro
#INT_FILES = $(BUI_DIR)/divergenceintegrator.int 
#FP_FILES = $(BUI_DIR)/fixedpointiteration.fp $(BUI_DIR)/standardfpiter.fp $(BUI_DIR)/jacobianfreenewton.fp $(BUI_DIR)/fpprecondjacobianfreenewton.fp
#OP_FILES = $(BUI_DIR)/pseudojacobian.op
# List scalar problems.

# PRO_FILES += $(BUI_DIR)/advection.pro $(BUI_DIR)/buckleyleverett.pro $(BUI_DIR)/burgers.pro $(BUI_DIR)/kpp.pro

# List systems.

#PRO_FILES += $(BUI_DIR)/euler.pro $(BUI_DIR)/advection.pro

MAIN_FILES = subcell

# ifeq ($(MFEM_USE_MPI), YES)
   #MAIN_FILES += pimpcg
   #PAUX_FILES += $(BUI_DIR)/pdofs.o
   #PAUX_FILES += $(BUI_DIR)/ptools.o
   #PMET_FILES += $(BUI_DIR)/pfixedpointiteration.met
#endif

# Setup valgrind test.

PROBLEM = 0
VALGRIND-CONFIG = -tf 0.001 -dt 0.001 -p $(PROBLEM) -e 0

## Remember some makefile rules.

# List keywords that are not associated with files (by default, all are).

.PHONY: all library subcell clean valgrind-test grid-convergence-test style # pimpcg inzufügen hinter impcg

# Delete the default suffixes.

.SUFFIXES:

# Define suffixes.

.SUFFIXES: .c .cpp

# Replace the default implicit rule for *.cpp files.

%: %.cpp $(CONFIG_MK)
	$(MFEM_CXX) $(MFEM_FLAGS) $< -ggdb -o $@ $(AUX_FILES) $(MET_FILES) $(MFEM_LIBS)

%.o: ../$(AUX_DIR)/%.cpp
	$(MFEM_CXX) $(MFEM_FLAGS) -c $^ -ggdb -o $@

%.met: ../$(MET_DIR)/%.cpp
	$(MFEM_CXX) $(MFEM_FLAGS) -c $^ -ggdb -o $@

#%.pro: ../$(PRO_DIR)/%.cpp
#	$(MFEM_CXX) $(MFEM_FLAGS) -c $^ -ggdb -o $@

#%.int: ../$(INT_DIR)/%.cpp
#	$(MFEM_CXX) $(MFEM_FLAGS) -c $^ -ggdb -o $@

#%.fp: ../$(FP_DIR)/%.cpp
#	$(MFEM_CXX) $(MFEM_FLAGS) -c $^ -ggdb -o $@

#%.op: ../$(OP_DIR)/%.cpp
#	$(MFEM_CXX) $(MFEM_FLAGS) -c $^ -ggdb -o $@

all: $(MAIN_FILES)

library: $(MET_FILES) $(AUX_FILES) 

subcell: library
	$(MFEM_CXX) $(MFEM_FLAGS) subcell.cpp -ggdb -o $@ $(AUX_FILES) $(MET_FILES) $(MFEM_LIBS) 

clean:
	@rm -f subcell errors.txt *gf* grid* $(BUI_DIR)/*

valgrind-serial:
	valgrind ./subcell $(VALGRIND-CONFIG) -c 1 -m meshes/inline-4segment.mesh -r 4 -o 1
	valgrind ./subcell $(VALGRIND-CONFIG) -c 0 -m meshes/periodic-square.mesh -r 2 -o 2 -s 1
	valgrind ./subcell $(VALGRIND-CONFIG) -c 0 -m meshes/beam-quad.mesh -r 1 -o 3
	valgrind ./subcell $(VALGRIND-CONFIG) -c 1 -m meshes/wall-bdr-4hex.mesh -r 0 -o 2

# TODO: Does not work yet.

#valgrind-parallel:
#	# mpirun -np 2 valgrind --leak-check=yes ./par-dgstab -tf 0.1 -r 0
#	# valgrind --gen-suppressions=yes pdgstab $(VALGRIND-CONFIG)
#	# valgrind --suppressions=$(MPI_HOME)/share/openmpi/openmpi-valgrind.supp pdgstab $(VALGRIND-CONFIG)
#	valgrind --suppressions=./my-mpi-suppressions.supp pdgstab $(VALGRIND-CONFIG)

# Generate an error message if the MFEM library is not built and exit.

$(MFEM_LIB_FILE):
	$(error The MFEM library is not built)

ASTYLE = astyle --options=$(MFEM_DIR)/config/mfem.astylerc
FORMAT_FILES = $(foreach dir,$(DIRS),"$(dir)/*.?pp")
FORMAT_FILES += subcell.cpp

style:
	@if ! $(ASTYLE) $(FORMAT_FILES) | grep Formatted; then\
	   echo "No source files were changed.";\
	fi
