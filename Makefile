# -*- Makefile -*-
#
# Description:	Quake's Hercules toolchain top directory Makefile

# User preferences (-include means that the file is optional. It won't raise an error if it doesn't exist.)
-include user.mk

# Set working directory
WORKDIR = $(CURDIR)

# Include other makefiles
include systemdef.mk
include common.mk

# Define function to run `make` in a directory for a given target
# 
# param ${1}: dir to run `make` on
# param ${2}: target
# NOTE: use `MAKE` instead of `make` for recursive calls to avoid infinite recursion
# NOTE 2: `MAKE` is a built-in variable in GNU Make
# NOTE 3: To use a variable in makefile, use ${VAR} or $(VAR)
# NOTE 4: The empty line after the function definition is required to indicate the end of the shell command 
# and prevents the shell from interpreting the next line as a command.
define run_make_in_dir
	@echo Building TARGET=${2} in DIR=${1}
	${MAKE} -C ${1} SYSTEM=${SYSTEM} WORKDIR=${WORKDIR} ${2}

endef

# Define function to run `make` for a target in the visualization directory if ENABLE_VIZ is set
make_target_in_viz_dir = $(if $(ENABLE_VIZ),$(call run_make_in_dir,$(VIS_DIR),${1}))
# NOTE: The following line is a commented line in the original Makefile that is not properly attached to any target
#	$(MAKE) -C $(VIS_DIR)     SYSTEM=$(SYSTEM) WORKDIR=$(WORKDIR)

.PHONY: all clean cleanall etree octor cvm forward

# Build all targets (default target since it is the first target in the file)
all:	etree octor cvm forward
	$(call make_target_in_viz_dir,all)

# Build the etree target
# NOTE: The `-C` option is used to change the working directory to the directory specified by the first argument
etree:
	$(MAKE) -C $(ETREE_DIR)   SYSTEM=$(SYSTEM) WORKDIR=$(WORKDIR)

# Build the octor target
octor:
	$(MAKE) -C $(OCTOR_DIR)   SYSTEM=$(SYSTEM) WORKDIR=$(WORKDIR)

# Build the cvm target
cvm:
	$(MAKE) -C $(CVM_DIR)     SYSTEM=$(SYSTEM) WORKDIR=$(WORKDIR)

# Build the forward target
forward:
	$(MAKE) -C $(FORWARD_DIR) SYSTEM=$(SYSTEM) WORKDIR=$(WORKDIR)

# Clean targets
clean:
	$(call make_target_in_viz_dir,clean)
	$(MAKE) -C $(ETREE_DIR)   WORKDIR=$(WORKDIR) clean
	$(MAKE) -C $(OCTOR_DIR)   WORKDIR=$(WORKDIR) clean
	$(MAKE) -C $(CVM_DIR)     WORKDIR=$(WORKDIR) clean
	$(MAKE) -C $(FORWARD_DIR) WORKDIR=$(WORKDIR) clean

# Define list of directories
MY_DIRS := $(ETREE_DIR) $(OCTOR_DIR) $(CVM_DIR) $(FORWARD_DIR)

# Define function to run `make` for target `run_make_in_dir` in each directory in the directory list
#
# param ${1}: dir list
# param ${2}: target
run_make_for_dirs = $(foreach MY_DIR,${1},$(call run_make_in_dir,${MY_DIR},${2}))

# Clean all targets
cleanall: clean
	@echo MY_DIRS="${MY_DIRS}"
	$(call run_make_for_dirs,${MY_DIRS},cleanall)

# $Id: Makefile,v 1.11 2010/07/13 20:34:21 rtaborda Exp $
# Update on 2023/04/14 by Clifford Yen