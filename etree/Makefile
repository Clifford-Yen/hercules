# NOTE: When called from the `Makefile` in the parent folder, `WORKDIR` should be given as an argument to the `make` command.
# Don't know why we need to define `WORKDIR` again here.
ifndef WORKDIR
	WORKDIR = $(PWD)/..
endif

-include $(WORKDIR)/systemdef.mk
-include $(WORKDIR)/common.mk

#
# Object modules of the library
# NOTE: `.o` files are compiled from source files (e.g. `.c` files) and contain machine code that can be linked together to form an executable or a library.
OBJECTS = dlink.o code.o  buffer.o schema.o xplatform.o btree.o etree.o wrapper.o
# NOTE: `.a` files are libraries that contain one or more object files (`.o` files).
TARGET = libetree.a

# The default target (the first target in the file)
all: $(TARGET)

# This target is used to build a static library file named `libetree.a` from the object files specified in the `OBJECTS` variable
# NOTE: The `ar` command is used to create, modify, and extract from archives. 
# The `r` option is used to replace or add files to the archive. 
# The `c` option is used to create the archive if it does not exist. 
# The `s` option is used to add an index to the archive.
# NOTE 2: The `ranlib` command is used to generate an index to the contents of an archive.
# In older versions of Unix, the `ranlib` command was necessary to create the index, but in modern Unix-like systems, the `ar` command includes the ability to create the index with the `s` option.
libetree.a: $(OBJECTS)
	ar rcs libetree.a $(OBJECTS) 
	ranlib libetree.a

.PHONY: clean
# This target is used to remove the object files and backup files (files ending with `~`) from the directory
# NOTE: The `-f` option is used to force the removal of files without prompting the user for confirmation. It also ignores the non-existence of files and arguments.
clean:
	rm -f $(OBJECTS) *~ 

.PHONY: cleanll
cleanall:
	rm -f $(TARGET) $(OBJECTS) *~	

