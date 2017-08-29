#######################################
# The list of libraries we are building seperated by spaces.
# The 'lib_' indicates that these build products will be installed
# in the $(libdir) directory. For example /usr/lib
lib_LTLIBRARIES = libpgchem.la

#######################################
# Build information for each library

# Sources for pgsparql
libpgchem_la_SOURCES = pgchem.c

# Linker options libTestProgram
libpgchem_la_LDFLAGS = 

# Compiler options. Here we are adding the include directory
# to be searched for headers included in the source code.
libpgchem_la_CPPFLAGS = $(POSTGRESQL_CPPFLAGS)
