dnl Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT(pgchem, 0.1)
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_AUX_DIR([build])
AC_PROG_CC

AC_CANONICAL_SYSTEM

AC_CONFIG_MACRO_DIR([m4])

dnl Initialize automake
AM_INIT_AUTOMAKE

dnl this allows us specify individual liking flags for each target
AM_PROG_CC_C_O 

dnl Initialize Libtool
LT_INIT

dnl check for postgresql
AX_LIB_POSTGRESQL(9.2.0)

libdir=$POSTGRESQL_LIBDIR

AC_CONFIG_FILES(Makefile src/Makefile)
AC_OUTPUT