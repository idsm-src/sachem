dnl Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT(sachem, 2.5.5)
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_AUX_DIR([build])
: ${CFLAGS=""}
: ${CXXFLAGS=""}
AC_PROG_CC

AC_CANONICAL_SYSTEM

AC_CONFIG_MACRO_DIR([m4])

dnl Initialize automake
AM_INIT_AUTOMAKE([subdir-objects tar-pax])

dnl this allows us specify individual liking flags for each target
AM_PROG_CC_C_O 

dnl Initialize Libtool
LT_INIT

dnl check for ant
AC_PATH_PROG([ANT], [ant], [])
if test "no$ANT" = no ; then
    AC_MSG_ERROR([ant not found])
fi

dnl check for PostgreSQL
AX_LIB_POSTGRESQL(9.6.0)

dnl check for JDK8
AX_JDK8


AC_CONFIG_FILES(Makefile src/Makefile java/Makefile extension/Makefile)
AC_OUTPUT
