dnl Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT(sachem, 0.1)
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_AUX_DIR([build])
AC_PROG_CC
AC_PROG_CXX

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
AX_LIB_POSTGRESQL(9.5.0)

dnl check for JDK8
AX_JDK8


dnl process third party libraries
AC_MSG_NOTICE([third party libraries])
if test ! -e build; then
    mkdir build
fi


dnl process clownfish
AC_MSG_NOTICE([apache-clownfish])
if test ! -e build/clownfish ; then
    cp -r $srcdir/lib/clownfish build
    chmod -R u+w build/clownfish
fi

cd build/clownfish/runtime/c
./configure || exit $?
CLOWNFISH_CPPFLAGS=
CLOWNFISH_LDFLAGS="-L$(pwd)"
cd -

AC_SUBST([CLOWNFISH_CPPFLAGS])
AC_SUBST([CLOWNFISH_LDFLAGS])


dnl process lucy
AC_MSG_NOTICE([apache-lucy])
if test ! -e build/lucy ; then
    cp -r $srcdir/lib/lucy build
    chmod -R u+w build/lucy
fi

cd build/lucy/c
./configure || exit $?
sed -i 's|cfc|../../clownfish/compiler/c/cfc --include=../../clownfish/runtime/c/autogen/share/clownfish/include|' Makefile
sed -i 's|-l clownfish|-L ../../clownfish/runtime/c -l clownfish|' Makefile
LUCY_CPPFLAGS="-I$(pwd)/autogen/include"
LUCY_LDFLAGS="-L$(pwd)"
cd -

AC_SUBST([LUCY_CPPFLAGS])
AC_SUBST([LUCY_LDFLAGS])


dnl process rdkit
AC_PATH_PROG([CMAKE], [cmake], [])
if test "no$CMAKE" = no ; then
    AC_MSG_ERROR([cmake not found])
fi

AC_MSG_NOTICE([IOCB RDKit])
if test ! -e build/rdkit ; then
    cp -r $srcdir/lib/rdkit build
    chmod -R u+w build/rdkit
fi

cd build/rdkit
$CMAKE -DCMAKE_CXX_FLAGS=-std=gnu++11 || exit $?
RDKIT_CPPFLAGS="-I$(pwd)/Code"
RDKIT_LDFLAGS="-L$(pwd)/lib"
cd -

AC_SUBST([RDKIT_CPPFLAGS])
AC_SUBST([RDKIT_LDFLAGS])


AC_CONFIG_FILES(Makefile src/Makefile src/fpsearch/Makefile java/Makefile lib/Makefile extension/Makefile)
AC_OUTPUT
