AC_DEFUN([AX_JDK],
[
    AC_ARG_WITH([jdk],
        AS_HELP_STRING([--with-jdk=@<:@ARG@:>@],
            [set the path to your JDK directory @<:@default=${JAVA_HOME}@:>@]
        ),
        [
        if test "$withval" == "yes"; then
            JDK_HOME=${JAVA_HOME}
        elif test "$withval" == "no"; then
            JDK_HOME=
        else
            JDK_HOME=${withval}
        fi
        ],
        [
            JDK_HOME=${JAVA_HOME}
        ]
    )


    AC_MSG_CHECKING([for JDK])

    if test x"${JDK_HOME}" != "x"; then
        AC_MSG_RESULT([${JDK_HOME}])
    else
        AC_MSG_RESULT([no])
        AC_MSG_ERROR([JDK was not found])
    fi


    AC_CHECK_FILE([${JDK_HOME}/jre/lib/amd64/server/libjvm.so],
        [JDK_LIB=${JDK_HOME}/jre/lib/amd64/server],
        [AC_CHECK_FILE([${JDK_HOME}/lib/server/libjvm.so],
            [JDK_LIB=${JDK_HOME}/lib/server],
            [AC_MSG_ERROR([JDK was not found])])
        ])

    AC_CHECK_FILE([${JDK_HOME}/include/jni.h],
        [JDK_INCLUDE=${JDK_HOME}/include],
        [AC_MSG_ERROR([JDK was not found])])


    AC_MSG_CHECKING([whether jni.h supports 21])
    
    if grep "^#define JNI_VERSION_21  0x00150000$" ${JDK_INCLUDE}/jni.h > /dev/null; then
        AC_MSG_RESULT([yes])
    else
        AC_MSG_RESULT([no])
        AC_MSG_ERROR([JDK was not found])
    fi
    

    JDK_CPPFLAGS="-I${JDK_INCLUDE} -I${JDK_INCLUDE}/linux"
    JDK_LDFLAGS="-Wl,-rpath,${JDK_LIB} -L${JDK_LIB} -ljvm"

    AC_SUBST([JDK_HOME])
    AC_SUBST([JDK_LDFLAGS])
    AC_SUBST([JDK_CPPFLAGS])
])
