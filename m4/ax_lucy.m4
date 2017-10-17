AC_DEFUN([AX_LUCY],
[
    AC_ARG_WITH([lucy],
        AC_HELP_STRING([--with-lucy=@<:@ARG@:>@],
            [set the prefix to your Lucy and Clownfish installation]
        ),
        [
            LUCY_CPPFLAGS="-I${withval}/include"
            LUCY_LDFLAGS="-L${withval}/lib"
        ],
        [
            LUCY_CPPFLAGS=
            LUCY_LDFLAGS=
        ])

    AC_SUBST([LUCY_CPPFLAGS])
    AC_SUBST([LUCY_LDFLAGS])
])
