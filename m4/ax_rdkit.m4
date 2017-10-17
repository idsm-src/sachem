AC_DEFUN([AX_RDKIT],
[
    AC_ARG_WITH([rdkit],
        AC_HELP_STRING([--with-lucy=@<:@ARG@:>@],
            [set the path to your local RDKit installation]
        ),
        [
            RDKIT_CPPFLAGS="-I${withval}/Code"
            RDKIT_LDFLAGS="-L${withval}/lib"
        ],
        [
            RDKIT_CPPFLAGS=
            RDKIT_LDFLAGS=
        ])

    AC_SUBST([RDKIT_CPPFLAGS])
    AC_SUBST([RDKIT_LDFLAGS])
])
