#
#
# AC_BLOSC
#
#
#
dnl @synopsis AC_BLOSC
dnl
dnl This macro test if BLOSC is to be used.
dnl Use in C code:
dnl     #ifdef BLOSC
dnl     #include "blosc.h"
dnl     #endif
dnl
dnl @version 2.0
dnl @author Zhenhuan (Steve) Gong
dnl @author Norbert Podhorszki
dnl
AC_DEFUN([AC_BLOSC],[

AC_MSG_NOTICE([=== checking for BLOSC ===])

AM_CONDITIONAL(HAVE_BLOSC,true)

AC_ARG_WITH(blosc,
        [  --with-blosc=DIR      Location of BLOSC library],
        [:], [with_blosc=no])

if test "x$with_blosc" == "xno"; then

   AM_CONDITIONAL(HAVE_BLOSC,false)

else

    save_CPPFLAGS="$CPPFLAGS"
    save_LIBS="$LIBS"
    save_LDFLAGS="$LDFLAGS"

    if test "x$with_blosc" == "xyes"; then
        dnl No path given
        BLOSC_CPPFLAGS=""
        BLOSC_LDFLAGS=""
        BLOSC_LIBS="-lblosc"
    else
        dnl Path given, first try path/lib64
        BLOSC_CPPFLAGS="-I$withval/include"
        BLOSC_LDFLAGS="-L$withval/lib64"
        BLOSC_LIBS="-lblosc"
    fi

    LIBS="$LIBS $BLOSC_LIBS"
    LDFLAGS="$LDFLAGS $BLOSC_LDFLAGS"
    CPPFLAGS="$CPPFLAGS $BLOSC_CPPFLAGS"

    dnl Find header file first
    AC_CHECK_HEADERS(blosc.h,
              ,
              [AM_CONDITIONAL(HAVE_BLOSC,false)])

    if test -z "${HAVE_BLOSC_TRUE}"; then
        dnl Try to link an example now
        AC_MSG_CHECKING([if blosc code can be linked with $BLOSC_LDFLAGS])
        AC_TRY_LINK(
            [#include <stdlib.h>
             #include "blosc.h"],
            [float in[50], out[50];
             int in_len = 50 * sizeof(float), out_len = in_len;
             int level = 5;
             blosc_init();		
             int zerr = blosc_compress (level, 1, sizeof(float), in_len, in, out, out_len);
             blosc_destroy();
             return (zerr <= 0);],
            [AC_MSG_RESULT(yes)],
            [AM_CONDITIONAL(HAVE_BLOSC,false)
             AC_MSG_RESULT(no)
            ])
            
        dnl If linking above failed, one reason might be that we looked in lib64/
        dnl instead of lib/
    fi

    LIBS="$save_LIBS"
    LDFLAGS="$save_LDFLAGS"
    CPPFLAGS="$save_CPPFLAGS"

    AC_SUBST(BLOSC_LIBS)
    AC_SUBST(BLOSC_LDFLAGS)
    AC_SUBST(BLOSC_CPPFLAGS)

    # Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
    if test -z "${HAVE_BLOSC_TRUE}"; then
            ifelse([$1],,[AC_DEFINE(HAVE_BLOSC,1,[Define if you have BLOSC.])],[$1])
            :
    else
            $2
            :
    fi
fi
])dnl AC_BLOSC
