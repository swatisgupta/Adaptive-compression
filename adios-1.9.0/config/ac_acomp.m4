#
#
# AC_ZLIB
#
#
#
dnl @synopsis AC_ACOMP
dnl
dnl This macro test if ACOMP is to be used.
dnl Use in C code:
dnl     #ifdef ACOMP
dnl     #include "acomp.h"
dnl     #endif
dnl
dnl @version 2.0
dnl @author Swati Singhal
dnl
AC_DEFUN([AC_ACOMP],[

AC_MSG_NOTICE([=== checking for ACOMP ===])

AM_CONDITIONAL(HAVE_ACOMP,true)

AC_ARG_WITH(acomp,
        [  --with-acomp=DIR      Location of ACOMP library],
        [:], [with_acomp=no])

if test "x$with_acomp" == "xno"; then

   AM_CONDITIONAL(HAVE_ACOMP,false)

else

    save_CPPFLAGS="$CPPFLAGS"
    save_LIBS="$LIBS"
    save_LDFLAGS="$LDFLAGS"

    if test "x$with_acomp" == "xyes"; then
        dnl No path given
        ACOMP_CPPFLAGS="-I/usr/include -I/homes/ssinghal/Builds/include" 
        ACOMP_LDFLAGS=""
        ACOMP_LIBS="-L/homes/ssinghal/Builds/lib -lbz2 -llzo2 -lacomp -lbz2 -llzo2 -L/usr/lib64 -lz -lm"
    else
        dnl Path given, first try path/lib64
        ACOMP_CPPFLAGS="-I$withval/include -I/homes/ssinghal/Builds/include"
        ACOMP_LDFLAGS=""
        ACOMP_LIBS="-L/homes/ssinghal/Builds/lib -lbz2 -llzo2 -lacomp -lbz2 -llzo2 -L$withval/lib64 -lz -lm"
    fi

    LIBS="$LIBS $ACOMP_LIBS"
    LDFLAGS="$LDFLAGS $ACOMP_LDFLAGS"
    CPPFLAGS="$CPPFLAGS $ACOMP_CPPFLAGS"

    dnl Find header file first
    AC_CHECK_HEADERS(acomp.h,
              ,
              [AM_CONDITIONAL(HAVE_ACOMP,false)])

    if test -z "${HAVE_ACOMP_TRUE}"; then
        dnl Try to link an example now
        AC_MSG_CHECKING([if acomp code can be linked with $ACOMP_LDFLAGS])
        AC_TRY_LINK(
            [#include <stdlib.h>
             #include "acomp.h"],
            [unsigned int res;
	     res = acomp_get_max_metadata_size();
	     return (res != 31);],	           
             [AC_MSG_RESULT(yes)],
             [AM_CONDITIONAL(HAVE_ACOMP,false)
             AC_MSG_RESULT(no)
            ])
            
        dnl If linking above failed, one reason might be that we looked in lib64/
        dnl instead of lib/
        if test -z "${HAVE_ACOMP_FALSE}"; then
            if test "x$with_lustre" != "xyes"; then
               COMP_LDFLAGS="-L$withval/lib"
               LDFLAGS="$LDFLAGS $ACOMP_LDFLAGS"
               AC_MSG_CHECKING([if acomp code can be linked with $ACOMP_LDFLAGS])
               AC_TRY_LINK(
               [#include <stdlib.h>
               #include "acomp.h"],
               [unsigned int res;
	       res = acomp_get_max_metadata_size();
	       return (res != 31);],	            
               [AC_MSG_RESULT(yes)],
               [AM_CONDITIONAL(HAVE_ACOMP,false)
               AC_MSG_RESULT(no)
              ])
            fi
        fi
    fi

    LIBS="$save_LIBS"
    LDFLAGS="$save_LDFLAGS"
    CPPFLAGS="$save_CPPFLAGS"

    AC_SUBST(ACOMP_LIBS)
    AC_SUBST(ACOMP_LDFLAGS)
    AC_SUBST(ACOMP_CPPFLAGS)

    # Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
    if test -z "${HAVE_ACOMP_TRUE}"; then
            ifelse([$1],,[AC_DEFINE(HAVE_ACOMP,1,[Define if you have ACOMP.])],[$1])
            :
    else
            $2
            :
    fi
fi
])dnl AC_ACOMP
