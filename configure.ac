# *- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([ACOMP], [1.0], [swati@cs.umd.edu])

AC_CONFIG_SRCDIR([src/acomp.c])
AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_HEADERS([config.h])

AM_INIT_AUTOMAKE

AC_CONFIG_MACRO_DIRS([m4])

LT_INIT

# Checks for programs.
AC_PROG_CC

# Checks for libraries.

# Checks for header files.
AC_CHECK_HEADERS([float.h limits.h stdlib.h string.h sys/time.h])

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.
AC_FUNC_MALLOC
AC_CHECK_FUNCS([gettimeofday pow sqrt])

# Create a Makefile
AC_CONFIG_FILES([Makefile src/Makefile])

AC_OUTPUT
