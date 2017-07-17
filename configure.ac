AC_PREREQ([2.63])
AC_INIT([Evolution], [1.0], [w.cunningham@northeastern.edu])
AM_INIT_AUTOMAKE

if [[ -e /home/$USER/.evolution ]] ; then
  source /home/$USER/.evolution
fi

if [[[ "$EVOLUTION_HOME_DIR" == "" ]]] ; then
  echo "export EVOLUTION_HOME_DIR=$PWD" >> $HOME/.bashrc
  source $HOME/.bashrc
fi

AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_SRCDIR([src/causet.cpp])
AC_CONFIG_HEADERS([inc/config.h])

AC_PROG_CC
AC_PROG_CXX
AC_PROG_RANLIB

AC_LANG_PUSH([C++])

AX_COMPILER_VERSION
AX_COMPARE_VERSION([$ax_cv_cxx_compiler_version], [ge], [4.8.1],, AC_MSG_ERROR([GNU compiler must be at least version 4.8.1]))
AX_CXX_COMPILE_STDCXX([11])

AX_FASTMATH([$FASTMATH_HOME], [mandatory])
AX_PATH_GSL([1.13],, AC_MSG_ERROR([Could not find the required version of GSL.]))
BOOST_REQUIRE([1.55.0])
AX_PLATFORM([$PLATFORM])
AC_OPENMP

AC_CHECK_LIB([m], [cos],, AC_MSG_ERROR([Could not find library "libm"]))
AC_CHECK_LIB([gslcblas], [cblas_dgemm],, AC_MSG_ERROR([Could not find library "libgslcblas"]))
AC_CHECK_LIB([stdc++], [main],, AC_MSG_ERROR([Could not find library "libstdc++"]))

AC_CHECK_FILES([src/CuResources.cpp src/causet.cpp src/subroutines.cpp],, AC_MSG_ERROR([Could not find Evolution source files.]))
AC_CHECK_FILES([inc/constants.h inc/CuResources.h inc/operations.h inc/causet.h inc/subroutines.h],, AC_MSG_ERROR([Could not find Evolution header files.]))

#Add AC_CHECK_HEADERS later

AC_FUNC_MALLOC
AC_TYPE_SIZE_T
AC_TYPE_INT64_T
AC_TYPE_UINT64_T
AC_C_INLINE
AC_C_CONST
AC_C_VOLATILE

CXXFLAGS="$CXXFLAGS -O3 -g -Wall -x c++ -DBOOST_NO_FENV_H -mpopcnt"
LDFLAGS="-lstdc++ -lm -lgsl -lgslcblas"
AX_GCC_X86_CPU_SUPPORTS([avx2], [AVXFLAGS="-mavx2 -DAVX2_ENABLED"], AC_MSG_RESULT([AVX2 features will be disabled.]))
AX_GCC_X86_CPU_SUPPORTS([popcnt],, AC_MSG_ERROR([CPU does not support the popcnt instruction.]))
AX_GCC_ARCHFLAG([no],, AC_MSG_ERROR([CPU architecture could not be identified.]))
AC_SUBST([ARCHFLAGS], [$ax_cv_gcc_archflag])

BINARY_NAME=Evolution$PLATFORM
COMPILE_STR="$CXX -c $CXXFLAGS $FASTMATH_FLAGS $AVXFLAGS $ARCHFLAGS $OPENMP_CXXFLAGS"
LINK_STR="$CXX src/*.o $OPENMP_CXXFLAGS $LDFLAGS $FASTMATH_LDFLAGS -o bin/$BINARY_NAME"

AC_SUBST([COMPILE_STR], [$COMPILE_STR])
AC_SUBST([LINK_STR], [$LINK_STR])
AC_SUBST([BINARY_NAME], [$BINARY_NAME])

AC_CONFIG_FILES([Makefile src/Makefile inc/Makefile printcolor/Makefile])
AC_OUTPUT
