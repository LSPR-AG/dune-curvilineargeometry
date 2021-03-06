dnl -*- autoconf -*-
# Macros needed to find dune-curvilineargeometry and dependent libraries.  They are called by
# the macros in ${top_src_dir}/dependencies.m4, which is generated by
# "dunecontrol autogen"

# Additional checks needed to build dune-curvilineargeometry
# This macro should be invoked by every module which depends on dune-curvilineargeometry, as
# well as by dune-curvilineargeometry itself
AC_DEFUN([DUNE_CURVILINEARGEOMETRY_CHECKS])

# Additional checks needed to find dune-curvilineargeometry
# This macro should be invoked by every module which depends on dune-curvilineargeometry, but
# not by dune-curvilineargeometry itself
AC_DEFUN([DUNE_CURVILINEARGEOMETRY_CHECK_MODULE],
[
  DUNE_CHECK_MODULES([dune-curvilineargeometry],[curvilineargeometry/curvilineargeometry.hh])
])
