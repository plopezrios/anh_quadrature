#!/bin/bash

# Select Fortran and C compilers based on hostname.
ccopts=""
case "$(hostname 2>/dev/null)" in
allogin*)
  . /usr/share/modules/init/bash
  module load ifort mpi.intel
  cc=mpiicc
  f90=mpiifort ;;
vortex*)
  . /etc/profile.d/modules.sh
  module load intel impi
  cc=mpiicc
  f90=mpiifort ;;
pc*.tcm.*)
  cc=gcc
  f90=gfortran ;;
*)
  cc=mpicc
  f90=mpif90 ;;
esac

# Choose compiler options.
case "$f90" in
mpif90*|mpgfortran*|gfortran*)
  if [ "$1" = -d ] ; then
    opts="-Wall -Wextra -fimplicit-none -O0 -fbounds-check -g -pg -pedantic\
       -fbacktrace -fcray-pointer"
  else
    opts="-O2 -fprotect-parens -march=native -fcray-pointer"
  fi ;;
mpifort*|mpiifort*|ifort*)
  if [ "$1" = -d ] ; then
    opts="-check all -check noarg_temp_created -fp-stack-check -g -O0\
       -implicitnone -std95 -traceback -warn all,nounused -debug all -ftrapuv\
       -assume buffered_io"
  else
    opts="-O3 -no-prec-div -no-prec-sqrt -funroll-loops -no-fp-port -ip\
       -complex-limited-range -assume protect_parens -assume buffered_io"
  fi ;;
mpnagfor*|nagfor*)
  if [ "$1" = -d ] ; then
    opts="-g -C=all -colour -f95 -strict95 -u -v -gline -nan\
       -no_underflow_warning -w=uda -O0"
  else
    opts="-O4 -Ounsafe"
  fi ;;
*)
  echo "Compiler '$f90' not configured."
  exit 1 ;;
esac

# Compile.
$f90 $opts -o anh_quadrature anh_quadrature.f90 -llapack 2>&1
rm -f *.mod *.o
