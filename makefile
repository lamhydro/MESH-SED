# ======================================================================
#
#						Makefile for SA_MESH
#
# ======================================================================

# ======================================================================
# Include definition file
# ======================================================================
include makefile.def

# ======================================================================
# Declaring variables
# ======================================================================
# The Compiler
# Ensure to disable the MPI stub if using an MPI compiler.

# SERIAL compilers
#   FC=gfortran for GNU compiler (Cygwin or MinGW)
#   FC=ifort for Intel compiler (Intel Visual Compiler or plato.usask.ca)
# ALSO UPDATE FTN90PP PRECOMPILER FOR SVS.
FC=gfortran
#FC=ifort

# MPI compilers for parallel computing
#   FC=mpifort for OpenMPI wrapper with either GNU or Intel compiler (Cygwin or plato.usask.ca)
# COMMENT MPI STUB IF USING AN MPI COMPILER.
#FC=mpifort

# Flags for compiling, profiling, and debugging - comment as necessary
#   -O2: Default optimization.
#   -O3 -ffast-math: faster optimization (for GCC/gfortran only).
#   -g: For debugging.
LFLAG=-c -g -O2 -Wall
#LFLAG=-c -Wall
#LFLAG=-c -O2
#LFLAG=-c -O3 -ffast-math
#LFLAG=-c -g


# ======================================================================
# Build SA_MESH executable and print message
# ======================================================================
all: ${OBJECTS}
	$(FC) -o sa_mesh_sed  $(OBJECTS)
#	$(FC) -o mpi_sa_mesh  $(OBJECTS)

#static: ${OBJECTS}
# For MinGW only (the Cygwin library cannot be statically linked to the binary):
#	$(FC) -o sa_mesh_static -static  $(OBJECTS)


# ======================================================================
# General rules
# ======================================================================
%.o: %.f
	$(FC) $(LFLAG) $<
%.o: %.F90
	$(FC) $(LFLAG) $<
%.o: %.f90
	$(FC) $(LFLAG) $<
%.o: %.for
	$(FC) $(LFLAG) $<

# ======================================================================
# Target to create dependencies by the program makedepf90
# ======================================================================
depend .depend:
#makedepf90 -o sa_mesh_sed $(SOURCES) > .depend
	makedepf90 -o sa_mesh_sed $(VPATH)*.f90 > .depend

# ======================================================================
# Cleaning object files
# ======================================================================
clean:
# 'rm' for Cygwin, 'del' for MinGW - comment as necessary
	rm *.mod *.o
#	del *.mod *.o

# ======================================================================
# Cleaning everything including the previously built executable
# ======================================================================
veryclean:
# 'rm' for Cygwin, 'del' for MinGW - comment as necessary
	rm *.mod *.o sa_mesh_sed
#	rm *.mod *.o mpi_sa_mesh
#	del *.mod *.o sa_mesh.exe
