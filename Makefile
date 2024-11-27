#Makefile for Basic 1D Electrostatic PIC code

# Makefile gfortran compiler with MacOS X

#OpenMP
#MPFC = gfortran -fopenmp
#MPCC = gcc -fopenmp
#NoOpenMP
#MPFC = gfortran
#MPCC = gcc

#FC90 = gfortran
#CC = gcc

#PFC = f2py --fcompiler=gnu95
#PCC = f2py --compiler=unix

#OPTS90 = -O3
#OPTS90 = -O3 -fdefault-real-8 -fdefault-double-8
#OPTS90 = -O3 -fcheck=bounds -fdefault-real-8 -Wall -std=f95

#CCOPTS = -O3 -std=c99
#CCOPTS = -O3 -Wall -std=c99

#Python
#PYOPTS90 = $(OPTS90) -fPIC
#PYCCOPTS = $(CCOPTS) -fPIC

#LEGACY = -std=legacy

# Makefile Intel compiler with Linux

# OpenMP
#MPFC = ifort -qopenmp
#MPCC = icc -qopenmp
#NoOpenMP
#MPFC = ifort
#MPCC = icc

#FC90 = ifort
#CC = gcc

#OPTS90 = -O3
#OPTS90 = -O3 -r8
#OPTS90 = -O3 -CB -r8 -warn all -std90

#CCOPTS = -O3 -std=c99
#CCOPTS = -O3 -no-vec -Wall -std=c99

#FF03 = -Tf

#Python
#PFC = f2py --fcompiler=ifort
#PCC = f2py --compiler=unix

#LOPTS = -liomp5

#Python
#PYOPTS90 = $(OPTS90) -fPIC
#PYCCOPTS = $(CCOPTS) -fPIC

# Makefile gfortran compiler with Linux

# OpenMP
MPFC = gfortran -fopenmp
MPCC = gcc -fopenmp
#NoOpenMP
#MPFC = gfortran
#MPCC = gcc
FC90 = gfortran
CC = gcc

FC90 = gfortran
CC = gcc

#Python
PFC = f2py --fcompiler=gnu95
PCC = f2py --compiler=unix

#OPTS90 = -O3
OPTS90 = -O3 -fdefault-real-8 -fdefault-double-8
#OPTS90 = -O3 -fcheck=all -fdefault-real-8 -fdefault-double-8
#OPTS90 = -O3 -fbounds-check -fdefault-real-8 -Wall -std=f95

CCOPTS = -O3 -std=c99
#CCOPTS = -O3 -Wall -std=c99

LOPTS = -lgomp

#Python
PYOPTS90 = $(OPTS90) -fPIC
PYCCOPTS = $(CCOPTS) -fPIC

LEGACY = -std=legacy

# Makefile PGI compiler with Linux

# OpenMP
#MPFC = pgf90
#MPCC = gcc -fopenmp

#FC90 = pgf90
#CC = gcc

#OPTS90 = -O3
#OPTS90 = -O3 -r8
#OPTS90 = -O3 -Mbounds -r8 -Mstandard

#CCOPTS = -O3 -std=c99
#CCOPTS = -O3 -Wall -std=c99

# Makefile Nag compiler with Linux

#FC90 = nagfor
#CC = gcc

#OPTS90 = -O3
#OPTS90 = -O3 -default_kinds:r=64
#OPTS90 = -O3 -default_kinds:r=64 -C=array -nan -w=all -f95

#CCOPTS = -O3 -std=c99
#CCOPTS = -O3 -Wall -std=c99


# Makefile Cray compiler with Linux

#FC90 = ftn
#CC = cc

#OPTS90 = -O 3
#OPTS90 = -O 3 -s real64
#OPTS90 = -O 3 -R b -s real64 -en

#CCOPTS = -O 3 -h c99 -h conform

# Linkage rules

all : espic1.py

python : espic1.py

fortran : espic1

# Fortran

espic1 : espic1.o push1.o diag1.o dtimer.o
	$(MPFC) $(OPTS90) $(LOPTS) -o espic1 espic1.o push1.o push1_h.o \
	wpush1.o diag1.o diag1_h.o wdiag1.o omplib.o dtimer.o

# Compilation rules

dtimer.o : dtimer.c
	$(CC) $(CCOPTS) -c dtimer.c

#OPENMP
omplib.o : omplib.f90
	$(MPFC) $(OPTS90) -o omplib.o -c omplib.f90

push1.o : push1.f
	$(MPFC) $(OPTS90) -o push1.o -c push1.f

push1_h.o : push1_h.f90
	$(FC90) $(OPTS90) -o push1_h.o -c push1_h.f90

wpush1.o : wpush1.f90 push1_h.o
	$(FC90) $(OPTS90) -o wpush1.o -c wpush1.f90

dwpush1.o : dwpush1.f90 push1_h.o
	$(FC90) $(OPTS90) -o dwpush1.o -c dwpush1.f90

diag1.o : diag1.f
	$(FC90) $(OPTS90) -o diag1.o -c diag1.f

diag1_h.o : diag1_h.f90
	$(FC90) $(OPTS90) -o diag1_h.o -c diag1_h.f90

wdiag1.o : wdiag1.f90 diag1_h.o
	$(FC90) $(OPTS90) -o wdiag1.o -c wdiag1.f90

dwdiag1.o : dwdiag1.f90 diag1_h.o
	$(FC90) $(OPTS90) -o dwdiag1.o -c dwdiag1.f90

espic1.o : espic1.f90 wpush1.o wdiag1.o omplib.o
	$(MPFC) $(OPTS90) -o espic1.o -c espic1.f90

# Python

espic1.py : pic1lib.so fomplib.so dtimer.so

pypush1.o : push1.f
	$(MPFC) $(PYOPTS90) -o pypush1.o -c push1.f

pydiag1.o : diag1.f
	$(FC90) $(PYOPTS90) -o pydiag1.o -c diag1.f

pic1lib.so : dwpush1.o pypush1.o dwdiag1.o pydiag1.o
	$(PFC) --f90flags="$(OPTS90)" $(LOPTS) -m pic1lib \
        -c dwpush1.f90 push1_h.o pypush1.o dwdiag1.f90 diag1_h.o \
        pydiag1.o

fomplib.so : omplib.f90
	$(PFC) --f90flags="$(OPTS90)" $(LOPTS) -m fomplib -c omplib.f90

dtimer.so : dtimer_hpy.f90 dtimer.c
	$(PFC) -m dtimer -c dtimer_hpy.f90 dtimer.c

clean :
	rm -f *.o *.mod *.pyf

clobber: clean
	rm -f espic1 *.so

list:
	echo fortran python

