FC      = gfortran
FCFLAGS = -O2 -Wall -Wno-unused-dummy-argument -Wno-unused-function -ffree-line-length-none -std=f2008 -fall-intrinsics
# FCFLAGS = -O2 -ffpe-trap=invalid,zero,overflow -fbacktrace -g
# FCFLAGS = -fopenmp -static-libgfortran -static-libgcc

DEST    = .

TARGET  = case0.dll
OBJS    = problem.o util.o

LIBPATH = ../orbitlib
INCPATH = ../orbitlib
LIBNAME = orbit
LIBOPT  = $(addprefix -L,$(LIBPATH)) $(addprefix -l,$(LIBNAME))
INCOPT  = $(addprefix -I,$(INCPATH))
# LIBBLAS = -llapack -lrefblas


################################################################################

.SUFFIXES: .f08
.PHONY: clean

all: $(TARGET)
test: test.exe

$(TARGET): $(OBJS)
	$(FC) $(FCFLAGS) -shared -o $@ $^ $(LIBOPT) $(LIBBLAS)

%.o: %.f08
	$(FC) $(FCFLAGS) -c -o $@ $< $(INCOPT)

clean:
	rm -f *.o *.mod *.smod

test.exe: test.o
	$(FC) -o $@ $<

################################################################################

problem.o: util.o
