FC      = gfortran
FCFLAGS = -O2 -Wall -Wno-unused-dummy-argument -Wno-unused-function
# FCFLAGS = -O2 -ffpe-trap=invalid,zero,overflow -fbacktrace -g
# FCFLAGS = -fopenmp -static-libgfortran -static-libgcc

DEST    = .

TARLIBD = liborbit.dll
TARLIBS = liborbit.a
OBJS    = orbit_debri.o orbit_base.o orbit_func.o

LIBPATH = ../fec/lib
INCPATH = ../fec/include
LIBNAME = fec
LIBOPT  = $(addprefix -L,$(LIBPATH)) $(addprefix -l,$(LIBNAME))
INCOPT  = $(addprefix -I,$(INCPATH))


################################################################################

.SUFFIXES: .f08
.PHONY: clean

all: static
dy: dynamic
st: static

dynamic: $(TARLIBD)
static: $(TARLIBS)

$(TARLIBD): $(OBJS)
	$(FC) $(FCFLAGS) -shared -o $@ $^ $(LIBOPT) -L$(LIBPATH) $(LIBBLAS)
	# cp $@ $(DEST)/$@

$(TARLIBS): $(OBJS)
	ar cr $@ $^

%.o: %.f08
	$(FC) $(FCFLAGS) -c -o $@ $< $(INCOPT)

%.mod: %.f08 %.o
	@:

%.smod: %.f08 %.o
	@:

clean:
	rm -f *.o *.mod *.smod

################################################################################

orbit_base.o: orbit_func.o
orbit_debri.o: orbit_func.o orbit_base.o
orbit_optimizer.o: orbit_func.o orbit_base.o
