FC      = gfortran
FCFLAGS = -O2 -Wall -Wno-unused-dummy-argument -Wno-unused-function
# FCFLAGS = -O2 -ffpe-trap=invalid,zero,overflow -fbacktrace -g
# FCFLAGS = -fopenmp -static-libgfortran -static-libgcc

DEST    = .

TARGET  = a.exe
OBJS    = main.o orbit_optimizer.o
OBJS    = main1_tsp.o
OBJS    = main2_op.o
OBJS    = main4.o class_prob1.o class_prob2.o class_opt1.o class_opt2.o
# OBJS    = plot_debris.o

LIBPATH = ../lib/fec/lib ../lib/orbit
INCPATH = ../lib/fec/include ../lib/orbit
LIBNAME = fec orbit
LIBOPT  = $(addprefix -L,$(LIBPATH)) $(addprefix -l,$(LIBNAME))
INCOPT  = $(addprefix -I,$(INCPATH))
LIBBLAS = -llapack -lrefblas


################################################################################

.SUFFIXES: .f08
.PHONY: clean

all: $(TARGET)
test: test.exe

$(TARGET): $(OBJS)
	$(FC) $(FCFLAGS) -o $@ $^ $(LIBOPT) $(LIBBLAS)

%.o: %.f08
	$(FC) $(FCFLAGS) -c -o $@ $< $(INCOPT)

clean:
	rm -f *.o *.mod *.smod

test.exe: test.o
	$(FC) -o $@ $<

################################################################################

main.o: orbit_optimizer.o
main3.o: class_prob1.o class_opt1.o class_opt2.o class_opt3.o
main4.o: class_prob1.o class_prob2.o class_opt1.o class_opt2.o

class_prob2.o: class_opt1.o
class_opt1.o: class_prob1.o
class_opt2.o: class_prob2.o
