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

LIBPATH = C:\yam\Studies\Programs\Fortran\FEC\src\lib ..\liborbit
LIBNAME = fec orbit
LIBOPT  = $(addprefix -L,$(LIBPATH)) $(addprefix -l,$(LIBNAME))
INCOPT  = $(addprefix -I,$(LIBPATH))
LIBBLAS = -llapack -lrefblas
INCLUDE = ..\include


################################################################################

.SUFFIXES: .f08
.PHONY: clean

all: $(TARGET)

$(TARGET): $(OBJS)
	$(FC) $(FCFLAGS) -o $@ $^ $(LIBOPT) $(LIBBLAS)

%.o: %.f08
	$(FC) $(FCFLAGS) -c -o $@ $< $(INCOPT)

clean:
	rm -f *.o *.mod *.smod

################################################################################

main.o: orbit_optimizer.o
main3.o: class_prob1.o class_opt1.o class_opt2.o class_opt3.o
main4.o: class_prob1.o class_prob2.o class_opt1.o class_opt2.o

class_prob2.o: class_opt1.o
class_opt1.o: class_prob1.o
class_opt2.o: class_prob2.o
