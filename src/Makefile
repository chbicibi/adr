FC      = gfortran
FCFLAGS = -O2 -fopenmp
TARGET  = ..\adr.exe
OBJS    = main.o base.o base_func.o base_indiv.o base_group.o         \
					ga_main.o ga_base.o ga_unit.o solver_orbit.o solver_debri.o \
          io_path.o io_input.o io_output.o

TARTEST = ..\test.exe
OBJTEST = xtest.o base.o base_func.o base_indiv.o solver_orbit.o solver_debri.o \
					io_path.o io_input.o io_output.o  xmod.o

.SUFFIXES: .f08
.PHONY:    clean

all: $(TARGET)

$(TARGET): $(OBJS)
	$(FC) $(FCFLAGS) -o $@ $^

%.o: %.f08
	$(FC) $(FCFLAGS) -c $<

%.mod: %.f08 %.o
	@:

clean:
	del *.o

main.o: base.mod base_indiv.mod base_group.mod ga_main.mod ga_base.mod io_input.mod

base.o:       io_input.mod io_output.mod solver_debri.mod
base_indiv.o: base_func.mod solver_debri.mod solver_orbit.mod
base_group.o: base_func.mod solver_debri.mod base_indiv.mod

ga_main.o:               base_indiv.mod base_group.mod ga_base.mod ga_unit.mod
ga_base.o: base_func.mod
ga_unit.o: base_func.mod base_indiv.mod base_group.mod ga_base.mod io_output.mod

solver_orbit.o: base_func.mod io_output.mod
solver_debri.o: base_func.mod solver_orbit.mod io_path.mod io_input.mod
class_orbit.o: solver_orbit.mod solver_debri.mod io_output.mod

io_input.o:  io_path.mod
io_output.o: io_path.mod base_func.mod# base_indiv.mod solver_orbit.mod solver_debri.mod


test: $(TARTEST)
$(TARTEST): $(OBJTEST)
	$(FC) $(FCFLAGS) -o $@ $^
xtest.o: base_func.mod io_path.mod  base.mod solver_orbit.mod solver_debri.mod io_output.mod xmod.mod

xtest1.o: base_func.mod io_path.mod base.mod solver_orbit.mod solver_debri.mod io_output.mod

xmod.o: base_func.mod solver_orbit.mod solver_debri.mod io_output.mod
