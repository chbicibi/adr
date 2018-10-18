FC      = gfortran
FCFLAGS = -Wall -O2

TARGET  = ..\main.exe
SRC     = main.f08 sub1.f08 sub2.f08
OBJS    = $(SRC:.f08=.o)

.SUFFIXES: .f08

all: $(TARGET)

$(TARGET): $(OBJS)
	$(FC) $(FCFLAGS) -o $@ $^

%.o: %.f08
	$(FC) $(FCFLAGS) -c $<

%.mod: %.f08 %.o
	@:

main.o: sub1.mod sub2.mod
