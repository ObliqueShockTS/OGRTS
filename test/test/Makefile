# Find all source files, create a list of corresponding object files
SRCS1=$(wildcard *.f90)
OBJS1=$(patsubst %.f90,%.o,$(SRCS1))

SRCS2=$(wildcard *.for)
OBJS2=$(patsubst %.for,%.o,$(SRCS2))

SRCS3=$(wildcard *.f)
OBJS3=$(patsubst %.f,%.o,$(SRCS3))

SRCS=$(SRCS1) $(SRCS2) $(SRCS3)

OBJS=$(OBJS1) $(OBJS2) $(OBJS3)

# Ditto for mods (They will be in both lists)
MODS=$(wildcard bspline*.f90)
MOD_OBJS=$(patsubst %.f90,%.o,$(MODS))

# Compiler/Linker settings
FC = gfortran
FLFLAGS = -g
FCFLAGS = -g -c #-Wall -Wextra -Wconversion -Og -pedantic -fcheck=bounds -fmax-errors=5
PROGRAM = a.out
PRG_OBJ = $(PROGRAM).o

# Clean the suffixes
#.SUFFIXES:

# Set the suffixes we are interested in
#.SUFFIXES: .f90 .o

# make without parameters will make first target found.
default : $(PROGRAM)

# Compiler steps for all objects
$(OBJS1) : %.o : %.f90
	$(FC) $(FCFLAGS) -o $@ $<

$(OBJS2) : %.o : %.for
	$(FC) $(FCFLAGS) -o $@ $<

$(OBJS3) : %.o : %.f
	$(FC) $(FCFLAGS) -o $@ $<

# Linker
$(PROGRAM) : $(OBJS)
	$(FC) $(FLFLAGS) -o $@ $^

debug:
	@echo "SRCS = $(SRCS)"
	@echo "OBJS = $(OBJS)"
	@echo "MODS = $(MODS)"
	@echo "MOD_OBJS = $(MOD_OBJS)"
	@echo "PROGRAM = $(PROGRAM)"
	@echo "PRG_OBJ = $(PRG_OBJ)"

clean:
	rm -rf $(OBJS) $(PROGRAM) $(patsubst %.o,%.mod,$(MOD_OBJS))

.PHONY: debug default clean

# Dependencies

# Main program depends on all modules
$(PRG_OBJ) : $(MOD_OBJS)

# Blocks and allocations depends on shared
cellcenter.o bspline_module.o bspline_sub_module.o bspline_oo_module.o : bspline_kinds_module.o

bspline_oo_module.o bspline_module.o : bspline_sub_module.o

bspline_module.o : bspline_oo_module.o
