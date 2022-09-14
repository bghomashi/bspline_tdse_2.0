IDIR =  inc/
SDIR =  src/
ODIR =	obj/
DDIR =	dep/
#csrc =  $(wildcard src/*.c)

ALLCCSRC = $(shell find $(SDIR) -name "*.cpp")
PRGS = $(notdir $(patsubst %_main.cpp,%,$(filter %main.cpp,$(ALLCCSRC))))

csrc =  $(shell find $(SDIR) -name "*.c")
ccsrc = $(filter-out %main.cpp,$(ALLCCSRC))

OBJ = 	$(subst $(SDIR), $(ODIR), $(csrc:.c=.o)) \
		$(subst $(SDIR), $(ODIR), $(ccsrc:.cpp=.o))
DEP =	$(subst $(ODIR), $(DDIR), $(OBJ:.o=.d))





PRG_SUFFIX = .out
BINS = $(patsubst %,%$(PRG_SUFFIX),$(PRGS))


$(info Executables="$(BINS)")

#preprocessor flags
EXTRAFLAGS = -g -O0 
#CPPFLAGS =  -DGLEW_STATIC -Iinc -IGL -I../freetype2 -I/usr/local/include/freetype2
CPPFLAGS =      -D__STDCPP_WANT_MATH_SPEC_FUNCS__ \
				-DCL_TARGET_OPENCL_VERSION=300 \
				-I$(PETSC_DIR)/include -I$(PETSC_DIR)/$(PETSC_ARCH)/include \
				-I$(SLEPC_DIR)/include -I$(SLEPC_DIR)/$(PETSC_ARCH)/include \
				-Iinc -I/usr/include \
				-Isrc -std=c++14
#linker flags (directories)
#LDFLAGS
#linker libraries
LDLIBS =-lm -lOpenCL -pthread -lgsl -lgslcblas -lhdf5 -L$(SLEPC_DIR)/$(PETSC_ARCH)/lib -L$(PETSC_DIR)/$(PETSC_ARCH)/lib -lpetsc -lslepc
#LDLIBS =-lm -lGL -lGLU -lglfw -lX11 -ldl -lfreetype -pthread -llapacke -lblas

CC = $(PETSC_DIR)/$(PETSC_ARCH)/bin/mpicxx
CPP = $(PETSC_DIR)/$(PETSC_ARCH)/bin/mpicxx
CXX = $(PETSC_DIR)/$(PETSC_ARCH)/bin/mpicxx

all: $(BINS)

#bspline_tdse.out: $(OBJ)
#	$(CXX) $(EXTRAFLAGS) -o $@ $^ $(LDLIBS)

%.out: $(ODIR)%_main.o $(OBJ)
	$(CXX) $(EXTRAFLAGS) -o $@ $^ $(LDLIBS)

-include $(DEP)							#include all dep files in the makefile

.PHONY: clean
clean:
	rm -f $(OBJ) bspline_tdse.out
	rm -f $(DEP)


%/ : 
	@mkdir -p $@

.SECONDEXPANSION:
$(DDIR)%.d: $(SDIR)%.cpp | $$(dir $$@)
	@$(CPP) $(CPPFLAGS) $(EXTRAFLAGS) $< -MM -MT $(subst $(DDIR), $(ODIR), $(@:.d=.o)) >$@


$(ODIR)%.o : $(SDIR)%.c | $$(dir $$@)
	$(CC) $(CPPFLAGS) $(EXTRAFLAGS) -c $< -o $@

$(ODIR)%.o : ../%.c | $$(dir $$@)
	$(CC) $(CPPFLAGS) -c $< -o $@

$(ODIR)%.o : $(SDIR)%.cpp | $$(dir $$@)
	$(CXX) $(CPPFLAGS) $(EXTRAFLAGS) -c $< -o $@


