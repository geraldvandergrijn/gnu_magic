
# Makefile for compiling the HELIOSAT cloudindex code
# on OSX and Linux
# Reto Stockli, MeteoSwiss December 2008

# Use this flag to allow for debugging the code
DEBUG = FALSE
# Use this flag to compile the converter with NetCDF support
# You need to define the environment variable NETCDF_ROOT which
# points to the NETCDF_ROOT/lib and NETCDF_ROOT/include directories
#orig HAS_NETCDF = TRUE
HAS_NETCDF = FALSE
#-----Macros---------------------------------

# Evaluate the Operating System
OS = $(shell uname -s)
# Evaluate the Processor type
MACHINE = $(shell uname -m)

# Check if we are on a supported system
ifeq ($(findstring $(OS), Darwin Linux),)
  $(error This makefile does only support Darwin or Linux at this point. Exiting.)
endif

# Compiler/OS specific options
ifeq ($(OS),Darwin)
# GCC/G++ 4.3.1 on OSX
ifeq ($(HAS_NETCDF),TRUE)
    LIBS = -L$(NETCDF_ROOT)/lib -lnetcdf_c++ -lnetcdf
    INCLUDES = -I. -I$(NETCDF_ROOT)/include -DHAS_NETCDF
else
    LIBS = 
    INCLUDES = -I.
endif
  CFLAGS = -c
  CXXFLAGS = -c
  LDFLAGS =
ifndef CC 
    CC = gcc-4
endif
ifndef CXX
  CXX = g++-4
endif
ifndef LD
  LD = g++-4
endif
endif

ifeq ($(OS),Linux)
# GCC/G++ 4.2 on Linux 32/64 bit
ifeq ($(HAS_NETCDF),TRUE)
  LIBS = -L$(NETCDF_ROOT)/lib -lnetcdf_c++ -lnetcdf
  INCLUDES = -I. -I$(NETCDF_ROOT)/include -DHAS_NETCDF
else
    LIBS = 
    INCLUDES = -I.
endif
  CFLAGS = -c
  CXXFLAGS = -c
  LDFLAGS =
  CC = gcc
  CXX = g++
  LD = g++
endif

# Optimization options
ifeq ($(DEBUG),TRUE)
  CFLAGS += -g2 -ggdb
  CXXFLAGS += -g2 -ggdb
else
  CFLAGS += -O2
  CXXFLAGS += -O2
endif

#-----Suffix Rules---------------------------
# set up C suffixes and relationship between .c and .o files

.SUFFIXES: .cpp .cc

.cpp.o:
	$(CXX) $(CXXFLAGS) $(INCLUDES) $<

.cc.o:
	$(CXX) $(CXXFLAGS) $(INCLUDES) $<

.c :
	$(LD) $(LDFLAGS) $< -o $@ $(LIBS)

#-----File Dependencies----------------------

SRC_CPP = OpenMTP_convert.cpp OpenMTP.cpp OpenMTP_ascii_header.cpp OpenMTP_binary_header.cpp OpenMTP_image.cpp OpenMTP_image_line.cpp OpenMTP_machine.cpp

OBJ = $(SRC_CPP:.cpp=.o) $(SRC_CC:.cc=.o)

OpenMTP_convert: $(OBJ)
	$(LD) $(LDFLAGS) -o $@ $(OBJ) $(LIBS)

#-----Other stuff----------------------------
#depend:
#	makedepend $(CFLAGS) -Y $(SRC)

clean:
	rm -f $(OBJ)