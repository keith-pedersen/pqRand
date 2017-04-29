# Set PQR_DIR to build pqrRand into your executable (using the %.x rule)
PQR_DIR = .
PQR_INC = $(PQR_DIR)/include

CXX = gcc
MARCH = native
STD = c++11 
# GCC flags, including many useful warnings
STABILITY_FLAGS = -pedantic-errors -fno-common -mfpmath=sse -mieee-fp #sse flag to avoid weird x87 registers (see https://gcc.gnu.org/wiki/FloatingPointMath)
STABILITY_WARNINGS = -Wall -Wextra -W -Wconversion -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -Wmissing-declarations -Wredundant-decls -Wmissing-field-initializers -Wlogical-op -Wunsafe-loop-optimizations -Wwrite-strings -Wundef -Wfloat-equal
PERFORMANCE_FLAGS = -O3 -march=$(MARCH) -msse4 -mavx2 -Winline -Wdisabled-optimization -Wpadded
BUILD_LIB_FLAGS = -fPIC
# 
CXXFLAGS = -std=$(STD) $(STABILITY_WARNINGS) $(PERFORMANCE_FLAGS) $(BUILD_LIB_FLAGS)

# The directory structure of pqRand
INCLUDE = ./include
SOURCE = ./source
EXAMPLES = ./examples

# external dependencies of pqRand
PQR_DEPENDENCIES = -lstdc++ -lm

INC_FLAGS = -I $(INCLUDE)
INC_FLAGS_EXTERN = -I $(PQR_INC)

LIB_FLAGS = $(PQR_DEPENDENCIES)
LIB_FLAGS_EXTERN = $(LIB_FLAGS) -L $(PQR_DIR) -lpqr



FILENAMES = pqRand pqrDistributions
OBJS = $(addsuffix .o, $(addprefix $(SOURCE)/, $(FILENAMES)))

all : libpqr.so $(EXAMPLES)/pqRandExample.x $(EXAMPLES)/random_device_tester.x \
	$(EXAMPLES)/xorshift1024star_tester.x $(EXAMPLES)/xorshift1024star_last3bits_tester.x 

libpqr.so: $(OBJS)
	$(CXX) $(CXXFLAGS) -shared $(OBJS) $(LIBFLAGS) -o $@

%.x : %.cpp
	$(CXX) $(CXXFLAGS) $(INC_FLAGS_EXTERN) $(LIB_FLAGS_EXTERN) $*.cpp -o $@
	
%.o : %.cpp 
	$(CXX) $(CXXFLAGS) $(INC_FLAGS) $(LIBFLAGS) $*.cpp -c -o $*.o
	
.PHONY: clean

clean:
	rm -f $(SOURCE)/*.o
	rm -f $(EXAMPLES)/*.x
