# (c) 2007 The Board of Trustees of the University of Illinois.

# Default language wide options

include ../../../Makefile.conf

OBJS=main.o 
BIN = bfs

APP_LDFLAGS=-lm -lstdc++
APP_CFLAGS += -w -ffast-math
APP_CXXFLAGS += -ffast-math

APP_CUDACFLAGS += -w -Xptxas -dlcm=cg --use_fast_math --ptxas-options=-v

# To use slow (but precise) math routines
#ifeq ($(slowmath), 1)
#	APP_CFLAGS += -fp-model precise
#else
#endif

# To use the MKL library
# ifeq ($(mkl), 1)
# 	APP_LDFLAGS += -L/opt/intel/mkl/9.1.023/lib/32 \
# 		-lmkl -lvml -lguide -lpthread
# 	INCLUDES += -I/opt/intel/mkl/9.1.023/include
# 	APP_CFLAGS += -DMKL
# endif


default: $(BIN)

$(BIN) : $(OBJS) parboil_cuda.o
	$(CUDACC) $^ -o $@ $(CUDALDFLAGS)

%.o : %.cu
	$(CUDACC) $< $(CUDACFLAGS) -w -c -o $@

%.o : %.c
	$(CC) $(CFLAGS) -c $< -o $@

parboil_cuda.o: $(PARBOIL_ROOT)/common/src/parboil_cuda.c
	$(CC) $(CFLAGS) -c $< -o $@

%.o : %.cc
	$(CXX) $(CXXFLAGS) -c $< -o $@

%.o : %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean :
	rm -rf *.o $(BIN)


