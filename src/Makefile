PKG_CC=clang++
unamestr := $(shell uname)
ifeq ($(unamestr),Linux)
	CFLAGS=-g -O3 -DNDEBUG -fpic -pipe -march=x86-64 -ffast-math -fopenmp -stdlib=libc++
endif
ifeq ($(unamestr),Darwin)
	CFLAGS=-g -O3 -DNDEBUG -fpic -pipe -march=x86-64 -ffast-math -stdlib=libc++
endif
#CFLAGS=-O2 -march=x86-64 -ffast-math -DHAVE_INLINE -DGSL_RANGE_CHECK_OFF -g -pg -fopenmp
LINKER=-L${CURDIR}/lib/ -Wl,-rpath,${CURDIR}/lib/

INCLUDE=-isystem${CURDIR}/ -isystem${CURDIR}/include/ -isystem${CURDIR}/include/R/ -w
#LIBS=-lstdc++ -lgfortran -lgsl -lgslcblas -lpthread -lm -static
LIBS=-lgsl -lgslcblas -lm 


all: install myHeader.o BicMix.so clean

BicMix.so: BicMix_c.o
	$(PKG_CC) -shared BicMix_c.o -o BicMix.so $(LINKER) $(LIBS)
BicMix_c.o:  BicMix_c.cpp myHeader.cpp
	$(PKG_CC) $(CFLAGS) $(INCLUDE) -c  BicMix_c.cpp
myHeader.o: myHeader.cpp
	$(PKG_CC) $(CFLAGS) $(INCLUDE) -c  myHeader.cpp
install:
	./install.sh
clean:
	rm *.o

