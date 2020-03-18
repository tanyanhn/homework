CC := g++

CPPFLAGS = -O0 -g -D_GLIBCXX_DEBUG --std=c++11

GSLLIB = $(shell gsl-config --libs)

.PHONY: clean all
.SUFFIXES: .o .cpp .ex

.cpp.o:
	$(CC) -o $@ -c $< $(CPPFLAGS)
.o.ex:
	$(CC) -o $@ $^ $(CPPFLAGS)

all: bisection_method.ex Newton_method.ex secant_method.ex

bisection_method.ex: bisection_method.o

Newton_method.ex: Newton_method.o
	$(CC) $^ -o $@ $(GSLLIB)

secant_method.ex: secant_method.o

bisection_method.o: bisection_method.cpp

Newton_method.o: Newton_method.cpp
#	$(CC) $< -c $(GSLLIB)

secant_method.o: secant_method.cpp

clean:
	-rm -f *.o *.ex core*