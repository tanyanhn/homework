CC := g++

CPPFLAGS = -O0 -g -D_GLIBCXX_DEBUG --std=c++11
MATLAB := matlab -nodesktop -nosplash -logfile `date +%Y_%m_%d-%H_%M_%S`.log -r

GSLLIB = -lgsl -lgslcblas -lm

LAPACKELIB = -llapacke -llapack -lblas

.PHONY: clean all
.SUFFIXES: .o .cpp .ex

.cpp.o:
	$(CC) -o $@ -c $< $(CPPFLAGS) $(LAPACKELIB) $(GSLLIB)
.o.ex:
	$(CC) -o $@ $^ $(CPPFLAGS) $(LAPACKELIB) $(GSLLIB)
.cpp.ex:
	$(CC) -o $@ $< $(CPPFLAGS) $(LAPACKELIB) $(GSLLIB)

all: main.ex

main.ex: main.cpp gsl_function.h CubicSpline.h
	$(CC) -o $@ $< $(CPPFLAGS) $(LAPACKELIB) $(GSLLIB)
# main.ex: main.o

run: main.ex p.m
	./main.ex output1 output2 output3 output4 output5 output6
	$(MATLAB) p

bsp1: bsp1.m plotTruncatedPowerFunc2Bsplines.m
	$(MATLAB) bsp1

bsp2: bsp2.m plotTruncatedPowerFunc2Bsplines.m
	$(MATLAB) bsp2

clear:
	-rm main.ex *.o core* *.log

clean:
	-rm -f main.ex *.output *.o *.h.gch *.obj core* test.ex *.log
	-rm -rf ./doxygen/
