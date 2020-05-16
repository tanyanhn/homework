CC := g++
CPPFLAGS = -O0 -g -D_GLIBCXX_DEBUG --std=c++11
MATLAB := matlab -nodesktop -nosplash -logfile Y.log -r
GSLLIB = -lgsl -lgslcblas -lm
LAPACKELIB = -llapacke -llapack -lblas

.PHONY: clean clear test run all doxygen
.SUFFIXES: .o .cpp .ex


.cpp.o:
	$(CC) -o $@ -c $< $(CPPFLAGS) $(LAPACKELIB) $(GSLLIB)
.o.ex:
	$(CC) -o $@ $^ $(CPPFLAGS) $(LAPACKELIB) $(GSLLIB)
.cpp.ex:
	$(CC) -o $@ $< $(CPPFLAGS) $(LAPACKELIB) $(GSLLIB)

all: FPN.ex machinecompute.ex

FPN.ex: FPN.cpp

machinecompute.ex: machinecompute.cpp

#main : main.o
#	$(CC) -o $@ $<
main.o: main.cpp


machinecompute: machinecompute.ex
	./machinecompute.ex
	$(MATLAB) fgh

FPN: FPN.ex
	./FPN.ex
	$(MATLAB) F

doxygen: doxygen.cfg
	doxygen doxygen.cfg

clear:
	-rm main.ex *.o core* *.log *.ex

clean:
	-rm -f main.ex *.output *.o *.h.gch *.obj core* test.ex *.log *.ex
	-rm -rf ./doxygen/
