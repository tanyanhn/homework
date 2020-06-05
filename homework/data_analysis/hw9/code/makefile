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

all: main.ex



#main : main.o
#	$(CC) -o $@ $<

main.ex: main.o normalequations.o QR_factorization.o

main.o: main.cpp Discrete_least_square.h normalequations.h

normalequations.o: normalequations.cpp normalequations.h Discrete_least_square.h

QR_factorization.o: QR_factorization.cpp QR_factorization.h Discrete_least_square.h


run: main.ex
	./main.ex
	$(MATLAB) calnorm


doxygen: doxygen.cfg
	doxygen doxygen.cfg

clear:
	-rm main.ex *.o core* *.log *.ex

clean:
	-rm -f main.ex *.output *.o *.h.gch *.obj core* test.ex *.log *.ex
	-rm -rf ./doxygen/
