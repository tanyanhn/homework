CC := g++
CPPFLAGS = -O0 -g -D_GLIBCXX_DEBUG --std=c++11
MATLAB := matlab -nodesktop -nosplash -logfile `date +%Y_%m_%d-%H_%M_%S`.log -r

.PHONY: clean clear test run all doxygen print
.SUFFIXES: .o .cpp .ex


.cpp.o:
	$(CC) -o $@ -c $< $(CPPFLAGS) -g

.o.ex:
	$(CC) -o $@ $^ $(CPPFLAGS) -g

all: main.ex

main.ex: main.o function.o

#main : main.o
#	$(CC) -o $@ $<
main.o: interpolation.cpp function.h
	$(CC) -o $@ -c $< $(CPPFLAGS) -g

function.o: function.cpp function.h

run: main.ex
	./main.ex

print: main.ex
	./main.ex
	$(MATLAB) p

doxygen: doxygen.cfg
	doxygen doxygen.cfg

clear:
	-rm main.ex *.o core* *.log

clean:
	-rm -f main.ex *.output *.o *.h.gch *.obj core* test.ex *.log
	-rm -rf ./doxygen/


