CC := g++
CPPFLAGS = -O0 -g -D_GLIBCXX_DEBUG --std=c++11
MATLAB := matlab -nodesktop -nosplash -logfile Y.log -r

.PHONY: clean clear test run all doxygen
.SUFFIXES: .o .cpp .ex


.cpp.o:
	$(CC) -o $@ -c $< $(CPPFLAGS)
.o.ex:
	$(CC) -o $@ $^ $(CPPFLAGS)



run:
	cd code/ \
&& make run

bsp1:
	make -C code/ bsp1

bsp2:
	make -C code/ bsp2

#doxygen: doxygen.cfg
#	doxygen doxygen.cfg

clear:
	cd code/  \
&& make clear
#-rm main.ex *.o core* *.log

clean:
	cd code/  \
&& make clean
