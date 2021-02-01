CC := g++
CPPFLAGS = -O0 -g -D_GLIBCXX_DEBUG --std=c++11
MATLAB := matlab -nodesktop -nosplash -logfile Y.log -r

.PHONY: clean clear test run all doxygen run4
.SUFFIXES: .o .cpp .ex




run:
	cd src/ \
	&& make run

run4:
	cd src/ \
	&& make run4

test:
	make -C test/ test


clear:
	cd src/  \
	&& make clear 
	make -C test/  clear 

clean:
	cd src/  \
	&& make clean
	make -C test/ clean
