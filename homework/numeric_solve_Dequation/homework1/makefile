CC := g++
CPPFLAGS = -O0 -g -D_GLIBCXX_DEBUG --std=c++11
MATLAB := matlab -nodesktop -nosplash -logfile Y.log -r

.PHONY: clean clear test run all doxygen
.SUFFIXES: .o .cpp .ex




run:
	cd exc2/ \
	&& make run

run4:
	cd exc2/ \
	&& make run4


test:
	make -C exc2/ test

anwser:
	make -C tex/ anwser

math:
	make -C tex/ math

design:
	make -C tex/ design

clear:
	cd exc2/  \
	&& make clear
	make -C tex/ clear

clean:
	cd exc2/  \
	&& make clean
	make -C tex/ clean
