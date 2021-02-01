CC := g++
CPPFLAGS = -O0 -g -D_GLIBCXX_DEBUG --std=c++11
MATLAB := matlab -nodesktop -nosplash -logfile Y.log -r

.PHONY: clean clear test run all doxygen
.SUFFIXES: .o .cpp .ex




run:
	make -C test/ test
	cd src/ \
	&& make run

test:
	make -C test/ test

story:
	make -C tex/ story

#doxygen: doxygen.cfg
#	doxygen doxygen.cfg

clear:
	cd src/ \
	&& make clear
	make -C test/ clear
	make -C output/ clear
	make -C tex/ clear
#-rm main.ex *.o core* *.lo

clean:
	cd src/ \
	&& make clean
	make -C test/ clean
	make -C output/ clean
	make -C tex/ clean
