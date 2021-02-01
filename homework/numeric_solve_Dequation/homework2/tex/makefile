CC := g++
TEX := xelatex
CPPFLAGS = -O0 -g -D_GLIBCXX_DEBUG -std=c++17
MATLAB := matlab -nodesktop -nosplash -logfile Y.log -r
GSLLIB = -lgsl -lgslcblas -lm
LAPACKELIB = -llapacke -llapack -lblas


story:
	$(TEX) story.tex



doxygen: doxygen.cfg
	doxygen doxygen.cfg

clear:
	-rm $(TARGET) *.o core* *.log *.ex -f *.d *.aux *.synctex.gz

clean:
	-rm -f $(TARGET) *.output *.o *.h.gch *.obj core* test.ex *.log *.ex *.d  *.aux *.synctex.gz
	-rm -rf ./doxygen/
