CC := g++
CPPFLAGS = -O0 -g -D_GLIBCXX_DEBUG -std=c++17
MATLAB := matlab -nodesktop -nosplash -logfile Y.log -r
GSLLIB = -lgsl -lgslcblas -lm
LAPACKELIB = -llapacke -llapack -lblas

VPATH = ../src
TESTFLAG =

TARGET := test.ex 

.PHONY: clean clear test run all doxygen
.SUFFIXES: .o .cpp .ex


.cpp.o:
	$(CC) -o $@ -c $< $(CPPFLAGS) $(LAPACKELIB) $(GSLLIB)
.hpp.o:
	$(CC) -o $@ -c $< $(CPPFLAGS) $(LAPACKELIB) $(GSLLIB)
.o.ex:
	$(CC) -o $@ $^ $(CPPFLAGS) $(LAPACKELIB) $(GSLLIB)
.cpp.ex:
	$(CC) -o $@ $< $(CPPFLAGS) $(LAPACKELIB) $(GSLLIB)


all:$(TARGET)

SRCS = test_Inverse.cpp test_ThreeBody.cpp ThreeBody.cpp TimeIntegrator.cpp\
	catch_with_main.cpp test_ABmethod.cpp norm.cpp Inverse.cpp test_AMmethod.cpp test_BDmethod.cpp
OBJS = $(patsubst %cpp,%o,$(SRCS))


$(TARGET): $(OBJS)
	$(CXX) $(CFLAGS) -o $@ $^ $(CPPFLAGS) $(LAPACKELIB)
#	$(CC) -o $@ $^ $(CPPFLAGS) $(LAPACKELIB) $(GSLLIB)


%.o:%.cpp
	$(CXX) $(CFLAGS) -c $< $(CPPFLAGS)  $(LAPACKELIB)

%.d:%.cpp
	@set -e; rm -f $@; $(CC) -MM $< $(INCLUDEFLAGS) > $@.$$$$; \
	sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
	rm -f $@.$$$$

-include $(OBJS:.o=.d)

test:$(TARGET)
	./$(TARGET) $(TESTFLAG)

doxygen: doxygen.cfg
	doxygen doxygen.cfg

clear:
	-rm $(TARGET) *.o core* *.log *.ex -f *.d

clean:
	-rm -f $(TARGET) *.output *.o *.h.gch *.obj core* test.ex *.log *.ex *.d
	-rm -rf ./doxygen/
