CC := g++
CPPFLAGS = -O0 -g -D_GLIBCXX_DEBUG -std=c++17
MATLAB := matlab -nodesktop -nosplash -logfile Y.log -r
GSLLIB = -lgsl -lgslcblas -lm
LAPACKELIB = -llapacke -llapack -lblas

VPATH = ../src : ../input

ARGUMENT = ../input/allinput.txt
ARGUMENT4 = ../input/allinput4.txt

TARGET := main.ex 

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

SRCS += $(wildcard ./*.cpp)


OBJS = $(patsubst %cpp,%o,$(SRCS))


$(TARGET): $(OBJS)
	$(CXX) $(CFLAGS) -o $@ $^ $(CPPFLAGS) $(LAPACKELIB)
#	$(CC) -o $@ $^ $(CPPFLAGS) $(LAPACKELIB) $(GSLLIB)


%.o:%.cpp
	$(CXX) $(CFLAGS) -c $< $(CPPFLAGS) $(LAPACKELIB)

%.d:%.cpp
	@set -e; rm -f $@; $(CC) -MM $< $(INCLUDEFLAGS) > $@.$$$$; \
	sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
	rm -f $@.$$$$


-include $(OBJS:.o=.d)

test:$(TARGET)
	./$(TARGET)

run:$(TARGET)
	./$(TARGET) $(ARGUMENT)

run4:$(TARGET)
	./$(TARGET) $(ARGUMENT4)

doxygen: doxygen.cfg
	doxygen doxygen.cfg

clear:
	-rm $(TARGET) *.o core* *.log *.ex -f *.d.* *.d

clean:
	-rm -f $(TARGET) *.output *.o *.h.gch *.obj core* test.ex *.log *.ex *.d.* *.d
	-rm -rf ./doxygen/
