

# e3ga makefile for OS X (run make)

all: test


e3ga.cpp: e3ga_cpp.xml
	g25 e3ga_cpp.xml


e3ga_parse_mv.cpp: e3ga.cpp

e3ga_test_main.cpp: e3ga.cpp



#.cpp.o: 

%.o: %.cpp 
	$(CXX) $(CXXFLAGS) -c -I/usr/local/include $< -o $@

test: \
	e3ga_parse_mv.o \
	e3ga.o e3ga_test_main.o
	$(CXX) $(CXXFLAGS) \
	e3ga.o e3ga_parse_mv.o \
	e3ga_test_main.o \
	-lm \
	-o test
	


clean:
	rm -f test *.o *.g e3ga.cpp e3ga.h doxyfile \
	e3gaParser.cpp e3gaLexer.cpp \
	e3ga_parse_mv.cpp e3ga_test_main.cpp \
		 

