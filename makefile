IDIR=${CURDIR}/src

	
program=novaCTF

SOURCES := $(shell find $(IDIR) -name '*.cpp')

CXX = g++ -O3 -s -DNDEBUG
CXXFLAGS = -std=c++0x
LDFLAGS = -L/usr/lib64 

OBJECTS = $(SOURCES:.cpp=.o)

all: ${program}

build: ${program}

debug: CXX = g++ -g -W -Wall -Werror 
debug: LDFLAGS += -g
debug: ${program}

 
${program}: CXXFLAGS += $(foreach d, $(includepath), -I$d)
${program}: LDFLAGS += $(foreach d, $(libpath), -L$d)
${program}: $(OBJECTS)
	$(CXX) -o $@ $^ $(LDFLAGS) -lfftw3 -lfftw3f

clean:
	rm -f src/*.o ${program} src/*.d

.cpp.o:
	$(CXX) -MD -MP $(CXXFLAGS) -o $@ -c $< $(LDFLAGS)
