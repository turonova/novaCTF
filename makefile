IDIR=${CURDIR}/src

	
program=novaCTF

SOURCES := $(shell find $(IDIR) -name '*.cpp')

CXX = g++ -O3 -s -DNDEBUG
#CXX = g++ -g -W -Wall -Werror
CXXFLAGS = -std=c++0x -I$(EBROOTFFTW)/include
LDFLAGS = -lfftw3 -lfftw3f -L/usr/lib64 -L$(EBROOTFFTW)/lib 
#LDFLAGS = -lfftw3 -lfftw3f -L/usr/lib64 -L/g/software/linux/pack/fftw-3.3.4/lib -L/g/software/linux/pack/fftw-3.1.2/lib64  
#LDFLAGS += -g

OBJECTS = $(SOURCES:.cpp=.o)
#OBJECTS +=  $(EBROOTFFTW)/lib/libfftw3.a $(EBROOTFFTW)/lib/libfftw3f.a

all: ${program}

build: ${program}

debug: CXX = g++ -g -W -Wall -Werror 
debug: LDFLAGS += -g
debug: ${program}

 
${program}: CXXFLAGS += $(foreach d, $(includepath), -I$d)
${program}: LDFLAGS += $(foreach d, $(libpath), -L$d)
${program}: $(OBJECTS)
	$(CXX) $(LDFLAGS) -o $@ $^

clean:
	rm -f src/*.o ${program} src/*.d

.cpp.o:
	$(CXX) -MD -MP $(CXXFLAGS) -o $@ -c $< $(LDFLAGS)
