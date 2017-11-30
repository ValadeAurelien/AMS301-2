compiler = mpic++
flags = -I. -I./eigen/ -g -O0

headers = $(wildcard *.hpp)
sources = $(wildcard *.cpp) $(experimental/wildcard *.cpp)
objects = $(sources:.cpp=.o)

executables: solver

%.o: %.cpp $(headers)
	$(compiler) -c -o $@ $< $(flags)

solver: $(objects)
	$(compiler) -o $@ $^ $(flags)

clean:
	rm -f solver *.o
