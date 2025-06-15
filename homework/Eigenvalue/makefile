CXX = g++
# -O3 for optimization in timing measurements
CXXFLAGS = -std=c++11 -Wall -Wextra -O3
# Define the path to your include files
INCLUDE_DIR = ../includes
# Add the include directory to the compiler flags
CPPFLAGS = -I$(INCLUDE_DIR)

.PHONY: all clean out test hydrogen


all: hydrogen test
	
test: test.exe
	./test.exe

hydrogen: hydrogen.exe
	./hydrogen.exe

test.exe: test.o matrix.o
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -o $@ $^


test.o: test.cpp $(INCLUDE_DIR)/matrix.h
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $<

hydrogen.exe: hydrogen.o matrix.o
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -o $@ $^

hydrogen.o: hydrogen.cpp $(INCLUDE_DIR)/matrix.h
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $<

# Rule to build the 'matrix.o' object from the matrix implementation file
matrix.o: $(INCLUDE_DIR)/matrix.cpp $(INCLUDE_DIR)/matrix.h
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $< -o $@


out: test.exe hydrogen.exe
	./test.exe > test_output.txt
	./hydrogen.exe > hydrogen_output.txt

# Clean up compiled objects and executables
clean:
	rm -f *.o main *.txt
