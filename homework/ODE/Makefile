CXX = g++
# -O3 for optimization in timing measurements
CXXFLAGS = -std=c++17 -Wall -Wextra -O3
# Define the path to your include files
INCLUDE_DIR = ../includes
# Add the include directory to the compiler flags
CPPFLAGS = -I$(INCLUDE_DIR)

.PHONY: all clean out run_test_timing

all: main 
	./main

# print the output to a file
out: main
	./main > out.txt

# Rule to build the 'main' executable (for general functionality testing, if any)
main: main.o $(INCLUDE_DIR)/matrix.o $(INCLUDE_DIR)/ODE.o
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -o $@ $^

main.o: main.cpp $(INCLUDE_DIR)/matrix.h
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $<

# Rule to build the 'matrix.o' object from the matrix implementation file
$(INCLUDE_DIR)/matrix.o: $(INCLUDE_DIR)/matrix.cpp $(INCLUDE_DIR)/matrix.h
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $< -o $@

$(INCLUDE_DIR)/ODE.o: $(INCLUDE_DIR)/ODE.cpp $(INCLUDE_DIR)/ODE.h
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $< -o $@


plot.exponential: exponential.txt
	echo ' \
		set terminal svg background rgb "white"; \
		set output "plots/exponential.svg"; \
		set key left ;\
		set xlabel "x";\
		set ylabel "y";\
		plot  \
			"exponential.txt" using 1:2 with lines title "Exponential", \
	' | gnuplot





clean:
	rm -f *.o main out out.txt $(INCLUDE_DIR)/*.o
