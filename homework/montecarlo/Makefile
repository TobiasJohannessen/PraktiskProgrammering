CXX = g++
# -O3 for optimization in timing measurements
CXXFLAGS = -std=c++17 -Wall -Wextra -O3
# Define the path to your include files
INCLUDE_DIR = ../includes
# Add the include directory to the compiler flags
CPPFLAGS = -I$(INCLUDE_DIR)

.PHONY: all clean out 

all: main plot.errors.circle

# Rule to build the 'main' executable (for general functionality testing, if any)
main: main.o 
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -o $@ $^
	./main > out.txt


main.o: main.cpp 
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $<




plot.errors.circle:
	echo ' \
		set terminal svg size 600,600 background rgb "white"; \
		set title "Monte Carlo Circle Area Estimation Errors (Logarithmic)"; \
		set output "plots/circle_errors, log.svg"; \
		set key left ;\
		set xlabel "Sample Size";\
		set ylabel "Error";\
		set logscale x; \
		set logscale y; \
		set xrange [90:1000000]; \
		f(x) = a/x**0.5; \
		fit f(x) "data/circle_errors.txt" using 1:3 via a; \
		plot  \
			"data/circle_errors.txt" using 1:3 with lines title "Error", \
			f(x) with lines title sprintf("Fitted O(1/\sqrt(N)) (a=%.2e)", a); \
	' | gnuplot
	echo ' \
		set terminal svg size 600,600 background rgb "white"; \
		set output "plots/circle_errors.svg"; \
		set key left ;\
		set xlabel "Sample Size";\
		set ylabel "Error";\
		set xrange [0:1000000]; \
		f(x) = a/x**0.5; \
		fit f(x) "data/circle_errors.txt" using 1:3 via a; \
		plot  \
			"data/circle_errors.txt" using 1:3 with lines title "Error", \
			f(x) with lines title sprintf("Fitted O(1/sqrt(N)) (a=%.2e)", a); \
	' | gnuplot
clean:
	rm -f *.o main data/*.txt $(INCLUDE_DIR)/*.o plots/* plots/*/*.svg
