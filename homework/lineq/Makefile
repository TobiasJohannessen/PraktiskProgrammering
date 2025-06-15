CXX = g++
CXXFLAGS = -std=c++11 -Wall -Wextra


main: main.o matrix.o
	$(CXX) $(CXXFLAGS) -o main main.o matrix.o
	./main

main.o: main.cpp ../includes/matrix.h
	$(CXX) $(CXXFLAGS) -c main.cpp

matrix.o: ../includes/matrix.cpp ../includes/matrix.h
	$(CXX) $(CXXFLAGS) -c ../includes/matrix.cpp

clean:
	rm -f *.o out.txt main

out: main
	./main > out.txt
