CXX = g++
CXXFLAGS = -std=c++11 -Wall -Wextra

out.txt: main
	./main > out.txt

main: main.o vec.o
	c++ -o main main.o vec.o

main.o: main.cpp vec.h
	c++ -c main.cpp

vec.o: vec.cpp vec.h
	c++ -c vec.o vec.cpp
clean:
	rm -f main.o out.txt vec.o main