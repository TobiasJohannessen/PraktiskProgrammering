out.txt: main.o main.cpp vec.h
	./main.o > out.txt

main.o: main.cpp vec.h
	c++ main.cpp -o main.o

clean:
	rm -f main.o out.txt