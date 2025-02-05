

Out.txt: hello.o
	$(CXX) hello.o -o Out.txt

hello.o: hello.cpp
	$(CXX) -c hello.cpp

clean:                           # a phoney target, no dependencies
	rm -f Out.txt hello.o  # remove secondary files


#test:
#	echo $(CXX)
