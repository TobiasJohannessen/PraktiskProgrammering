sfuns.o : sfuns.cpp sfuns.h
	$(CXX) $(CXXFLAGS) -c sfuns.o sfuns.cpp #-c stands for compile only

vec.o: vec.cpp vec.h sfuns.h
	c++ -c vec.o vec.cpp 

clean :
	$(RM)  *.o 
