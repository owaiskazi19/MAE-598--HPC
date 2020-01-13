SOURCES = rk4.cc pendulum.cc broyden.cc
OBJECTS = rk4.o pendulum.o broyden.o
LIBS = -lm
CXX = g++
CXXFLAGS = -std=c++11 -O3 -c -Wall

broyden: $(OBJECTS)
	$(CXX) $(OBJECTS) $(LIBS) -o broyden

broyden.o: rk4.h pendulum.h
rk4.o: rk4.h
pendulum.o: rk4.h pendulum.h

clean:
	@rm -rf $(OBJECTS) a.out core broyden
