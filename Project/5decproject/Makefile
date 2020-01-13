SOURCES = rk4.cc pendulum.cc wada.cc
OBJECTS = rk4.o pendulum.o wada.o
LIBS = -lm
CXX = mpic++
CXXFLAGS = -std=c++11 -O3 -c -Wall

wada: $(OBJECTS)
	$(CXX) $(OBJECTS) $(LIBS) -o wada

wada.o: rk4.h pendulum.h
rk4.o: rk4.h
pendulum.o: rk4.h pendulum.h

clean:
	@rm -rf $(OBJECTS) a.out core wada
