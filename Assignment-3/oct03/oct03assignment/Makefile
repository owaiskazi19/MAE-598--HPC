#Adds source files
#Adds objects file which are created
#Library needed	
#Defines compiler
#Flags needed		
#Creates executable file of name 'wada'
#Adds dependecies
SOURCES = wada.cc pendulum.cc rk4.cc rk4.h pendulum.h 
OBJECTS = wada.o pendulum.o rk4.o				
LIBS = -lm												
CXX = g++ 								
CXXFLAGS = -std=c++11 -O -c -Wall 									

wada: $(OBJECTS) 
	$(CXX) $(OBJECTS) $(LIBS) -o wada

wada.o: rk4.h pendulum.h								
rk4.o: rk4.h
pendulum.o: rk4.h pendulum.h

clean:													
		@rm -rf $(OBJECTS) a.out core wada