#Adds source files
SOURCES = wada.cc pendulum.cc rk4.cc rk4.h pendulum.h 
#Adds objects file which are created
OBJECTS = wada.o pendulum.o rk4.o		
#Library needed			
LIBS = -lm		
#Defines compiler										
CC = g++ -std=c++11			
#Flags needed							
OPTIONS = -O -c -Wall									

#Creates executable file of name 'wada'
wada: $(OBJECTS) 
			$(CC) $(OBJECTS) $(LIBS) -o wada

#Adds dependecies
wada.o: rk4.h pendulum.h								
rk4.o: rk4.h
pendulum.o: rk4.h pendulum.h

clean:													
		@rm -rf $(OBJECTS) a.out core wada