file = file.o file.h
common = common.h 
deTruss = deTruss.o deTruss.h
calPek = calPek.o calPek.h
myG = myG.o myG.h

APXTarget: APX.o $(file) $(common) $(deTruss) $(calPek) $(myG)
	g++ -rdynamic -o APX.out APX.o file.o deTruss.o calPek.o myG.o
	
APX.o: APX.cpp file.h $(common) deTruss.h calPek.h myG.h
	g++ -rdynamic -c APX.cpp

file.o: file.cpp $(common) file.h myG.h
	g++ -rdynamic -c file.cpp
	
deTruss.o: deTruss.cpp $(common) deTruss.h myG.h calPek.h file.h
	g++ -rdynamic -c deTruss.cpp
	
calPek.o: calPek.cpp calPek.h deTruss.h $(common) myG.h
	g++ -rdynamic -c calPek.cpp
	
myG.o: myG.cpp myG.h $(common)
	g++ -rdynamic -c myG.cpp
	
clean:
	rm *.o