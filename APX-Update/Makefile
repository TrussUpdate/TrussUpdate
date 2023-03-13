file = file.o file.h
common = common.h 
deTruss = deTruss.o deTruss.h
calPek = calPek.o calPek.h
myG = myG.o myG.h

total: APXUpdateTarget randomTarget
	echo -e "\ndone!\n"

randomTarget: randomP.o $(file) $(common) $(myG)
	g++ -rdynamic -o randomP.out randomP.o file.o myG.o

randomP.o: randomP.cpp file.h $(common) myG.h
	g++ -rdynamic -c randomP.cpp
	
APXUpdateTarget: APXUpdate.o $(file) $(common) $(deTruss) $(calPek) $(myG)
	g++ -rdynamic -o APXUpdate.out APXUpdate.o file.o deTruss.o calPek.o myG.o

APXUpdate.o: APXUpdate.cpp file.h $(common) deTruss.h calPek.h  myG.h
	g++ -rdynamic -c APXUpdate.cpp

file.o: file.cpp $(common) file.h myG.h deTruss.h
	g++ -rdynamic -c file.cpp
	
deTruss.o: deTruss.cpp $(common) deTruss.h myG.h calPek.h file.h
	g++ -rdynamic -c deTruss.cpp
	
calPek.o: calPek.cpp calPek.h deTruss.h $(common) myG.h
	g++ -rdynamic -c calPek.cpp
	
myG.o: myG.cpp myG.h $(common)
	g++ -rdynamic -c myG.cpp
	
clean:
	rm *.o