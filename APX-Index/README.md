# README #

## File Organization ##

* **APX.cpp**: main function, a sample program
* **calPek.h, calPek.cpp**: implementation of building/updating support vector
* **deTruss.h, deTruss.cpp**: implementation of peeling algorithms to calculate probabilistic trussness
* **myG.h, myG.cpp**: a class for graph management
* **file.h, file.cpp**: implementation of file I/O
* **common.h**: some global definition

## How to Use the Code? ##

* First, use "make" to get the executable file "APX.out".
* Then, run "./APX.out dataFile valueA valueB resultPath" to construct an index for the graph before update.
* "dataFile" is the graph data file. "valueA" is the first parameter of algorithm APX, represents the minimum probabilistic trussness. "valueB" is the second parameter of algorithm APX, represents the minimum difference between probabilistic trussnesses. "resultPath" is the output folder and it must exist.
* Example:
*       ./APX.out Fruit-Fly.txt 0.001 0.001 ./index

## Data Format ##

* The file "dataFile" containing the graph information is of the following format:

    u_1,v_1,p_1
    u_2,v_2,p_2
    ...

where (u_i, v_i) is the edge in the graph and p_i is the existing probability.
