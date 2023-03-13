# README #

## File Organization ##

* **APXUpdate.cpp**: main function, a sample program
* **calPek.h, calPek.cpp**: implementation of building/updating support vector
* **deTruss.h, deTruss.cpp**: implementation of peeling algorithms to calculate probabilistic trussness
* **myG.h, myG.cpp**: a class for graph management
* **file.h, file.cpp**: implementation of file I/O
* **common.h**: some global definition
* **randomP.cpp**: a program to generate queries that represent the dynamic updates of the graph

## How to Use the Code? ##

* First, use "make" to get the executable file "randomP.out" and "APXUpdate.out".
* Then, run "./randomP.out indexPath queryPath sampleNumber type batchNumber " to generate queries.
* "indexPath" is the index of the graph. "queryPath" is the output folder and must exist. "sampleNumber" is the number of edges sampled. "type" is either 2 or -2 and represents existing probability increase/decrease. "batchNumber" is the number of queries generated.
* Example:
*       ./randomP.out ./index ./query 100 -2 1
* Next, using the index and the query file as input, you can update the graph and get the new index. Run "./APXUpdate.out indexPath queryFile resultPath".
* "indexPath" is the index of the graph. "queryFile" is the query file. "resultPath" is the output folder and must exist, but this parameter is optional. 
* Example:
*       ./APXUpdate.out ./index query/query_-20_100_1.txt
* Or:
*       ./APXUpdate.out ./index query/query_-20_100_1.txt ./result

## Data Format ##

* The file "queryFile" containing the update information is of the following format:

    u_1,v_1,p_1
    u_2,v_2,p_2
    ...

where (u_i, v_i) is the edge in the graph and p_i is the updated existing probability.
