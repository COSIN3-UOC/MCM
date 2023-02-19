#  MCM: Microscopic congestion model

The code here reproduces partially the experiments published for the Microscopic Congestion Model developed by Solé-Ribalta et al, Royal Society Open Science 2016 (https://doi.org/10.1098/rsos.160098)
        
## System Requirements 	
	
- Matlab R2018a or above
- C99 and C++ compiler 
- Boost C++ library (https://www.boost.org/)

## Structure and setup of the MCM packet

The structure is very simple, a main file that runs the experiments "toy_directed_weighted_nonHomStarts.m" and a directory “mcm” with auxiliary files. Computationally costly functions are implemented in C and C++ and need to be compiled before running the toy example. 

There are three files that need to be compiled: cSPCongestion_statMem_dir_weighted_local.c, sp_computeMCMParams_dir_weighted_dynObsBwComp_C_aux.cpp, SPDataPathDegeneration_staticMemory.c withing the Matlab interface. One needs to make sure that Matlab is correctly configured to compile C and C++ files, usually one does that using command "mex -setup c++" but see the documentation [here](https://es.mathworks.com/help/matlab/matlab_external/changing-default-compiler.html) if required. Compilation line is indicated within the C and C++ files (first line).

File SPDataPathDegeneration_staticMemory.c defines the structures that will be used to run the Montecarlo simulations (ABM). To compile use:

```sh
mex -largeArrayDims SPDataPathDegeneration_staticMemory.c
```
### Configuration of the C file
Definition "#define MAX_COLA 2000" within the C file needs to be adjusted to be at least the number of nodes of the network. That is, if the size of network used is 2000 the MAX_COLA must have a value larger than 2000. We recall the user to recompile the file once modified.

File cSPCongestion_statMem_dir_weighted_local is the main file that runs the ABM simulation. NOTE: some but withing the code, in some runs it basically does a segmentation fault.

```sh
mex -largeArrayDims cSPCongestion_statMem_dir_weighted_local.c
```
### Configuration of the C file
Definition "#define MAX_COLA 501000" withing this C sets the maximum number of packets that will be stored within each node queue. This value depends on the simulation time and the degree of congestion of the system. As set the value should be enough for long simulation times. We recall the user to recompile the file once modified.


File sp_computeMCMParams_dir_weighted_dynObsBwComp_C_aux an auxiliar file that compute the theoretical values of the model. This file is used by SP_computeEtaGivenRho_dir_weigh_nonHomStart, the main function to obtain the theoretical predictions of the model.

```sh
mex -largeArrayDims sp_computeMCMParams_dir_weighted_dynObsBwComp_C_aux.cpp -I/usr/local/Cellar/boost/1.80.0/include/
```

### Configuration of the C file
Definitions "#define MAX_COLA 2000" and "#define MAX_NODES 2000" within the C file needs to be adjusted to be at least the number of nodes of the network. That is, if the size of network used is 2000 the MAX_COLA must have a value larger than 2000. We recall the user to recompile the file once modified.

File SPEdgeNodeBetweennessC_BrandesWeightedTestEdgeBW_local_fheap computes the betweeness of the network. The function is used to compute the theoretical critical injection rate. 

```sh
mex -largeArrayDims SPEdgeNodeBetweennessC_BrandesWeightedTestEdgeBW_local_fheap.cpp -I/usr/local/Cellar/boost/1.80.0/include/
```

### Configuration of the C file
Definitions "#define MAX_COLA 2000" and "#define MAX_NODES 2000" within the C file needs to be adjusted to be at least the number of nodes of the network. That is, if the size of network used is 2000 the MAX_COLA must have a value larger than 2000. We recall the user to recompile the file once modified.

See that the theoretical model requires to use some boost libraries, so one need to link to the version installed in the computer.

**Additionally, one can do a batch compilation of the C files by using the function "batchCompileCFiles.m" within the mcm directory.**

## Use examples: 

Once configured and compiled, running the toy example very simple. Just run the file "toy_directed_weighted_nonHomStarts.m". The file runs a toy example to obtain a plot like the one in Figure 2 (panel c) of [1]. Similar plots for the microscopic variables could be obtained working out the different outputs of the given by functions cSPCongestion_statMem_dir_weighted_local and SP_computeEtaGivenRho_dir_weigh_nonHomStart. The toy example is configured to run the analysis on a weighted non-symmetric transportation network of 50 nodes. 

The code should output something like the following plot.
<img src="https://github.com/COSIN3-UOC/MCM/blob/main/outputExample.png?raw=true" width="600">


# Citations
Solé-Ribalta, A., Gómez, S., & Arenas, A. (2016). A model to identify urban traffic congestion hotspots in complex networks. Royal Society open science, 3(10), 160098.

Solé-Ribalta, A., Gómez, S., & Arenas, A. (2018). Decongestion of urban areas with hotspot pricing. Networks and Spatial Economics, 18, 33-50.
