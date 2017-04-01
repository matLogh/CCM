CCM, Cross-correlation Correction Method

Matus Balogh, Slovak Academy of Sciences
matus.balogh@savba.sk


Application of discrete cross-correlation to correct time-dependent detector output 

Prerequisities: ROOT (see https://root.cern.ch/building-root)
As an author I want to apologize for the coding itself, I am way too far from being professional programmer.

HOW TO USE:
Please use only public functions listed in Process.h. Example of using CCM is in main.cpp. In function StartCCM specify number of threads you want to use in computation, using more than your CPU have may be counter-productive.
For now there is no MakeFile, to compile use:
"g++ -std=c++0x `root-config --libs --cflags` -c CheckCCM.cpp Process.cpp Cross_correlation.cpp main.cpp -o CCM"
and to run:
"./CCM"
Note, that using a huge matrices (TH2D) is not recomended since the code is making its partial copy of TH2D into C-style array. Thus larger matrix means larger memory usage. 

The code was tested on i5-6200U CPU with 8 GB of RAM on 250 GB SSD drive. Computation time for matrix of 4500 time bins and 2 ROIs is around 16 seconds. 

FOR BUG REPORT OR ENHANCEMENT REQUEST PLEASE DO NOT HESITATE TO WRITE ME AN EMAIL.
