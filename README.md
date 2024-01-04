# E-WMMSE
This is the code implementation for the E-WMMSE algorithm.
CLick here for the original paper link:  
[A Novel Extrapolation Technique to Accelerate WMMSE](https://ieeexplore.ieee.org/document/10096806)  
# Code Introduction
WMMSE.m : The main function for the [WMMSE algorithm](https://ieeexplore.ieee.org/document/5756489/).  
E-WMMSE.m : The main function for the [E-WMMSE algorithm](https://ieeexplore.ieee.org/document/10096806).    
find_U.m : The function for finding the U in each iteration.  
find_W.m : The function for finding the W in each iteration.  
find_V.m : The function for finding the V in each iteration.   
sumrate.m : The function for computing the sum rate.  
Test_E_WMMSE.m : This is a function used to test E-WMMSE performance, enter the required parameters and the function would return the number of iterations, running time and sum rate information  
Test_WMMSE.m : This is a function used to test WMMSE performance, enter the required parameters and the function would return the number of iterations, running time and sum rate information  

Test.m : This script is used to access the performance gap between the two algorithms WMMSE and E-WMMSE. The indicators include running time, number of iterations and final sum rate. Currently this script only supports the simulation scenario of a single base station.  
figs : The folder that stores the results in different scenario configurations.  
# Result
Run Test.m in matlab and get the following figures, one for running time and the other for sum rate:  
![Running time comparison](/figs/result1.png)  
![Sum rate comparison](/figs/result2.png)  

# Computer specs：
CPU : 13600K (5.3 GHz, 6 Performance-cores, 8 Efficient-cores)  
Motherboard : ASUS PRIME Z790-P
DRAM : 64G DDR5 6000MHz (KINGBANK)  
Disk : 2T (SHPP41-2000GM)  
GPU : NVIDIA Geforce RTX 4070    
OS : Windows 11 Pro (23H2)  
MATLAB Version : R2023a  

