## Test environments
* local OS X install, R 3.6.3
* ubuntu 18.04.5 (on travis-ci), R 3.6.3
* win-builder (devel)

Hello, 
 
 - The only modification was to remove the base package from Imports in the description file
 
There was 1 NOTE:
installed size is 19.3Mb
     sub-directories of 1Mb or more:
     libs 18.9Mb

It seems that this is a common note among packages using C++ through the Rcpp library.

