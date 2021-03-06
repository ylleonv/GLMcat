## Test environments
* local OS X install, R 3.6.3
* ubuntu 18.04.5 (on travis-ci), R 3.6.3
* win-builder (devel)

Hello, I hereby submit version 0.2.0 of GLMcat
 
 - The only modification was an update to the design matrix of the function DiscreteCM()

- I also fixed the bug detected by r-patched-solaris-x86 and the version of the package.
 
There was 1 NOTE:
installed size is 19.3Mb
     sub-directories of 1Mb or more:
     libs 18.9Mb

It seems that this is a common note among packages using C++ through the Rcpp library.
