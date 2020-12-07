## Test environments
* local OS X install, R 3.6.3
* ubuntu 18.04.5 (on travis-ci), R 3.6.3
* win-builder (devel)

## R CMD check results
WARNING:
    ‘gs+qpdf’ made some significant size reductions:
       compacted ‘GLMcat_tutorial.pdf’ from 387Kb to 95Kb
    consider running tools::compactPDF(gs_quality = "ebook") on these files

NOTE:
    - installed size is 31.7Mb
      sub-directories of 1Mb or more:
      libs  31.1Mb
    - Package in Depends/Imports which should probably only be in LinkingTo: 'BH'
    
There were no ERRORs

## Downstream dependencies
