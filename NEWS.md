# periodicDNA 0.3.2 (dev)  

* IMPORTANT: 
  - Implemented data-raw for reproducibility

* MINOR: 
  - Changed xlim of norm. distr. plot in plotPeriodicityResults()

# periodicDNA 0.3.1 (2020-05-05)  

* IMPORTANT: 
  - rollmean(k=3) is now applied *before* normalisation *as well*, 
      on the raw distribution vector
  - plotPeriodicityResults() output returns one single plot (with cowplot)
  - getPeriodicityTrack() now returns the Rle
  - Improved plotting functions -now show shuffled for plotPeriodicityResults()
  - Added ggplot2 theming

* MINOR:
    * Changed many variable names (all to snake_case)
    * sampleGRanges is now full-fledged function 
        (GRanges, DNAStringSet, character and BSgenome methods)
    * sampleGenome is an alias for sampleGRanges.character
    * Added sacCer3 to getPeriodicity BSgenomes
    * Added DNAString method for getPeriodicity
    * Added a vignette describing the internal steps
    * Clarified user-level functions in README
    * Added ce11_TSSs data
    * Renamed generateperiodicitytrack as getPeriodicityTrack
    * Renamed variables in getFPI and getPeriodicity
    * Created a utility char2BSgenome()

# periodicDNA 0.3.0 (2020-05-03)  

* Added tests
* Added getFPI function
* cleaned-up functions
* cleaned-up function dependencies
* Added toy data
* Added vignette

# periodicDNA 0.2.1 (2020-03-04)  

* Added Travis build check
* Simplified README.md
