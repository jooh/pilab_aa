This code provides fast, reproducible pattern information analysis of functional
imaging by combining analysis code from [pilab](https://github.com/jooh/pilab)
with distributed batch processing and provenance functionality from [automatic
analysis](https://github.com/rhodricusack/automaticanalysis) 5 (AA).

# dependencies
* MATLAB R2013A (newer versions will *generally* work but you may bump into
  problems during the visualisation modules since we don't yet support the new
  handle graphics system introduced in R2014B).
* [automatic analysis](https://github.com/rhodricusack/automaticanalysis) (see
  also [my fork](https://github.com/jooh/automaticanalysis) for the
  bleeding edge version)
* [pilab](https://github.com/jooh/pilab)

# installation
The easiest way to use these modules in AA is to use addpath to add this
directory to your MATLABPATH in your startup.m file. AA completely resets your
path during cluster-based execution so simply adding the code to your current
session may not work.

# examples
Example usage can be found in the analysis code for my [face-space fMRI
project](https://github.com/jooh/facedistid_analysis.git).
