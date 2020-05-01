# bayesian-position-decoding
This repository applies bayesian population decoding of the electrophysiological activity of place cells in order to estimate position. It was based on the approach used by the Frank Lab for bayesian population decoding of sorted place cell data, and some of these functions can be found in the archive folder. 

The main script you will use is the `script-bayesian-population-decoding.m` file. This analysis has four main steps:
  1. Get place cell activity (spike rate across position normalized by occupancy) using all of the data when the animal is running above a specific speed.
  2. Get the training data (spike rate across position bins) using a subsection of the data (the training set, can be defined in multiple ways)
  3. Get the testing data (spike rate across temporal bins) using the other subsection of the data (the testing set, can be defined in multiple ways)
  4. Decode position during theta and plot the estimated vs. actual position and decoding error. 
  
The second script to use is `script_rippleseqdecoding.m` found in the `replaydecoding` subfolder. This script can be used to decode replay sequences during sharp-wave ripple events using the same principles as described above. Note that this script is much less refined and debugged. 

NOTE - this is a first pass analysis and further refinement is likely needed
