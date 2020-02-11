# bayesian-position-decoding
This repository applies bayesian population decoding of the electrophysiological activity of place cells in order to estimate position. It was based on the approach used by the Frank Lab for bayesian population decoding of sorted place cell data, and some of these functions can be found in the archive folder. 

The main script you will use is the `script-bayesian-population-decoding.m` file. This analysis has four main steps:
  1. Get place cell activity (spike rate across position normalized by occupancy) using all of the data when the animal is running above a specific speed.
  2. Get the training data (spike rate across position bins) using a subsection of the data (the training set, can be defined in multiple ways)
  3. Get the testing data (spike rate across temporal bins) using the other subsection of the data (the testing set, can be defined in multiple ways)
  4. Decode position during theta and plot the estimated vs. actual position and decoding error. There is also code to decode replay sequences during sharp-wave ripple events that can be added in the future.
