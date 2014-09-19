HomeoNorm
=========

normalize a set of homeologous read counts from different libraries

Simon Renny-Byfield, Iowa State University September 2014

A function to normalize libraries with homeologous counts.
The script assumes that each individual library has been 
through polyCat and been split into At,Dt and N counts.
Each sequence of three columns will be treated as a single library
and the normalization number for each will be calculated based 
on the total counts (At, Dt and N) from the library (now in three seperate 
columns). This approach ensures the ratios betwee At, Dt and N are maintained
within a sample, even after normalization between samples has
occured.

Input is a table of read counts, with At, Dt, and N columns for each library.
