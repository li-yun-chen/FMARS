# FMARS
This program uses MARS to improve Frechet regression. Y is an element of a metric space, and X is an element of a Euclidean space.

In Demo1, Y is a density function. In Demo2, Y is SPD.

FMARS.R is the main program, and to speed up the computation, the program uses parallel computing. The forward process of MARS in this file is written in C++. The corresponding R version code is in MARS.R.
