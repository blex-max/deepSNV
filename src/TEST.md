The functions in bam2r-pileup.cpp have been exposed for testing, and a test suite has been provided using Catch2 (v3).
The test executable can be compiled via make using `Makefile.test`. A helper bash script, `maketest.sh`, is provided
which simply runs make clean, compiles the test executable, and then runs the test executable. The required dependencies
are htslib and Catch2. By default the makefile will use pkg-config to attempt to locate these on your system, however
these arguments can be overriden at the command line when calling make, like `argname=some/path/`, so you can easily
provide your own include paths. Note that Catch2 is NOT required for installing deepSNV as an R library, only for
running these tests.
