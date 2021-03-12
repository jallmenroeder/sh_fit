# BRDF fitting for Linearly Transformed Spherical Harmonics (LTSH)

## THIS IS A WORK IN PROGRESS REPOSITORY

Implementation of the BRDF approximation process for LTSH. This Code is reimplementation of the code used in my [Bachelor's Thesis](http://www.jallmenroeder.de/2020/11/19/linearly-transformed-spherical-harmonics/) with improved sampling. The original python implementation used code not intended for publication, so I decided to do an improved reimplementation in C++. Currently, the implementation has issues that make the resulting fits less desirable than the original fits from the python implementation. The LTC fitting implementation contains a lot of NaNs and the LTSH quality is often worse than that from the python fitting. However the code shows how we generally approached the problem and might help people doing further research with LTSH. I hope to find the time to fix the implementation and make it a true improvement of the python version. 
## Dependencies
The code is built using CMake and C++14. It depends on the [GLM library](https://github.com/g-truc/glm) and [GSL](https://www.gnu.org/software/gsl/doc/html/index.html) for the minimization algorithms and SH evaluation.
## How to use
TODO
