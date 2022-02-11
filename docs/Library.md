# Construction of a library function
This script contains a library that specifies the candidate basis functions to be used to construct ```Theta(x)```. The library consists of two variants
 1) Polynomial Basis up to the 3<sup> rd</sup> order.
 2) sine, cosine and tangent Trigonometric Basis.

# Motivation

The completeness of these bases is justified for our set of problems, as most real dynamical systems’ vector
field can be approximated efficiently either in the polynomial or the Fourier basis. 

# Module
In ```src/Library.jl```, we define three library types:
1) Polynomial 
2) Trigonometric
3) Polytrigonometric  

The polynomial basis function ```basis(X :: Array,o :: PolynomialBasis)``` returns a matrix containing all the: 
1) Linear (x, y, z)
2) Quadratic (x<sup>2</sup>, xy, y<sup>2</sup>, yz, z<sup>2</sup>, xz) and
3) cubic (x<sup>3</sup>, x<sup>2</sup>y, x<sup>2</sup>z, xy<sup>2</sup>, xz<sup>2</sup>, xyz,y<sup>3</sup>, y<sup>2</sup>z, yz<sup>2</sup>, z<sup>3</sup>) terms of the coordinates (x,y) or (x,y,z).

The trigonometric basis function ```basis(X :: Array,o :: TrigBasis)``` contains the trigonometric functions
1) ```cos(x)```
2) ```sin(x)```
3) ```tan(x)``` applied to the coordinates
On the other hand, the poly-trigonometric basis function ```basis(X :: Array,o :: PolyTrigBasis)``` is a superposition of both the polynomial and the trigonometric bases.
# Remarks
one is encouraged to try a set of functions with a lower order first before one uses the entire basis set that is at your disposal. Unfortunately, this part of the implementation is still manual, in that the user is asked to choose the basis.