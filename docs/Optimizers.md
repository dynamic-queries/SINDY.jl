# Motivation
This script is implemented after we have obtained the temporal derivative $\dot x$  and the library function $\Theta$, to look for a set of sparse regression solvers to determine $\Xi$. It would seem trivial to solve this regression problem. But it turns out a least square regression solution is not expected owing to two arguments.
1) Physical laws are rarely dense in polynomials
2) Least Squares does not generalize to datasets with large number of parameters as is the case when solving PDEs.

While one can write down Lagrangians in less than ideal coordinates that are so dense they can fill up a page, it is an outlier. We are interested in the general class of problems that would be governed by classical mechanics.

# STLSQ Algorithm 

We use the cookie-cutter recipe that the author’s provide in their paper - Sequentially threshold least squares (STLSQ). The algorithm is as follows.
```
1) Define a threshold - lambda.
2) Perform least square optimization on theta and v.
3) Truncate intermediate result using the threshold lambda.
4) Repeat 2.
```

# Remarks
This optimizer fails in those scenarios when the condition number of the basis matrix is very large, so one needs algorithms that use relaxation methods that realize this like the [SR3](https://ieeexplore.ieee.org/document/8573778/) algorithm. Our optimizer is not optimized, as a result it might be more efficient to use a library implementation for the future releases of our implementation.
