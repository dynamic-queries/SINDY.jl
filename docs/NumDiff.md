### Numerical Differentiation of Discrete Data.

- Given a sequence of Discrete Data (a time series for instance)
- Our objective is to compute the derivatives obtained from these samples.
- We shall focus initially on first order derivatives and build our way to higher order ones.


#### Algorithm 1 : No noise

- Assuming that you data has no noise the first order derivative is simply given by the Euler approximation.
- <img src="https://render.githubusercontent.com/render/math?math= v_t = \frac{x_{t+1} - x_t}{δ}">
- Implementation of this function is embarissingly trivial but we do it for completeness.
