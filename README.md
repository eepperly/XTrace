# XTrace

This repository contains code for the paper _XTrace: Making the most of every sample in stochastic trace estimation_ by [Ethan N. Epperly](https://www.ethanepperly.com), [Joel A. Tropp](https://tropp.caltech.edu), and [Robert J. Webber](https://rwebber.people.caltech.edu).

This paper developes new algorithms for estimating the _trace_

$$
\text{trace}(\boldsymbol{A}) = \sum_{i=1}^n a_{ii}
$$

of an $n\times n$ matrix $\boldsymbol{A}$. We assume we only have access to the matrix $A$ through the matrix–vector product operation $\boldsymbol{x} \mapsto \boldsymbol{Ax}$. We have two algorithms for trace estimation: XTrace, appropriate for an arbitrary square matrix, and XNysTrace, appropriate for real symmetric (or complex Hermitian) positive definite matrices.

As an example of what one can do with such an algorithm, we can use XTrace to estimate the number of triangles in a large network. The number of triangles in a network with adjacency matrix $\boldsymbol{M}$ is $\text{trace}(\boldsymbol{M}^3/6)$. We can use XTrace to estimate this quantity by applying XTrace to $\boldsymbol{A} = \boldsymbol{M}^3/6$. The matrix $\boldsymbol{A}$ is _implicitly defined_ in terms of the matrix $\boldsymbol{M}$, and directly forming the matrix $\boldsymbol{A}$ to sum its diagonal entries is prohibitively expense for large, sparse networks. However, we can compute matrix–vector products with $\boldsymbol{A} = \boldsymbol{M}^3/6$ by repeated multiplications with $\boldsymbol{M}$, i.e.,

$$
\boldsymbol{Ax} = \frac{1}{6} \boldsymbol{M}(\boldsymbol{M}(\boldsymbol{Mx})).
$$

This allows us to estimate the trace of $\boldsymbol{A}$, equivalently the number of triangles in the network, without ever forming $\boldsymbol{A}$ explicitly:

```
>> load('as-Skitter.mat'); M = Problem.A; % See https://sparse.tamu.edu/SNAP/as-Skitter
>> fprintf('Approximately %e triangles in this network\n', xtrace(@(X) M*(M*(M*X))/6, 100, size(A,1)))
Approximately 2.911425e+07 triangles in this network
```

We also have an algorithm, XDiag, which estimates the _diagonal_ of a matrix $\boldsymbol{A}$ using matrix–vector products with both $\boldsymbol{A}$ and $\boldsymbol{A}^\top$.

## Summary of code

This repository contains the following folders:
- Our best implementations of the XTrace, XNysTrace, and XDiag stochastic trace and diagonal estimators are in the `code/` folder.
- Implementations of other trace and diagonal estimators, retooled to use a similar interface to our estimators, are in the `existing_estimators/` folder.
- The `tests/` folder contains code to reproduce the figures in the paper.
- The `figs/` folder holds the output of figures created by the scripts in the `tests/` folder.

## Interface

All trace estimators use the interface

```
[trace_est, err_est] = estimator(A, m[, n, test_vector_type, adjoint_function])
```

The outputs are an estimate `trace_est` of `trace(A)` and an estimate `err_est` of the error `abs(trace_est - trace(A))`. The error estimator is only computed for XTrace and XNysTrace and is outputted as `NaN` for other methods; our diagonal estimators simply omit this output entirely. The input matrix `A` can is an `n` by `n` numeric array or a function handle which computes the action of the matrix `@(X) A*X` on a rectanglar matrix `X` of size `n` by `k`. The implementations of XTrace (`xtrace`), XNysTrace (`xnystrace`), and XDiag (`xdiag`) are designed to work with either real or complex inputs. The input `m` allocates the total number of matrix–vector products to be used for trace estimation. The following optional arguments can be specified in any order:

- `n`: the dimension of the square input matrix `A`. This argument _must_ be specified if the matrix `A` is inputted as a function handle. If `A` is a numeric array, this argument is ignored if provided.
- `test_vector_type`: character array or string for type of test vector to be used; see below.
- `adjoint_function`: XDiag requires multiplications with both `A` and its (conjugate) transpose `A'`. If `A` is specified by a function handle, then this optional argument allows the user to specify a function handle which computes the action of the adjoint of the matrix `@(X) A'*X`. If this argument is not specified, it is assumed that the input matrix is Hermitian and the provided function handle implementing `@(X) A*X` is used. **XDiag will produce incorrect results for non-Hermitian matrices `A` specified by function handle unless this argument is set**; if `A` is Hermitian or specified as a numeric array, this argument can be ignored.

We also provide implementations `xtrace_tol` and `xnystrace_tol` of XTrace and XNysTrace which run until the error estimate `err_est` satisfies `err_est <= abs(trace_est) * reltol + abstol`; `reltol` and `abstol` specify relative and absolute tolerances. For these adaptive implementations, the interface is

```
[trace_est, err_est, m] = estimator_tol(A, abstol, reltol[, n, test_vector_type])
```

These scripts output the total number of matrix–vector products used as the final output argument `m`. We recommend that one set the tolerances conservatively with respect to the desired accuracy, as the error estimate for XTrace and XNysTrace tends to underestimate the true error and is subject to random fluctuations.

## Running the tests

This repository contains code to reproduce figures from the paper. Scripts to do this are as follows:

- Figures 1 and 3: `tests/basic_tests.m`
- Figures 2 and 5a: `tests/tfim_partition.m`
- Figure 4: `tests/test_vectors.m`
- Figure 5b: `tests/tfim_energy.m`
- Figure 6: `tests/networks.m`
- Figure SM1: `tests/compare_adaptive_hutchpp.m`

In order to run these some of these scripts, you will need to download:

1. The [following code](https://www.mathworks.com/matlabcentral/fileexchange/29576-matrix-exponential-times-a-vector) and store it in a folder `expcode` at root level,
2. and the matrix [`yeast.mat` from the SuiteSparse matrix collection](https://sparse.tamu.edu/Pajek/yeast) and store it in `tests/`.

## Test vectors

We recommend that users use the default test vectors for best performance for XTrace, XNysTrace, and XDiag. For testing purposes and for other edge cases, we provide the following options for test vectors:

- **Improved `'improved'`(recommended):** This implements the normalization approach described in section 2.3 of the paper. This can lead to substantial reductions in the error for matrices with flat portions in their spectrum, leading us to recommend this approach as the default for XTrace and XNYsTrace. For algorithms other than XTrace and XNysTrace, this defaults to the `'sphere'` option, described below.
- **Random signs `'signs'`:** Test vectors are composed of independent uniform random signs $\pm 1$. This is a sensible default for most problems, though `'sphere'` is slightly better in the worst case. This is the only supported option for `XDiag` currently.
- **Uniform on the sphere `'sphere'`:** Test vectors are uniformly random vectors drawn from the sphere of radius `sqrt(n)`. This is the worst-case optimal distribution of test vectors for the Girard–Hutchinson trace estimator and is a natural choice.
- **Gaussian `'gaussian'`:** Test vectors are composed of independent standard normal entries. This is here for testing purposes, but should not be used in practice; use `'signs'` or `'sphere'` instead.
- **Orthogonal `'orth'`:** Test vectors are chosen to be the orthogonalized, $\sqrt{n}$-length normalized columns of a Gaussian random matrix. For XTrace and XNysTrace, this performs terribly in practice and should not be used.

Prepending a letter `c` in front of the test vector name gives complex-valued versions (e.g., `cgaussian` uses [standard complex Gaussian random variables](https://en.wikipedia.org/wiki/Complex_normal_distribution) in place of real Gaussians and `csigns` uses uniform values from the complex uniform circle in place of uniform $\pm1$'s).
