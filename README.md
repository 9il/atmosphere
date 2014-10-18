atmosphere
=============
Statistical package.

Testing in progress. NOT READY FOR PRODUCTION!

# Numeric methods

## General nonparametric algorithms

Each algorithm solvs optimization problem f(p) -> min, 
where p is discrete probability distribution with k elements.

```
f(p) -> min
f = u(Wp),
W - matrix(n rows, k columns),
u - convex function.
```

## Separates normal variance mean mixtures.
```
///normal variance mean mixtures
f_sample(p, alpha) -> max
f_sample = u(W_sample(alpha)p),
W_sample(alpha) - matrix(n rows, k columns) of posterior probabilities,
u(ω) =  Σ_j log(ω_j).
```

#Instalation
Use [dub package manager](https://github.com/D-Programming-Language/dub) for instalation.
## Dependencies
1. [D Simple Matrix](https://github.com/9il/simple_matrix)
2. [D cblas header](https://github.com/9il/cblas)
3. BLAS library

##C interface
C header file cab be found [here](https://github.com/9il/atmosphere_gm/tree/master/include).
