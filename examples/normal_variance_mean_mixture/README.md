
## [output.csv](data/output.csv)

### Column name prefixes
1. GLM - fixed alpha, gradient descent
2. СLM - fixed alpha, coordinate descent (all coordinates per iteration)
3. NVMME - floating alpha, EM
4. NVMMG - floating alpha, EM and gradient descent
5. NVMMG - floating alpha, EM and coordinate descent (all coordinates per iteration)

### Column name suffixes
1. _iter - count of iterations (max 1000)
2. _time - time of computation
3. _lh - log2 likelihood
4. _alpha - computed alpha

## Testing 
[Testing engine for Atmosphere GM](https://github.com/9il/atmosphere_gm_test) Statistical package.
