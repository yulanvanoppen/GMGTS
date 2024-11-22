# GMGTS: Gradient Matching Global Two-Stage
_MATLAB®  code package for (ODE-based) nonlinear mixed-effects (NLME) model inference using GMGTS as presented in [1]_

&nbsp;

## Quickstart
Install MATLAB® R2021a with IQM Tools Pro Version 1.2.2 (02.01.2017) and Monolix 2021R2.

Run `example.m` to generate data from an example ODE-based ME model and recover the random effects distribution using GMGTS.

Run the `inference_*.m` files to reproduce the results in [1].

The folder `FP_data/` contains the FP maturation microscopy data used for inference.

All subclasses and auxiliary files related to the Generator and Estimator classes are contained in `Generator_Estimator/`.

The `Generator` class generates heterogeneous single-cell time series data from a specified NLME model.

The `Estimator` class infers underlying ODE parameter (random effect) distributions from single-cell time series using the GMGTS [1] and Global Two-Stage (GTS) [2] methods.

&nbsp;


## Usage
### `GMGTS()` utility function
`GMGTS` attempts to recover the distribution parameters `b` and `D` of a mixed-effects model
```
    dx_i(t)/dt = g(t, x_i(t)) beta_i + h(t, x_i(t))
    beta_i ~ (log)Normal(b, D)                          (i = 1, ..., N)
```
from measurements `(t_j, y_ij)`, where `y_ij` are vectors of observed components of `x_i(t_j)` perturbed by measurement noise.

`GMGTS` builds on the GTS framework, using gradient matching to obtain cell-specific estimates after smoothing the measurements.
 
`out = GMGTS(model_file, data, ...)` infers random effect distributions of the system specified in the IQM `model_file` (instructions for setting up the model file are given below) from measurements given in `data`. Here `data` is either a `TxLxN`-dimensional array with measurements at `T` time points for `L` observables and `N` cells (individuals), or a `1x1 struct` with its measurements stored in a field named `y` or `traces`. The `estimates` are returned as a struct containing the inferred random effect mean `b` and covariance matrix `D`, individual estimates `beta`, and predicted states. Additional arguments are passed to the System and Estimator constructors, see the details below.

`out = GMGTS(model_file, data, t, ...)` assumes which the measurements were taken at time points `t`. If `data` is a struct, `t` is ignored and assumed to be a field of `data`.

`out = GMGTS(model_file, data, t, observed, ...)` specifies the indices of the observables with respect to the system determined by `model_file` through `observed`. If `data` is a struct, `observed` is ignored and assumed to be a field of `data`.

`out = GMGTS(model_file, data, t, observed, init, ...)` integrates the ODE system from the initial values given in `init` to make state predictions. If `data` is a struct, `init` is ignored and assumed to be a field of `data`.

`out = GMGTS(_, 'Plot', false, ...)` disables plots with parameter estimates, the inferred random effects distribution, model predictions, and any smoothed measurements (enabled by default).

`[out, estimator] = GMGTS(model_file, data, ...)` also returns the instantiated `Estimator` object.

See `System` and `Estimator` for a description of additional input arguments.

&nbsp;


### `System` class ODE functions and simulations
The `System` class handles ODE integration, evaluation of (functions derived from) the system's right-hand side, and corresponding partial derivatives as part of the Estimator and Generator classes.

`system = System(model_file)` instantiates a `System` object by processing an IQM tools `model_file`. The specified ODE system should be linear in parameters, otherwise, unexpected results will follow. The model file should contain the model name, the differential equations and initial conditions for each state, and nominal values for the parameters (model reactions and functions are not yet supported). The `model_file` argument should be the path to a `txt` file. Its contents should be structured as follows:
```
********** MODEL NAME
Simple model
********** MODEL STATES
d/dt(A) = k1 - k2*A
d/dt(B) = k2*A
A(0) = 1
B(0) = 0
********** MODEL PARAMETERS
k1 = 0.1
k2 = 0.5
```

`system = System(model_file, 'FixedParameters', names)` treats the parameters specified in `names` as constant, fixing them at their corresponding nominal values prescribed in `model_file`.

`system = System(model_file, 'FixedParameters', names, 'FixedValues', values)` alternatively fixes the parameters in `names` at the given `values`. The arguments `names` and `values` should have equal numbers of elements.

&nbsp;


### Output

### Examples

&nbsp;

## References
[1] van Oppen, Yulan B. and Milias-Argeitis, Andreas (2024). Gradient matching accelerates mixed-effects inference for biochemical networks. _Manuscript submitted for publication_.

[2] Davidian, Marie. (2017) _Nonlinear models for repeated measurement data._ London: Routledge.

&nbsp;

## DISCLAIMER
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
