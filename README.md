# GMGTS: Gradient Matching Global Two-Stage
_MATLAB®  code package for (ODE-based) nonlinear mixed-effects (NLME) model inference using GMGTS as presented in [1]_

&nbsp;

## Quickstart
Install MATLAB® R2021a with IQM Tools Pro Version 1.2.2 (02.01.2017) and Monolix 2021R2.

Run the `inference_*.m` files to reproduce the results in [1].

The folder `FP_data/` contains the FP maturation microscopy data used for inference.

All subclasses and auxiliary files related to the Generator and Estimator classes are contained in `Generator_Estimator/`.

The `Generator` class generates heterogeneous single-cell time series data from a specified NLME model.

The `Estimator` class infers underlying ODE parameter (random effect) distributions from single-cell time series using the GMGTS [1] and Global Two-Stage (GTS) [2] methods.

&nbsp;

## Usage
### Syntax
`estimates = GMGTS(model_file, data, ...)` infers random effect distributions of the system specified in the IQM `model_file` (instructions for setting up the model file are given below) from measurements given in `data`. Here `data` is either a `TxLxN`-dimensional array with measurements at `T` time points for `L` observables and `N` cells (individuals), or a `1x1 struct` with its measurements stored in a field named `y` or `traces`. The `estimates` are returned as a struct containing the inferred random effect mean `b` and covariance matrix `D`, individual estimates `beta`, and predicted states. Additional arguments are passed to the System and Estimator constructors, see the details below. &nbsp;

`estimates = GMGTS(model_file, data, t, ...)` assumes which the measurements were taken at time points `t`. If `data` is a struct, `t` is ignored and assumed to be a field of `data`. &nbsp;

`estimates = GMGTS(model_file, data, t, observed, ...)` specifies the indices of the observables with respect to the system determined by `model_file` through `observed`. If `data` is a struct, `observed` is ignored and assumed to be a field of `data`. &nbsp;

`estimates = GMGTS(model_file, data, t, observed, init, ...)` integrates the ODE system from the initial values given in `init` to make state predictions. If `data` is a struct, `init` is ignored and assumed to be a field of `data`. &nbsp;

`[estimates, estimator] = GMGTS(model_file, data, ...)` also returns the instantiated `Estimator` object.&nbsp;

### Input arguments

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
