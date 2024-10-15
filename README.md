# GMGTS: Gradient Matching Global Two Stage
_MATLAB®  code package (ODE-based) nonlinear mixed-effects (NLME) model inference using GMGTS presented in [1]_

&nbsp;

## Quickstart
Install MATLAB® R2022a with IQM Tools Pro Version 1.2.2 (02.01.2017) and Monolix 2021R2.

Run the `inference_*.m` files to reproduce the results in [1].

The folder `FP_data/` contains the FP maturation microscopy data used for inference.

All subclasses and auxiliary files related to the Generator and Estimator classes are contained in `Generator_Estimator/`.

The `Generator` class is generates heterogeneous single-cell time series data from a specified NLME model.

The `Estimator` class infers underlying ODE parameter (random effect) distributions from single-cell time series using the GMGTS [1] and Global Two Stage (GTS) [2] methods.

&nbsp;

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
