# dfmux_calc

### Introduction

The goal of this code is to estimate the readout noise of a full DfMUX system, including a SQUID, wiring harness, parasitics. The bolometer's affect on the readout circuit is included, as well as its Johnson noise. However, phonon and photon noise are not treated in `dfmux_calc`, as they are not considered "readout noise."

### Usage

Please see `LiteBIRD_SSAA_parameters.ipynb` for an example, as well as a brief explanation of how to set up SPICE, which is an optional dependency.

### Revision History

 - Aug 2022: first version committed to https://github.com/megan-russell/dfmux_calc
 - Mar 2024: simplified version forked to https://github.com/tijmen/dfmux_calc

### Citations
This code is heavily based on the following (which are referred to in comments when relevant)
Joshua Montgomery's PhD Thesis - https://escholarship.mcgill.ca/concern/theses/1c18dm29r
2021 paper SPT-3G assessment of DfMUX noise https://arxiv.org/abs/2103.16017
Joshua Montgomery's Masters Thesis - https://www.proquest.com/docview/2514727898?pq-origsite=gscholar&fromopenview=true