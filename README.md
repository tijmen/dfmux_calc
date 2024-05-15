# dfmux_calc

Calculations related to the Simons Array / SPT-3G / LiteBIRD DfMux readout system for TES bolometers.

 - `dfmux_calc.py` is the main entry point. The goal of this code is to estimate the readout noise of a full DfMUX system, including a SQUID, wiring harness, and some parasitics. The bolometer's effect on the readout circuit is included, as well as its Johnson noise. However, phonon and photon noise are not treated in `dfmux_calc`, as they are not considered "readout noise."
 - `current_sharing.py` is needed by dfmux_calc if SPICE is enabled. 
 - `LiteBIRD_SSAA_parameters.ipynb` shows an example of how to use `dfmux_calc` for calculation of the LiteBIRD SQUID electrical performance requirements. It also has a brief explanation of how to install SPICE, which is an optional dependency.
 - `mntes.py` calculates the properties of a TES in the DfMux readout system including nonlinearity. 

### Revision History

 - Aug 2022: first version committed to https://github.com/megan-russell/dfmux_calc
 - Mar 2024: forked and simplified, hosted at https://github.com/tijmen/dfmux_calc

### Citations

 - Joshua Montgomery's PhD Thesis - https://escholarship.mcgill.ca/concern/theses/1c18dm29r
 - 2021 SPT-3G paper https://arxiv.org/abs/2103.16017
 - Joshua Montgomery's Masters Thesis - https://www.proquest.com/docview/2514727898?pq-origsite=gscholar&fromopenview=true