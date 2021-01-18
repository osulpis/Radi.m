# RADI v20

1-D Reaction-Advection-Diffusion-Irrigation (RADI) Diagenetic Sediment Module.

Source code by O. Sulpis and M. Wilhelmus. Optimised for MATLAB by [M.P. Humphreys](https://mvdh.xyz).

Requirements for RADI: `calc_pco2()` (in this repo) and [CO2SYS](https://github.com/jamesorr/CO2SYS-MATLAB/blob/master/src/CO2SYS.m).

Requirements for IC_W29: `gsw_rho` from the [Gibbs-SeaWater Oceanographic Toolbox](http://www.teos-10.org/software.htm).

## Instructions

```matlab
IC_W2 % set environmental conditions
rerun = 0; % for a fresh start; 1 for a re-run
RADI % run the model
```

```matlab
delta_phi = -phiBeta.*(phi0 - phiInf).*exp(-phiBeta*depths);
% delta_phi(1) = 0; % don't do this
delta_phiS = -delta_phi;
delta_tort2 = -2*delta_phi./phi; % not used in Julia
delta_D_bio = -2*depths.*D_bio/lambda_b^2; % not used in Julia
```
