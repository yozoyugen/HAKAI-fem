In some cases, parameter adjustment may be required.


[crash-tube-80-350-solid.inp]
Matlab:
HAKAI.m around line 1506
->  k = young * S / Lmax * 10.0;
(default: k = young * S / Lmax * 1.0;)

Julia:
HAKAI_j.jl around line 3263
->  k = young * S / Lmax * 10.0
(default: k = young * S / Lmax * 1.0)
