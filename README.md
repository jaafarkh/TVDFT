# TVDFT
This repository is dedicated to provide c++ implementation of Time Variant Discrete Fourier Transform with Non-Orthogonality Compensation for Order Tracking Analysis.
The software generates a simulated signal with close and cross orders and analyze it twice
with non-orthog. compensation enabled and with non-orthog. compensation disabled.
You can grab the data and plot it in Matlab or Excel.
The project runs on VS 2017 or later.
Here is the results.
when OCM is disabled, data is very correlated:
https://github.com/jaafarkh/TVDFT/blob/master/TVDFT/ocm-disabled.png

when OCM is enabled, data looks fine:
https://github.com/jaafarkh/TVDFT/blob/master/TVDFT/ocm-enabled.png

please cite as:
Alsalaet, J., Najim, S., and Ali, A. (August 27, 2014). "Order Tracking Analysis Using Generalized Fourier Transform With Nonorthogonal Basis.
" ASME. J. Vib. Acoust. December 2014; 136(6): 061002. https://doi.org/10.1115/1.4028269

# Project Dependency
none

# License
Please read https://github.com/jaafarkh/spectcorr/blob/add-license-1/LICENSE
