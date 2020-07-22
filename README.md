# Matlab Codes to reproduce some of the results of [Viens and Denolle (2019, BSSA)](https://pubs.geoscienceworld.org/ssa/bssa/article/571631/long-period-ground-motions-from-past-and-virtual)

## Description:
* In [Viens and Denolle (2019, BSSA)](https://pubs.geoscienceworld.org/ssa/bssa/article/571631/long-period-ground-motions-from-past-and-virtual), we use the ambient seismic field recorded by offshore and onshore seismic stations in Japan to simulate the ground motions of subduction earthquakes. This repository contains the Matlab functions and codes to reproduce some results of the paper (e.g., waveforms for the NAGH and KNHH stations in Figure 8 and Figures 4a and 4c.).

* The **"Codes"** folder contains three Matlab codes: 
  - **RS_spectra.m**: function to compute the response spectra.
  - **func_readsac.m**: function to read sac files. 
  - **Large_earthquake_simulation.m**: Main code and simulates the 2004 Mw 7.2 Off the Kii peninsula earthquake using a simple elliptic-like source model. The code produces three plots: 
    1) The comparison between the simulated and observed waveforms. 
    2) The amplitude of the moment rate functions over the fault plane as well as the total moment rate function. 
    3) The comparison between the observed and simulated response spectra. 
  
* The **Data** folder contains: 
  - The **GFs** folder: Z-Z IRFs, which have been calibrated to absolute amplitude using a Mw 5.5 event.  
  - The **Earthquake** folder: Recorded vertical waveforms from the Mw 7.2 earthquake. 
  - The **Meta_data.mat** file contains the latitude and longitude of the stations. 
  - All the waveforms have a 10 Hz sampling rate and have been bandpass filtered between 4 and 10 s using a four-pole two-pass Butterworth filter. 
  
## Additional info:
* The Python codes used to compute impulse response functions with the deconvolution method (see [Viens et al. (2017)](https://academic.oup.com/gji/article/210/1/210/3747441) for details) are available at: https://github.com/lviens/2017_GJI
