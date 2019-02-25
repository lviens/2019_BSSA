# Codes to reproduce some of the results of Viens and Denolle (Submitted to BSSA)

* The Matlab codes in the "Codes" folder reproduce the waveforms for the NAGH and KNHH stations presented in Figure 8 as well as Figures 4a and 4c.

* The **Codes** folder contains three Matlab codes: 
  - The **RS_spectra.m** code is a function to compute response spectra.
  - The **func_readsac.m** code is a function to read sac files. 
  - The **Large_earthquake_simulation.m** is the main code and simulates the 2004 Mw 7.2 Off the Kii peninsula earthquake using a simple elliptic-like source model. The code produces three plots: 
    1) The comparison between the simulated and observed waveforms. 
    2) The amplitude of the moment rate functions over the fault plane as well as the total moment rate function. 
    3) The comparison between the observed and simulated response spectra. 
  
* The **Data** folder contains: 
  - The **GFs** folder: Z-Z IRFs, which have been calibrated to absolute amplitude using a Mw 5.5 event.  
  - The **Earthquake** folder: Recorded vertical waveforms from the Mw 7.2 earthquake. 
  - The **Meta_data.mat** file contains the latitude and longitude of the stations. 
  - All the waveforms have a 10 Hz sampling rate and have been bandpass filtered between 4 and 10 s using a four-pole and two-pass Butterworth filter. 
  
* The Python codes to compute impulse response functions (IRFs) are available at: https://github.com/lviens/2017_GJI
