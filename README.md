# Codes to reproduce some of the results of Viens and Denolle (Submitted to BSSA)

- The Matlab codes in the "Codes" folder reproduce the waveforms for the NAGH and KNHH stations presented in Figure 8 as well as Figures 4a and 4c. <br/>

- The Codes folder contains three Matlab codes:  <br/>
 - The RS_spectra.m code is a function to compute response spectra <br/>
 - The func_readsac.m code is a function to read sac files. <br/>
 - The Large_earthquake_simulation.m is the main code and simulates the 2004 Mw 7.2 Off the Kii peninsula earthquake using a simple elliptic-like source model. The code produces three plots: <br/>
  1) The comparison between the simulated and observed waveforms. <br/>
  2) The amplitude of the moment rate functions over the fault plane as well as the total moment rate function. <br/>
  3) The comparison between the observed and simulated response spectra. <br/>
  
- The Data folder contains:  <br/>
 - the Z-Z IRFs (GFs folder), which have been calibrated to absolute amplitude using a Mw 5.5 event.  <br/>
 - The Earthquake folder contains the vertical waveforms from the Mw 7.2 earthquake. All the waveforms have a 10 Hz sampling rate and have been bandpass filtered between 4 and 10 s using a four-pole and two-pass Butterworth filter. <br/>
 - The Meta_data.mat file contains the latitude and longitude of the stations. <br/>
  
  
- The Python codes to compute impulse response functions (IRFs) are available at: https://github.com/lviens/2017_GJI <br/>
