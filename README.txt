PyALTIWAVES - Python-based ALong-Track Inventory of WAVelet based EddieS
------------------------------------------------------------------------

	--> Automatic census of eddies from along-track altimetry data.

	
This library has been base on the work of Dussurget et al., 2011 in the Bay of Biscay and updated with the work accomplished in the Golf of Mexico by Le Hénaff et al, 2012.
Wavelet analysis is applied along-track to detect gaussian-shaped signals such as Oceanic Eddies. Once detected, their properties can thus be derived over a given area or period. Properties include several definitions of the spatial length-scales (eg. decorrelation scale and size of the core of a solid-body rotating eddy), amplitude, rotation sense, relative vorticity, etc.


Useful references:
++++++++++++++++++
  * Dussurget, R., F. Birol, et R.A. Morrow. in prep. « Constructing fine-scale multi-mission altimeter maps for regional and coastal applications »

  * Le Henaff, M., V. H. Kourafalou, R. Dussurget, R. Lumpkin. 2013. « Cyclonic activity in the eastern Gulf of Mexico: characterization from along-track altimetry and in situ drifter trajectories ». Progress in Oceanography, Available online 14 September 2013, ISSN 0079-6611, http://dx.doi.org/10.1016/j.pocean.2013.08.002. (http://www.sciencedirect.com/science/article/pii/S0079661113001626)

  * Dussurget, R, F Birol, R.A. Morrow, et P. De Mey. 2011. « Fine Resolution Altimetry Data for a Regional Application in the Bay of Biscay ». Marine Geodesy 2 (34): 1‑30. doi:10.1080/01490419.2011.584835.

  * Lilly, J.M., P.B. Rhines, F. Schott, K. Lavender, J. Lazier, U. Send, et E. D’Asaro. 2003. « Observations of the Labrador Sea eddy field ». Progress In Oceanography 59 (1) (octobre): 75‑176. doi:10.1016/j.pocean.2003.08.013.

  * Torrence, C., et G.P. Compo. 1998. « A Practical Guide to Wavelet Analysis ». Bulletin of the American Meteorological Society 79 (1) (janvier 1): 61‑78.


Installation :
++++++++++++++
Refer to the `Download <https://code.google.com/p/py-altiwaves/wiki/Download>`_ page.

Special notes:
++++++++++++++

Author :
========
Renaud DUSSURGET, renaud.dussurget

Licensing & copyright:
======================
`GNU Lesser General Public License <http://www.gnu.org/licenses/>`_. Copyright 2012 Renaud Dussurget::
   
   PyAltiWAVES is free software: you can redistribute it and/or modify it under
   the terms of the GNU Lesser General Public License as published by the Free
   Software Foundation, either version 3 of the License, or (at your option)
   any later version.
   PyAltiWAVES is distributed in the hope that it will be useful, but WITHOUT
   ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
   FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License
   for more details.
   
   You should have received a copy of the GNU Lesser General Public License along
   with PyAltiWAVES.  If not, see <http://www.gnu.org/licenses/>.



History :
=========
   * Adapted from CTOH/LEGOS IDL ALTI_WAVES library (RD).
   * November 2012 : Converted to Python
      - Added computation of eddy core scale and relative vorticity
      - Improved 2D eddy detection
      - Modified decorrelation scales (asses the scale on both sides of eddy) : still require some work when computing the scales.
      - Added computing temporal variability of parameters.
      - Added a test case using simulated red noise data.
      - Tested on DUACS MERSEA NRT & DT products over Med : spatial and temporal variability of eddy properties.
      - Validation over GoM : good agreement over long time series, improved consistency. 
   * February 2013 :
      - Adapted to 20Hz altimetry data analysis
      - Added spectral analysis functions (compute mean spectrum and perdiodogram).
      - Added rankine eddy fitting to data : better estimation of relative vorticity.
