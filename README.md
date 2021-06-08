# SpawnSeis hydrophone measurements
This is the desciption of sound alalysis for an experiment where we wanted to test if cod changed behaviour or escaped the area when exposed to sound from a double airgun. This analysis focus on analyzing the recorded sound level in a small bay where cod was monitored by telemetry while it was exposed to different sound treatments for 3 hour intervalls for 5 days. The treatmets were:
"Seismic" with airgun shots every 10 seconds from a vessel driving back and forth in front of the cod bay.
"Boat control" wiht the vessel without airgun driving bac and forth in fron to of the cod bay.
"Silent control" with the vessel leaving the area.
The treatments were performed in random order as described in the document treatments.csv.

The same experiment was performed twice, first in 2020 and then in 2021. This is the analysis for both years.

## Hydrophones
 The document SpawnSeisHydrophoneMetaData.csv describes serialnumber, frequenc range, and sensitiity for the hydrophones used in the experiment. For this analysis we use only the hydrophones which was placed inside the cod bay denoted H1 H2 H5 H6 and H8. The hyrophones was of type Naxys Ethernet Hydrophone model 02345, and they were using battery which had to be loaded during the experiment. For this reason there are not sound recordings for all treatments, but for several examples of each type.
 
Each period of recording for each hydrohone is called a deployment. Parameters like start time, stop time and position for each deployment is described in the document SpawnSeisHydrophoneDeploymentMetaData.csv. Two of the hydrophones was placed close to the bottom, one at the entry of the bay, and one at the inner part of the bay. These are refered to as bottom moored hydrophones. Two hydrohones was also placed in the center of the bay at two different depht, both hanging below each other under a bouy. These are refered to as vertical hyrophone aray. 

The hyrdophones logged data for 22 seconds with 8 seconds brake between each recording. This means that all pulses is not recorded, b

# Calibration
The sensitivity given for the Naxys Ethernet Hydrophones is taken from the datasheet. This gives an approximate value. The hydrophones was calibrated before and after the experiment by a B&K Type 4229 calibrator and a coupler designed for the Naxys hydrophones. Rms value of the signal from the calibrator with the coupler is 148.5 dB re uPa^2 or 26.6 Pa, refered to as "knownPa". The path to the recordings of the calibration signal is given in the document SpawnSeisHydrophoneDeploymentMetaData.csv. For some cases the calibration was only made before the experiment.
The calibration was done in the linear domain by multiplying the raw data with a calibration factor, cal, estimated as cal=knownPa/rms(calibration-data).


## Filtering
The Naxys hydrophones had a frequency range from 5 Hz to 300 kHz according to the datasheet. For the hyrdrophone in the inner part of the bay the signal to noise ratio for the airgun signals was very poor. This was mainly because of low frequency ambient noise with frequencies below 5 kHz. Even if this was outside the defined frequency range of the hyrdrophone it was still detected. Since the calibration for the hydrophone is not valid under 5 Hz, and since we get noise in this range a filter was applied. A bandpass filter from 5-10000 Hz was used. The upper frequency of the band was set at a rather high value so it would not affect the steep rize time and the peak of the airgun pulses. The higher frequencies was regarded as unwanted since they can not be heard by cod, and they do not originate from any of our sound treatments. 

The filter used was a 6th-order Butterworth bandpass filter from 5-10000 Hz.

## Metrics
For the airgun pulses, each peak was detected by a function "findpeaks" in matlab. The minimum spacing between the detected peaks was set to 8 seconds since there was 10 seconds between each airgun pulse.
Then a range around the detected peaks (3 seconds before, and two seconds after) was sent to analyses. The peaks that did not have this range around the them, if the peak was close to the end of the file, was disregarded. The rest of the peaks was analyzed further.

The same method was applied to the Boat contro and the silent control where the findpeaks function found the highest peak in the ambient noise.

For each of the selected 5 second interwall, 1 second around the highest peak (0.3 sec before and 0.7 sec after peak) was used for further analysis. The ambient noise for a one second period selected 2 seconds before the peak was also analyzed. 
The following metrics was estimated for the 1 second sequences for peak and for noise:

pospeakpressure: maximum peak found in the 1 second long sequence for peak. max(1sec_data)
negpeakpressure : minimum peak found in the 1 second long sequence for peak. min(1sec_data)
pospeakpressureN:maximum peak found in the 1 second long sequence for noise. max(1sec_data)
negpeakpressureN: minumum peak found in the 1 second long sequence for noise. min(1sec_data)
Ex: The exposure level which is a measure of the energy in the signal in linear domain. Ex=dt*sum(1sec_data.^2)
ExN : The exposure level which is a measure of the energy in the noise in linear domain. Ex=dt*sum(1sec_data.^2)
SEL : The sound exposure level which is a mesure of the energy in the singal in dB: 10*log10(Ex/1e-12);
SELN :The sound exposure level which is a mesure of the energy in the noise in dB: 10*log10(Ex/1e-12);

For selected pulses it is also possible to do a frequency analyses:
The 1 second signal sequence was then tapered by a tukeywindow to set the start and stop values of the signal seqence to zero.
Then a fast fourier transform was used to estimate the frequency spectrum of the signal.
The energy spectral density was estimated as 

ESD=((abs(P1)).^2)*df) (Energy spectral density, also called "frequency integrated sound exposure spectral density (ISO)

where P1 fast fourier transform of the pressure,





