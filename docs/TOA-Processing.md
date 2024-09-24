# TOA Processing 

## Read IQ data from ring buffer
https://dspsr.sourceforge.net/manuals/dspsr/shm.shtml

## Dedisperse pulse proile
Electromagnetic waves propagating through the inter- stellar medium (ISM) experience phase dispersion that effects a frequency-dependent group velocity.
```bash
dspsr -P accel.txt pointing.dat
```

## Fold Pulse Profile
A new survey candidate has very few parameters: the spin period, trial dispersion measure, and perhaps an acceleration. These can be input to dspsr using a simple text file; for example: 
```
SOURCE: SC0001_04151
PERIOD: 1.362819517 s
DM: 50.6
ACC: 42.70 (m/s/s)
RA: 09:10:18
DEC: -72:07:35
```

 This file is input to dspsr as a phase predictor using the -P command-line option; e.g. if the above example were saved as accel.txt and the raw data were in pointing.dat: 

## Matched Filter Search for Pulses
To perform matched filtering, save the template pulse profile in a file (e.g. template.std) and correlate each pulse with the template as in the following example: 
```
# disable log messages
init verbose --

# load the template
init load template template.std

# form a total intensity copy of the raw data
push
fscrunch
pscrunch

# cross correlate the profile with the template
profile correlate template

# keep only profiles with a peak greater than 15 sigma
test $all:max > ( $off:avg + 15*$off:rms )

# restore original file
pop
```

## Estimate Arrival Times 
Can be performed using two different methods of arrival time estimation: scalar template matching using only the total intensity and matrix template matching using all four Stokes parameters
```bash
# Scalar template matching mode 
cd $PSRCHIVE_DATA/mem/pulsar
pat -F -s ../std/standard.FF *.zz > uncal_unsmooth_stm.tim

# Matrix template matching mode 
pat -Fpc -s ../std/standard.FF.sm *.calib > cal_smooth_mtm.tim

```

```bash
# Run tempo2 to evaluate the arrival times; e.g.
tempo2 -f ../pulsar.par uncal_unsmooth_stm.tim
```

In python, this calculation would be:
```python
import psrchive

arrtim = psrchive.ArrivalTime()
arrtim.set_shift_estimator('PGS')        # Set algorithm (see 'pat -A' help)
arrtim.set_format('Tempo2')              # set TOA format
arrtim.set_format_flags('IPTA')          # set some TOA flags
arrtim.set_attributes(['chan','subint']) # More TOA flags

# Load template profile
std = psrchive.Archive_load('J1713+0747.Rcvr1_2.GUPPI.9y.x.sum.sm')
std.pscrunch()
arrtim.set_standard(std)

# Load observation profiles
obs = psrchive.Archive_load('guppi_55616_J1713+0747_0009.12y.x.ff')
obs.pscrunch()
arrtim.set_observation(obs)

# Result is a tuple of TOA strings:
toas = arrtim.get_toas()
```

## Fit Timing Residuals 
TEMPO2 typically assumes that the position of the observatory is well known, and processes TOA measurements to refine the pulse profile model for a known pulsar. However, if we assume that the pulsar timing model is well-formed we can use the incoming TOA to refine the barycentric coordinates of the observer. If the observer is stationary this can be worked out directly from the timing-residial, but if the observer is dynamically positioned, the processing must be integrated with a local velocity model and kalman filtering algorithm. Additionally, due to the limitations imposed by the pulse phase estimation problem, we assume that the positional error is within the half-phase of the pulsar. 

Process multiple pulsars using TEMPO2
```bash
tempo2 -f psr1.par psr1.tim -f psr2.par psr2.tim -f psr3.par psr3.tim 
```
This will output clock correctional files and also the fitted pulsar model, including the pre- and post-fit RMS residuals

`-f x.par y.tim` specifies the parameter (.par) and arrival time (.tim) files to use for subsequent
processing. If only a .tim file is present (without -f option) then the parameter file will be assumed
to be y.par.

`-fit x`, turn on fitting for parameter ’x’ (this command-line option can be repeated for multiple
parameters)

`-residuals` outputs the residuals to a file called “residuals.dat”.

## TEMPO2 Plugin 
"""
sudo apt update
sudo apt install libgsl-dev
sudo apt install proj-bin
sudo apt install libgdal-dev 
"""
To build, navigate to the plugins/ directory and run the following
```
sudo g++ -I/usr/local/tempo2/include -fPIC -shared -o ${TEMPO2}/plugins/pulsar_positioning_Linux_plug.t2 plugin/pulsar_positioning_plug.C -ltiff -lgeotiff -lgdal -lgsl
```

To execute the positioning script, run the following
```
tempo2 -gr pulsar_positioning
```

```
sudo valgrind --leak-check=full          --show-leak-kinds=all          --track-origins=yes          --verbose          --log-file=valgrind-out.txt               ./../../../../usr/local/
tempo2/bin/tempo2 -gr pulsar_positioning -f rpnav/simulate/inputs/navigate.par rpnav/simulate/inputs/navigate.tim -observer NAVIGATE
```