# Dependencies

The RPNAV package has the following dependencies:
    - PSRCAT
    - TEMPO2
    - PSARCHIVE 
    - DSPSR
    - PSRDADA

## PSRCAT 
PSRCAT is a library for interacting with the ANTF pulsar catalogue of all known pulsars, maintained by the Australian National Telescope Facility (ANTF). The catalogue includes all published rotation-powered pulsars and their unique profile, including those detected only at high energies. The pulsar profiles consist of their physical properies and their uncertainties. [1]. This package is a dependency of most, if not all, of the following packages and is used in RPNAV to identify the observed pulsar from its reference characteristics, and provide position information that is used in the navigation algorithms.[1] The python package psrqpy will also be used, which provides a tool for interacting with the ANTF pulsar catalogue. [2] 

## Install PSRCAT
See `Install Dependencies` section.

# TEMPO2

TEMPO2 is a widely-used library that is used in high precision timing of radio pulsars and the analysis of pulse arrival time estimates. In the RPNAV package, TEMPO2 is used to calculate the pulse timing residuals, given an known observed pulsar (in the form of a `.par` file) and and an arrival time file format (in the form of a  `.tim` file). The software can be used to calculate the barycentric arrival times, adjust for propagation delays and frequency-dependent parameters, form the timing residuals and carry out the weighted least-squares fit on the pulsar model. The resultant precise pulsar parameters and measured TOAs are used to calculate the timing delta and associated position correction. [3] [4]


## Install TEMPO2
See `Install Dependencies` section.

## PSRCHIVE 
PSRCHIVE is a library for the analysis of pulsar astronomical data and timing. In the RPNAV package, PARCHIVE is used to match the folded pulsar profile to a standard template, thereby identifying the observed profile.[5]
PARCHIVE also has a python interface that can be used for pulse profile matching and TOA calculation. 

### Install PSRCHIVE
Firstly, make sure you have the lastest version of GNU autotools installed (instructions can be found [here](https://psrchive.sourceforge.net/third/autotools/)). 

The install and build instructions for the package can be found (here)[https://psrchive.sourceforge.net/current/build.shtml]. 
```bash
git submodule add git://git.code.sf.net/p/psrchive/code psrchive
```

After downloading the source code from the Git repo, it can be compiled and installed with the following commands: 

```bash
cd psrchive
./bootstrap
./configure
```

### Install Dependencies 
The PSRCHIVE has inbuilt process to install the package dependencies listed above, most notable TEMPO2, PSRCAT, FFTW (a library for fast fourier transforms), and CFITSIO (a library for handling the Flexible Image Transport System data format). The installation instructions can be found [here](https://psrchive.sourceforge.net/third/install.shtml). 

```bash
./configure
./packages/fftw.csh
./packages/cfitsio.csh
./packages/pgplot.csh
# set the PGPLOT_DIR and PGPLOT_FONT environment variables
./packages/tempo2.csh
# set the TEMPO2 environment variable
./packages/psrcat.csh
# set the PSRCAT_FILE environment variable
./configure
```

To persist the environment variables configuring the location of the packages, it is recommended that you add the following to `.bashrc`:

```bash
# Pulsar definitions
export LD_LIBRARY_PATH=/usr/local/lib
export TEMPO2=/usr/local/tempo2
export PGPLOT_DIR=/usr/local/pgplot
export PGPLOT_FONT=$PGPLOT_DIR/grfont.dat
export PSRCAT_FILE=/usr/local/psrcat/psrcat.db
export PSRCAT_RUNDIR=/usr/local/bin/psrcat
```

To make the tempo2 scripts and ptasimulate scripts readily accessible, further add the following to `.bashrc`
```
alias ptasimulate="/home/naverill/external/ptasimulate/ptaSimulate"
alias tempo2="/usr/local/tempo2/bin/tempo2"
```

## DSPSR
DSPSR3 enables real-time phase-coherent dispersion removal and pulse folding for radio pulsar astronomical data.  The DSPSR software processes continuous streams of radio pulsar astronomical data, producing integrated statistics such as the phase-resolved average polarisation of the pulsar signal

- computation of phase-resolved averages (folding) using either a polynomial approximation to the pulsar phase model, a constant period, or acceleration search parameters
- simultaneous folding of multiple pulsars, such as globular cluster or double pulsars

 [6]

### Install DSPSR
The install and build instructions can be found (here)[https://dspsr.sourceforge.net/current/build.shtml]. 
```bash
git clone git://git.code.sf.net/p/dspsr/code dspsr
```

## PSRDADA

PSRDADA is library that supports Data Acquisition and Distributed Analysis (DADA) systems and is commonly used to handle baseband recording and processing instrumentation for pulsar astronomy.

### Install PSRDADA
The official install instructions can be found (here)[https://psrdada.sourceforge.net/download.shtml]
```
git clone git://git.code.sf.net/p/psrdada/code psrdada
```

## References
[1] Manchester, R. N., Hobbs, G. B., Teoh, A. & Hobbs, M., Astron. J., 129, 1993-2006 (2005) (astro-ph/0412641)
[2] Pitkin, (2018). psrqpy: a python interface for querying the ATNF pulsar catalogue. Journal of Open Source Software, 3(22), 538
[3] (TEMPO2 User Manual)[https://www.jb.man.ac.uk/research/pulsar/Resources/tempo2_manual.pdf]
[4] (TEMPO2 examples)[https://www.jb.man.ac.uk/~pulsar/Resources/tempo2_examples_ver1.pdf]
[5] van Straten, Willem, Paul Demorest, and Stefan Osłowski. "Pulsar data analysis with PSRCHIVE." arXiv preprint arXiv:1205.6276 (2012). (https://arxiv.org/pdf/1205.6276)
[6] van Straten, W., & Bailes, M. (2011). DSPSR: Digital Signal Processing Software for Pulsar Astronomy. Publications of the Astronomical Society of Australia, 28(1), 1–14. doi:10.1071/AS10021 