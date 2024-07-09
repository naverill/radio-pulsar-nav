# Dependencies

The RPNAV package has the following dependencies:
    - PSRCAT
    - PSARCHIVE 
    - DSPSR
    - TEMPO2

## Install PSRCHIVE
Firstly, make sure you have the lastest version of GNU autotools installed (instructions can be found [here](https://psrchive.sourceforge.net/third/autotools/)). 

The install and build instructions for the package can be found (here)[https://psrchive.sourceforge.net/current/build.shtml]. 
```
git submodule add git://git.code.sf.net/p/psrchive/code psrchive
```

After downloading the source code from the Git repo, it can be compiled and installed with the following commands: 

```
cd psrchive
./bootstrap
./configure
```

The PSRCHIVE has inbuilt process to install the package dependencies listed above, most notable TEMPO2 and PSRCAT. The installation instructions can be found [here](https://psrchive.sourceforge.net/third/install.shtml). To persist the environment variables configuring the location of the packages, it is recommended that you add the following to `.bashrc`:

```
# Pulsar definitions
LD_LIBRARY_PATH=/usr/local/lib
TEMPO2=/usr/local/tempo2
PGPLOT_DIR=/usr/local/pgplot
PGPLOT_FONT=$PGPLOT_DIR/grfont.dat
PSRCAT_FILE=/usr/local/psrcat/psrcat.db
```

it is recommended 

## Install DSPSR
The install and build instructions can be found (here)[https://dspsr.sourceforge.net/current/build.shtml]. 
```
git clone git://git.code.sf.net/p/dspsr/code dspsr
```
