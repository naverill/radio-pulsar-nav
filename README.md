# Radio Pulsar Navigation Package

## Dependencies

### Install psrcat

The module can be downloaded (here)[https://www.atnf.csiro.au/research/pulsar/psrcat/download.html].
```
> cd psrcat
> source makeit
```

### Install GNUplot
```
brew install gnuplot
```

### Install PSRCHIVE
The install and build instructions can be found (here)[https://psrchive.sourceforge.net/current/build.shtml]
```
git submodule add git://git.code.sf.net/p/psrchive/code psrchive
```

After downloading the source code from the Git repo, it can be compiled and installed with the following commands: 

```
cd psrchive
./bootstrap
./configure
```

### Install DSPSR
The install and build instructions can be found (here)[https://dspsr.sourceforge.net/current/build.shtml]
```
git clone git://git.code.sf.net/p/dspsr/code dspsr
```

### Install TEMPO2
```
git clone https://bitbucket.org/psrsoft/tempo2.git
```

# Build
Build by making a build directory (i.e. build/), run cmake in that dir, and then use 
make to build the desired target.

```
> mkdir build && cd build
> cmake .. -DCMAKE_BUILD_TYPE=[Debug | Coverage | Release]
> make
> ./simulate    # Run simulator
```
