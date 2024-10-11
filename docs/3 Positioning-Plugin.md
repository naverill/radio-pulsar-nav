# Pulsar Positioning (PPOS) Plugin 

The PPOS Plugin is a TEMPO2 extension designed to localise an observer from observation of two or more pulsars. The plugin assumes that if the observer is located on the Earth's surface, observation of at least two pulsars and an Earth model is sufficient to resolve a position in 3D space. The algorithm uses a locally stored digital elevation model from the JAXA Advanced Land Observing Satellite (ALOS) dataset to resolve geodetic (latitude/longitude) coordinates into geocentric positions and fits the timing residuals to determine the true location of the observer. 

## Set Up
### Install
Download the repository from the following link and install according to the instructions in `1 Dependencies.md`.
```sh
git clone https://bitbucket.org/psrnav/tempo2/src/master/
cd tempo2
git checkout feature/positioning
```


### Dependencies
The package relies on a number of external plugins for usage.  
```sh
sudo apt update     
sudo apt install libgsl-dev       # GNU Scientific Library for gradient descent algorithm
sudo apt install proj-bin         # PROJ for coordinate transformation
sudo apt install libgdal-dev      # GDAL for TIFF processing 
```

###  Load DEM Files
The JAXA ALOS Global Digital Surface Model dataset is a collection of GEOTiff files with a file resolution of one degree latitude/longitude and a pixel resolution of one arcsecond (approx 30m). Each pixel stores an average elevation in metres that indicates the height of that point above sea level. The files can be downloaded from [JAXA Earth Observation Research Center](https://www.eorc.jaxa.jp/ALOS/en/aw3d30/data/) website. The code will look for the DSM files in the local install of tempo2 in the following folder:

```sh
mkdir ${TEMPO2}/map_data
```
where 
```sh
{$TEMPO} is the location of `TEMPO2` (usualy /usr/local/tempo2).
```

## Build
To build the plugin, navigate to the plugins/ directory of the respository and run the following:
```sh
sudo g++ -I/usr/local/tempo2/include -fPIC -shared -o {$TEMPO2}/plugins/pulsar_positioning_{$LOGIN_ARCH}_plug.t2 plugin/pulsar_positioning_plug.C -lgdal -lgsl -lproj
```
where 
```sh
{$LOGIN_ARCH} is the result of `uname` (usualy Linux).
```

This will create the plugin's `.t2` file in the TEMPO2 install folder.


## Run 
To execute the positioning script, run the following
```sh
tempo2 -gr pulsar_positioning -f ${/path/to/parfile} ${/path/to/timfile} -f ... 
```

The script has the following configuration parameters
```
-f   Define the absolute paths to the .par and .tim files [required]
        Form:
            -f {/path/to/.par} {/path/to/.tim}
-a   Define the algorithm to use for fitting [required]
        Values:
            grid: Grid search algorithm
            grde: Gradient Descent algorithm
-i   Define the maximum number of fit iterations for TEMPO2 
        Default: 2
-e   Define the method to check the fit error 
        Values:
            rms: Root mean squared error
            chi: reduced chi-squared error  
        Default: rms
-o   Define the observer code
-l   Define the starting search coordinates
        Form:
            -l {long (deg)} {lat (deg)}
-n  Fill NODATA errors with 0 elevation
-n  Randomly generate initial starting position
```

```sh
tempo2 -gr pulsar_positioning -f J0835-4510.tdb.par J0835-4510.tim -f J1939+2134.tdb.par J1939+2134.tim -l 150 -32 -a grde -e chi
```

```
archive time 
pat 
```