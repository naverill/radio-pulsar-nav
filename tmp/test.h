/**
 * @brief Pulsar Positioning (PPOS) Plugin
 * 
 * @details The PPOS Plugin is a TEMPO2 extension designed to localise an observer from 
 * observation of two or more pulsars. The plugin assumes that if the observer is 
 * located on the Earth's surface, observation of at least two pulsars and an Earth 
 * model is sufficient to resolve a position in 3D space. The algorithm uses a locally 
 * stored digital elevation model from the JAXA ALOS dataset to resolve geodetic 
 * (latitude/longitude) coordinates into geocentric positions and fits the timing 
 * residuals to determine the true location of the observer. 
 */
#ifndef PULSAR_POSITIONING_PLUG_H
#define PULSAR_POSITIONING_PLUG_H

/**********************************************************
 *  Includes
 *********************************************************/
#include "geokeys.h"
#include "tempo2.h"

/**********************************************************
 *  Macro Definitions
 *********************************************************/
#define MAX_OBSERVER_LEN            20     /**< @brief Maximum allowable length for observer code names */
#define MAX_GRID_SIZE               3600   /**< @brief Maximum allowable pixel size for GEOTiff in x and y directions */
#define MAX_GSL_ITER                500    /**< @brief Maximum iterations for gradient descent algorithm */
#define MAX_GRADIENT_DESCENT_ERROR  100000 /**< @brief Arbitrarily high value to set gradient descent iteration that fail to resolve a position */
#define MAX_FILENAME_LEN            100    /**< @brief Maximum length for absolute paths to files */
#define GDAL_NODATA_VALUE           -9999  /**< @brief Value in GEOTiff indicating no value was recorded */
#define EARTH_RADIUS                6371   /**< @brief Average radius of Earth */
#define EARTH_RADIUS_FLATENNING     0.003353   /**< @brief Error in Earth radius*/

/**********************************************************
 *  Structure Definitions
 *********************************************************/
/**
 * @brief Status Enum for Pulsar Positioning plugin
 */
typedef enum PPOS_Status_t
{
    PPOS_SUCCESS                = 0,   /**< @brief Operation was successful */
    PPOS_ERROR_GDAL             = -1,  /**< @brief  Error occurred with GDAL API */
    PPOS_ERROR_PROJ             = -2,  /**< @brief  Error occurred with PROJ API */
    PPOS_ERROR_FOPEN            = -3,  /**< @brief  Error occurred opening file */
    PPOS_INVALID_OBSERVER       = -4,  /**< @brief  TEMPO2 observatory code was invalid */
    PPOS_INVALID_GRID_SIZE      = -5,  /**< @brief  Elevation model grid size is invalid */
    PPOS_INVALID_PTR            = -7,  /**< @brief  Null pointer encountered */
    PPOS_ERROR_GSL              = -8,  /**< @brief  Error occurred in GSL library during gradient descent algorithm */
    PPOS_WARNING_OUTOFBOUNDS    = -9,  /**< @brief  Pixel coordiinates exceed model grid size  */
    PPOS_ERROR_INVALID_POSITION = -10, /**< @brief  Invalid geodetic position encountered  */
    PPOS_WARNING_NODATA         = -11, /**< @brief  No data value encountered in elevation model */
    PPOS_ERROR_NOT_IMPLEMENTED  = -12, /**< @brief  Operation has not been implemented */
    PPOS_INVALID_PARAMS         = -13, /**< @brief  Invalid fit paramaters */
} PPOS_Status_t;

/**
 * @brief Enum for type of fit algorithm
 */
typedef enum PPOS_FitAlgorithm
{
    PPOS_FitAlgorithm_UNKNOWN           = 0,    /**< @brief Unknown fit algorithm */
    PPOS_FitAlgorithm_GRID              = 1,    /**< @brief Grid search fit algorithm */
    PPOS_FitAlgorithm_GRADIENT_DESCENT  = 2,    /**< @brief Gradient descent algorithm */
} PPOS_FitAlgorithm_t;

/**
 * @brief Enum for checking the fit of the timing residuals
 */
typedef enum PPOS_FitCriterion
{
    PPOS_FitCriterion_RMS               = 1,    /**< @brief Root Mean Square (RMS) timing residuals  */
    PPOS_FitCriterion_CHISQUARED        = 2,    /**< @brief Reduced Chi-Squared value */
} PPOS_FitCriterion_t;

/**
 * @brief Geotransformation parameters
 * 
 * @details A geotransform is an affine transformation from the image coordinate 
 *  space (row, column), also known as (pixel, line) to the georeferenced coordinate 
 *  space (projected or geographic coordinates
 */
typedef struct PPOS_GTransform
{
    double xCoord;      /**< @brief x-coordinate of the upper-left corner of the upper-left pixel. */
    double wePixelRes;  /**< @brief w-e pixel resolution / pixel width. */
    double rowRot;      /**< @brief row rotation (typically zero) */
    double yCoord;      /**< @brief y-coordinate of the upper-left corner of the upper-left pixel. */
    double colRot;      /**< @brief column rotation (typically zero) */
    double nsPixelRes;  /**< @brief n-s pixel resolution / pixel height (negative value for a north-up image) */
} PPOS_GTransform_t;


/**
 * @brief Elevation Model Configuration
 */
typedef struct PPOS_ElevationModel
{
    char *filePath;/**< @brief Path to folder with all GEOTiff files */
    char *fileName;/**< @brief Name of elevation model GEOTiff file */
    PJ* P;                          /**< @brief PROJ  projection/transformation information */
    GDALDatasetH  hDataset;         /**< @brief GDAL Dataset */
    GDALRasterBandH hBand;          /**< @brief Elevation Raster band */
    PPOS_GTransform_t gTransform;   /**< @brief Coordinate transformation object */
    size_t xSize;                   /**< @brief Size of GEOTiff X-axis */
    size_t ySize;                   /**< @brief Size of GEOTiff Y-axis */
    bool initialised;               /**< @brief Model has been initialised */
} PPOS_ElevationModel_t;

/**
 * @brief Fit Algorithm Configuration
 */
typedef struct PPOS_FitParams
{
    int argc;                               /**< @brief Number of input arguments */
    char **argv;                            /**< @brief Input arguments */
    pulsar *psrPtr;                         /**< @brief Array of pulsars to use in fitting routine */
    int npsr;                               /**< @brief Number of pulsars to fit */
    PPOS_ElevationModel_t *elevModelPtr;    /**< @brief Pointer to elevation model configuration */
    PPOS_FitAlgorithm_t fitAlg;             /**< @brief Position fitting algorithm  */
    PPOS_FitCriterion_t fitCriteria;        /**< @brief Method of calculating error of position fit based on timing residuals  */
    float *pafScanline;                     /**< @brief Stored Y-axis GEOTiff (saved for memory efficiency) */
    int yLine;                              /**< @brief Current Y-axis value */
    int maxiter;                            /**< @brief Number of fitting iterations for TEMPO2 */
    bool nodataFill;                        /**< @brief Fill NODATA errors with zero elevation */
} PPOS_FitParams_t;


/**********************************************************
 *  Internal Function Definitions
 *********************************************************/
/**
 * @brief Main function to configure and run the fit algorithm
 * 
 * @param[in]   paramsPtr  Fit algorithm configuration
 * @param[in]   observer   TEMPO2 code for the observatory
 * 
 * @returns Result of the operation
 * @retval  PPOS_SUCCESS            #PPOS_SUCCESS
 * @retval  PPOS_INVALID_OBSERVER   #PPOS_INVALID_OBSERVER   
 */
PPOS_Status_t PPOS_processPosition(PPOS_FitParams_t *paramsPtr, char *observer, PJ_LPZT startPos);

/**
 * @brief Run the grid search Algorithm find the observer position
 * 
 * @details Iterate over all grid cells in the input TIFF file and fit residuals
 * 
 * @param[in]   paramsPtr  Fit algorithm configuration
 * 
 * @returns Result of the operation
 * @retval  PPOS_SUCCESS        #PPOS_SUCCESS
 * @retval  PPOS_ERROR_FOPEN    #PPOS_ERROR_FOPEN    
 */
PPOS_Status_t PPOS_fitGrid(PPOS_FitParams_t *params, PJ_LPZT startPos);

/**
 * @brief Run the gradient descent search algorithm to find the observer position
 * 
 * @details TODO
 * 
 * @param[in]   paramsPtr  Fit algorithm configuration
 * @param[in]   initPosPtr  Fit algorithm configuration
 * 
 * @returns Result of the operation
 * @retval  PPOS_SUCCESS        #PPOS_SUCCESS
 * @retval  PPOS_ERROR_GSL  #PPOS_ERROR_GSL      
 */
PPOS_Status_t PPOS_fitPosition(PPOS_FitParams_t *paramsPtr, PJ_LPZT startPos);

/**
 * @brief Run the gradient descent search algorithm to find the observer position
 * 
 * @details TODO
 * 
 * @param[in]   v           Position vector
 * @param[in]   params      Fit algorithm configuration
 * 
 * @returns Result of the operation
 * @retval  PPOS_SUCCESS        #PPOS_SUCCESS
 * @retval  PPOS_ERROR_GSL  #PPOS_ERROR_GSL      
 */
double PPOS_getPositionError(const gsl_vector *v, void *params);

/**
 * @brief Fit the timing residuals based on the observer position using TEMPO2
 * 
 * @details TODO
 * 
 * @param[in]   paramsPtr   Fit algorithm configuration
 * @param[in]   posPtr      Observer position (in geodetic coordinate frame WS)
 * 
 * @returns Result of the operation
 * @retval  PPOS_SUCCESS        #PPOS_SUCCESS
 * @retval  PPOS_ERROR_GSL  #PPOS_ERROR_GSL      
 */
PPOS_Status_t PPOS_callFit(PPOS_FitParams_t *paramsPtr, PJ_XYZT *posPtr);

/**
 * @brief Calculate the fit
 * 
 * @details TODO
 * 
 * @param[in]   psrPtr      List of pulsars
 * @param[in]   npsr        Number of pulsars
 * 
 * @returns Calculate residual error
 */
long double PPOS_calculateError(PPOS_FitParams_t *paramsPtr);

/**
 * @brief Calculate the RMS timing residual post-fit
 * 
 * @details TODO
 * 
 * @param[in]   psrPtr      List of pulsars
 * @param[in]   npsr        Number of pulsars
 * 
 * @returns Calculate residual error
 */
long double PPOS_calculateRMSResiduals(pulsar *psrPtr, int npsr);

/**
 * @brief Calculate the Reduced Chi-Squared value post-fit
 * 
 * @details Chi-Squared calculation from 
 * Purver, Mark Benedict. High-precision pulsar timing: The stability of integrated pulse 
 * profiles and their representation by analytic templates. Diss. The University of Manchester 
 * (United Kingdom), 2011. [Equation 2.2]
 * 
 * @param[in]   psrPtr      List of pulsars
 * @param[in]   npsr        Number of pulsars
 * 
 * @returns Calculate residual error
 */
long double PPOS_calculateChiSquared(pulsar *psrPtr, int npsr);

/**
 * @brief Load the elevation model from an ALOS geotiff file
 * 
 * @details TIFF contains a digital elevation model over the range of 1 degree lat/long
 * 
 * @param[in]   filePath        Path to the elevation model file
 * @param[in]   fileName        Name of the elevation model file
 * @param[in]   elevModelPtr    Pointer to the elevation model object
 * 
 * @returns Result of the operation
 * @retval  PPOS_SUCCESS               #PPOS_SUCCESS
 * @retval  PPOS_ERROR_PROJ            #PPOS_ERROR_PROJ
 * @retval  PPOS_ERROR_GDAL            #PPOS_ERROR_GDAL
 * @retval  PPOS_INVALID_GRID_SIZE     #PPOS_INVALID_GRID_SIZE
 */
PPOS_Status_t PPOS_loadElevationModel(PPOS_ElevationModel_t *elevModelPtr);

/**
 * @brief Convert Latitude/Longitude coordinate to geocentric coordinates 
 * 
 * @details The function gets the elevation from the digital model and converts the geodetic 
 *  position into geocentric coordinates
 * 
 * @param[in]   elevModelPtr    Pointer to the elevation model object
 * @param[in]   geodet          Geodetic coordinate
 * @param[out]  geocenPtr       Pointer to the geocentric coordinate
 * 
 * @returns Result of the operation
 * @retval  PPOS_SUCCESS        #PPOS_SUCCESS
 * @retval  PPOS_ERROR_PROJ     #PPOS_ERROR_PROJ
 */
PPOS_Status_t PPOS_convertLatLongToGeocentric(PPOS_ElevationModel_t *elevModelPtr, PJ_COORD geodet, PJ_COORD *geocenPtr);

/**
 * @brief Convert pixel coordinates in DEM to Latitude/Longitude coordinate 
 * 
 * @param[in]   elevModelPtr    Pointer to the elevation model object
 * @param[in]   xPixel          X pixel coordinate
 * @param[in]   yPixel          Y pixel coordinate
 * @param[out]  geodetPtr       Pointer toGeodetic coordinate
 * 
 * @returns Result of the operation
 * @retval  PPOS_SUCCESS        #PPOS_SUCCESS
 * @retval  PPOS_ERROR_PROJ     #PPOS_ERROR_PROJ
 */
PPOS_Status_t PPOS_convertPixelIndexToLatLong(PPOS_ElevationModel_t *elevModelPtr, int16_t xPixel, int16_t yPixel, PJ_LPZT *geodet);

/**
 * @brief Convert Latitude/Longitude coordinate to pixel coordinates in DEM 
 * 
 * @param[in]   elevModelPtr    Pointer to the elevation model object
 * @param[in]   geodet          Geodetic coordinate
 * @param[in]   xPixelPtr       Pointer to X pixel coordinate
 * @param[in]   yPixelPtr       Pointer to Y pixel coordinate
 * 
 * @returns Result of the operation
 * @retval  PPOS_SUCCESS                #PPOS_SUCCESS
 * @retval  PPOS_WARNING_OUTOFBOUNDS    #PPOS_WARNING_OUTOFBOUNDS
 */
PPOS_Status_t PPOS_convertLatLongToPixelIndex(PPOS_ElevationModel_t *elevModelPtr, PJ_LPZT geodet, int16_t *xPixelPtr, int16_t *yPixelPtr);

/**
 * @brief Convert geodetic coordinate to geocentric 
 * 
 * @param[in]   P               Pointer to PROJ transformation object               
 * @param[in]   geodet          Geodetic coordinate
 * @param[out]  geocenPtr       Pointer to the geocentric coordinate
 * 
 * @returns Result of the operation
 * @retval  PPOS_SUCCESS     #PPOS_SUCCESS
 */
PPOS_Status_t PPOS_convertGeodeticToGeocentric(PJ* P, PJ_LPZT geodet, PJ_XYZT *geocent);

/**
 * @brief Set current position of observer 
 * 
 * @brief Configures the observer position for TEMPO2 fitting routine 
 * 
 * @param[in]   coordPtr       Pointer to Geocentric observer position               
 * 
 * @returns Result of the operation
 * @retval  PPOS_SUCCESS        #PPOS_SUCCESS
 * @retval  PPOS_INVALID_PTR    #PPOS_INVALID_PTR
 */
PPOS_Status_t PPOS_setObserverLocation(PJ_XYZT *coordPtr);

/**
 * @brief Load a new Y-band from the elevation model
 * 
 * @param[in]   elevModel       Pointer to the elevation model               
 * @param[out]  pafScanline     Array of elevation values               
 * @param[out]  yLine           Y coordinate of row                
 * 
 * @returns Result of the operation
 * @retval  PPOS_SUCCESS     #PPOS_SUCCESS
 * @retval  PPOS_ERROR_GDAL  #PPOS_ERROR_GDAL
 */
PPOS_Status_t PPOS_loadCoordinates(PPOS_ElevationModel_t *elevModel, float *pafScanline, int yLine);

/**
 * @brief Load a new elevation model from an adjacent grid
 * 
 * @details Called when the lat/long value exceeds the bounds of the current grid. The function uses the 
 *  out-of-bounds coordinate to calculate the new model and loads it into the  elevation model pointer. Assumes 
 * the elevation model naming convention follows the ALOS naming format 
 * 
 * @param[in, out]  elevModel       Pointer to the elevation model               
 * @param[out]      pafScanline     Array of elevation values               
 * @param[out]      yLine           Y coordinate of row                
 * 
 * @returns Result of the operation
 * @retval  PPOS_SUCCESS              #PPOS_SUCCESS
 * @retval  PPOS_ERROR_GDAL           #PPOS_ERROR_GDAL
 * @retval  PPOS_WARNING_OUTOFBOUNDS  #PPOS_WARNING_OUTOFBOUNDS
 */
PPOS_Status_t PPOS_loadAdjacentGrid(PPOS_ElevationModel_t *elevModelPtr, int16_t xPixel, int16_t yPixel);

/**
 * @brief Load a new elevation model
 * 
 * @details Called when the lat/long value exceeds the bounds of the current grid. The function uses the 
 *  out-of-bounds coordinate to calculate the new model and loads it into the  elevation model pointer. Assumes 
 * the elevation model naming convention follows the ALOS naming format 
 * 
 * @param[in, out]  elevModel       Pointer to the elevation model               
 * @param[out]      pafScanline     Array of elevation values               
 * @param[out]      yLine           Y coordinate of row                
 * 
 * @returns Result of the operation
 * @retval  PPOS_SUCCESS              #PPOS_SUCCESS
 * @retval  PPOS_ERROR_GDAL           #PPOS_ERROR_GDAL
 * @retval  PPOS_WARNING_OUTOFBOUNDS  #PPOS_WARNING_OUTOFBOUNDS
 */
PPOS_Status_t PPOS_loadGrid(PPOS_ElevationModel_t *elevModelPtr, PJ_LPZT geodet);

/**
 * @brief TODO
 */
PPOS_Status_t PPOS_estimateSolutions(PPOS_FitParams_t *paramsPtr);

/**
 * @brief Opens the results file
 * 
 * @param[in]  fileName       Name of results file                             
 * 
 * @returns Result of the operation
 * @retval  PPOS_SUCCESS              #PPOS_SUCCESS
 * @retval  PPOS_ERROR_GDAL           #PPOS_ERROR_GDAL
 */
PPOS_Status_t PPOS_openResults(char *fileName);

/**
 * @brief Write one row of residual errors to results file
 * 
 * @param[in]  paramsPtr    Fit algorithm configuration
 * @param[in]  fileName     Name of results file                             
 * @param[in]  residErrMap  Residual error array                             
 * 
 * @returns Result of the operation
 * @retval  PPOS_SUCCESS              #PPOS_SUCCESS
 * @retval  PPOS_ERROR_FOPEN          #PPOS_ERROR_FOPEN
 */
PPOS_Status_t PPOS_writeResults(PPOS_FitParams_t *paramsPtr, char * fileName, long double residErrMap[MAX_GRID_SIZE]);

/**
 * @brief Display help
 * 
 * @param None
 * 
 * @return None 
 * */
void help();

#endif /** PULSAR_POSITIONING_PLUG_H */
