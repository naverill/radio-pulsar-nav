/* Plugin to perform positioning from pulsar TOAs */

/**********************************************************
 *  Includes
 *********************************************************/
/** Standard C packages */
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

/** External packages  */
#include <proj.h>
#include <gdal.h>
#include <gdal_priv.h>
#include <gsl/gsl_multimin.h>

/** Tempo2 packages */
#include "t2fit.h"
#include "pulsar_positioning_plug.h"

using namespace std;

/**********************************************************
 *  Global Parameters Definitions
 *********************************************************/
observatory *observerLocPtr;

/**********************************************************
 *  External Function(s)
 *********************************************************/
/* The main function called from the TEMPO2 package is 'graphicalInterface' */
/* Therefore this function is required in all plugins                       */
extern "C" int graphicalInterface(int argc,char *argv[],pulsar *psrPtr,int *npsrPtr) 
{
    PPOS_Status_t status = PPOS_SUCCESS;
    PPOS_FitParams_t params;
    PPOS_ElevationModel_t elevModel;
    char parFile[MAX_PSR_VAL][MAX_FILELEN];
    char timFile[MAX_PSR_VAL][MAX_FILELEN];
    char filePath[] = "/usr/local/tempo2/map_data";
    char fileName[] = " ";

    /** Configurable fit parameters */
    int maxiter = 2;
    char observer[] = "NAVIGATE";
    params.fitAlg = PPOS_FitAlgorithm_UNKNOWN;
    params.fitCriteria = PPOS_FitCriterion_RMS;
    params.nodataFill = false;
    PJ_LPZT startPos;
    startPos.lam = 0;
    startPos.phi = 0;
    *npsrPtr = 0;  /* For a graphical interface that only shows results for one pulsar */

    printf("Graphical Interface: pulsarPositioning\n");
    printf("Author:              N. Averill\n");
    printf("Version:             1.0\n");
    printf(" --- type 'h' for help information\n");

    /* Obtain all parameters from the command line */
    for (int arg = 2; arg < argc; arg++)
    {
        /** Load par/tim files */
        if (strcmp(argv[arg],"-f")==0)
        {
            strcpy(parFile[*npsrPtr],argv[++arg]); 
            strcpy(timFile[*npsrPtr],argv[++arg]);
            (*npsrPtr)++;
        }
        /** Set TEMPO2 iteration count */
        else if (strcmp(argv[arg],"-i")==0) 
        {
	        sscanf(argv[++arg],"%d",&maxiter);
        }
        /** Set observer name */
        else if (strcmp(argv[arg],"-o")==0) 
        {
	        sscanf(argv[++arg],"%s", observer);
        }
        /** Read long/lat of start position */
        else if (strcmp(argv[arg],"-l")==0) 
        {
	        sscanf(argv[++arg],"%lf", &startPos.lam);
	        sscanf(argv[++arg],"%lf", &startPos.phi);
        }
        /** Fill nodata with zero elevation */
        else if (strcmp(argv[arg],"-n") == 0) 
        {
	        params.nodataFill = true;
        }
        /** Set fit error estimator  */
        else if (strcmp(argv[arg],"-e")==0) 
        {
            if (strcmp(argv[arg],"rms") == 0)
            {
                params.fitCriteria = PPOS_FitCriterion_RMS;
            }
            else if (strcmp(argv[arg],"chi") == 0)
            {
                params.fitCriteria = PPOS_FitCriterion_CHISQUARED;
            }
            arg++;
        }
        else if (strcmp(argv[arg],"-a")==0) 
        {
            if (strcmp(argv[arg + 1],"grid") == 0)
            {
                params.fitAlg = PPOS_FitAlgorithm_GRID;
            }
            else if (strcmp(argv[arg + 1], "grde") == 0)
            {
                params.fitAlg = PPOS_FitAlgorithm_GRADIENT_DESCENT;
            }
            arg++;
        }
    }

    if (0 >= *npsrPtr)
    {
        status = PPOS_INVALID_PARAMS;
        printf("[%d:%s] ERROR must input a .par and .tim file\n", __LINE__, __func__);
        help();
    }

    if (PPOS_SUCCESS == status)
    {
        /* Load the parameters */
        readParfile(psrPtr, parFile, timFile,*npsrPtr); 

        /* Load the arrival times */
        readTimfile(psrPtr,timFile,*npsrPtr); 

        GDALAllRegister();

        /** Load params into Fit Config */
        params.argc = argc;
        params.argv = argv;
        params.maxiter = maxiter;
        params.npsr = *npsrPtr;
        params.psrPtr = psrPtr;
        params.elevModelPtr = &elevModel;

        params.elevModelPtr->fileName = (char*) malloc(sizeof(char) * MAX_FILENAME_LEN);
        if (NULL == params.elevModelPtr->fileName)
        {
            status = PPOS_INVALID_PTR;
        }
    }

    if (PPOS_SUCCESS == status)
    {
        params.elevModelPtr->filePath = (char*) malloc(sizeof(char) * MAX_FILENAME_LEN);
        if (NULL == params.elevModelPtr->fileName)
        {
            status = PPOS_INVALID_PTR;
        }
        else
        {
            sprintf(params.elevModelPtr->filePath, "%s", filePath);
        }
    }

    /** Run positioning */
    if (PPOS_SUCCESS == status)
    {
        status = PPOS_processPosition(&params, observer, startPos);
        if (PPOS_SUCCESS != status)
        {
            printf("[%d:%s] ERROR %d processing observer position\n", __LINE__, __func__, status);
        }
    }

    if (NULL != params.elevModelPtr->fileName)
    {
        free(params.elevModelPtr->fileName);
    }
    if (NULL != params.elevModelPtr->filePath)
    {
        free(params.elevModelPtr->filePath);
    }
    return status;
}

/**********************************************************
 *  Internal Function(s)
 *********************************************************/
/**
 * Description in header file
 */
PPOS_Status_t PPOS_processPosition(PPOS_FitParams_t *paramsPtr, char *observer, PJ_LPZT startPos)
{
    PPOS_Status_t status = PPOS_SUCCESS;
    observatory posTrue;    

    /** Load true observatory position */
    if (PPOS_SUCCESS == status)
    {
        observerLocPtr = getObservatory(observer);
        if (NULL == observerLocPtr)
        {
            printf("[%d:%s] ERROR: Invalid observatory name %s\n", __LINE__, __func__, observer);
            status = PPOS_INVALID_OBSERVER;
        }
        else
        {
            memcpy(&posTrue, observerLocPtr, sizeof(*observerLocPtr));
        }
    }

    if (PPOS_SUCCESS == status)
    {   
        paramsPtr->yLine = -1;

        /** Fit for entire grid and calculate statistics */
        if (PPOS_FitAlgorithm_GRID == paramsPtr->fitAlg)
        {
            status = PPOS_fitGrid(paramsPtr, startPos);
        }
        /** Fit for entire grid and calculate statistics */
        else if (PPOS_FitAlgorithm_GRADIENT_DESCENT == paramsPtr->fitAlg)
        {
            /** Set initial guess to top left corner of elevation model */
            status = PPOS_fitPosition(paramsPtr, startPos);
        }
        else
        {
            status = PPOS_ERROR_NOT_IMPLEMENTED;
        }
    }

    if (NULL != paramsPtr->elevModelPtr->hDataset)
    {
        GDALClose(paramsPtr->elevModelPtr->hDataset);
    }
    return status;
}

/**
 * Description in header file
 */
PPOS_Status_t PPOS_fitGrid(PPOS_FitParams_t *paramsPtr, PJ_LPZT startPos)
{
    PPOS_Status_t status = PPOS_SUCCESS;
    char fileName[MAX_FILENAME_LEN];
    /** Create results map */
    long double residErrMap[MAX_GRID_SIZE];

    /** Format output path name */
    char *fExt;
    size_t index;
    fExt = strchr(paramsPtr->elevModelPtr->fileName, '.');
    if (NULL == fExt)
    {
        status =  PPOS_ERROR_FOPEN;
    }

    /** Load grid */
    if (PPOS_SUCCESS == status)
    {
        status = PPOS_loadGrid(paramsPtr->elevModelPtr, startPos);
        if (PPOS_SUCCESS != status)
        {
            printf("[%d:%s] ERROR %d loading elevation model grid\n", __LINE__, __func__, status);
        }
    }

    /** Find index of file extension */
    if (PPOS_SUCCESS == status)
    {
        index = (fExt - paramsPtr->elevModelPtr->fileName);
        sprintf(fileName, "%-.*s_rmse.csv", (int)index, paramsPtr->elevModelPtr->fileName);

        PPOS_openResults(fileName);
    }

    /** Fit for Earth model grid */
    for (int yi = 0; yi < paramsPtr->elevModelPtr->ySize; yi++)
    {
        memset(residErrMap, 0, sizeof(residErrMap));
        paramsPtr->pafScanline = (float *) CPLMalloc(sizeof(float) * paramsPtr->elevModelPtr->xSize);

        status = PPOS_loadCoordinates(paramsPtr->elevModelPtr, paramsPtr->pafScanline, yi);
        if (PPOS_SUCCESS != status)
        {
            CPLFree(paramsPtr->pafScanline);
            continue;
        }

        for (int xi = 0; xi < paramsPtr->elevModelPtr->xSize; xi++)
        {
            PJ_LPZT geodet;
            status = PPOS_convertPixelIndexToLatLong(paramsPtr->elevModelPtr, xi, yi, &geodet);
            geodet.z = paramsPtr->pafScanline[xi]; /* Set elevation in meters */
            PJ_XYZT geocent;
            status = PPOS_convertGeodeticToGeocentric(paramsPtr->elevModelPtr->P, geodet, &geocent);

            if (PPOS_SUCCESS == status)
            {
                status = PPOS_callFit(paramsPtr, &geocent);
            }

            if (PPOS_SUCCESS == status)
            {
                residErrMap[xi] = PPOS_calculateError(paramsPtr);
            }

        }
        CPLFree(paramsPtr->pafScanline);

        status = PPOS_writeResults(paramsPtr, fileName, residErrMap);
    }
    return status;
}

/**
 * Description in header file
 */
PPOS_Status_t PPOS_fitPosition(PPOS_FitParams_t *paramsPtr, PJ_LPZT startPos)
{
    PPOS_Status_t status = PPOS_SUCCESS;
    const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
    gsl_multimin_fminimizer *s = NULL;
    gsl_vector *stepSize;
    gsl_vector *latlong;
    gsl_multimin_function minex_func;
    int nDim = 2;
    size_t iter = 0;
    int result;
    double size;

    paramsPtr->pafScanline = (float *) CPLMalloc(sizeof(float) * paramsPtr->elevModelPtr->xSize);

    if (PPOS_SUCCESS == status)
    {
        /* Starting point */
        latlong = gsl_vector_alloc (nDim);
        gsl_vector_set (latlong, 0, startPos.lam); /* Longitude in degrees */
        gsl_vector_set (latlong, 1, startPos.phi); /* Latitude in degrees */

        /* Set initial step sizes to 1 */
        stepSize = gsl_vector_alloc (nDim);
        gsl_vector_set_all (stepSize, 1);

        /* Initialize method and iterate */
        minex_func.n = nDim;
        minex_func.f = PPOS_getPositionError;
        minex_func.params = paramsPtr;

        s = gsl_multimin_fminimizer_alloc (T, 2);
        gsl_multimin_fminimizer_set (s, &minex_func, latlong, stepSize);

        do
        {
            iter++;
            result = gsl_multimin_fminimizer_iterate(s);

            if (result)
            {
                break;
            }

            size = gsl_multimin_fminimizer_size (s);
            result = gsl_multimin_test_size (size, 1e-2);

            if (GSL_SUCCESS == result)
            {
                printf ("[%d:%s] converged to minimum at\n", __LINE__, __func__);
            }

            printf ("[%d:%s] %ld %10.3e %10.3e f() = %7.3f size = %.3f\n",
                    __LINE__, __func__,
                    iter,
                    gsl_vector_get (s->x, 0),
                    gsl_vector_get (s->x, 1),
                    s->fval, size);
        }
        while (GSL_CONTINUE == result && iter < MAX_GSL_ITER);

        if (GSL_SUCCESS != result)
        {
            status = PPOS_ERROR_GSL;
        }

        CPLFree(paramsPtr->pafScanline);
        gsl_vector_free(latlong);
        gsl_vector_free(stepSize);
        gsl_multimin_fminimizer_free (s);
    }

    return status;
}

/**
 * Description in header file
 */
double PPOS_getPositionError(const gsl_vector *v, void *paramsPtr)
{
    PPOS_Status_t status = PPOS_SUCCESS;
    long double error = MAX_GRADIENT_DESCENT_ERROR;

    PJ_LPZT geodet;
    PJ_XYZT geocen;
    int16_t xPixel;
    int16_t yPixel;

    geodet.lam = gsl_vector_get(v, 0);  /* Longitude in degrees */
    geodet.phi = gsl_vector_get(v, 1);  /* Latitude in degrees */

    PPOS_FitParams_t *p = (PPOS_FitParams_t *)paramsPtr;

    /** Attempt to load model */
    if (!p->elevModelPtr->initialised)
    {
        status = PPOS_loadGrid(p->elevModelPtr, geodet);
    }
    /** Model is loaded, check position coordinates and load new model if necessary */
    else
    {
        status = PPOS_convertLatLongToPixelIndex(p->elevModelPtr, geodet, &xPixel, &yPixel);

        /** Position exceeds current earth model, load new */
        if (PPOS_WARNING_OUTOFBOUNDS == status)
        {
            status = PPOS_loadAdjacentGrid(p->elevModelPtr, xPixel, yPixel);
            /**
             * If model is loaded, update pixel values 
             * if we can't load new model, set error to arbitrarily large value */
            if (PPOS_SUCCESS == status)
            {
                /** Recalculate pixel positions */
                p->yLine = -1;
                status = PPOS_convertLatLongToPixelIndex(p->elevModelPtr, geodet, &xPixel, &yPixel);
            }
        }
    }

    /** Y coordinate is different from last iteration, overwrite scanline load new elevation data */
    if ((PPOS_SUCCESS == status) && (yPixel != p->yLine))
    {
        CPLFree(p->pafScanline);
        p->pafScanline = (float *) CPLMalloc(sizeof(float) * p->elevModelPtr->xSize);
        status = PPOS_loadCoordinates(p->elevModelPtr, p->pafScanline, yPixel);
    }

    /** Transform coordinates */
    if (PPOS_SUCCESS == status)
    {
        /** Valid elevation value */
        if (GDAL_NODATA_VALUE != p->pafScanline[xPixel])
        {
            geodet.z = p->pafScanline[xPixel]; /* Height in meters */ 
        }
        /** Assume zero elevation */
        else if (p->nodataFill)
        {
            geodet.z = 0;
        }
        else
        {
            /** do nothing */
            status = PPOS_WARNING_NODATA;
        }
    }
    /** Only fill data if parameter is set and position is valid geodetic */
    else if (PPOS_ERROR_INVALID_POSITION != status && p->nodataFill)
    {
        geodet.z = 0;
        status = PPOS_SUCCESS;
    }

    if (PPOS_SUCCESS == status)
    {
        status = PPOS_convertGeodeticToGeocentric(p->elevModelPtr->P, geodet, &geocen);
    }

    if (PPOS_SUCCESS == status)
    {
        status = PPOS_callFit(p, &geocen);
    }

    if (PPOS_SUCCESS == status)
    {
        error = PPOS_calculateError(p);
    }

    return error;
}

/**
 * Description in header file
 */
PPOS_Status_t PPOS_callFit(PPOS_FitParams_t *paramsPtr, PJ_XYZT *posPtr)
{
    PPOS_Status_t status = PPOS_SUCCESS;
    double globalParameter;

    status = PPOS_setObserverLocation(posPtr);

    if (PPOS_SUCCESS == status)
    {
        preProcess(paramsPtr->psrPtr, paramsPtr->npsr, paramsPtr->argc, paramsPtr->argv);

        /* Do two iterations for pre- and post-fit residuals*/
        for (int i = 0; i < paramsPtr->maxiter; i++)
        {
            /* Form the barycentric arrival times */
            formBatsAll(paramsPtr->psrPtr, paramsPtr->npsr);
            /* Form the residuals */
            formResiduals(paramsPtr->psrPtr, paramsPtr->npsr,1);    
            /* Do the fitting */
            if (0 == i) 
            {
                t2Fit(paramsPtr->psrPtr, paramsPtr->npsr, 0);   
            }
            /* Display the output */
            else 
            {
                // textOutput(psrPtr, npsr, globalParameter,0,0,0,"");  
            }
        }

        calcRMS(paramsPtr->psrPtr, paramsPtr->npsr);
    }

    return status;
}


long double PPOS_calculateError(PPOS_FitParams_t *paramsPtr)
{
    long double residErr = 0;

    if (PPOS_FitCriterion_RMS == paramsPtr->fitCriteria)
    {
        residErr = PPOS_calculateRMSResiduals(paramsPtr->psrPtr, paramsPtr->npsr);
    }
    else if (PPOS_FitCriterion_CHISQUARED == paramsPtr->fitCriteria) 
    {
        residErr = PPOS_calculateChiSquared(paramsPtr->psrPtr, paramsPtr->npsr);
    }
    else
    {
        residErr = INFINITY;
    }
    return residErr;
}

/**
 * Description in header file
 */
long double PPOS_calculateRMSResiduals(pulsar *psrPtr, int npsr)
{
    long double residErr = 0;
    for (int psr = 0; psr < npsr; psr++)
    {
        for (int obs = 0; obs < psrPtr[psr].nobs; obs++)
        {
            residErr += pow(psrPtr[psr].obsn[obs].residual / (long double)(psrPtr[psr].obsn[obs].toaErr * 1e-6), 2); 
        }
    }
    return sqrt(residErr);
}

/**
 * Description in header file
 */
long double PPOS_calculateChiSquared(pulsar *psrPtr, int npsr)
{
    long double chiSq = 0;

    for (int psr = 0; psr < npsr; psr++)
    {
        chiSq += psrPtr[psr].fitChisq/(double)psrPtr[psr].fitNfree;
    }
    return chiSq / npsr;
}

/**
 * Description in header file
 */
PPOS_Status_t PPOS_convertLatLongToGeocentric(PPOS_ElevationModel_t *elevModelPtr, PJ_LPZT geodet, PJ_XYZT *geocenPtr)
{
    PPOS_Status_t status = PPOS_SUCCESS;
    int16_t xPixel;
    int16_t yPixel;
    PPOS_GTransform_t gTransform = elevModelPtr->gTransform;
    float *pafScanline = (float *) CPLMalloc(sizeof(float) * elevModelPtr->xSize);

    status = PPOS_convertLatLongToPixelIndex(elevModelPtr, geodet, &xPixel, &yPixel);

    /** Load elevation data */
    if (PPOS_SUCCESS == status)
    {
        status = PPOS_loadCoordinates(elevModelPtr, pafScanline, yPixel);
    }

    /** Transform coordinates */
    if (PPOS_SUCCESS == status)
    {
        geodet.z = pafScanline[xPixel]; /* Height in meters */
        status = PPOS_convertGeodeticToGeocentric(elevModelPtr->P, geodet, geocenPtr);
    }

    CPLFree(pafScanline);

    return status;
}

/**
 * Description in header file
 */
PPOS_Status_t PPOS_convertPixelIndexToLatLong(PPOS_ElevationModel_t *elevModelPtr, int16_t xPixel, int16_t yPixel, PJ_LPZT *geodet)
{
    PPOS_Status_t status = PPOS_SUCCESS;
    PPOS_GTransform_t gTransform = elevModelPtr->gTransform;

    geodet->lam = gTransform.xCoord + xPixel * gTransform.wePixelRes + yPixel * gTransform.rowRot;
    geodet->phi = gTransform.yCoord + xPixel * gTransform.colRot + yPixel * gTransform.nsPixelRes;

    return status;
}

/**
 * Description in header file
 */
PPOS_Status_t PPOS_convertLatLongToPixelIndex(PPOS_ElevationModel_t *elevModelPtr, PJ_LPZT geodet, int16_t *xPixelPtr, int16_t *yPixelPtr)
{
    PPOS_Status_t status = PPOS_SUCCESS;
    PPOS_GTransform_t gTransform = elevModelPtr->gTransform;

    double lat = geodet.phi;  /* Latitude in degrees */
    double lon = geodet.lam;  /* Longitude in degrees */

    /** Convert LAT/LONG to pixel coordinates */
    int16_t xPixel = (int)(((lon - gTransform.xCoord) / gTransform.wePixelRes));
    if (xPixel >= elevModelPtr->xSize || xPixel < 0)
    {
        status = PPOS_WARNING_OUTOFBOUNDS;
    }

    int16_t yPixel = (int)(((lat - gTransform.yCoord) / gTransform.nsPixelRes));
    if (yPixel >= elevModelPtr->ySize || yPixel < 0)
    {
        status = PPOS_WARNING_OUTOFBOUNDS;
    }

    *xPixelPtr = xPixel;
    *yPixelPtr = yPixel;
    return status;
}

/**
 * Description in header file
 */
PPOS_Status_t PPOS_convertGeodeticToGeocentric(PJ* P, PJ_LPZT geodet, PJ_XYZT *geocent)
{
    PPOS_Status_t status = PPOS_SUCCESS;
    
    /** Define transform coordinates */
    PJ_COORD tgeodet;
    PJ_COORD tgeocent;
    tgeodet.lpzt = geodet;
    /** @todo add in current time */
    tgeodet.lpzt.t = HUGE_VAL;    /* Time (important only for time-dependent projections) */

    /** Transform to geocentric */
    tgeocent = proj_trans(P, PJ_FWD, tgeodet);

    *geocent = tgeocent.xyzt;
    return status;
}

/**
 * Description in header file
 */
PPOS_Status_t PPOS_setObserverLocation(PJ_XYZT *coordPtr)
{
    PPOS_Status_t status = PPOS_SUCCESS;

    if (NULL == coordPtr || NULL == observerLocPtr)
    {
        status = PPOS_INVALID_PTR;
    }

    if (PPOS_SUCCESS == status)
    {
        observerLocPtr->x = coordPtr->x;
        observerLocPtr->y = coordPtr->y;
        observerLocPtr->z = coordPtr->z;
    }
    return status;
}

/**
 * Description in header file
 */
PPOS_Status_t PPOS_loadElevationModel(PPOS_ElevationModel_t *elevModelPtr)
{
    PPOS_Status_t status = PPOS_SUCCESS;
    CPLErr result;
    char absPath[MAX_FILENAME_LEN];

    /** Create coordinate transformation object */
    elevModelPtr->P = proj_create_crs_to_crs(PJ_DEFAULT_CTX,
        "+proj=latlong +datum=WGS84",
        "+proj=geocent +datum=WGS84", NULL);

    if (NULL == elevModelPtr->P)
    {
        status = PPOS_ERROR_PROJ;
    }

    /** Load ALOS GTIFF */
    if (PPOS_SUCCESS == status)
    {
        sprintf(absPath, "%s/%s", elevModelPtr->filePath, elevModelPtr->fileName);
        elevModelPtr->hDataset = GDALOpen(absPath, GA_ReadOnly);
        if (NULL == elevModelPtr->hDataset) 
        {
            status = PPOS_ERROR_GDAL;
            const char *msg = CPLGetLastErrorMsg();
            printf("[%d:%s] ERROR: Unable to open GDAL file %s\n", __LINE__, __func__, absPath);
            printf("[%d:%s] ERROR: %s\n", __LINE__, __func__, msg);
        }
    }

    /** Get pixel tramsformation parameters */
    if (PPOS_SUCCESS == status)
    {
        double gTransform[6];
        result = GDALGetGeoTransform(elevModelPtr->hDataset, &gTransform[0]);
        if (CE_None != result && CE_Debug != result) 
        {
            status = PPOS_ERROR_GDAL;
            const char *msg = CPLGetLastErrorMsg();
            printf("[%d:%s] ERROR: Unable to get GeoTransform\n", __LINE__, __func__);
            printf("[%d:%s] ERROR: %s\n", __LINE__, __func__, msg);
        }
        else
        {
            elevModelPtr->gTransform = {gTransform[0], gTransform[1], gTransform[2], 
                gTransform[3], gTransform[4], gTransform[5]};   
        }
    }

    /** Get Elevation raster band  */
    if (PPOS_SUCCESS == status)
    {
        elevModelPtr->hBand = GDALGetRasterBand(elevModelPtr->hDataset, 1);
        if (NULL == elevModelPtr->hBand) 
        {
            status = PPOS_ERROR_GDAL;
            const char *msg = CPLGetLastErrorMsg();
            printf("[%d:%s] ERROR: Unable to get GDAL raster band\n", __LINE__, __func__);
            printf("[%d:%s] ERROR: %s\n", __LINE__, __func__, msg);
        }
    }

    /** Get grid size */
    if (PPOS_SUCCESS == status)
    {
        elevModelPtr->xSize = GDALGetRasterBandXSize(elevModelPtr->hBand);
        elevModelPtr->ySize = GDALGetRasterBandYSize(elevModelPtr->hBand);

        if (elevModelPtr->xSize > MAX_GRID_SIZE || elevModelPtr->ySize > MAX_GRID_SIZE)
        {
            status = PPOS_INVALID_GRID_SIZE;
        }
    }
    return status;
}

/**
 * Description in header file
 */
PPOS_Status_t PPOS_loadCoordinates(PPOS_ElevationModel_t *elevModelPtr, float *pafScanline, int yLine)
{
    PPOS_Status_t status = PPOS_SUCCESS;
    CPLErr result;
    
    /** load new raster latitude */
    result = GDALRasterIO(elevModelPtr->hBand, GF_Read, 0, yLine, elevModelPtr->xSize, 1,
        pafScanline, elevModelPtr->xSize, 1, GDT_Float32, 0, 0 );
    if (CE_None != result && CE_Debug != result) 
    {
        status = PPOS_ERROR_GDAL;
        const char *msg = CPLGetLastErrorMsg();
        printf("[%d:%s] ERROR: Unable to read GDAL raster yPixel = %d \n", __LINE__, __func__, yLine);
        printf("[%d:%s] ERROR: %s\n", __LINE__, __func__, msg);
    }

    return status;
}

/**
 * Description in header file
 */
PPOS_Status_t PPOS_loadGrid(PPOS_ElevationModel_t *elevModelPtr, PJ_LPZT geodet)
{
    PPOS_Status_t status = PPOS_SUCCESS;
    GDALDatasetH oldhDataset = elevModelPtr->hDataset;
    char fname[MAX_FILENAME_LEN];
    char nS;
    char eW;

    if (geodet.phi >= 0 && geodet.phi <= 90)
    {
        nS = 'N';
    }
    else if (geodet.phi < 0 && geodet.phi >= -90)
    {
        nS = 'S';
    }
    else
    {
        status = PPOS_ERROR_INVALID_POSITION;
    }

    if (geodet.lam >= 0 && geodet.lam <= 180)
    {
        eW = 'E';
    }
    else if (geodet.lam < 0 && geodet.lam > -180)
    {
        eW = 'W';
    }
    else
    {
        status = PPOS_ERROR_INVALID_POSITION;
    }

    if (PPOS_SUCCESS == status)
    {
        /** ALOS naming is based on Geographic coordinates at lower left corner of lower left pixel */
        sprintf(fname, "ALPSMLC30_%c%03d%c%03d_DSM.tif", nS, (int)floor(abs(geodet.phi - 1)), eW, (int)floor(abs(geodet.lam)));
        strcpy(elevModelPtr->fileName, fname);

        elevModelPtr->hDataset = NULL;
        status = PPOS_loadElevationModel(elevModelPtr);
        if (PPOS_SUCCESS != status)
        {
            printf("[%d:%s] WARNING: failed to load elevation model %s RC = %d\n", __LINE__, __func__, fname, status);
        }
    }
    
    if (PPOS_SUCCESS == status)
    {
        elevModelPtr->initialised = true;
        
        if (NULL != oldhDataset)
        {
            GDALClose(oldhDataset);
        }
    }

    return status;
}

/**
 * Description in header file
 */
PPOS_Status_t PPOS_loadAdjacentGrid(PPOS_ElevationModel_t *elevModelPtr, int16_t xPixel, int16_t yPixel)
{
    PPOS_Status_t status = PPOS_SUCCESS;
    PJ_LPZT geodet;

    geodet.lam = elevModelPtr->gTransform.xCoord; /** Longitude in degrees */
    geodet.phi = elevModelPtr->gTransform.yCoord; /** Latitude in degrees */

    /** Pixel exceeds east boundary */
    if (xPixel < 0)
    {
        geodet.lam -= floor(1 + abs(xPixel) / (int)elevModelPtr->xSize);
    }
    /** Pixel exceeds west boundary */
    else if (xPixel >= elevModelPtr->xSize)
    {
        geodet.lam += floor(xPixel / (int)elevModelPtr->xSize);
    }

    /** Pixel exceeds north boundary */
    if (yPixel < 0)
    {
        geodet.phi += floor(1 + abs(yPixel) / (int)elevModelPtr->ySize);
    } 
    /** Pixel exceeds south boundary */
    else if (yPixel >= elevModelPtr->ySize)
    {
        geodet.phi -=floor(yPixel / (int)elevModelPtr->ySize);
    }

    status = PPOS_loadGrid(elevModelPtr, geodet);

    return status;
}

/**
 * Description in header file
 */
PPOS_Status_t PPOS_estimateSolutions(PPOS_FitParams_t *paramsPtr)
{
    PPOS_Status_t status = PPOS_SUCCESS;

    /** Normal vector in direction from centre of Earth to observer (this is zenith in NEU coordinate system)*/
    double eNorm[3] = {cos(M_PI)*cos(0), sin(M_PI)*sin(0), sin(0)}; 
    double rho[paramsPtr->npsr];    /** Distance from centre of the Earth to pulsar plane */
    double radI[paramsPtr->npsr];   /** Radius of intersection with the Earth */
    double pNorm[paramsPtr->npsr][3];

    double centreDiffNorm[paramsPtr->npsr][paramsPtr->npsr];    /** Difference between  */

    pulsar *psrPtr = paramsPtr->psrPtr;

    int16_t nIntersect = 0;
    for (int pj = 0; pj < paramsPtr->npsr; pj++)
    {
        /** Calculate circle of intersection characteristics for pulsar */
        rho[pj] = -EARTH_RADIUS * dotproduct(eNorm, psrPtr[pj].posPulsar);
        radI[pj] = EARTH_RADIUS * sqrt(1 - pow(dotproduct(eNorm, psrPtr[pj].posPulsar), 2));
        pNorm[pj][0] = psrPtr[pj].posPulsar[0];
        pNorm[pj][1] = psrPtr[pj].posPulsar[1];
        pNorm[pj][2] = psrPtr[pj].posPulsar[2];

        pulsar pB = psrPtr[pj];
        for (int pi = 0; pi < paramsPtr->npsr; pi++)
        {
            pulsar pA = psrPtr[pi];

            if (pj >= pi)
            {
                continue;
            }
            double eA = dotproduct(eNorm, pA.posPulsar);
            double aNorm[3] = {eA * pNorm[pj][0], eA * pNorm[pj][1], eA * pNorm[pj][2]};
            double eB = dotproduct(eNorm, pA.posPulsar);
            double bNorm[3] = {eB * pNorm[pi][0], eB * pNorm[pi][1], eB * pNorm[pi][2]};

            double cDiff[3] = {bNorm[0] - aNorm[0], bNorm[1] - aNorm[1], bNorm[2] - aNorm[2]};

            double diffNorm = EARTH_RADIUS * sqrt(pow(cDiff[0], 2) + pow(cDiff[1], 2) + pow(cDiff[2], 2));

            centreDiffNorm[pj][pi] = diffNorm;
            if (diffNorm < (radI[pj] + radI[pj]))
            {
                nIntersect += 2;
            }
            else if (diffNorm == (radI[pj] + radI[pj]))
            {
                nIntersect++;
            }
            else
            {
                /** No intersection */
            }
        }
    }
    return status;
}

/**
 * Description in header file
 */
PPOS_Status_t PPOS_openResults(char *fileName)
{
    PPOS_Status_t status = PPOS_SUCCESS;

    /** Display Statistucs */
    FILE *fpt;
    fpt = fopen(fileName, "w+");

    if (NULL == fpt)
    {
        printf("[%d:%s] ERROR: Unable to open results file %s\n", __LINE__, __func__, fileName);
        status = PPOS_ERROR_FOPEN;
    }

    /** Close file */
    fclose(fpt);

    return status;
}

/**
 * Description in header file
 */
PPOS_Status_t PPOS_writeResults(PPOS_FitParams_t *paramsPtr, char *fileName, long double residErrMap[MAX_GRID_SIZE])
{
    PPOS_Status_t status = PPOS_SUCCESS;

    /** Display Statistucs */
    FILE *fpt;
    fpt = fopen(fileName, "a+");

    if (NULL == fpt)
    {
        printf("[%d:%s] ERROR: Unable to open results file %s\n", __LINE__, __func__, fileName);
        status = PPOS_ERROR_FOPEN;
    }

    /** Read results map into file */
    if (PPOS_SUCCESS == status)
    {
        for (int xi = 0; xi < paramsPtr->elevModelPtr->xSize; xi++)
        {
            fprintf(fpt, "%Lf,", residErrMap[xi]);
        }
        fprintf(fpt, "\n");
    }

    /** Close file */
    fclose(fpt);
    return status;
}

/**
 * Description in header file
 */
void help()
{
  /* This function should contain usage information about the plugin which should (in general) be accessed */
  /* by the user pressing 'h'                                                                        */
  printf("\n\n");
  printf("This plugin determines the positioning of a pulsar observer \n\n");
  printf("usage:\n");
  printf("\t-f   Define the absolute paths to the .par and .tim files [required]: -f {/path/to/.par} {/path/to/.tim}\n");
  printf("\t-a   Define the algorithm to use for fitting, [VALUES = (grid, grde)] [required]\n");
  printf("\t-i   Define the maximum number of fit iterations for TEMPO2 [DEFAULT = 2]\n");
  printf("\t-e   Define the method to check the fit error [VALUES = (rms, chi)] [DEFAULT = rms] \n");
  printf("\t-o   Define the observer code\n");
  printf("\t-l   Define the initial geodetic position: -l {long (deg)} {lat (deg)}\n");
  printf("\t-n   Fill NODATA errors with 0 elevation\n");
}
