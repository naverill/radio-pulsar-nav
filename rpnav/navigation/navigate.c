/*    Tempo2 Plugin to localise an observer in 3D space 
 *
 *
 *
 *    Copyright (C) 2006,2007,2008,2009, George Hobbs, Russel Edwards
 *    This file is part of TEMPO2.
 *
 *    TEMPO2 is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *    TEMPO2 is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *    You should have received a copy of the GNU General Public License
 *    along with TEMPO2.  If not, see <http://www.gnu.org/licenses/>.
 */

/*
 *    If you use TEMPO2 then please acknowledge it by citing
 *    Hobbs, Edwards & Manchester (2006) MNRAS, Vol 369, Issue 2,
 *    pp. 655-672 (bibtex: 2006MNRAS.369..655H)
 *    or Edwards, Hobbs & Manchester (2006) MNRAS, VOl 372, Issue 4,
 *    pp. 1549-1574 (bibtex: 2006MNRAS.372.1549E) when discussing the
 *    timing model.
 */
#include <stdio.h>
#include <math.h>
#include <tempo2.h>


void help() /* Display help */
{
    printf("-----------------------------------------\n");
    printf("-h\t\t help menu\n");
    printf("-t\t\t path to .tim file\n");
    printf("-p\t\t path to .par file\n");
    printf("-x\t\t Initial estimate for observer in x axis (WGS84)\n");
    printf("-y\t\t Initial estimate for observer in y axis (WGS84)\n");
    printf("-z\t\t Initial estimate for observer in z axis (WGS84)\n");
    printf("-xu\t\t Uncertainty for observer in x axis\n");
    printf("-yu\t\t Uncertainty for observer in y axis\n");
    printf("-zu\t\t Uncertainty for observer in z axis\n");
    //printf("-c\t\t Reference coordinate system for position estimate. Valid options\n");
    //printf("are WGS84");
    printf("-----------------------------------------\n");
}

extern "C" int tempoOutput(int argc,char *argv[],pulsar *psr,int npsr){
    char parFile[MAX_PSR][MAX_FILELEN];
    char timFile[MAX_PSR][MAX_FILELEN];
    char resFile[MAX_FILELEN];
    float x, y, z, xu, yu, zu;

    printf("Radio Pulsar Navigation\n");
    printf("Author:              N. Averill\n");
    printf("Version:             1.0\n");

    int i;
    if (arc < 2){
        printf("Invalid number of arguments: input -h for help\n");
    }
    for (i=2; i<argc, i++){
        if !(strcmp(argv[i], "-h")){
            help()
        }
        else if !(strcmp(argv[i], "-p")){
            strcpy(parFile[*npsr],argv[i+1]);
            (*npsr)++;
        }
        else if !(strcmp(argv[i], "-t")){
            strcpy(timFile[*npsr],argv[i+2]);
        }
        else if !(strcmp(argv[i], "-x")){
            strcpy(timFile[*npsr],argv[i+2]);
        }
        else if !(strcmp(argv[i], "-y")){
            strcpy(timFile[*npsr],argv[i+2]);
        }
        else if !(strcmp(argv[i], "-z")){
            strcpy(timFile[*npsr],argv[i+2]);
        }
        else if !(strcmp(argv[i], "-xu")){
            strcpy(timFile[*npsr],argv[i+2]);
        }
        else if !(strcmp(argv[i], "-yu")){
            strcpy(timFile[*npsr],argv[i+2]);
        }
        else if !(strcmp(argv[i], "-zu")){
            strcpy(timFile[*npsr],argv[i+2]);
        }
    }

    readParfile(psr,parFile,timFile,*npsr);
    readTimfile(psr,timFile,*npsr);

}



