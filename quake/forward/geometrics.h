/* -*- C -*- */

/* @copyright_notice_start
 *
 * This file is part of the CMU Hercules ground motion simulator developed
 * by the CMU Quake project.
 *
 * Copyright (C) Carnegie Mellon University. All rights reserved.
 *
 * This program is covered by the terms described in the 'LICENSE.txt' file
 * included with this software package.
 *
 * This program comes WITHOUT ANY WARRANTY; without even the implied warranty
 * of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * 'LICENSE.txt' file for more details.
 *
 *  @copyright_notice_end
 */
#ifndef GEOMETRICS_H
#define GEOMETRICS_H

#ifdef PROJ
	#include <proj.h>
#endif
#include "psolve.h"


#define PI 3.14159265358979323846264338327

/*
 * This struct was also in psolve.h, and was necessary to make it visible
 * to the nonlinear module
 *
 * typedef struct vector3D_t{
 *
 *     double x [ 3 ];
 *
 * }vector3D_t;
 */



vector3D_t compute_global_coords( vector3D_t origin, vector3D_t local,
				  double dip, double rake ,double strike );


vector3D_t compute_centroid( vector3D_t* p );


double compute_area(vector3D_t* p);


int compute_1D_grid( double cellSize, int numberOfCells, int pointsInCell,
		     double minimumEdgeTriangle, double* grid );


vector3D_t compute_domain_coords( vector3D_t point, double azimuth );

// Define a structure to hold the UTM coordinates
typedef struct {
    double e; // Easting
    double n; // Northing
} UTMCoordinates_t;

typedef struct {
    int zone;
    char hemisphere;
	UTMCoordinates_t refPoint;
	#ifdef PROJ
		PJ *P;
	#endif
} UTMZone_t;

vector3D_t
compute_domain_coords_linearinterp(
	double lon , double lat, double* longcorner, double *latcorner,
	double domainlengthetha, double domainlengthcsi, UTMZone_t* utmZone);

UTMZone_t getUTMZone(double latitude, double longitude);
#ifdef PROJ
	UTMCoordinates_t convertLatLonToUTM(double latitude, double longitude, PJ *P);
	void convertUTMToLatLon(UTMCoordinates_t utmCoordinates, PJ *P, double *latitude, double *longitude);
#endif

#endif /* GEOMETRICS_H */
