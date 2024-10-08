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

#include <inttypes.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "geometrics.h"

/*
 * Transforms a point from the local coordinates of the fault to the
 * global coordinates.
 *
 */
vector3D_t
compute_global_coords( vector3D_t origin, vector3D_t local, double dip,
		       double rake, double strike )
{
    double x, y, z;
    double d, l, p;
    int i;

    vector3D_t xGlobal;

    d = dip    * PI / 180;
    l = rake   * PI / 180;
    p = strike * PI / 180;

    x = local.x[ 0 ];
    y = local.x[ 1 ];
    z = local.x[ 2 ];

    xGlobal.x[ 0 ] =
	( cos( p ) * cos( l ) + sin( p ) * cos( d ) * sin( l )) * x -
	(-cos( p ) * sin( l ) + sin( p ) * cos( d ) * cos( l )) * y -
	(-sin( p ) * sin( d ) ) * z;

    xGlobal.x[ 1 ] =
	( sin( p ) * cos( l ) - cos( p ) * cos( d ) * sin( l )) * x -
	(-sin( p ) * sin( l ) - cos( p ) * cos( d ) * cos( l )) * y -
	( cos( p ) * sin( d ) ) * z;

    xGlobal.x[ 2 ] =
	- sin( d ) * sin( l ) * x
	+ sin( d ) * cos( l ) * y 
	+ cos( d ) * z;

    for (i = 0; i < 3; i++) {
	xGlobal.x[ i ] = xGlobal.x[ i ] + origin.x[ i ];
    }

    return xGlobal;
}


vector3D_t compute_centroid( vector3D_t* p )
{
    vector3D_t centroid;
    int i, j;

    for (i = 0; i < 3; i++) {
	centroid.x[ i ] = 0;

	for (j = 0; j < 3; j++) {
	    centroid.x[ i ] += p[ j ].x[ i ];
	}

	centroid.x[ i ] = centroid.x[ i ] / 3;
    }

    return centroid;
}


double compute_area( vector3D_t* p )
{
    double area;
    double x1, x2, x3, y1, y2, y3;

    x1 = p[ 0 ].x[ 0 ];
    y1 = p[ 0 ].x[ 1 ];
    x2 = p[ 1 ].x[ 0 ];
    y2 = p[ 1 ].x[ 1 ];
    x3 = p[ 2 ].x[ 0 ];
    y3 = p[ 2 ].x[ 1 ];

    area = .5 * ( x1 * ( y2 - y3 ) - x2 * ( y1 - y3 ) + x3 * ( y1 - y2 ) );

    return area;
}


int compute_1D_grid( double cellSize, int numberOfCells, int pointsInCell,
		     double minimumEdgeTriangle, double* grid )
{
    int i, j, k;

    k = 0;

    for (i = 0; i < numberOfCells; i++) {
	for (j = 0; j < pointsInCell; j++) {
	    grid [ k ] = i * cellSize + j * minimumEdgeTriangle;
	    k = k + 1;
	}
    }

    grid [ k ] = i * cellSize;
    k = k + 1;

    return k;
}


/**
 * Rotates the x domain and y domain components of the given point if the
 * rectangular prism is not alinged with the north.
 *
 * \param azimuth in degrees
 */
vector3D_t compute_domain_coords( vector3D_t point, double azimuth )
{
    int	       iRow, iCol;
    double     transformation[3][3];
    vector3D_t pointTransformed;

    for (iRow = 0; iRow < 3; iRow++) {
	for (iCol = 0; iCol < 3; iCol++) {
	    transformation[iRow][iCol] = 0;
	}
    }

    azimuth = PI * azimuth / 180;

    transformation[ 0 ][ 0 ] =  cos( azimuth );
    transformation[ 1 ][ 0 ] = -sin( azimuth );
    transformation[ 0 ][ 1 ] =  sin( azimuth );
    transformation[ 1 ][ 1 ] =  cos( azimuth );
    transformation[ 2 ][ 2 ] =  1;

    pointTransformed.x[0] = 0;
    pointTransformed.x[1] = 0;
    pointTransformed.x[2] = 0;

    for (iRow = 0; iRow < 3; iRow++) {
	for (iCol = 0; iCol < 3; iCol++) {
	    pointTransformed.x[iRow]
		+= transformation[iRow][iCol] * point.x[iCol];
	}
    }

    return pointTransformed;
}

// NOTE: For both approaches, the function take the same input arguments to 
// maintain compatibility. As a result, some of the arguments are not going to 
// be used in either case.
#ifdef PROJ

    UTMZone_t getUTMZone(double latitude, double longitude) {
        UTMZone_t utmZone;
        utmZone.zone = floor((longitude + 180) / 6) + 1;
        utmZone.hemisphere = (latitude >= 0) ? 'N' : 'S';
        // Save PROJ context and PJ for later use
        PJ_CONTEXT *C;
        PJ *P;
        char projString[64]; // Buffer to hold the PROJ string
        C = proj_context_create();
        // Construct the PROJ string dynamically based on zone and hemisphere
        snprintf(projString, sizeof(projString), "+proj=utm +zone=%d +%c +ellps=WGS84 +datum=WGS84 +units=m +no_defs", utmZone.zone, utmZone.hemisphere == 'N' ? 'n' : 's');
        // Create a PJ for the transformation using the constructed PROJ string
        P = proj_create(C, projString);
        if (P == NULL) {
            printf("Projection creation failed\n");
            // Terminate the program
            exit(1);
        }
        utmZone.P = P;
        return utmZone;
    }

    // Function to convert latitude and longitude to UTM coordinates
    UTMCoordinates_t convertLatLonToUTM(double latitude, double longitude, PJ *P) {
        if (latitude < -80.0 || latitude > 84.0) {
            printf("Error: Invalid latitude value. Latitude must be between -80 and 84 degrees.\n");
            return (UTMCoordinates_t){-1, -1}; // Return error coordinates
        }
        PJ_COORD p;
        // Assign latitude and longitude to the PJ_COORD
        // NOTE: the third argument in proj_coord() is the elevation, while the fourth 
        // argument is the time. Since we don't need them, we set them to 0.
        p = proj_coord(proj_torad(longitude), proj_torad(latitude), 0, 0);
        // Transform the coordinates
        p = proj_trans(P, PJ_FWD, p);
        // Return the transformed coordinates
        return (UTMCoordinates_t){p.enu.e, p.enu.n};
    }

    void convertUTMToLatLon(UTMCoordinates_t utmCoordinates, PJ *P, double *latitude, double *longitude) {
        PJ_COORD p;
        p = proj_coord(utmCoordinates.e, utmCoordinates.n, 0, 0);
        p = proj_trans(P, PJ_INV, p);
        *longitude = proj_todeg(p.lp.lam);
        *latitude = proj_todeg(p.lp.phi);
    }

    /**
     * Computes the domain coordinates of a point based on UTM projection.
     */
    vector3D_t compute_domain_coords_linearinterp(double lon, double lat,
        double* longcorner, double* latcorner, 
        double domainlengthetha, double domainlengthcsi, 
        UTMZone_t* utmZone)
    {
        // NOTE: longcorner, latcorner, domainlengthetha, domainlengthcsi are not used in this implementation
        vector3D_t domainCoords;
        // NOTE: utmZone->refPoint is the origin in UTM coordinates
        UTMCoordinates_t utmCoordinates = convertLatLonToUTM(lat, lon, utmZone->P);
        if ( (utmZone->refPoint.e == -1 && utmZone->refPoint.n == -1) || (utmCoordinates.e == -1 && utmCoordinates.n == -1) ) {
            printf("Conversion from latitude/longitude to UTM failed\n");
            return (vector3D_t){-1, -1, -1}; // Return error coordinates
        }
        domainCoords.x[0] = utmCoordinates.n - utmZone->refPoint.n;
        domainCoords.x[1] = utmCoordinates.e - utmZone->refPoint.e;
        domainCoords.x[2] = 0;
        return domainCoords;
    }

#else

/**
 * Computes the domain coordinates of a point based on corners of the global
 * coordinates (long lat depth).
 */
vector3D_t
compute_domain_coords_linearinterp( double  lon , double lat,
    double* longcorner, double* latcorner,
    double  domainlengthetha, double  domainlengthcsi,
    UTMZone_t* utmZone)
{
    // NOTE: utmZone is not used in this implementation
    int i;
    double Ax,Ay,Bx,By,Cx,Cy,Dx,Dy;
    double res;
    double Xi[4],Yi[4],XN[2];
    double M[2][2],F[2],DXN[2];
    double X,Y;
    vector3D_t domainCoords;

    X=lat;
    Y=lon;

    for (i = 0; i < 4; i++) {
	Xi[i] = latcorner[i];
	Yi[i] = longcorner[i];
    }

    Ax = 4*X - (Xi[0] + Xi[1] + Xi[2] + Xi[3]);
    Ay = 4*Y - (Yi[0] + Yi[1] + Yi[2] + Yi[3]);

    Bx = -Xi[0] + Xi[1] + Xi[2] - Xi[3];
    By = -Yi[0] + Yi[1] + Yi[2] - Yi[3];

    Cx = -Xi[0] - Xi[1] + Xi[2] + Xi[3];
    Cy = -Yi[0] - Yi[1] + Yi[2] + Yi[3];

    Dx = Xi[0] - Xi[1] + Xi[2] - Xi[3];
    Dy = Yi[0] - Yi[1] + Yi[2] - Yi[3];

    /* initial values for csi and etha*/
    XN[0] = 0;
    XN[1] = 0;

    res = 1e10;

    while (res > 1e-6) {
	M[0][0] = Bx + Dx * XN[1];
	M[0][1] = Cx + Dx * XN[0];
	M[1][0] = By + Dy * XN[1];
	M[1][1] = Cy + Dy * XN[0];

	F[0] = -Ax + Bx * XN[0] + Cx * XN[1] + Dx * XN[0] * XN[1];
	F[1] = -Ay + By * XN[0] + Cy * XN[1] + Dy * XN[0] * XN[1];

	DXN[0] = -(F[0]*M[1][1] - F[1]*M[0][1])
	     / (M[0][0]*M[1][1] - M[1][0]*M[0][1]);

	DXN[1] = -(F[1]*M[0][0] - F[0]*M[1][0])
	     / (M[0][0]*M[1][1] - M[1][0]*M[0][1]);

	res = fabs( F[0] ) + fabs( F[1] );

	XN[0] = XN[0] + DXN[0];
	XN[1] = XN[1] + DXN[1];
    }

    domainCoords.x[0] = .5 * (XN[0] + 1) * domainlengthcsi;
    domainCoords.x[1] = .5 * (XN[1] + 1) * domainlengthetha;
    domainCoords.x[2] = 0;

    return domainCoords;
}

#endif