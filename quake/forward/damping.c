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
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "psolve.h"
#include "quake_util.h"
#include "damping.h"
#include "stiffness.h"

//#include "quake_util.h"
#include "util.h"
#include "timers.h"
#include "cvm.h"

#include "nonlinear.h"
#include "topography.h"


/* -------------------------------------------------------------------------- */
/*                             Global Variables                               */
/* -------------------------------------------------------------------------- */

static int32_t  myLinearElementsCount;
static int32_t *myLinearElementsMapper;


/* -------------------------------------------------------------------------- */
/*                      Initialization of parameters for
 *                   Nonlinear and Topography compatibility                   */
/* -------------------------------------------------------------------------- */

/**
 * Counts the number of nonlinear and topo elements in my local mesh
 */
void trad_elements_count(int32_t myID, mesh_t *myMesh) {

    int32_t eindex;
    int32_t count = 0;

    for (eindex = 0; eindex < myMesh->lenum; eindex++) {

/*        if ( ( isThisElementNonLinear(myMesh, eindex) == NO ) &&
             ( BelongstoTopography(myMesh, eindex)    == NO )  ) {
            count++;
        }*/
        if ( ( isThisElementNonLinear(myMesh, eindex) == NO ) &&
             ( IsDampingElement(myMesh, eindex)    == YES )  ) {
            count++;
        }

    }

    if ( count > myMesh-> lenum ) {
        fprintf(stderr,"Thread %d: damping_linear_elements_count: "
                "more elements than expected\n", myID);
        MPI_Abort(MPI_COMM_WORLD, ERROR);
        exit(1);
    }

    myLinearElementsCount = count;

    return;
}

/**
 * Re-counts and stores the nonlinear and topo element's indices to a static local array
 * that will serve as mapping tool to the local mesh elements table.
 */
void trad_elements_mapping(int32_t myID, mesh_t *myMesh) {

    int32_t eindex;
    int32_t count = 0;

    XMALLOC_VAR_N(myLinearElementsMapper, int32_t, myLinearElementsCount);

    for (eindex = 0; eindex < myMesh->lenum; eindex++) {

        if ( ( isThisElementNonLinear(myMesh, eindex) == NO ) &&
             ( IsDampingElement(myMesh, eindex)    == YES ) ) {
            myLinearElementsMapper[count] = eindex;
            count++;
        }


/*        if ( ( isThisElementNonLinear(myMesh, eindex) == NO ) &&
             ( BelongstoTopography(myMesh, eindex)    == NO ) ) {
            myLinearElementsMapper[count] = eindex;
            count++;
        }*/

    }

    if ( count != myLinearElementsCount ) {
        fprintf(stderr,"Thread %d: linear_elements_mapping: "
                "more elements than the count\n", myID);
        MPI_Abort(MPI_COMM_WORLD, ERROR);
        exit(1);
    }

    return;
}

void damp_init(int32_t myID, mesh_t *myMesh) {

    trad_elements_count(myID, myMesh);
    trad_elements_mapping(myID, myMesh);
}

/* -------------------------------------------------------------------------- */
/*                       Damping Contribution Methods                         */
/* -------------------------------------------------------------------------- */

void damping_addforce(mesh_t *myMesh, mysolver_t *mySolver, fmatrix_t (*theK1)[8], fmatrix_t (*theK2)[8]){

    fvector_t localForce[8];
    int       i,j;
    int32_t   eindex;
    int32_t   lin_eindex;

    fvector_t deltaDisp[8];

    /* loop on the number of elements */
    /* loop on the number of elements */
    for (lin_eindex = 0; lin_eindex < myLinearElementsCount; lin_eindex++)
    {
        elem_t *elemp;
        e_t    *ep;

        eindex = myLinearElementsMapper[lin_eindex];
        elemp = &myMesh->elemTable[eindex];
        ep = &mySolver->eTable[eindex];

        /* Step 1. calculate the force due to the element stiffness */
        memset(localForce, 0, 8 * sizeof(fvector_t));

        /* compute the diff between disp(tm1) and disp(tm2) */

        /* the_E1_timer -= MPI_Wtime();*/
        for (i = 0; i < 8; i++) {
            fvector_t *tm1Disp, *tm2Disp;
            int32_t    lnid;

            lnid = elemp->lnid[i];

            tm1Disp = mySolver->tm1 + lnid;
            tm2Disp = mySolver->tm2 + lnid;

            deltaDisp[i].f[0] = tm1Disp->f[0] - tm2Disp->f[0];
            deltaDisp[i].f[1] = tm1Disp->f[1] - tm2Disp->f[1];
            deltaDisp[i].f[2] = tm1Disp->f[2] - tm2Disp->f[2];
        }

        if(vector_is_zero( deltaDisp ) != 0) {
            /*the_E3_timer += MPI_Wtime();  */

            for (i = 0; i < 8; i++)
            {
                fvector_t* toForce;
                toForce = &localForce[i];

                for (j = 0; j < 8; j++)
                {
                    fvector_t *myDeltaDisp;

                    /* contribution by ( - b * deltaT * Ke * ( Ut - Ut-1 ) ) */
                    /* But if myDeltaDisp is zero avoids multiplications     */
                    myDeltaDisp = &deltaDisp[j];

                    MultAddMatVec(&theK1[i][j], myDeltaDisp, -ep->c3, toForce);
                    MultAddMatVec(&theK2[i][j], myDeltaDisp, -ep->c4, toForce);

                }
            }
        }
        for (i = 0; i < 8; i++) {
            int32_t lnid;
            fvector_t *nodalForce;

            lnid = elemp->lnid[i];

            nodalForce = mySolver->force + lnid;
            nodalForce->f[0] += localForce[i].f[0];
            nodalForce->f[1] += localForce[i].f[1];
            nodalForce->f[2] += localForce[i].f[2];
        }

    } /* for all the elements */

    return;
}

/**
 * Introduce BKT Model: Compute and add the force due to the element
 *              damping.
 */

void calc_conv_optimized(mesh_t *myMesh, mysolver_t *mySolver, double theFreq, double theDeltaT, double theDeltaTSquared, damping_type_t typeOfDamping){

    int32_t eindex;
    int i;

    double rmax = 2. * M_PI * theFreq * theDeltaT;
    int32_t   lin_eindex;

    double cdt;

    double g0, g02, cg0, eg0, g0k, g02k, cg0k, eg0k;
    double g1, g12, cg1, eg1, g1k, g12k, cg1k, eg1k;
    double g2, g22, cg2, eg2, g2k, g22k, cg2k, eg2k;

    if (typeOfDamping == BKT) {
        cdt = 2. * M_PI * theFreq * theDeltaT;
    } else {
        cdt = theDeltaT;
    }

    for (lin_eindex = 0; lin_eindex < myLinearElementsCount; lin_eindex++)
    {

        elem_t *elemp;
        edata_t *edata;

        eindex = myLinearElementsMapper[lin_eindex];
        elemp = &myMesh->elemTable[eindex];
        edata = (edata_t *)elemp->data;

        // SHEAR RELATED CONVOLUTION
        g0  = cdt * edata->g0_shear;
        g02 = g0 / 2.;
        cg0 = g02 * ( 1. - g0 );
        eg0 = exp( -g0 );

        g1  = cdt * edata->g1_shear;
        g12 = g1 / 2.;
        cg1 = g12 * ( 1. - g1 );
        eg1 = exp( -g1 );

        g0k  = cdt * edata->g0_kappa;
        g02k = g0k / 2.;
        cg0k = g02k * ( 1. - g0k );
        eg0k = exp( -g0k );

        g1k  = cdt * edata->g1_kappa;
        g12k = g1k / 2.;
        cg1k = g12k * ( 1. - g1k );
        eg1k = exp( -g1k );

        if (typeOfDamping >= BKT3) {
            g2  = cdt * edata->g2_shear;
            g22 = g2 / 2.;
            cg2 = g22 * ( 1. - g2 );
            eg2 = exp( -g2 );

            g2k  = cdt * edata->g2_kappa;
            g22k = g2k / 2.;
            cg2k = g22k * ( 1. - g2k );
            eg2k = exp( -g2k );
        }

        for(i = 0; i < 8; i++)
        {
            int32_t     lnid, cindex;
            fvector_t  *tm1Disp, *tm2Disp;
            fvector_t  *f0_tm1, *f1_tm1, *f2_tm1;

            lnid = elemp->lnid[i];

            /* cindex is the index of the node in the convolution vector */
            cindex = eindex * 8 + i;

            tm1Disp = mySolver->tm1 + lnid;
            tm2Disp = mySolver->tm2 + lnid;

            // SHEAR RELATED CONVOLUTION
            if ( (edata->g0_shear != 0) && (edata->g1_shear != 0) ) {
                f0_tm1 = mySolver->conv_shear_1 + cindex;
                f1_tm1 = mySolver->conv_shear_2 + cindex;

                f0_tm1->f[0] = cg0 * tm1Disp->f[0] + g02 * tm2Disp->f[0] + eg0 * f0_tm1->f[0];
                f0_tm1->f[1] = cg0 * tm1Disp->f[1] + g02 * tm2Disp->f[1] + eg0 * f0_tm1->f[1];
                f0_tm1->f[2] = cg0 * tm1Disp->f[2] + g02 * tm2Disp->f[2] + eg0 * f0_tm1->f[2];

                f1_tm1->f[0] = cg1 * tm1Disp->f[0] + g12 * tm2Disp->f[0] + eg1 * f1_tm1->f[0];
                f1_tm1->f[1] = cg1 * tm1Disp->f[1] + g12 * tm2Disp->f[1] + eg1 * f1_tm1->f[1];
                f1_tm1->f[2] = cg1 * tm1Disp->f[2] + g12 * tm2Disp->f[2] + eg1 * f1_tm1->f[2];

                if (typeOfDamping >= BKT3) {
                    f2_tm1 = mySolver->conv_shear_3 + cindex;
                    f2_tm1->f[0] = cg2 * tm1Disp->f[0] + g22 * tm2Disp->f[0] + eg2 * f2_tm1->f[0];
                    f2_tm1->f[1] = cg2 * tm1Disp->f[1] + g22 * tm2Disp->f[1] + eg2 * f2_tm1->f[1];
                    f2_tm1->f[2] = cg2 * tm1Disp->f[2] + g22 * tm2Disp->f[2] + eg2 * f2_tm1->f[2];
                }
            }

            // DILATATION RELATED CONVOLUTION
            if ( (edata->g0_kappa != 0) && (edata->g1_kappa != 0) ) {
                f0_tm1 = mySolver->conv_kappa_1 + cindex;
                f1_tm1 = mySolver->conv_kappa_2 + cindex;

                f0_tm1->f[0] = cg0k * tm1Disp->f[0] + g02k * tm2Disp->f[0] + eg0k * f0_tm1->f[0];
                f0_tm1->f[1] = cg0k * tm1Disp->f[1] + g02k * tm2Disp->f[1] + eg0k * f0_tm1->f[1];
                f0_tm1->f[2] = cg0k * tm1Disp->f[2] + g02k * tm2Disp->f[2] + eg0k * f0_tm1->f[2];

                f1_tm1->f[0] = cg1k * tm1Disp->f[0] + g12k * tm2Disp->f[0] + eg1k * f1_tm1->f[0];
                f1_tm1->f[1] = cg1k * tm1Disp->f[1] + g12k * tm2Disp->f[1] + eg1k * f1_tm1->f[1];
                f1_tm1->f[2] = cg1k * tm1Disp->f[2] + g12k * tm2Disp->f[2] + eg1k * f1_tm1->f[2];

                if (typeOfDamping >= BKT3) {
                    f2_tm1 = mySolver->conv_kappa_3 + cindex;
                    f2_tm1->f[0] = cg2k * tm1Disp->f[0] + g22k * tm2Disp->f[0] + eg2k * f2_tm1->f[0];
                    f2_tm1->f[1] = cg2k * tm1Disp->f[1] + g22k * tm2Disp->f[1] + eg2k * f2_tm1->f[1];
                    f2_tm1->f[2] = cg2k * tm1Disp->f[2] + g22k * tm2Disp->f[2] + eg2k * f2_tm1->f[2];
                }
            }
        } // For local nodes (0:7)

    } // For all elements

    return;

}

void calc_conv(mesh_t *myMesh, mysolver_t *mySolver, double theFreq, double theDeltaT, double theDeltaTSquared, damping_type_t typeOfDamping){

    int32_t eindex;
    int i;

    double rmax = 2. * M_PI * theFreq * theDeltaT;
    int32_t   lin_eindex;

    double cdt;

    if (typeOfDamping == BKT) {
        cdt = 2. * M_PI * theFreq * theDeltaT;
    } else {
        cdt = theDeltaT;
    }

    for (lin_eindex = 0; lin_eindex < myLinearElementsCount; lin_eindex++)
    {

        elem_t *elemp;
        edata_t *edata;

        eindex = myLinearElementsMapper[lin_eindex];
        elemp = &myMesh->elemTable[eindex];
        edata = (edata_t *)elemp->data;

        double g0, g02, cg0, eg0;
        double g1, g12, cg1, eg1;
        double g2, g22, cg2, eg2;

        // SHEAR RELATED CONVOLUTION

        if ( (edata->g0_shear != 0) && (edata->g1_shear != 0) ) {

            g0  = cdt * edata->g0_shear;
            g02 = g0 / 2.;
            cg0 = g02 * ( 1. - g0 );
            eg0 = exp( -g0 );

            g1  = cdt * edata->g1_shear;
            g12 = g1 / 2.;
            cg1 = g12 * ( 1. - g1 );
            eg1 = exp( -g1 );

            if (typeOfDamping >= BKT3) {
                g2  = cdt * edata->g2_shear;
                g22 = g2 / 2.;
                cg2 = g22 * ( 1. - g2 );
                eg2 = exp( -g2 );
            }

            for(i = 0; i < 8; i++)
            {
                int32_t     lnid, cindex;
                fvector_t  *tm1Disp, *tm2Disp;
                fvector_t  *f0_tm1, *f1_tm1, *f2_tm1;

                lnid = elemp->lnid[i];

                /* cindex is the index of the node in the convolution vector */
                cindex = eindex * 8 + i;

                tm1Disp = mySolver->tm1 + lnid;
                tm2Disp = mySolver->tm2 + lnid;

                f0_tm1 = mySolver->conv_shear_1 + cindex;
                f1_tm1 = mySolver->conv_shear_2 + cindex;

                f0_tm1->f[0] = cg0 * tm1Disp->f[0] + g02 * tm2Disp->f[0] + eg0 * f0_tm1->f[0];
                f0_tm1->f[1] = cg0 * tm1Disp->f[1] + g02 * tm2Disp->f[1] + eg0 * f0_tm1->f[1];
                f0_tm1->f[2] = cg0 * tm1Disp->f[2] + g02 * tm2Disp->f[2] + eg0 * f0_tm1->f[2];

                f1_tm1->f[0] = cg1 * tm1Disp->f[0] + g12 * tm2Disp->f[0] + eg1 * f1_tm1->f[0];
                f1_tm1->f[1] = cg1 * tm1Disp->f[1] + g12 * tm2Disp->f[1] + eg1 * f1_tm1->f[1];
                f1_tm1->f[2] = cg1 * tm1Disp->f[2] + g12 * tm2Disp->f[2] + eg1 * f1_tm1->f[2];

                if (typeOfDamping >= BKT3) {
                    f2_tm1 = mySolver->conv_shear_3 + cindex;
                    f2_tm1->f[0] = cg2 * tm1Disp->f[0] + g22 * tm2Disp->f[0] + eg2 * f2_tm1->f[0];
                    f2_tm1->f[1] = cg2 * tm1Disp->f[1] + g22 * tm2Disp->f[1] + eg2 * f2_tm1->f[1];
                    f2_tm1->f[2] = cg2 * tm1Disp->f[2] + g22 * tm2Disp->f[2] + eg2 * f2_tm1->f[2];
                }

            } // For local nodes (0:7)

        } // end if null coefficients

        // DILATATION RELATED CONVOLUTION

        if ( (edata->g0_kappa != 0) && (edata->g1_kappa != 0) ) {

            g0  = cdt * edata->g0_kappa;
            g02 = g0 / 2.;
            cg0 = g02 * ( 1. - g0 );
            eg0 = exp( -g0 );

            g1  = cdt * edata->g1_kappa;
            g12 = g1 / 2.;
            cg1 = g12 * ( 1. - g1 );
            eg1 = exp( -g1 );

            if (typeOfDamping >= BKT3) {
                g2  = cdt * edata->g2_kappa;
                g22 = g2 / 2.;
                cg2 = g22 * ( 1. - g2 );
                eg2 = exp( -g2 );
            }

            for(i = 0; i < 8; i++)
            {
                int32_t     lnid, cindex;
                fvector_t  *tm1Disp, *tm2Disp;
                fvector_t  *f0_tm1, *f1_tm1, *f2_tm1;

                lnid = elemp->lnid[i];

                /* cindex is the index of the node in the convolution vector */
                cindex = eindex * 8 + i;

                tm1Disp = mySolver->tm1 + lnid;
                tm2Disp = mySolver->tm2 + lnid;

                f0_tm1 = mySolver->conv_kappa_1 + cindex;
                f1_tm1 = mySolver->conv_kappa_2 + cindex;

                f0_tm1->f[0] = cg0 * tm1Disp->f[0] + g02 * tm2Disp->f[0] + eg0 * f0_tm1->f[0];
                f0_tm1->f[1] = cg0 * tm1Disp->f[1] + g02 * tm2Disp->f[1] + eg0 * f0_tm1->f[1];
                f0_tm1->f[2] = cg0 * tm1Disp->f[2] + g02 * tm2Disp->f[2] + eg0 * f0_tm1->f[2];

                f1_tm1->f[0] = cg1 * tm1Disp->f[0] + g12 * tm2Disp->f[0] + eg1 * f1_tm1->f[0];
                f1_tm1->f[1] = cg1 * tm1Disp->f[1] + g12 * tm2Disp->f[1] + eg1 * f1_tm1->f[1];
                f1_tm1->f[2] = cg1 * tm1Disp->f[2] + g12 * tm2Disp->f[2] + eg1 * f1_tm1->f[2];

                if (typeOfDamping >= BKT3) {
                    f2_tm1 = mySolver->conv_kappa_3 + cindex;
                    f2_tm1->f[0] = cg2 * tm1Disp->f[0] + g22 * tm2Disp->f[0] + eg2 * f2_tm1->f[0];
                    f2_tm1->f[1] = cg2 * tm1Disp->f[1] + g22 * tm2Disp->f[1] + eg2 * f2_tm1->f[1];
                    f2_tm1->f[2] = cg2 * tm1Disp->f[2] + g22 * tm2Disp->f[2] + eg2 * f2_tm1->f[2];
                }

            } // For local nodes (0:7)

        } // end if null coefficients

    } // For all elements

    return;

}


/**
 * new_damping: Compute and add the force due to the element
 *              damping.
 */
void constant_Q_addforce_optimized(mesh_t *myMesh, mysolver_t *mySolver, double theFreq, double theDeltaT, double theDeltaTSquared, damping_type_t typeOfDamping)
{
    /* \todo use mu_and_lamda to compute first,second and third coefficients */

    int i;
    fvector_t localForce[8];
    int32_t   eindex;
    fvector_t damping_vector_shear[8], damping_vector_kappa[8];

    int32_t   lin_eindex;
    double rmax;

    if (typeOfDamping == BKT) {
        rmax = 2. * M_PI * theFreq * theDeltaT;
    } else {
        rmax = theDeltaT;
    }

    /* theAddForceETime -= MPI_Wtime(); */

    /* loop on the number of elements */
    for (lin_eindex = 0; lin_eindex < myLinearElementsCount; lin_eindex++)
    {
        elem_t *elemp;
        e_t    *ep;
        edata_t *edata;

        double a0_shear, a1_shear, a2_shear,
        a0_kappa, a1_kappa, a2_kappa,
        b_shear, b_kappa,
        csum, csumk;

        eindex = myLinearElementsMapper[lin_eindex];
        elemp = &myMesh->elemTable[eindex];
        edata = (edata_t *)elemp->data;
        ep = &mySolver->eTable[eindex];

        a0_shear = edata->a0_shear;
        a1_shear = edata->a1_shear;
        a2_shear = edata->a2_shear;
        b_shear  = edata->b_shear;

        a0_kappa   = edata->a0_kappa;
        a1_kappa   = edata->a1_kappa;
        a2_kappa   = edata->a2_kappa;
        b_kappa    = edata->b_kappa;

        csum = a0_shear + a1_shear + b_shear;
        csumk = a0_kappa + a1_kappa + b_kappa;

        if ( typeOfDamping >= BKT3 ) {
            csum  += a2_shear;
            csumk += a2_kappa;
        }

        double coef_shear = b_shear / rmax;
        double coef_kappa = b_kappa / rmax;

        for (i = 0; i < 8; i++) {

            fvector_t *tm1Disp, *tm2Disp, *f0_tm1, *f1_tm1, *f2_tm1;
            int32_t    lnid, cindex;

            cindex = eindex * 8 + i;
            lnid = elemp->lnid[i];

            tm1Disp = mySolver->tm1 + lnid;
            tm2Disp = mySolver->tm2 + lnid;

            if ( csum != 0 ) {
                f0_tm1  = mySolver->conv_shear_1 + cindex;
                f1_tm1  = mySolver->conv_shear_2 + cindex;

                damping_vector_shear[i].f[0] = coef_shear * (tm1Disp->f[0] - tm2Disp->f[0])
                                                                                     - (a0_shear * f0_tm1->f[0] + a1_shear * f1_tm1->f[0])
                                                                                     + tm1Disp->f[0];

                damping_vector_shear[i].f[1] = coef_shear * (tm1Disp->f[1] - tm2Disp->f[1])
                                                                                     - (a0_shear * f0_tm1->f[1] + a1_shear * f1_tm1->f[1])
                                                                                     + tm1Disp->f[1];

                damping_vector_shear[i].f[2] = coef_shear * (tm1Disp->f[2] - tm2Disp->f[2])
                                                                                     - (a0_shear * f0_tm1->f[2] + a1_shear * f1_tm1->f[2])
                                                                                     + tm1Disp->f[2];

                if ( typeOfDamping >= BKT3 ) {
                    f2_tm1  = mySolver->conv_shear_3 + cindex;
                    damping_vector_shear[i].f[0] -= a2_shear * f2_tm1->f[0];
                    damping_vector_shear[i].f[1] -= a2_shear * f2_tm1->f[1];
                    damping_vector_shear[i].f[2] -= a2_shear * f2_tm1->f[2];
                }
            } else {
                damping_vector_shear[i].f[0] = tm1Disp->f[0];
                damping_vector_shear[i].f[1] = tm1Disp->f[1];
                damping_vector_shear[i].f[2] = tm1Disp->f[2];
            }

            if ( csumk != 0 ) {
                f0_tm1  = mySolver->conv_kappa_1 + cindex;
                f1_tm1  = mySolver->conv_kappa_2 + cindex;

                damping_vector_kappa[i].f[0] = coef_kappa * (tm1Disp->f[0] - tm2Disp->f[0])
                                                                                     - (a0_kappa * f0_tm1->f[0] + a1_kappa * f1_tm1->f[0])
                                                                                     + tm1Disp->f[0];

                damping_vector_kappa[i].f[1] = coef_kappa * (tm1Disp->f[1] - tm2Disp->f[1])
                                                                                     - (a0_kappa * f0_tm1->f[1] + a1_kappa * f1_tm1->f[1])
                                                                                     + tm1Disp->f[1];

                damping_vector_kappa[i].f[2] = coef_kappa * (tm1Disp->f[2] - tm2Disp->f[2])
                                                                                     - (a0_kappa * f0_tm1->f[2] + a1_kappa * f1_tm1->f[2])
                                                                                     + tm1Disp->f[2];

                if ( typeOfDamping >= BKT3 ) {
                    f2_tm1  = mySolver->conv_kappa_3 + cindex;
                    damping_vector_kappa[i].f[0] -= a2_kappa * f2_tm1->f[0];
                    damping_vector_kappa[i].f[1] -= a2_kappa * f2_tm1->f[1];
                    damping_vector_kappa[i].f[2] -= a2_kappa * f2_tm1->f[2];
                }
            } else {

                damping_vector_kappa[i].f[0] = tm1Disp->f[0];
                damping_vector_kappa[i].f[1] = tm1Disp->f[1];
                damping_vector_kappa[i].f[2] = tm1Disp->f[2];

            }
        } // end for nodes in the element

        double kappa = -0.5625 * (ep->c2 + 2. / 3. * ep->c1);
        double mu = -0.5625 * ep->c1;

        double atu[24] = {0.0};
        double firstVec[24] = {0.0};

        memset(localForce, 0, 8 * sizeof(fvector_t));

        if(vector_is_zero( damping_vector_shear ) != 0) {

            aTransposeU( damping_vector_shear, atu );
            firstVector_mu( atu, firstVec, mu);

        }

        if(vector_is_zero( damping_vector_kappa ) != 0) {

            aTransposeU( damping_vector_kappa, atu );
            firstVector_kappa( atu, firstVec, kappa);

        }

        au( localForce, firstVec );

        for (i = 0; i < 8; i++) {
            int32_t lnid;
            fvector_t *nodalForce;

            lnid = elemp->lnid[i];

            nodalForce = mySolver->force + lnid;
            nodalForce->f[0] += localForce[i].f[0];
            nodalForce->f[1] += localForce[i].f[1];
            nodalForce->f[2] += localForce[i].f[2];
        }

    } /* for all the elements */

    return;

}

void constant_Q_addforce(mesh_t *myMesh, mysolver_t *mySolver, double theFreq, double theDeltaT, double theDeltaTSquared, damping_type_t typeOfDamping)
{
    /* \todo use mu_and_lamda to compute first,second and third coefficients */

    int i;
    fvector_t localForce[8];
    int32_t   eindex;
    fvector_t damping_vector_shear[8], damping_vector_kappa[8];

    int32_t   lin_eindex;
    double rmax;

    if (typeOfDamping == BKT) {
        rmax = 2. * M_PI * theFreq * theDeltaT;
    } else {
        rmax = theDeltaT;
    }

    /* theAddForceETime -= MPI_Wtime(); */

    /* loop on the number of elements */
    for (lin_eindex = 0; lin_eindex < myLinearElementsCount; lin_eindex++)
    {
        elem_t *elemp;
        e_t    *ep;
        edata_t *edata;

        double a0_shear, a1_shear, a2_shear,
               a0_kappa, a1_kappa, a2_kappa,
               b_shear, b_kappa,
               csum;

        eindex = myLinearElementsMapper[lin_eindex];
        elemp = &myMesh->elemTable[eindex];
        edata = (edata_t *)elemp->data;
        ep = &mySolver->eTable[eindex];

        // SHEAR CONTRIBUTION

        a0_shear = edata->a0_shear;
        a1_shear = edata->a1_shear;
        a2_shear = edata->a2_shear;
        b_shear  = edata->b_shear;

        csum = a0_shear + a1_shear + b_shear;
        if ( typeOfDamping >= BKT3 ) {
            csum += a2_shear;
        }

        if ( csum != 0 ) {

            double coef_shear = b_shear / rmax;

            for (i = 0; i < 8; i++) {

                fvector_t *tm1Disp, *tm2Disp, *f0_tm1, *f1_tm1, *f2_tm1;
                int32_t    lnid, cindex;

                cindex = eindex * 8 + i;
                lnid = elemp->lnid[i];

                tm1Disp = mySolver->tm1 + lnid;
                tm2Disp = mySolver->tm2 + lnid;
                f0_tm1  = mySolver->conv_shear_1 + cindex;
                f1_tm1  = mySolver->conv_shear_2 + cindex;

                damping_vector_shear[i].f[0] = coef_shear * (tm1Disp->f[0] - tm2Disp->f[0])
                                             - (a0_shear * f0_tm1->f[0] + a1_shear * f1_tm1->f[0])
                                             + tm1Disp->f[0];

                damping_vector_shear[i].f[1] = coef_shear * (tm1Disp->f[1] - tm2Disp->f[1])
                                             - (a0_shear * f0_tm1->f[1] + a1_shear * f1_tm1->f[1])
                                             + tm1Disp->f[1];

                damping_vector_shear[i].f[2] = coef_shear * (tm1Disp->f[2] - tm2Disp->f[2])
                                             - (a0_shear * f0_tm1->f[2] + a1_shear * f1_tm1->f[2])
                                             + tm1Disp->f[2];

                if ( typeOfDamping >= BKT3 ) {
                    f2_tm1  = mySolver->conv_shear_3 + cindex;
                    damping_vector_shear[i].f[0] -= a2_shear * f2_tm1->f[0];
                    damping_vector_shear[i].f[1] -= a2_shear * f2_tm1->f[1];
                    damping_vector_shear[i].f[2] -= a2_shear * f2_tm1->f[2];
                }

            } // end for nodes in the element

        } else {

            for (i = 0; i < 8; i++) {

                fvector_t *tm1Disp, *tm2Disp;
                int32_t    lnid;

                lnid = elemp->lnid[i];
                tm1Disp = mySolver->tm1 + lnid;
                tm2Disp = mySolver->tm2 + lnid;

                damping_vector_shear[i].f[0] = tm1Disp->f[0];
                damping_vector_shear[i].f[1] = tm1Disp->f[1];
                damping_vector_shear[i].f[2] = tm1Disp->f[2];

            } // end for nodes in the element

        } // end if for coefficients

        // DILATATION CONTRIBUTION

        a0_kappa   = edata->a0_kappa;
        a1_kappa   = edata->a1_kappa;
        a2_kappa   = edata->a2_kappa;
        b_kappa    = edata->b_kappa;

        csum = a0_kappa + a1_kappa + b_kappa;
        if ( typeOfDamping >= BKT3 ) {
            csum += a2_kappa;
        }

        if ( csum != 0 ) {

            double coef_kappa = b_kappa / rmax;

            for (i = 0; i < 8; i++) {

                fvector_t *tm1Disp, *tm2Disp, *f0_tm1, *f1_tm1, *f2_tm1;
                int32_t    lnid, cindex;

                cindex = eindex * 8 + i;

                lnid = elemp->lnid[i];

                tm1Disp = mySolver->tm1 + lnid;
                tm2Disp = mySolver->tm2 + lnid;

                f0_tm1  = mySolver->conv_kappa_1 + cindex;
                f1_tm1  = mySolver->conv_kappa_2 + cindex;

                damping_vector_kappa[i].f[0] = coef_kappa * (tm1Disp->f[0] - tm2Disp->f[0])
                                             - (a0_kappa * f0_tm1->f[0] + a1_kappa * f1_tm1->f[0])
                                             + tm1Disp->f[0];

                damping_vector_kappa[i].f[1] = coef_kappa * (tm1Disp->f[1] - tm2Disp->f[1])
                                             - (a0_kappa * f0_tm1->f[1] + a1_kappa * f1_tm1->f[1])
                                             + tm1Disp->f[1];

                damping_vector_kappa[i].f[2] = coef_kappa * (tm1Disp->f[2] - tm2Disp->f[2])
                                             - (a0_kappa * f0_tm1->f[2] + a1_kappa * f1_tm1->f[2])
                                             + tm1Disp->f[2];

                if ( typeOfDamping >= BKT3 ) {
                    f2_tm1  = mySolver->conv_kappa_3 + cindex;
                    damping_vector_kappa[i].f[0] -= a2_kappa * f2_tm1->f[0];
                    damping_vector_kappa[i].f[1] -= a2_kappa * f2_tm1->f[1];
                    damping_vector_kappa[i].f[2] -= a2_kappa * f2_tm1->f[2];
                }

            } // end for nodes in the element

        } else {

            for (i = 0; i < 8; i++) {

                fvector_t *tm1Disp, *tm2Disp;
                int32_t    lnid;

                lnid = elemp->lnid[i];
                tm1Disp = mySolver->tm1 + lnid;
                tm2Disp = mySolver->tm2 + lnid;

                damping_vector_kappa[i].f[0] = tm1Disp->f[0];
                damping_vector_kappa[i].f[1] = tm1Disp->f[1];
                damping_vector_kappa[i].f[2] = tm1Disp->f[2];

            } // end for nodes in the element

        } // end if for coefficients

        double kappa = -0.5625 * (ep->c2 + 2. / 3. * ep->c1);
        double mu = -0.5625 * ep->c1;

        double atu[24];
        double firstVec[24];

        if ( edata->topoBkt == 0 ) {

            for(i = 0; i<24; i++)
                firstVec[i] = 0.;

            memset(localForce, 0, 8 * sizeof(fvector_t));

            if(vector_is_zero( damping_vector_shear ) != 0) {

                aTransposeU( damping_vector_shear, atu );
                firstVector_mu( atu, firstVec, mu);

            }

            if(vector_is_zero( damping_vector_kappa ) != 0) {

                aTransposeU( damping_vector_kappa, atu );
                firstVector_kappa( atu, firstVec, kappa);

            }

            au( localForce, firstVec );

            for (i = 0; i < 8; i++) {
                int32_t lnid;
                fvector_t *nodalForce;

                lnid = elemp->lnid[i];

                nodalForce = mySolver->force + lnid;
                nodalForce->f[0] += localForce[i].f[0];
                nodalForce->f[1] += localForce[i].f[1];
                nodalForce->f[2] += localForce[i].f[2];
            }

        } else {

            fvector_t localForce[8];
            memset( localForce, 0, 8 * sizeof(fvector_t) );

            mu    = edata->rho * edata->Vs * edata->Vs;
            kappa = edata->rho * ( edata->Vp * edata->Vp - 4./3. * edata->Vs * edata->Vs );

            if( vector_is_zero( damping_vector_shear ) != 0 ||
                vector_is_zero( damping_vector_kappa ) != 0) {

                TetraForcesBKT(  damping_vector_shear, damping_vector_kappa, localForce,
                        edata->edgesize,
                        mu, kappa,
                        edata->topo_eindex );
            }

            /* Loop over the 8 element nodes:
             * Add the contribution calculated above to the node
             * forces carried from the source and stiffness.
             */

            for (i = 0; i < 8; i++) {

                int32_t    lnid;
                fvector_t *nodalForce;

                lnid = elemp->lnid[i];

                nodalForce = mySolver->force + lnid;

                nodalForce->f[0] -= localForce[i].f[0] * theDeltaTSquared;
                nodalForce->f[1] -= localForce[i].f[1] * theDeltaTSquared;
                nodalForce->f[2] -= localForce[i].f[2] * theDeltaTSquared;

            } /* element nodes */

        }

    } /* for all the elements */

    return;

}


void constant_Q_addforce_simplifiedConv_old(mesh_t *myMesh, mysolver_t *mySolver, double theFreq, double theDeltaT, double theDeltaTSquared, damping_type_t typeOfDamping)
{
    /* \todo use mu_and_lamda to compute first,second and third coefficients */

    int i;
    fvector_t localForce[8];
    int32_t   eindex;
    fvector_t damping_vector_shear[8], damping_vector_kappa[8];

    int32_t   lin_eindex;
    double rmax, dt=theDeltaT;


    if (typeOfDamping == BKT) {
        rmax = 2. * M_PI * theFreq * theDeltaT;
    } else {
        rmax = theDeltaT;
    }

    /* theAddForceETime -= MPI_Wtime(); */

    /* loop on the number of elements */
    for (lin_eindex = 0; lin_eindex < myLinearElementsCount; lin_eindex++)
    {
        elem_t *elemp;
        e_t    *ep;
        edata_t *edata;

        double a0_shear, a1_shear, a2_shear,
        a0_kappa, a1_kappa, a2_kappa,
        g0_shear, g1_shear, g2_shear,
        g0_kappa, g1_kappa, g2_kappa,
        b_shear, b_kappa,
        csum, csumk, phi_shear, psi_shear, phi_kappa, psi_kappa;

        eindex = myLinearElementsMapper[lin_eindex];
        elemp = &myMesh->elemTable[eindex];
        edata = (edata_t *)elemp->data;
        ep = &mySolver->eTable[eindex];

        a0_shear = edata->a0_shear;
        a1_shear = edata->a1_shear;
        a2_shear = edata->a2_shear;
        b_shear  = edata->b_shear;

        g0_shear = edata->g0_shear;
        g1_shear = edata->g1_shear;
        g2_shear = edata->g2_shear;

        a0_kappa   = edata->a0_kappa;
        a1_kappa   = edata->a1_kappa;
        a2_kappa   = edata->a2_kappa;
        b_kappa    = edata->b_kappa;

        g0_kappa   = edata->g0_kappa;
        g1_kappa   = edata->g1_kappa;
        g2_kappa   = edata->g2_kappa;

        double coef_shear = b_shear / rmax;
        double coef_kappa = b_kappa / rmax;

        for (i = 0; i < 8; i++) {

            fvector_t *tm1Disp, *tm2Disp;
            int32_t    lnid, cindex;

            cindex = eindex * 8 + i;
            lnid = elemp->lnid[i];

            tm1Disp = mySolver->tm1 + lnid;
            tm2Disp = mySolver->tm2 + lnid;


            damping_vector_shear[i].f[0] = coef_shear * (tm1Disp->f[0] - tm2Disp->f[0])
                                           - (    a0_shear * g0_shear * ( dt/2.0 * ( (1.0 - g0_shear*dt) * tm1Disp->f[0] + tm2Disp->f[0] ) + exp(-g0_shear*dt) * dt / 2.0 * ( 1.0 - g0_shear * dt ) * tm2Disp->f[0] )
                                                + a1_shear * g1_shear * ( dt/2.0 * ( (1.0 - g1_shear*dt) * tm1Disp->f[0] + tm2Disp->f[0] ) + exp(-g1_shear*dt) * dt / 2.0 * ( 1.0 - g1_shear * dt ) * tm2Disp->f[0] ) )
                                           + tm1Disp->f[0];
            damping_vector_shear[i].f[1] = coef_shear * (tm1Disp->f[1] - tm2Disp->f[1])
                                           - (  a0_shear * g0_shear * ( dt/2.0 * ( (1.0 - g0_shear*dt) * tm1Disp->f[1] + tm2Disp->f[1] ) + exp(-g0_shear*dt) * dt / 2.0 * ( 1.0 - g0_shear * dt ) * tm2Disp->f[1] )
                                              + a1_shear * g1_shear * ( dt/2.0 * ( (1.0 - g1_shear*dt) * tm1Disp->f[1] + tm2Disp->f[1] ) + exp(-g1_shear*dt) * dt / 2.0 * ( 1.0 - g1_shear * dt ) * tm2Disp->f[1] ) )
                                           + tm1Disp->f[1];
            damping_vector_shear[i].f[2] = coef_shear * (tm1Disp->f[2] - tm2Disp->f[2])
                                           - (  a0_shear * g0_shear * ( dt/2.0 * ( (1.0 - g0_shear*dt) * tm1Disp->f[2] + tm2Disp->f[2] ) + exp(-g0_shear*dt) * dt / 2.0 * ( 1.0 - g0_shear * dt ) * tm2Disp->f[2] )
                                              + a1_shear * g1_shear * ( dt/2.0 * ( (1.0 - g1_shear*dt) * tm1Disp->f[2] + tm2Disp->f[2] ) + exp(-g1_shear*dt) * dt / 2.0 * ( 1.0 - g1_shear * dt ) * tm2Disp->f[2] ) )
                                           + tm1Disp->f[2];


            damping_vector_kappa[i].f[0] = coef_kappa * (tm1Disp->f[0] - tm2Disp->f[0])
                                           - (  a0_kappa * g0_kappa * ( dt/2.0 * ( (1.0 - g0_kappa*dt) * tm1Disp->f[0] + tm2Disp->f[0] ) + exp(-g0_kappa*dt) * dt / 2.0 * ( 1.0 - g0_kappa * dt ) * tm2Disp->f[0] )
                                              + a1_kappa * g1_kappa * ( dt/2.0 * ( (1.0 - g1_kappa*dt) * tm1Disp->f[0] + tm2Disp->f[0] ) + exp(-g1_kappa*dt) * dt / 2.0 * ( 1.0 - g1_kappa * dt ) * tm2Disp->f[0] ) )
                                           + tm1Disp->f[0];
            damping_vector_kappa[i].f[1] = coef_kappa * (tm1Disp->f[1] - tm2Disp->f[1])
                                           - (  a0_kappa * g0_kappa * ( dt/2.0 * ( (1.0 - g0_kappa*dt) * tm1Disp->f[1] + tm2Disp->f[1] ) + exp(-g0_kappa*dt) * dt / 2.0 * ( 1.0 - g0_kappa * dt ) * tm2Disp->f[1] )
                                              + a1_kappa * g1_kappa * ( dt/2.0 * ( (1.0 - g1_kappa*dt) * tm1Disp->f[1] + tm2Disp->f[1] ) + exp(-g1_kappa*dt) * dt / 2.0 * ( 1.0 - g1_kappa * dt ) * tm2Disp->f[1] ) )
                                           + tm1Disp->f[1];
            damping_vector_kappa[i].f[2] = coef_kappa * (tm1Disp->f[2] - tm2Disp->f[2])
                                           - (  a0_kappa * g0_kappa * ( dt/2.0 * ( (1.0 - g0_kappa*dt) * tm1Disp->f[2] + tm2Disp->f[2] ) + exp(-g0_kappa*dt) * dt / 2.0 * ( 1.0 - g0_kappa * dt ) * tm2Disp->f[2] )
                                           + a1_kappa * g1_kappa * ( dt/2.0 * ( (1.0 - g1_kappa*dt) * tm1Disp->f[2] + tm2Disp->f[2] ) + exp(-g1_kappa*dt) * dt / 2.0 * ( 1.0 - g1_kappa * dt ) * tm2Disp->f[2] ) )
                                           + tm1Disp->f[2];

        } // end for nodes in the element

        double kappa = -0.5625 * (ep->c2 + 2. / 3. * ep->c1);
        double mu = -0.5625 * ep->c1;

        double atu[24] = {0.0};
        double firstVec[24] = {0.0};

        memset(localForce, 0, 8 * sizeof(fvector_t));

        if(vector_is_zero( damping_vector_shear ) != 0) {

            aTransposeU( damping_vector_shear, atu );
            firstVector_mu( atu, firstVec, mu);

        }

        if(vector_is_zero( damping_vector_kappa ) != 0) {

            aTransposeU( damping_vector_kappa, atu );
            firstVector_kappa( atu, firstVec, kappa);

        }

        au( localForce, firstVec );

        for (i = 0; i < 8; i++) {
            int32_t lnid;
            fvector_t *nodalForce;

            lnid = elemp->lnid[i];

            nodalForce = mySolver->force + lnid;
            nodalForce->f[0] += localForce[i].f[0];
            nodalForce->f[1] += localForce[i].f[1];
            nodalForce->f[2] += localForce[i].f[2];
        }

    } /* for all the elements */

    return;

}

void constant_Q_addforce_simplifiedConvolution(mesh_t *myMesh, mysolver_t *mySolver, double theFreq, double theDeltaT, double theDeltaTSquared, damping_type_t typeOfDamping)
{
    /* \todo use mu_and_lamda to compute first,second and third coefficients */

    int i;
    fvector_t localForce[8];
    int32_t   eindex;
    fvector_t damping_vector_shear[8], damping_vector_kappa[8];

    int32_t   lin_eindex;
    double    dt=theDeltaT;


    /* loop on the number of elements */
    for (lin_eindex = 0; lin_eindex < myLinearElementsCount; lin_eindex++)
    {
        elem_t *elemp;
        e_t    *ep;
        edata_t *edata;

        double a0_shear, a1_shear, a2_shear,
        a0_kappa, a1_kappa, a2_kappa,
        g0_shear, g1_shear, g2_shear,
        g0_kappa, g1_kappa, g2_kappa,
        b_shear, b_kappa,
        phi_shear, psi_shear, phi_kappa, psi_kappa;

        eindex = myLinearElementsMapper[lin_eindex];
        elemp = &myMesh->elemTable[eindex];
        edata = (edata_t *)elemp->data;
        ep = &mySolver->eTable[eindex];

        a0_shear = edata->a0_shear;
        a1_shear = edata->a1_shear;
        a2_shear = edata->a2_shear;
        b_shear  = edata->b_shear;

        g0_shear = edata->g0_shear;
        g1_shear = edata->g1_shear;
        g2_shear = edata->g2_shear;

        a0_kappa   = edata->a0_kappa;
        a1_kappa   = edata->a1_kappa;
        a2_kappa   = edata->a2_kappa;
        b_kappa    = edata->b_kappa;

        g0_kappa   = edata->g0_kappa;
        g1_kappa   = edata->g1_kappa;
        g2_kappa   = edata->g2_kappa;


        phi_shear = dt/2.0 * ( a0_shear * g0_shear * ( 1 - g0_shear * dt ) + a1_shear * g1_shear * ( 1 - g1_shear * dt ) );
        phi_kappa = dt/2.0 * ( a0_kappa * g0_kappa * ( 1 - g0_kappa * dt ) + a1_kappa * g1_kappa * ( 1 - g1_kappa * dt ) );

        psi_shear = dt/2.0 * ( a0_shear * g0_shear * ( 1 + exp(-g0_shear*dt)*( 1 - g0_shear * dt) )
                             + a1_shear * g1_shear * ( 1 + exp(-g1_shear*dt)*( 1 - g1_shear * dt) ) );

        psi_kappa = dt/2.0 * ( a0_kappa * g0_kappa * ( 1 + exp(-g0_kappa*dt) * ( 1 - g0_kappa * dt) )
                             + a1_kappa * g1_kappa * ( 1 + exp(-g1_kappa*dt) * ( 1 - g1_kappa * dt) ) );

        for (i = 0; i < 8; i++) {

            fvector_t *tm1Disp, *tm2Disp;
            int32_t    lnid, cindex;

            cindex = eindex * 8 + i;
            lnid = elemp->lnid[i];

            tm1Disp = mySolver->tm1 + lnid;
            tm2Disp = mySolver->tm2 + lnid;

            damping_vector_shear[i].f[0] = ( 1.0 + b_shear/dt - phi_shear ) * tm1Disp->f[0] - ( b_shear/dt + psi_shear ) * tm2Disp->f[0];
            damping_vector_shear[i].f[1] = ( 1.0 + b_shear/dt - phi_shear ) * tm1Disp->f[1] - ( b_shear/dt + psi_shear ) * tm2Disp->f[1];
            damping_vector_shear[i].f[2] = ( 1.0 + b_shear/dt - phi_shear ) * tm1Disp->f[2] - ( b_shear/dt + psi_shear ) * tm2Disp->f[2];

            damping_vector_kappa[i].f[0] = ( 1.0 + b_kappa/dt - phi_kappa ) * tm1Disp->f[0] - ( b_kappa/dt + psi_kappa ) * tm2Disp->f[0];
            damping_vector_kappa[i].f[1] = ( 1.0 + b_kappa/dt - phi_kappa ) * tm1Disp->f[1] - ( b_kappa/dt + psi_kappa ) * tm2Disp->f[1];
            damping_vector_kappa[i].f[2] = ( 1.0 + b_kappa/dt - phi_kappa ) * tm1Disp->f[2] - ( b_kappa/dt + psi_kappa ) * tm2Disp->f[2];

        } // end for nodes in the element

        double kappa = -0.5625 * (ep->c2 + 2. / 3. * ep->c1);
        double mu = -0.5625 * ep->c1;

        double atu[24] = {0.0};
        double firstVec[24] = {0.0};

        if(vector_is_zero( damping_vector_shear ) != 0) {

            aTransposeU( damping_vector_shear, atu );
            firstVector_mu( atu, firstVec, mu);

        }

        if(vector_is_zero( damping_vector_kappa ) != 0) {

            aTransposeU( damping_vector_kappa, atu );
            firstVector_kappa( atu, firstVec, kappa);

        }

        memset(localForce, 0, 8 * sizeof(fvector_t));

        au( localForce, firstVec );

        for (i = 0; i < 8; i++) {
            int32_t lnid;
            fvector_t *nodalForce;

            lnid = elemp->lnid[i];

            nodalForce = mySolver->force + lnid;
            nodalForce->f[0] += localForce[i].f[0];
            nodalForce->f[1] += localForce[i].f[1];
            nodalForce->f[2] += localForce[i].f[2];
        }

    } /* for all the elements */

    return;

}


void conv_and_bktForceCombined(mesh_t *myMesh, mysolver_t *mySolver, double theFreq, double theDeltaT, double theDeltaTSquared, damping_type_t typeOfDamping){

    int32_t eindex;
    int i;

    double rmax;
    int32_t   lin_eindex;

    double cdt;

    double g0, g02, cg0, eg0, g0k, g02k, cg0k, eg0k;
    double g1, g12, cg1, eg1, g1k, g12k, cg1k, eg1k;
    double g2, g22, cg2, eg2, g2k, g22k, cg2k, eg2k;

    fvector_t localForce[8];
    fvector_t damping_vector_shear[8], damping_vector_kappa[8];

    if (typeOfDamping == BKT) {
        cdt = 2. * M_PI * theFreq * theDeltaT;
    } else {
        cdt = theDeltaT;
    }

    for (lin_eindex = 0; lin_eindex < myLinearElementsCount; lin_eindex++)
    {

        elem_t *elemp;
        edata_t *edata;
        e_t    *ep;

        eindex = myLinearElementsMapper[lin_eindex];
        elemp = &myMesh->elemTable[eindex];
        edata = (edata_t *)elemp->data;

        // SHEAR RELATED CONVOLUTION
        g0  = cdt * edata->g0_shear;
        g02 = g0 / 2.;
        cg0 = g02 * ( 1. - g0 );
        eg0 = exp( -g0 );

        g1  = cdt * edata->g1_shear;
        g12 = g1 / 2.;
        cg1 = g12 * ( 1. - g1 );
        eg1 = exp( -g1 );

        g0k  = cdt * edata->g0_kappa;
        g02k = g0k / 2.;
        cg0k = g02k * ( 1. - g0k );
        eg0k = exp( -g0k );

        g1k  = cdt * edata->g1_kappa;
        g12k = g1k / 2.;
        cg1k = g12k * ( 1. - g1k );
        eg1k = exp( -g1k );

        if (typeOfDamping >= BKT3) {
            g2  = cdt * edata->g2_shear;
            g22 = g2 / 2.;
            cg2 = g22 * ( 1. - g2 );
            eg2 = exp( -g2 );

            g2k  = cdt * edata->g2_kappa;
            g22k = g2k / 2.;
            cg2k = g22k * ( 1. - g2k );
            eg2k = exp( -g2k );
        }

        for(i = 0; i < 8; i++)
        {
            int32_t     lnid, cindex;
            fvector_t  *tm1Disp, *tm2Disp;
            fvector_t  *f0_tm1, *f1_tm1, *f2_tm1;

            lnid = elemp->lnid[i];

            /* cindex is the index of the node in the convolution vector */
            cindex = eindex * 8 + i;

            tm1Disp = mySolver->tm1 + lnid;
            tm2Disp = mySolver->tm2 + lnid;

            // SHEAR RELATED CONVOLUTION
            if ( (edata->g0_shear != 0) && (edata->g1_shear != 0) ) {
                f0_tm1 = mySolver->conv_shear_1 + cindex;
                f1_tm1 = mySolver->conv_shear_2 + cindex;

                f0_tm1->f[0] = cg0 * tm1Disp->f[0] + g02 * tm2Disp->f[0] + eg0 * f0_tm1->f[0];
                f0_tm1->f[1] = cg0 * tm1Disp->f[1] + g02 * tm2Disp->f[1] + eg0 * f0_tm1->f[1];
                f0_tm1->f[2] = cg0 * tm1Disp->f[2] + g02 * tm2Disp->f[2] + eg0 * f0_tm1->f[2];

                f1_tm1->f[0] = cg1 * tm1Disp->f[0] + g12 * tm2Disp->f[0] + eg1 * f1_tm1->f[0];
                f1_tm1->f[1] = cg1 * tm1Disp->f[1] + g12 * tm2Disp->f[1] + eg1 * f1_tm1->f[1];
                f1_tm1->f[2] = cg1 * tm1Disp->f[2] + g12 * tm2Disp->f[2] + eg1 * f1_tm1->f[2];

                if (typeOfDamping >= BKT3) {
                    f2_tm1 = mySolver->conv_shear_3 + cindex;
                    f2_tm1->f[0] = cg2 * tm1Disp->f[0] + g22 * tm2Disp->f[0] + eg2 * f2_tm1->f[0];
                    f2_tm1->f[1] = cg2 * tm1Disp->f[1] + g22 * tm2Disp->f[1] + eg2 * f2_tm1->f[1];
                    f2_tm1->f[2] = cg2 * tm1Disp->f[2] + g22 * tm2Disp->f[2] + eg2 * f2_tm1->f[2];
                }
            }

            // DILATATION RELATED CONVOLUTION
            if ( (edata->g0_kappa != 0) && (edata->g1_kappa != 0) ) {
                f0_tm1 = mySolver->conv_kappa_1 + cindex;
                f1_tm1 = mySolver->conv_kappa_2 + cindex;

                f0_tm1->f[0] = cg0k * tm1Disp->f[0] + g02k * tm2Disp->f[0] + eg0k * f0_tm1->f[0];
                f0_tm1->f[1] = cg0k * tm1Disp->f[1] + g02k * tm2Disp->f[1] + eg0k * f0_tm1->f[1];
                f0_tm1->f[2] = cg0k * tm1Disp->f[2] + g02k * tm2Disp->f[2] + eg0k * f0_tm1->f[2];

                f1_tm1->f[0] = cg1k * tm1Disp->f[0] + g12k * tm2Disp->f[0] + eg1k * f1_tm1->f[0];
                f1_tm1->f[1] = cg1k * tm1Disp->f[1] + g12k * tm2Disp->f[1] + eg1k * f1_tm1->f[1];
                f1_tm1->f[2] = cg1k * tm1Disp->f[2] + g12k * tm2Disp->f[2] + eg1k * f1_tm1->f[2];

                if (typeOfDamping >= BKT3) {
                    f2_tm1 = mySolver->conv_kappa_3 + cindex;
                    f2_tm1->f[0] = cg2k * tm1Disp->f[0] + g22k * tm2Disp->f[0] + eg2k * f2_tm1->f[0];
                    f2_tm1->f[1] = cg2k * tm1Disp->f[1] + g22k * tm2Disp->f[1] + eg2k * f2_tm1->f[1];
                    f2_tm1->f[2] = cg2k * tm1Disp->f[2] + g22k * tm2Disp->f[2] + eg2k * f2_tm1->f[2];
                }
            }
        } // For local nodes (0:7)

        // =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
        // =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=* apply force =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
        // =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

        double a0_shear, a1_shear, a2_shear,
        a0_kappa, a1_kappa, a2_kappa,
        b_shear, b_kappa,
        csum, csumk;

        if (typeOfDamping == BKT) {
            rmax = 2. * M_PI * theFreq * theDeltaT;
        } else {
            rmax = theDeltaT;
        }

        //eindex = myLinearElementsMapper[lin_eindex];
        //elemp = &myMesh->elemTable[eindex];
        //edata = (edata_t *)elemp->data;
        ep = &mySolver->eTable[eindex];

        a0_shear = edata->a0_shear;
        a1_shear = edata->a1_shear;
        a2_shear = edata->a2_shear;
        b_shear  = edata->b_shear;

        a0_kappa   = edata->a0_kappa;
        a1_kappa   = edata->a1_kappa;
        a2_kappa   = edata->a2_kappa;
        b_kappa    = edata->b_kappa;

        csum = a0_shear + a1_shear + b_shear;
        csumk = a0_kappa + a1_kappa + b_kappa;

        if ( typeOfDamping >= BKT3 ) {
            csum  += a2_shear;
            csumk += a2_kappa;
        }

        double coef_shear = b_shear / rmax;
        double coef_kappa = b_kappa / rmax;

        for (i = 0; i < 8; i++) {

            fvector_t *tm1Disp, *tm2Disp, *f0_tm1, *f1_tm1, *f2_tm1;
            int32_t    lnid, cindex;

            cindex = eindex * 8 + i;
            lnid = elemp->lnid[i];

            tm1Disp = mySolver->tm1 + lnid;
            tm2Disp = mySolver->tm2 + lnid;

            if ( csum != 0 ) {
                f0_tm1  = mySolver->conv_shear_1 + cindex;
                f1_tm1  = mySolver->conv_shear_2 + cindex;

                damping_vector_shear[i].f[0] = coef_shear * (tm1Disp->f[0] - tm2Disp->f[0])
                                                                                     - (a0_shear * f0_tm1->f[0] + a1_shear * f1_tm1->f[0])
                                                                                     + tm1Disp->f[0];

                damping_vector_shear[i].f[1] = coef_shear * (tm1Disp->f[1] - tm2Disp->f[1])
                                                                                     - (a0_shear * f0_tm1->f[1] + a1_shear * f1_tm1->f[1])
                                                                                     + tm1Disp->f[1];

                damping_vector_shear[i].f[2] = coef_shear * (tm1Disp->f[2] - tm2Disp->f[2])
                                                                                     - (a0_shear * f0_tm1->f[2] + a1_shear * f1_tm1->f[2])
                                                                                     + tm1Disp->f[2];

                if ( typeOfDamping >= BKT3 ) {
                    f2_tm1  = mySolver->conv_shear_3 + cindex;
                    damping_vector_shear[i].f[0] -= a2_shear * f2_tm1->f[0];
                    damping_vector_shear[i].f[1] -= a2_shear * f2_tm1->f[1];
                    damping_vector_shear[i].f[2] -= a2_shear * f2_tm1->f[2];
                }
            } else {
                damping_vector_shear[i].f[0] = tm1Disp->f[0];
                damping_vector_shear[i].f[1] = tm1Disp->f[1];
                damping_vector_shear[i].f[2] = tm1Disp->f[2];
            }

            if ( csumk != 0 ) {
                f0_tm1  = mySolver->conv_kappa_1 + cindex;
                f1_tm1  = mySolver->conv_kappa_2 + cindex;

                damping_vector_kappa[i].f[0] = coef_kappa * (tm1Disp->f[0] - tm2Disp->f[0])
                                                                                     - (a0_kappa * f0_tm1->f[0] + a1_kappa * f1_tm1->f[0])
                                                                                     + tm1Disp->f[0];

                damping_vector_kappa[i].f[1] = coef_kappa * (tm1Disp->f[1] - tm2Disp->f[1])
                                                                                     - (a0_kappa * f0_tm1->f[1] + a1_kappa * f1_tm1->f[1])
                                                                                     + tm1Disp->f[1];

                damping_vector_kappa[i].f[2] = coef_kappa * (tm1Disp->f[2] - tm2Disp->f[2])
                                                                                     - (a0_kappa * f0_tm1->f[2] + a1_kappa * f1_tm1->f[2])
                                                                                     + tm1Disp->f[2];

                if ( typeOfDamping >= BKT3 ) {
                    f2_tm1  = mySolver->conv_kappa_3 + cindex;
                    damping_vector_kappa[i].f[0] -= a2_kappa * f2_tm1->f[0];
                    damping_vector_kappa[i].f[1] -= a2_kappa * f2_tm1->f[1];
                    damping_vector_kappa[i].f[2] -= a2_kappa * f2_tm1->f[2];
                }
            } else {

                damping_vector_kappa[i].f[0] = tm1Disp->f[0];
                damping_vector_kappa[i].f[1] = tm1Disp->f[1];
                damping_vector_kappa[i].f[2] = tm1Disp->f[2];

            }
        } // end for nodes in the element

        double kappa = -0.5625 * (ep->c2 + 2. / 3. * ep->c1);
        double mu = -0.5625 * ep->c1;

        double atu[24] = {0.0};
        double firstVec[24] = {0.0};

        memset(localForce, 0, 8 * sizeof(fvector_t));

        if(vector_is_zero( damping_vector_shear ) != 0) {

            aTransposeU( damping_vector_shear, atu );
            firstVector_mu( atu, firstVec, mu);

        }

        if(vector_is_zero( damping_vector_kappa ) != 0) {

            aTransposeU( damping_vector_kappa, atu );
            firstVector_kappa( atu, firstVec, kappa);

        }

        au( localForce, firstVec );

        for (i = 0; i < 8; i++) {
            int32_t lnid;
            fvector_t *nodalForce;

            lnid = elemp->lnid[i];

            nodalForce = mySolver->force + lnid;
            nodalForce->f[0] += localForce[i].f[0];
            nodalForce->f[1] += localForce[i].f[1];
            nodalForce->f[2] += localForce[i].f[2];
        }

    } // For all elements

    return;

}

void BKT_TU_transf( double theFreq, double theDeltaT, double theDeltaTSquared, damping_type_t typeOfDamping,
        fvector_t *convShear1, fvector_t *convKappa1,
        fvector_t *convShear2, fvector_t *convKappa2,
        fvector_t *convShear3, fvector_t *convKappa3,
        edata_t *edata,
        fvector_t  *tm1Disp, fvector_t *tm2Disp,
        fvector_t  *damping_vector_shear, fvector_t *damping_vector_kappa){

    // *damping_vector_shear, fvector_t *damping_vector_kappa: Transformed displacements according to the bkt approach
    // fvector_t  *tm1Disp, fvector_t *tm2Disp:                Current and previous displacements

    double rmax, cdt;

    double g0, g02, cg0, eg0, g0k, g02k, cg0k, eg0k;
    double g1, g12, cg1, eg1, g1k, g12k, cg1k, eg1k;
    double g2, g22, cg2, eg2, g2k, g22k, cg2k, eg2k;


    if (typeOfDamping == BKT) {
        cdt = 2. * M_PI * theFreq * theDeltaT;
    } else {
        cdt = theDeltaT;
    }

    // SHEAR RELATED CONVOLUTION
    g0  = cdt * edata->g0_shear;
    g02 = g0 / 2.;
    cg0 = g02 * ( 1. - g0 );
    eg0 = exp( -g0 );

    g1  = cdt * edata->g1_shear;
    g12 = g1 / 2.;
    cg1 = g12 * ( 1. - g1 );
    eg1 = exp( -g1 );

    g0k  = cdt * edata->g0_kappa;
    g02k = g0k / 2.;
    cg0k = g02k * ( 1. - g0k );
    eg0k = exp( -g0k );

    g1k  = cdt * edata->g1_kappa;
    g12k = g1k / 2.;
    cg1k = g12k * ( 1. - g1k );
    eg1k = exp( -g1k );

    if (typeOfDamping >= BKT3) {
        g2  = cdt * edata->g2_shear;
        g22 = g2 / 2.;
        cg2 = g22 * ( 1. - g2 );
        eg2 = exp( -g2 );

        g2k  = cdt * edata->g2_kappa;
        g22k = g2k / 2.;
        cg2k = g22k * ( 1. - g2k );
        eg2k = exp( -g2k );
    }

    fvector_t  *f0_tm1, *f1_tm1, *f2_tm1;

    // SHEAR RELATED CONVOLUTION
    if ( (edata->g0_shear != 0) && (edata->g1_shear != 0) ) {
        f0_tm1 = convShear1;
        f1_tm1 = convShear2;

        f0_tm1->f[0] = cg0 * tm1Disp->f[0] + g02 * tm2Disp->f[0] + eg0 * f0_tm1->f[0];
        f0_tm1->f[1] = cg0 * tm1Disp->f[1] + g02 * tm2Disp->f[1] + eg0 * f0_tm1->f[1];
        f0_tm1->f[2] = cg0 * tm1Disp->f[2] + g02 * tm2Disp->f[2] + eg0 * f0_tm1->f[2];

        f1_tm1->f[0] = cg1 * tm1Disp->f[0] + g12 * tm2Disp->f[0] + eg1 * f1_tm1->f[0];
        f1_tm1->f[1] = cg1 * tm1Disp->f[1] + g12 * tm2Disp->f[1] + eg1 * f1_tm1->f[1];
        f1_tm1->f[2] = cg1 * tm1Disp->f[2] + g12 * tm2Disp->f[2] + eg1 * f1_tm1->f[2];

        if (typeOfDamping >= BKT3) {
            f2_tm1 = convShear3;
            f2_tm1->f[0] = cg2 * tm1Disp->f[0] + g22 * tm2Disp->f[0] + eg2 * f2_tm1->f[0];
            f2_tm1->f[1] = cg2 * tm1Disp->f[1] + g22 * tm2Disp->f[1] + eg2 * f2_tm1->f[1];
            f2_tm1->f[2] = cg2 * tm1Disp->f[2] + g22 * tm2Disp->f[2] + eg2 * f2_tm1->f[2];
        }
    }

    // DILATATION RELATED CONVOLUTION
    if ( (edata->g0_kappa != 0) && (edata->g1_kappa != 0) ) {
        f0_tm1 = convKappa1;
        f1_tm1 = convKappa2;

        f0_tm1->f[0] = cg0k * tm1Disp->f[0] + g02k * tm2Disp->f[0] + eg0k * f0_tm1->f[0];
        f0_tm1->f[1] = cg0k * tm1Disp->f[1] + g02k * tm2Disp->f[1] + eg0k * f0_tm1->f[1];
        f0_tm1->f[2] = cg0k * tm1Disp->f[2] + g02k * tm2Disp->f[2] + eg0k * f0_tm1->f[2];

        f1_tm1->f[0] = cg1k * tm1Disp->f[0] + g12k * tm2Disp->f[0] + eg1k * f1_tm1->f[0];
        f1_tm1->f[1] = cg1k * tm1Disp->f[1] + g12k * tm2Disp->f[1] + eg1k * f1_tm1->f[1];
        f1_tm1->f[2] = cg1k * tm1Disp->f[2] + g12k * tm2Disp->f[2] + eg1k * f1_tm1->f[2];

        if (typeOfDamping >= BKT3) {
            f2_tm1 = convKappa3;
            f2_tm1->f[0] = cg2k * tm1Disp->f[0] + g22k * tm2Disp->f[0] + eg2k * f2_tm1->f[0];
            f2_tm1->f[1] = cg2k * tm1Disp->f[1] + g22k * tm2Disp->f[1] + eg2k * f2_tm1->f[1];
            f2_tm1->f[2] = cg2k * tm1Disp->f[2] + g22k * tm2Disp->f[2] + eg2k * f2_tm1->f[2];
        }
    }

    // =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
    // =*=*=*=*=*=*=*=*=*=*=*= get transformed displacements =*=*=*=*=*=*=*=*=*=*=*
    // =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

    double a0_shear, a1_shear, a2_shear,
    a0_kappa, a1_kappa, a2_kappa,
    b_shear, b_kappa,
    csum, csumk;

    if (typeOfDamping == BKT) {
        rmax = 2. * M_PI * theFreq * theDeltaT;
    } else {
        rmax = theDeltaT;
    }

    a0_shear = edata->a0_shear;
    a1_shear = edata->a1_shear;
    a2_shear = edata->a2_shear;
    b_shear  = edata->b_shear;

    a0_kappa   = edata->a0_kappa;
    a1_kappa   = edata->a1_kappa;
    a2_kappa   = edata->a2_kappa;
    b_kappa    = edata->b_kappa;

    csum = a0_shear + a1_shear + b_shear;
    csumk = a0_kappa + a1_kappa + b_kappa;

    if ( typeOfDamping >= BKT3 ) {
        csum  += a2_shear;
        csumk += a2_kappa;
    }

    double coef_shear = b_shear / rmax;
    double coef_kappa = b_kappa / rmax;

    if ( csum != 0 ) {
        f0_tm1  = convShear1;
        f1_tm1  = convShear2;

        damping_vector_shear->f[0] = coef_shear * (tm1Disp->f[0] - tm2Disp->f[0])
                                   - (a0_shear * f0_tm1->f[0] + a1_shear * f1_tm1->f[0])
                                   + tm1Disp->f[0];

        damping_vector_shear->f[1] = coef_shear * (tm1Disp->f[1] - tm2Disp->f[1])
                                   - (a0_shear * f0_tm1->f[1] + a1_shear * f1_tm1->f[1])
                                   + tm1Disp->f[1];

        damping_vector_shear->f[2] = coef_shear * (tm1Disp->f[2] - tm2Disp->f[2])
                                   - (a0_shear * f0_tm1->f[2] + a1_shear * f1_tm1->f[2])
                                   + tm1Disp->f[2];

        if ( typeOfDamping >= BKT3 ) {
            f2_tm1  = convShear3;
            damping_vector_shear->f[0] -= a2_shear * f2_tm1->f[0];
            damping_vector_shear->f[1] -= a2_shear * f2_tm1->f[1];
            damping_vector_shear->f[2] -= a2_shear * f2_tm1->f[2];
        }
    } else {
        damping_vector_shear->f[0] = tm1Disp->f[0];
        damping_vector_shear->f[1] = tm1Disp->f[1];
        damping_vector_shear->f[2] = tm1Disp->f[2];
    }

    if ( csumk != 0 ) {
        f0_tm1  = convKappa1;
        f1_tm1  = convKappa2;

        damping_vector_kappa->f[0] = coef_kappa * (tm1Disp->f[0] - tm2Disp->f[0])
                                   - (a0_kappa * f0_tm1->f[0] + a1_kappa * f1_tm1->f[0])
                                   + tm1Disp->f[0];

        damping_vector_kappa->f[1] = coef_kappa * (tm1Disp->f[1] - tm2Disp->f[1])
                                   - (a0_kappa * f0_tm1->f[1] + a1_kappa * f1_tm1->f[1])
                                   + tm1Disp->f[1];

        damping_vector_kappa->f[2] = coef_kappa * (tm1Disp->f[2] - tm2Disp->f[2])
                                   - (a0_kappa * f0_tm1->f[2] + a1_kappa * f1_tm1->f[2])
                                   + tm1Disp->f[2];

        if ( typeOfDamping >= BKT3 ) {
            f2_tm1  = convKappa3;
            damping_vector_kappa->f[0] -= a2_kappa * f2_tm1->f[0];
            damping_vector_kappa->f[1] -= a2_kappa * f2_tm1->f[1];
            damping_vector_kappa->f[2] -= a2_kappa * f2_tm1->f[2];
        }
    } else {

        damping_vector_kappa->f[0] = tm1Disp->f[0];
        damping_vector_kappa->f[1] = tm1Disp->f[1];
        damping_vector_kappa->f[2] = tm1Disp->f[2];

    }

    return;

}

