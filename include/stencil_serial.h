/* -*- Mode: C; c-basic-offset:4 ; indent-tabs-mode:nil ; -*- */
/*
 * See COPYRIGHT in top-level directory.
 */

#define _GNU_SOURCE

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <getopt.h>
#include <time.h>
#include <float.h>
#include <math.h>

struct timespec start, end;
double elapsed_time_ns;


#define NORTH 0
#define SOUTH 1
#define EAST  2
#define WEST  3

#define SEND 0
#define RECV 1

#define OLD 0
#define NEW 1

#define _x_ 0
#define _y_ 1

typedef unsigned int uint;



//? ======================= function prototypes =======================

int initialize (int, char**, int*, int*, int*, int*, int**, double*, double**, int*, int*, int*);

int memory_release ( double*, int* );

extern int inject_energy ( const int, const int, const int*, const double, const int[2], double*);

extern int update_plane (const int, const int[2], const double*, double* );


extern int get_total_energy( const int[2], const double*, double* );


//? ======================= function definition for inline functions =======================

inline int inject_energy ( 
        const int       periodic,
        const int       Nsources,
        const int      *Sources,
        const double    energy,
        const int       mysize[2],
        double         *plane
    ) {

    #define IDX( i, j ) ( (j)*(mysize[_x_]+2) + (i) )
    int lim = 2*Nsources; 
        for (int s = 0; s < lim; s+=2) { // optimized
            
            int x = Sources[s];
            int y = Sources[s+1];
            plane[IDX(x, y)] += energy;

            if ( periodic )
                {
                    if ( x == 1 )
                        plane[IDX(mysize[_x_]+1, y)] += energy;
                    if ( x == mysize[_x_] )
                        plane[IDX(0, y)] += energy;
                    if ( y == 1 )
                        plane[IDX(x, mysize[_y_]+1)] += energy;
                    if ( y == mysize[_y_] )
                        plane[IDX(x, 0)] += energy;
                }
        }
    #undef IDX
    
    return 0;
}

inline int update_plane (const int periodic, const int size[2], const double *old, double *new) {
    /*
    * calculate the new energy values
    * the old plane contains the current data, 
    * the new plane will store the updated data
    *
    * NOTE: in parallel, every MPI task will perform the
    *       calculation for its patch
    */

    register const int fxsize = size[_x_]+2;
    register const int xsize = size[_x_];
    register const int ysize = size[_y_];
    
   #define IDX( i, j ) ( (j)*fxsize + (i) )

    // HINT: you may attempt to
    //       (i)  manually unroll the loop
    //       (ii) ask the compiler to do it
    // for instance
    // #pragma GCC unroll 4
    //
    // HINT: in any case, this loop is a good candidate
    //       for openmp parallelization

    //? ······················ Optimized loop ······················
    double alpha = 0.6;
    double constant =  (1-alpha) / 4.0;
    double result;
    
    for (int j = 1; j <= ysize; j++)
        for ( int i = 1; i <= xsize; i++)
            {
                //
                // five-points stencil formula
                //

                // simpler stencil with no explicit diffusivity
                // always conserve the smoohed quantity
                // alpha here mimics how much "easily" the heat
                // travels
                
                result = old[ IDX(i,j) ] *alpha + (old[IDX(i-1, j)] + old[IDX(i+1, j)] + old[IDX(i, j-1)] + old[IDX(i, j+1)]) * constant;

                new[ IDX(i,j) ] = result;
                
            }

    if ( periodic )
        /*
         * propagate boundaries if they are periodic
         *
         * NOTE: when is that needed in distributed memory, if any?
         */
        {
            for ( int i = 1; i <= xsize; i++ )
                {
                    new[ IDX(i, 0) ] = new[ IDX(i, ysize) ];
                    new[ IDX(i, ysize+1) ] = new[ IDX(i, 1) ];
                }
            for ( int j = 1; j <= ysize; j++ )
                {
                    new[ IDX( 0, j) ] = new[ IDX(xsize, j) ];
                    new[ IDX( xsize+1, j) ] = new[ IDX(1, j) ];
                }
        }
    
    return 0;

   #undef IDX
}

inline int get_total_energy( const int size[2], const double *plane, double *energy ) {
    /*
    * NOTE: this routine a good candidate for openmp parallelization
    */

    register const int ysize = size[_y_]; // optimization
    register const int xsize = size[_x_]; // read once
    
    #define IDX( i, j ) ( (j)*(xsize+2) + (i) )

        #if defined(LONG_ACCURACY)    
            long double totenergy = 0;
        #else
            double totenergy = 0;    
        #endif

        // HINT: you may attempt to
        //       (i)  manually unroll the loop
        //       (ii) ask the compiler to do it
        // for instance

        #pragma GCC unroll 4
        for ( int j = 1; j <= ysize; j++ )
            for ( int i = 1; i <= xsize; i++ )
                totenergy += plane[ IDX(i, j) ];
    
    #undef IDX

    *energy = (double)totenergy;
    return 0;
}
                            
