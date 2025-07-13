#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <getopt.h>
#include <time.h>
#include <math.h>

#include <omp.h>
#include <mpi.h>

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

// tags for MPI communication
#define TAG_N 0
#define TAG_S 1
#define TAG_E 2
#define TAG_W 3


typedef unsigned int uint;

typedef uint    vec2_t[2];
typedef double *restrict buffers_t[4];

typedef struct {
    double   * restrict data;
    vec2_t     size;
} plane_t;

//? ======================= function prototypes =======================

extern int inject_energy (const int, const int, const vec2_t*, const double, plane_t*, const vec2_t);

extern int update_plane ( const int, const vec2_t, const plane_t*, plane_t*);

extern int get_total_energy(plane_t*, double*);

int initialize (MPI_Comm *, int, int, int, char**, vec2_t*, vec2_t*, int*, int*, int*, int*, int*, int*, vec2_t  **, double   *, plane_t  *, buffers_t * );

int memory_release (const int *, buffers_t *, plane_t *, vec2_t **);

int output_energy_stat (int, plane_t *, double, int, MPI_Comm *);


//? ======================= function definition for inline functions =======================

inline int inject_energy ( 
        const int       periodic,
        const int       Nsources,
        const vec2_t    *Sources,
        const double    energy,
        plane_t         *plane,
        const vec2_t    N
    ) {

    const uint xsize = plane->size[_x_]+2;
    double * restrict data = plane->data;
    
    #define IDX( i, j ) ( (j)*xsize + (i) )
        for (int s = 0; s < Nsources; s++) {
            
            uint x = Sources[s][_x_];
            uint y = Sources[s][_y_];
            
            
            data[ IDX(x,y) ] += energy;

            if ( periodic ) {
                if ( (N[_x_] == 1)  ) { 
                    // in this case there is only a column of tasks
                    // we proceed as in the serial version
                    data[IDX(0, y)]                    += data[IDX(plane->size[_x_]+1, y)]; // West from East
                    data[IDX(plane->size[_x_]+1, y)]   += data[IDX(1, y)];                 // East from West
                }
                
                if ( (N[_y_] == 1) ) {
                    // in this case there is only a row of tasks
                    // we proceed as in the serial version  
                    data[IDX(x, 0)]                    += data[IDX(x, plane->size[_y_]+1)]; // North from South
                    data[IDX(x, plane->size[_y_]+1)]   += data[IDX(x, 1)];                 // South from North
                }
            }                
        }
    #undef IDX
        
    return 0;
}

inline int update_plane ( 
        const int       periodic, 
        const vec2_t    N,         // the grid of MPI tasks
        const plane_t  *oldplane,
        plane_t        *newplane
    ) {

    const uint xsize = oldplane->size[_x_];
    const uint ysize = oldplane->size[_y_];

    const uint fxsize = xsize+2;
    
    #define IDX( i, j ) ( (j)*fxsize + (i) )
    
        // HINT: you may attempt to
        //       (i)  manually unroll the loop
        //       (ii) ask the compiler to do it
        // for instance
        // #pragma GCC unroll 4
        //
        // HINT: in any case, this loop is a good candidate
        //       for openmp parallelization

        double * restrict old = oldplane->data;
        double * restrict new = newplane->data;

        double alpha = 0.6;
        double constant =  (1-alpha) / 4.0;
        double result;
        
        for (uint j = 1; j <= ysize; j++) {
            for ( uint i = 1; i <= xsize; i++)
                {
                    // NOTE: (i-1,j), (i+1,j), (i,j-1) and (i,j+1) always exist even
                    //       if this patch is at some border without periodic conditions;
                    //       in that case it is assumed that the +-1 points are outside the
                    //       plate and always have a value of 0, i.e. they are an
                    //       "infinite sink" of heat
                    
                    // five-points stencil formula

                    result = old[ IDX(i,j) ] * alpha + (old[IDX(i-1, j)] + old[IDX(i+1, j)] + old[IDX(i, j-1)] + old[IDX(i, j+1)]) * constant;

                    new[ IDX(i,j) ] = result;
                }
        }

        if ( periodic ) {
            // if there is only a column of tasks, the periodicity on the X axis is local
            if ( N[_x_] == 1 ) {
                // copy the values of the first column to the right ghost column (xsize+1)
                // and the values of the last column to the left ghost column (0)
                for ( uint j = 1; j <= ysize; j++ ) {
                    new[ IDX( 0, j) ]       = new[ IDX(xsize, j) ];
                    new[ IDX( xsize+1, j) ] = new[ IDX(1, j) ];
                }
            }
    
            // if there is only a row of tasks, the periodicity on the Y axis is local
            if ( N[_y_] == 1 ) {
                // copy the values of the first row to the bottom ghost row (ysize+1)
                // and the values of the last row to the top ghost row (0)
                for ( uint i = 1; i <= xsize; i++ ) {
                    new[ IDX( i, 0 ) ]       = new[ IDX(i, ysize) ];
                    new[ IDX( i, ysize+1) ] = new[ IDX(i, 1) ];
                }
            }
        }
    
    #undef IDX
    return 0;
}



inline int get_total_energy( plane_t *plane, double  *energy ) {
    /*
    * NOTE: this routine a good candiadate for openmp
    *                   parallelization
    */

    const uint xsize = plane->size[_x_];
    const uint ysize = plane->size[_y_];
    const uint fsize = xsize+2;

    double * restrict data = plane->data;
    
    #define IDX( i, j ) ( (j)*fsize + (i) )

    #if defined(LONG_ACCURACY)    
        long double totenergy = 0;
    #else
        double totenergy = 0;    
    #endif

    // HINT: you may attempt to
    //       (i)  manually unroll the loop
    //       (ii) ask the compiler to do it
    // for instance
    // #pragma GCC unroll 4
    for ( uint j = 1; j <= ysize; j++ )
        for ( uint i = 1; i <= xsize; i++ )
            totenergy += data[ IDX(i, j) ];

    
    #undef IDX

    *energy = (double)totenergy;
    return 0;
}



