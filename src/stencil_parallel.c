#include "stencil_parallel.h"

int verbose = 0;
int seed = 0; // default is random itself

//? ---------------------------------------------------------------
//?                          Main function 
//? ---------------------------------------------------------------

int main(int argc, char **argv) {

	MPI_Comm myCOMM_WORLD;      // MPI communicator
	int  Rank, Ntasks;			// Rank: the rank of the calling process, Ntasks: the total number of MPI ranks

	int neighbours[4];       	// 0: North, 1: South, 2: East, 3: West
	int  Niterations;          	// number of iterations
	int  periodic;             	// periodic boundary condition
	vec2_t S, N;                // size of the plane and size of MPI tasks

	int      Nsources;          // number of heat sources
	int      Nsources_local;    // number of heat sources on the local task
	vec2_t  *Sources_local;     // coordinates of the heat sources on the local task
	double   energy_per_source; // energy per source

	plane_t   planes[2];   		// two planes for the two iterations (current and next)
	buffers_t buffers[2]; 		// two buffers for the two iterations (send and receive)
  
	int output_energy_stat_perstep; // whether to output energy statistics per step
	
	//? ························ initialize MPI envrionment ························
	{
		int level_obtained;
		
		// NOTE: change MPI_FUNNELED if appropriate
		//
		MPI_Init_thread( &argc, &argv, MPI_THREAD_FUNNELED, &level_obtained );
		if ( level_obtained < MPI_THREAD_FUNNELED ) {
			printf("MPI_thread level obtained is %d instead of %d\n", level_obtained, MPI_THREAD_FUNNELED );
			MPI_Finalize();
			exit(1); }
		
		MPI_Comm_rank(MPI_COMM_WORLD, &Rank);
		MPI_Comm_size(MPI_COMM_WORLD, &Ntasks);
		MPI_Comm_dup (MPI_COMM_WORLD, &myCOMM_WORLD);
	}
  
	//? ························ argument checking and setting ························


	int ret = initialize (&myCOMM_WORLD, Rank, Ntasks, argc, argv, &S, &N, &periodic, &output_energy_stat_perstep,neighbours,
					  &Niterations, &Nsources, &Nsources_local, &Sources_local, &energy_per_source, planes, buffers);

	if ( ret ) {
		printf("task %d is opting out with termination code %d\n", Rank, ret );
		
		MPI_Finalize();
		return 0;
    }

	//? -------------------------------- Core loop --------------------------------
	
	int current = OLD;
	double t1 = MPI_Wtime();   /* take wall-clock time */
	
	for (int iter = 0; iter < Niterations; ++iter) {
      
		MPI_Request reqs[8];
		int nreqs = 0;
		
		//? - - - - - - - - - - inject energy from new sources - - - - - - - - - - -

		inject_energy( periodic, Nsources_local, Sources_local, energy_per_source, &planes[current], N );

		//?- - - - - - - - - - - - - prepare the buffers - - - - - - - - - - - - - -

		uint ysize = planes[current].size[_y_];
		uint xsize = planes[current].size[_x_];

		uint xframe = xsize + 2;

		// Checking if buffers are allocated
		if (buffers[SEND][WEST] != NULL && buffers[SEND][EAST] != NULL) {
			for (uint i = 0; i < ysize; i++) {
				// WEST: first effective column (excluding frame)
				buffers[SEND][WEST][i] = planes[current].data[(i + 1) * xframe + 1];
				// EAST: last effective column (excluding frame)
				buffers[SEND][EAST][i] = planes[current].data[(i + 1) * xframe + xsize];
			}
		}

		buffers[SEND][NORTH] = &(planes[current].data[xframe + 1]); 		// the first effective row
		buffers[SEND][SOUTH] = &(planes[current].data[ysize * xframe + 1]); // the last effective row
		buffers[RECV][NORTH] = &(planes[current].data[1]);
		buffers[RECV][SOUTH] = &(planes[current].data[(ysize + 1) * xframe + 1]);
		
		//? - - - - - - - - - - - perform the halo communications - - - - - - - - - -
		
		//     (1) use Send / Recv
		//     (2) use Isend / Irecv
		//         --> can you overlap communication and compution in this way?

		if (neighbours[EAST] != MPI_PROC_NULL) {
			MPI_Isend(buffers[SEND][EAST], (int)ysize, MPI_DOUBLE, neighbours[EAST], TAG_E, myCOMM_WORLD, &reqs[nreqs++]);
			MPI_Irecv(buffers[RECV][EAST], (int)ysize, MPI_DOUBLE, neighbours[EAST], TAG_W, myCOMM_WORLD, &reqs[nreqs++]);
		}
		if (neighbours[WEST] != MPI_PROC_NULL) {
			MPI_Isend(buffers[SEND][WEST], (int)ysize, MPI_DOUBLE, neighbours[WEST], TAG_W, myCOMM_WORLD, &reqs[nreqs++]);
			MPI_Irecv(buffers[RECV][WEST], (int)ysize, MPI_DOUBLE, neighbours[WEST], TAG_E, myCOMM_WORLD, &reqs[nreqs++]);
		}
		if (neighbours[NORTH] != MPI_PROC_NULL) {
			MPI_Isend(buffers[SEND][NORTH], (int)xsize, MPI_DOUBLE, neighbours[NORTH], TAG_N, myCOMM_WORLD, &reqs[nreqs++]);
			MPI_Irecv(buffers[RECV][NORTH], (int)xsize, MPI_DOUBLE, neighbours[NORTH], TAG_S, myCOMM_WORLD, &reqs[nreqs++]);
		}
		if (neighbours[SOUTH] != MPI_PROC_NULL) {
			MPI_Isend(buffers[SEND][SOUTH], (int)xsize, MPI_DOUBLE, neighbours[SOUTH], TAG_S, myCOMM_WORLD, &reqs[nreqs++]);
			MPI_Irecv(buffers[RECV][SOUTH], (int)xsize, MPI_DOUBLE, neighbours[SOUTH], TAG_N, myCOMM_WORLD, &reqs[nreqs++]);
		}
		
		MPI_Waitall(nreqs, reqs, MPI_STATUSES_IGNORE);
		
		//? - - - - - - - - - - - - - - copy the haloes data - - - - - - - - - - - - - -

		if (neighbours[WEST] != MPI_PROC_NULL) {
			for (uint i = 0; i < ysize; i++) {
				planes[current].data[(i + 1) * xframe + 0] = buffers[RECV][WEST][i];
			}
		}
		if (neighbours[EAST] != MPI_PROC_NULL) {
			for (uint i = 0; i < ysize; i++) {
				planes[current].data[(i + 1) * xframe + (xsize + 1)] = buffers[RECV][EAST][i];
			}
		}

		//? - - - - - - - - - - - - - - - update grid points - - - - - - - - - - - - - -

		update_plane( periodic, N, &planes[current], &planes[!current] );


		//! ------------------------- Print the surface at each step -------------------------
		//!
		//! warning:   this is a very slow operation, and should be used only for debugging
		//!  		               since it looses the parallelization
		//!
		//! ----------------------------------------------------------------------------------
		if (verbose > 1) {
			// Gather local grid sizes and data to rank 0
			uint p_xsize = planes[!current].size[_x_];
			uint p_ysize = planes[!current].size[_y_];
			uint p_xframe = p_xsize + 2;
			size_t local_size = p_xsize * p_ysize;

			// Get the global grid decomposition
			int Grid_x, Grid_y;
			MPI_Comm_size(myCOMM_WORLD, &Ntasks);
			MPI_Comm_rank(myCOMM_WORLD, &Rank);
			Grid_x = (int)N[_x_];
			Grid_y = (int)N[_y_];

			// Special case: Ntasks == 1, just print the local grid directly
			if (Ntasks == 1) {
				printf("Merged Surface at step %d:\n\n", iter);
				for (uint j = 0; j < p_ysize; j++) {
					for (uint i = 0; i < p_xsize; i++) {
						printf("%8.4f ", planes[!current].data[(j+1) * p_xframe + (i+1)]);
					}
					printf("\n");
				}
				printf("\n");
			} else {
				// Each rank sends its local grid size and data to rank 0
				int *all_xsize = NULL, *all_ysize = NULL;
				int *displs = NULL, *recvcounts = NULL;

				if (Rank == 0) {
					all_xsize = (int*)malloc((size_t)Ntasks * sizeof(int));
					all_ysize = (int*)malloc((size_t)Ntasks * sizeof(int));
				}
				int my_sizes[2] = {(int)p_xsize, (int)p_ysize};
				MPI_Gather(&my_sizes[0], 1, MPI_INT, all_xsize, 1, MPI_INT, 0, myCOMM_WORLD);
				MPI_Gather(&my_sizes[1], 1, MPI_INT, all_ysize, 1, MPI_INT, 0, myCOMM_WORLD);

				// Gather all local data to rank 0
				double *localbuf = (double*)malloc(local_size * sizeof(double));
				for (uint j = 0; j < p_ysize; j++)
					for (uint i = 0; i < p_xsize; i++)
						localbuf[j * p_xsize + i] = planes[!current].data[(j+1) * p_xframe + (i+1)];

				// Compute displacements and receive counts for Gatherv
				if (Rank == 0) {
					recvcounts = (int*)malloc((size_t)Ntasks * sizeof(int));
					displs = (int*)malloc((size_t)Ntasks * sizeof(int));
					int offset = 0;
					for (int r = 0; r < Ntasks; r++) {
						recvcounts[r] = all_xsize[r] * all_ysize[r];
						displs[r] = offset;
						offset += recvcounts[r];
					}
				}
				double *gathered = NULL;
				if (Rank == 0) {
					int total = 0;
					for (int r = 0; r < Ntasks; r++) total += all_xsize[r] * all_ysize[r];
					gathered = (double*)malloc((size_t)total * sizeof(double));
				}
				MPI_Gatherv(localbuf, (int)local_size, MPI_DOUBLE,
				            gathered, recvcounts, displs, MPI_DOUBLE, 0, myCOMM_WORLD);
				free(localbuf);

				// Print the merged grid at rank 0
				if (Rank == 0) {
					// Compute the global grid size
					int global_x = 0, global_y = 0;
					for (int gx = 0; gx < Grid_x; gx++) global_x += all_xsize[gx];
					for (int gy = 0; gy < Grid_y; gy++) global_y += all_ysize[gy * Grid_x];

					// Build a map of where each rank's patch goes
					int *x_offsets = (int*)malloc((size_t)Grid_x * sizeof(int));
					int *y_offsets = (int*)malloc((size_t)Grid_y * sizeof(int));
					x_offsets[0] = 0;
					y_offsets[0] = 0;
					for (int gx = 1; gx < Grid_x; gx++)
						x_offsets[gx] = x_offsets[gx-1] + all_xsize[gx-1];
					for (int gy = 1; gy < Grid_y; gy++)
						y_offsets[gy] = y_offsets[gy-1] + all_ysize[(gy-1)*Grid_x];

					// Allocate the global grid
					double **global_grid = (double**)malloc((size_t)global_y * sizeof(double*));
					for (int j = 0; j < global_y; j++)
						global_grid[j] = (double*)malloc((size_t)global_x * sizeof(double));
					// Fill with zeros
					for (int j = 0; j < global_y; j++)
						for (int i = 0; i < global_x; i++)
							global_grid[j][i] = 0.0;

					// Place each rank's patch in the global grid
					int *patch_idx = (int*)malloc((size_t)Ntasks * sizeof(int));
					for (int r = 0; r < Ntasks; r++) patch_idx[r] = 0;
					for (int gy = 0; gy < Grid_y; gy++) {
						for (int gx = 0; gx < Grid_x; gx++) {
							int r = gy * Grid_x + gx;
							int px = x_offsets[gx];
							int py = y_offsets[gy];
							int xs = all_xsize[r];
							int ys = all_ysize[r];
							int base = displs[r];
							for (int j = 0; j < ys; j++) {
								for (int i = 0; i < xs; i++) {
									global_grid[py + j][px + i] = gathered[base + j * xs + i];
								}
							}
						}
					}

					// Enhanced pretty-printing of the merged surface with improved grid lines and intersections
					// Now using UTF-8 string literals for box-drawing characters to avoid multi-character constant warnings

					printf("Merged Surface at step %d:\n\n", iter);

					// Precompute vertical separator columns for each region
					int *vert_seps = (int*)malloc(((size_t)Grid_x-1) * sizeof(int));
					int sep_sum = 0;
					for (int gx = 0; gx < Grid_x-1; gx++) {
						sep_sum += all_xsize[gx];
						vert_seps[gx] = sep_sum;
					}

					// Precompute horizontal separator rows for each region
					int *horiz_seps = (int*)malloc(((size_t)Grid_y-1) * sizeof(int));
					sep_sum = 0;
					for (int gy = 0; gy < Grid_y-1; gy++) {
						sep_sum += all_ysize[gy * Grid_x];
						horiz_seps[gy] = sep_sum;
					}

					for (int j = 0; j < global_y; j++) {
						for (int i = 0; i < global_x; i++) {
							printf("%8.4f", global_grid[j][i]);
							// Add vertical separator between task regions, aligned
							for (int s = 0; s < Grid_x-1; s++) {
								if (i + 1 == vert_seps[s]) {
									printf(" │");
									break;
								}
							}
						}
						printf("\n");
						// Add horizontal separator between task regions, aligned
						for (int s = 0; s < Grid_y-1; s++) {
							if (j + 1 == horiz_seps[s]) {
								// Print left margin to align with numbers
								for (int i = 0; i < global_x; i++) {
									if (i == 0) printf("  ───────");
									else printf("────────");
									// printf("──────");
									// Add intersection or end
									for (int t = 0; t < Grid_x-1; t++) {
										if (i + 1 == vert_seps[t]) {
											printf("┼─");
											break;
										}
									}
								}
								printf("\n");
								break;
							}
						}
					}
					printf("\n");

					// Clean up
					for (int j = 0; j < global_y; j++) free(global_grid[j]);
					free(global_grid);
					free(x_offsets);
					free(y_offsets);
					free(patch_idx);
					free(all_xsize);
					free(all_ysize);
					free(recvcounts);
					free(displs);
					free(gathered);
					free(vert_seps);
					free(horiz_seps);
				}
				MPI_Barrier(myCOMM_WORLD);
			}
		}
		//! -------------------------------- End of printing ---------------------------------

		//? - - - - - - - - - - - - - - output energy statistics - - - - - - - - - - - - - - -
		if ( output_energy_stat_perstep )
			output_energy_stat ( iter, &planes[!current], (iter+1) * Nsources*energy_per_source, Rank, &myCOMM_WORLD );
	
		// swap plane indexes for the new iteration
		current = !current;
      
    }
  
	t1 = MPI_Wtime() - t1;

	if (verbose > 0) printf("Time taken: %f seconds\n", t1);

	output_energy_stat ( -1, &planes[current], Niterations * Nsources*energy_per_source, Rank, &myCOMM_WORLD );
	
	memory_release( buffers, planes, Sources_local );
	
	
	MPI_Finalize();
	return 0;
}


 
//? -------------------------------------------------------------------------
//?                              initialization
//? -------------------------------------------------------------------------

int simple_factorization(uint, int*, uint**);

int initialize_sources(int, MPI_Comm  *, uint[2], int, int*, vec2_t**, vec2_t, vec2_t);

int memory_allocate (const int*, vec2_t, buffers_t*, plane_t*);

int initialize ( 
	MPI_Comm *Comm,               // MPI communicator	
	int      Me,                  // the rank of the calling process
	int      Ntasks,              // the total number of MPI ranks
	int      argc,                // the argc from command line
	char   **argv,                // the argv from command line
	vec2_t  *S,                   // the size of the plane
	vec2_t  *N,                   // two-uint array defining the MPI tasks' grid
	int     *periodic,            // periodic-boundary tag
	int     *output_energy_stat,  // whether to output energy statistics
	int     *neighbours,          // four-int array that gives back the neighbours of the calling task
	int     *Niterations,         // how many iterations
	int     *Nsources,            // how many heat sources
	int     *Nsources_local,      // number of heat sources on the local task
	vec2_t **Sources_local,       // coordinates of the heat sources on the local task
	double  *energy_per_source,   // how much heat per source
	plane_t *planes,              // plane data
	buffers_t *buffers            // buffers for communication
	//int     *injection_frequency  // frequency of injection of energy
) {
	
	int ret;

	//? ················ set deffault values ················

	(*S)[_x_]         = 10000;
	(*S)[_y_]         = 10000;
	*periodic         = 0;
	*Nsources         = 4;
	*Nsources_local   = 0;
	*Sources_local    = NULL;
	*Niterations      = 1000;
	*energy_per_source = 1.0;
	*output_energy_stat = 0;

	if ( planes == NULL ) {
			// manage the situation
	}

	planes[OLD].size[0] = planes[OLD].size[1] = 0;
	planes[NEW].size[0] = planes[NEW].size[1] = 0;
  
  	for ( int i = 0; i < 4; i++ ) 
		neighbours[i] = MPI_PROC_NULL;


	// set random seed according to time
	seed = (int)time(NULL);
	
	//? ··················· process the commadn line ·······················
	
	while ( 1 ) {
		int opt;
		while((opt = getopt(argc, argv, ":hx:y:e:E:n:o:p:v:s:")) != -1) {
			switch( opt ) {
				case 'x': 
					(*S)[_x_] = (uint)atoi(optarg);
					break;
				
				case 'y': 
					(*S)[_y_] = (uint)atoi(optarg);
					break;
				
				case 'e': 
					*Nsources = atoi(optarg);
					break;
				
				case 'E': 
					*energy_per_source = atof(optarg);
					break;
				
				case 'n':
					*Niterations = atoi(optarg);
					break;
				
				case 'o': 
					*output_energy_stat = (atoi(optarg) > 0);
					break;

				case 'p': 
					*periodic = (atoi(optarg) > 0);
					break;

				case 'v': 
					verbose = atoi(optarg);
					break;

				case 's':
					seed = atoi(optarg);
					break;
				case 'h':
					if ( Me == 0 )
						printf( "\nvalid options are ( values btw [] are the default values ):\n"
								"-x    x size of the plate [10000]\n"
								"-y    y size of the plate [10000]\n"
								"-e    how many energy sources on the plate [4]\n"
								"-E    how many energy sources on the plate [1.0]\n"
								"-n    how many iterations [1000]\n"
								"-p    whether periodic boundaries applies  [0 = false]\n\n"
								"-o    output energy statistics [0 = false]\n"
								"-v    verbose level [0]\n"
								"-h    print this help message\n"
						);
				return 1;
				
				case ':':
					printf( "option -%c requires an argument\n", optopt);
					break;
				
				case '?':
					printf(" -------- help unavailable ----------\n");
					break;
			}
		}
		
		if ( opt == -1 )
		break;
	}

	//? ······················ Check the parameters ······················

	/*
	if ( freq <= 0 )
		*injection_frequency = 1;
	else if ( freq <= 0 ) {
		printf("Error: injection_frequency must be positive\n");
		return 1;
	} else {
		freq = (freq > 1.0 ? 1.0 : freq );
		*injection_frequency = freq * *Niterations;
	}
	*/

	if ( Ntasks <= 0 ) {
		printf("Error: number of tasks must be positive\n");
		return 1;
	} else if ( (uint)Ntasks > (*S)[_x_] * (*S)[_y_] ) {
		printf("Error: number of tasks must be less than or equal to the number of grid points\n");
		return 1;
	}
	
	if ( *Niterations <= 0 ) {
		printf("Error: number of iterations must be positive\n");
		return 1;
	} else if ( *Niterations > 1000000 ) {
		printf("Error: number of iterations must be less than 1000000\n");
		return 1;
	}
	
	if ( *Nsources <= 0 ) {
		printf("Error: number of sources must be positive\n");
		return 1;
	} else if ( (uint)*Nsources > (*S)[_x_] * (*S)[_y_] ) {
		printf("Error: number of sources must be less than or equal to the number of grid points\n");
		return 1;
	}
	
	if ( *energy_per_source <= 0 ) {
		printf("Error: energy per source must be positive\n");
		return 1;
	}
	
	if ( *periodic < 0 || *periodic > 1 ) {
		printf("Error: periodic boundary condition must be 0 or 1\n");
		return 1;
	}
	
	if ( *output_energy_stat < 0 || *output_energy_stat > 1 ) {
		printf("Error: output energy statistics must be 0 or 1\n");
		return 1;
	}
	
	/*
	* find a suitable domain decomposition
	* very simple algorithm, you may want to
	* substitute it with a better one
	*
	* the plane Sx x Sy will be solved with a grid
	* of Nx x Ny MPI tasks
	*/

	vec2_t Grid;
	// formfactor is the ratio of the x and y dimensions of the plane if x > y, or y / x if y > x
	double formfactor = ((*S)[_x_] >= (*S)[_y_] ? (double)(*S)[_x_]/(*S)[_y_] : (double)(*S)[_y_]/(*S)[_x_] );
	int    dimensions = 2 - (Ntasks <= ((int)formfactor+1) );
	
	if ( dimensions == 1 ) {
		// if dimensions is 1, means that we can split the plane in a single row or column
		if ( (*S)[_x_] >= (*S)[_y_] )
			Grid[_x_] = (uint)Ntasks, Grid[_y_] = 1;
		else
			Grid[_x_] = 1, Grid[_y_] = (uint)Ntasks;
    } else {
		// if dimensions is 2, means that we can split the plane in a grid of Nx x Ny tasks
		int   Nf;
		uint *factors = NULL;
		uint  first = 1;
		ret = simple_factorization( (uint)Ntasks, &Nf, &factors );
		if (ret != 0) {
			printf("Error: factorization failed\n");
			return 1;
		}
		
		for ( int i = 0; (i < Nf) && (((uint)Ntasks/first)/first > formfactor); i++ )
			first *= factors[i];
	
		uint factor1 = first;
		uint factor2 = (uint)Ntasks/first;

		if ( (*S)[_x_] >= (*S)[_y_] ) {
			// wide data, make grid wide
			Grid[_x_] = (factor1 > factor2) ? factor1 : factor2;
			Grid[_y_] = (factor1 > factor2) ? factor2 : factor1;
		} else {
			// tall or square data, make grid tall
			Grid[_y_] = (factor1 > factor2) ? factor1 : factor2;
			Grid[_x_] = (factor1 > factor2) ? factor2 : factor1;
		}
	}


	(*N)[_x_] = Grid[_x_];
	(*N)[_y_] = Grid[_y_];


	//? ·········· my cooridnates in the grid of processors ··········
	int X = Me % (int)Grid[_x_];
	int Y = Me / (int)Grid[_x_];

	//? ···················· find my neighbours ······················
	if ( *periodic ) {
		// Horizontal neighbours
		if (Grid[_x_] > 1 || *periodic) {
			neighbours[EAST]  = (int)(((uint)Y)*Grid[_x_] + (uint)(X + 1 ) % Grid[_x_]);
			neighbours[WEST]  = (int)(((uint)Y)*Grid[_x_] + (uint)(X - 1 + (int)Grid[_x_]) % Grid[_x_]);
		}
		// Vertical neighbours
		if (Grid[_y_] > 1 || *periodic) {
			neighbours[NORTH] = (int)(((uint)(Y - 1 + (int)Grid[_y_]) % Grid[_y_]) * Grid[_x_] + (uint)X);
			neighbours[SOUTH] = (int)(((uint)(Y + 1) % Grid[_y_]) * Grid[_x_] + (uint)X);
		}
	} else {
		// Horizontal neighbours
		if ( Grid[_x_] > 1 ) {  
			neighbours[EAST]  = ( X < (int)Grid[_x_]-1 ? Me+1 : MPI_PROC_NULL );
			neighbours[WEST]  = ( X > 0 ? Me-1 : MPI_PROC_NULL ); 
		}
		// Vertical neighbours
		if ( Grid[_y_] > 1 ) {
			neighbours[NORTH] = ( Y > 0 ? Me - (int)Grid[_x_]: MPI_PROC_NULL );
			neighbours[SOUTH] = ( Y < (int)Grid[_y_]-1 ? Me + (int)Grid[_x_] : MPI_PROC_NULL );
		}
	}


	//? ··················· the size of my patch ···················

	/*
	* every MPI task determines the size sx x sy of its own domain
	* REMIND: the computational domain will be embedded into a frame
	*         that is (sx+2) x (sy+2)
	*         the outern frame will be used for halo communication or
	*/

	vec2_t mysize;
	uint s = (*S)[_x_] / Grid[_x_];
	uint r = (*S)[_x_] % Grid[_x_];
	// the last r tasks will have an extra grid point
	mysize[_x_] = s + (X < (int)r);

	// same logic for the y direction
	s = (*S)[_y_] / Grid[_y_];
	r = (*S)[_y_] % Grid[_y_];
	mysize[_y_] = s + (Y < (int)r);

	//? ························ verbose output ···························

	if ( verbose > 0 ) {
		if ( Me == 0 ) {
			printf("Tasks are decomposed in a grid %d x %d\n\n",
				Grid[_x_], Grid[_y_] );
			fflush(stdout);
		}
		
		MPI_Barrier(*Comm);
      
		if (Me == 0) {
			printf("Neighbours:\n\n");
			printf("   Task   N     E     S     W\n");
			printf("  ============================\n");
			fflush(stdout);
		}
		MPI_Barrier(*Comm);
		for (int t = 0; t < Ntasks; t++) {
			if (t == Me) {
				printf("%5d %5d %5d %5d %5d\n",
					Me,
					neighbours[NORTH],
					neighbours[EAST],
					neighbours[SOUTH],
					neighbours[WEST]
				);
				fflush(stdout);
			}
			MPI_Barrier(*Comm);
		}
		if (Me == 0) {
			printf("\n");
			fflush(stdout);
		}
		MPI_Barrier(*Comm);
	}

	//? ···················· initialize the buffers ·······················
	for ( int b = 0; b < 2; b++ ) {
		for ( int d = 0; d < 4; d++ ) {
			buffers[b][d] = NULL;
		}
	}

	//? ···················· allocate the needed memory ···················
	ret = memory_allocate(neighbours, mysize, buffers, planes);

	if (ret != 0) {
		printf("Error: failed to allocate memory for the buffers\n");
		return 1;
	}
		
	//? ··················· allocae the heat sources ······················
	ret = initialize_sources(Me, Comm, mysize, *Nsources, Nsources_local, Sources_local, *S, *N);

	if (ret != 0) {
		printf("Error: failed to initialize the sources\n");
		return 1;
	}

	return 0;  
}

int simple_factorization( uint Ntasks, int *Nfactors, uint **factors ) {
	/*
	* rought factorization;
	* assumes that A is small, of the order of <~ 10^5 max,
	* since it represents the number of tasks
	*/
	
	int N = 0;
	uint f = 2;
	uint _A_ = Ntasks;
	
	uint temp_A = Ntasks;
	while ( f * f <= temp_A ) {
		while( temp_A % f == 0 ) {
			N++;
			temp_A /= f;
		}
		f++;
	}
	if (temp_A > 1) {
		N++;
	}

	*Nfactors = N;
	if (N == 0) {
		*factors = NULL;
		return 0;
	}

	uint *_factors_ = (uint*)malloc( (size_t)N * sizeof(uint) );
	if (_factors_ == NULL) return 1;

	N   = 0;
	f   = 2;
	_A_ = Ntasks;

	while ( f * f <= _A_ ) {
		while( _A_ % f == 0 ) {
			_factors_[N++] = f;
			_A_ /= f;
		}
		f++;
	}
	if (_A_ > 1) {
		_factors_[N++] = _A_;
	}

	*factors = _factors_;
	return 0;
}

int initialize_sources(
    int       Me,
    MPI_Comm *Comm,
    vec2_t    mysize,
    int       Nsources,
    int      *Nsources_local,
    vec2_t  **Sources_local,
    vec2_t    S,
    vec2_t    N
) {
    uint *global_sources_idx = (uint*)malloc((size_t)Nsources * sizeof(uint));
    if (global_sources_idx == NULL) {
        perror("Failed to allocate global_sources_idx");
        return 1;
    }

	//? - - - - - - - Generate the global sources - - - - - - - - - 
	// I decided to opt for this approach to ensure reproducibility
    if (Me == 0) {
        srand48(42); // Fixed seed for reproducibility
        for (int i = 0; i < Nsources; i++) {
            // Generate global coordinates (gx, gy) in the range [0, S-1]
            uint gx = (uint)(lrand48() % S[_x_]);
            uint gy = (uint)(lrand48() % S[_y_]);
            
            // Calculate 1D index based on 0
            global_sources_idx[i] = gy * S[_x_] + gx;
        }
    }

	//? - - - - - - - Broadcast the global sources - - - - - - - - - 
    MPI_Bcast(global_sources_idx, Nsources, MPI_UNSIGNED, 0, *Comm);

	//? - - - - - - - Determine the local sources - - - - - - - - - 
    int my_grid_x = Me % (int)N[_x_];
    int my_grid_y = Me / (int)N[_x_];
    
    uint global_offset_x = 0;
    for (int i = 0; i < my_grid_x; i++) {
        global_offset_x += S[_x_] / N[_x_] + (i < (int)(S[_x_] % N[_x_]));
    }
    uint global_offset_y = 0;
    for (int j = 0; j < my_grid_y; j++) {
        global_offset_y += S[_y_] / N[_y_] + (j < (int)(S[_y_] % N[_y_]));
    }

    int nlocal = 0;
    for (int i = 0; i < Nsources; i++) {
        uint g_idx = global_sources_idx[i];
        // Recover the global coordinates (base-0) from the index
        uint gx = g_idx % S[_x_];
        uint gy = g_idx / S[_x_];

        if (gx >= global_offset_x && gx < global_offset_x + mysize[_x_] &&
            gy >= global_offset_y && gy < global_offset_y + mysize[_y_]) {
            nlocal++;
        }
    }
    *Nsources_local = nlocal;

    if (nlocal > 0) {
        *Sources_local = (vec2_t*)malloc((size_t)nlocal * sizeof(vec2_t));
        if (*Sources_local == NULL) {
            perror("Failed to allocate memory for local sources");
            free(global_sources_idx);
            return 1;
        }

        int current_local_source = 0;
        for (int i = 0; i < Nsources; i++) {
            uint g_idx = global_sources_idx[i];
            uint gx = g_idx % S[_x_];
            uint gy = g_idx / S[_x_];

            if (gx >= global_offset_x && gx < global_offset_x + mysize[_x_] &&
                gy >= global_offset_y && gy < global_offset_y + mysize[_y_]) {
                
                // Convert from global (base-0) to local (base-1 for the framed grid)
                (*Sources_local)[current_local_source][_x_] = gx - global_offset_x + 1;
                (*Sources_local)[current_local_source][_y_] = gy - global_offset_y + 1;
                current_local_source++;
            }
        }
    } else {
		// if there are no local sources,
		// set the pointer to NULL
        *Sources_local = NULL; 
    }

    free(global_sources_idx);
    return 0;
}

int memory_allocate ( 
			const int *neighbours,  // four-int array that gives back the neighbours of the calling task
			vec2_t    mysize,       // size of the local task
			buffers_t *buffers_ptr, // buffers for communication
			plane_t   *planes_ptr   // plane data
		) {
    /*
     * This function allocates the memory buffers required for both computation 
	 * and communication in the stencil parallel application.
	 *
     *
     * For computation, it allocates two memory regions: one for the "OLD" plane (previous step) 
	 * and one for the "NEW" plane (to store updated results for the current step). 
	 *
	 * These are accessed via the planes_ptr array:
     *   planes_ptr[OLD] refers to the "OLD" region,
     *   planes_ptr[NEW] refers to the "NEW" region.
     *
	 *
     * For communication, it allocates send and receive buffers for each neighbor (if any)
	 * Each buffer can hold up to mysizex or mysizey double values, depending on the direction.
     */

  	if (planes_ptr == NULL ) {
		printf("Error: an invalid pointer has been passed\n");
		return 1;
	}

	if (buffers_ptr == NULL ) {
		printf("Error: an invalid pointer has been passed\n");
		return 1;
	}

	//? ············· allocate memory for data ··············
	// we allocate the space needed for the plane plus a contour frame
	// that will contains data form neighbouring MPI tasks
 
	const uint xsize = mysize[_x_];
	const uint ysize = mysize[_y_];

  	unsigned int frame_size = (xsize+2) * (ysize+2);

	planes_ptr[OLD].data = (double*)malloc(frame_size * sizeof(double));
	if ( planes_ptr[OLD].data == NULL ) {
		// manage the malloc fail
		printf("Error: failed to allocate memory for the old plane\n");
		return 1;
	}
	
	memset(planes_ptr[OLD].data, 0, frame_size * sizeof(double));
	
	planes_ptr[NEW].data = (double*)malloc(frame_size * sizeof(double));
	if ( planes_ptr[NEW].data == NULL ) {
		// manage the malloc fail
		printf("Error: failed to allocate memory for the new plane\n");
		return 1;
	}

	memset(planes_ptr[NEW].data, 0, frame_size * sizeof(double));

	planes_ptr[OLD].size[_x_] = xsize;
	planes_ptr[OLD].size[_y_] = ysize;
	planes_ptr[NEW].size[_x_] = xsize;
	planes_ptr[NEW].size[_y_] = ysize;
	
	//? ················· allocate buffers ·················
	// buffers for north and south communication are not really needed
	// in fact, they are already contiguous, just the first and last line of every rank's plane
	// you may just make some pointers pointing to the correct positions
	// or, if you prefer, just go on and allocate buffers also for north and south communications

	if (neighbours[EAST] != MPI_PROC_NULL || neighbours[WEST] != MPI_PROC_NULL) {
		buffers_ptr[SEND][EAST] = (double*)malloc(ysize * sizeof(double));
		buffers_ptr[RECV][EAST] = (double*)malloc(ysize * sizeof(double));
		if (buffers_ptr[SEND][EAST] == NULL || buffers_ptr[RECV][EAST] == NULL) {
			printf("Error: failed to allocate memory for the buffers\n");
			return 1;
		}
		buffers_ptr[SEND][WEST] = (double*)malloc(ysize * sizeof(double));
		buffers_ptr[RECV][WEST] = (double*)malloc(ysize * sizeof(double));
		if (buffers_ptr[SEND][WEST] == NULL || buffers_ptr[RECV][WEST] == NULL) {
			printf("Error: failed to allocate memory for the buffers\n");
			return 1;
		}
	}

  	return 0;
}

 int memory_release ( 
			buffers_t *buffers,
			plane_t   *planes,
			vec2_t    *Sources_local
		) {

	// free the planes
	if (planes != NULL) {
		if ( planes[OLD].data != NULL ) free (planes[OLD].data);
		if ( planes[NEW].data != NULL ) free (planes[NEW].data);
	}

	// free the buffers
	if (buffers != NULL) {
		if (buffers[SEND][EAST] != NULL) free(buffers[SEND][EAST]);
		if (buffers[RECV][EAST] != NULL) free(buffers[RECV][EAST]);

		if (buffers[SEND][WEST] != NULL) free(buffers[SEND][WEST]);
		if (buffers[RECV][WEST] != NULL) free(buffers[RECV][WEST]);
	}

	if (Sources_local != NULL && *Sources_local != NULL) free(*Sources_local);
      
  	return 0;
}

int output_energy_stat (int step, plane_t *plane, double budget, int Me, MPI_Comm *Comm) {

	double system_energy = 0;
	double tot_system_energy = 0;
	get_total_energy ( plane, &system_energy );
	
	MPI_Reduce ( &system_energy, &tot_system_energy, 1, MPI_DOUBLE, MPI_SUM, 0, *Comm );
  
	if ( Me == 0 ) {
		if ( step >= 0 ) {
			printf(" [ step %4d ] ", step ); fflush(stdout);
		}
      
		printf( "total injected energy is %g, "
			"system energy is %g "
			"( in avg %g per grid point)\n",
			budget,
			tot_system_energy,
			tot_system_energy / (plane->size[_x_]*plane->size[_y_]) 
		);
    }
  
  	return 0;
}
