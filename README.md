# Parallel Stencil Computation: A Hybrid MPI+OpenMP Heat Equation Solver

## Overview

This project presents a high-performance parallel implementation of a 2D heat equation solver using a five-point stencil method. The primary goal is to demonstrate and analyze the scalability of stencil computations on modern High-Performance Computing (HPC) architectures.

The implementation leverages a hybrid parallel programming model:
-   **MPI (Message Passing Interface)** is used for distributed-memory parallelism, decomposing the computational domain across multiple nodes of a cluster.
-   **OpenMP (Open Multi-Processing)** is used for shared-memory parallelism, parallelizing the computational loops within a single node to utilize its multi-core architecture.

The solver simulates the diffusion of heat over a 2D grid, updating each point's temperature based on its neighbors at each timestep. This pattern requires efficient communication of boundary data ("halos") between adjacent subdomains, which is a core challenge in parallel scientific computing.

### Project Structure

```
.
├── analysis/           # Python scripts for plotting results
├── assignment/         # Original assignment description and templates
├── data/               # CSV output files from scalability experiments
├── include/            # Header files (.h) for the source code
├── plots/              # Generated plots (speedup, efficiency, timing)
├── scripts/            # Batch and shell scripts for running experiments on SLURM
├── src/                # C source code (.c) for serial and parallel versions
├── LICENSE
├── Makefile
└── README.md
```

## Building the Project

The project is built using a standard `Makefile` which utilizes `mpicc` as the compiler wrapper, enabling both MPI and OpenMP functionalities.

To compile all executables:
```bash
make all
```

This command will generate two primary executables in the root directory:
-   `stencil_parallel`: The hybrid MPI+OpenMP implementation.
-   `stencil_serial`: A reference serial implementation for correctness checking and performance baseline.

To clean the build artifacts:
```bash
make clean
```

#### Build Configurations

The `Makefile` provides several build targets for different purposes. All executables are placed in the `build/` directory.

*   **Release (default):** For maximum performance.
    ```bash
    make all
    ```
    This builds `build/stencil_serial` and `build/stencil_parallel` with full optimizations (`-Ofast`, `-flto`).

*   **Debug:** For use with debuggers like GDB.
    ```bash
    make debug
    ```
    This builds unoptimized (`-O0`) versions with debug symbols (`-g`), named `build/serial_debug` and `build/stencil_parallel_debug`.

*   **Verbose Builds:** To enable compile-time logging for observing program behavior. These are optimized builds that also include debug symbols.
    -   For basic output (timing, neighbor info):
        ```bash
        make verbose1
        ```
    -   For detailed output (including grid visualization):
        ```bash
        make verbose2
        ```

*   **Verbose Debug:** The ultimate debugging build, combining no optimization, debug symbols, and maximum verbosity.
    ```bash
    make verbose_debug
    ```

For a full list of targets and options, run `make help`.

## Running the Code

### Serial Version
The serial code can be run directly from the command line:
```bash
./stencil_serial [options]
```

### Parallel Version (MPI + OpenMP)
The parallel version requires `mpirun` to launch multiple processes. The number of OpenMP threads per process is controlled via an environment variable.

```bash
# Set the number of threads for OpenMP
export OMP_NUM_THREADS=<N_THREADS>

# Run the MPI application
mpirun -np <N_MPI_TASKS> ./stencil_parallel [options]
```

#### Command-Line Options
All flags are optional. Default values are set within the code, optimized for runs on HPC systems.

-   `-x <width>`: Global width of the simulation grid.
-   `-y <height>`: Global height of the simulation grid.
-   `-n <iterations>`: Number of time steps to simulate.
-   `-e <n_sources>`: Number of heat sources.
-   `-E <energy>`: Energy value injected by each source per step.
-   `-p <0|1>`: Enable (1) or disable (0) periodic boundary conditions.
-   `-o <0|1>`: Enable (1) or disable (0) energy statistics output at each step.

> [!WARNING]
> The default grid dimensions are large and intended for HPC environments. For local testing on a personal computer, it is recommended to specify smaller dimensions (e.g., `-x 512 -y 512`).

### Running on an HPC Cluster (SLURM)

The `scripts/` directory contains scripts to automate scalability studies on SLURM-based clusters like Leonardo.

-   **`go_dcgp.sbatch`**: A template SLURM batch script for submitting a single job to the data-centric partition. You should modify the `#SBATCH` directives and the execution command as needed.
-   **`omp_scaling.sh`**, **`strong_scaling.sh`**, **`weak_scaling.sh`**: Helper shell scripts that automatically generate and submit a series of jobs to perform the complete scalability analysis.

## Scalability Analysis and Results

The performance of the implementation was evaluated on the **Leonardo supercomputer** (CINECA), specifically on its data-centric (DCGP) partition. Each DCGP node features:
-   **2x Intel(R) Xeon(R) Platinum 8480+ CPUs**
-   **112 physical cores** in total (56 per socket)
-   Hyper-Threading **disabled** (1 thread per core)
-   A complex **NUMA architecture** with 8 nodes.

The following analyses were conducted based on the requirements.

- #### OpenMP Thread Scaling

    An OpenMP thread scaling analysis was performed on a single node using one MPI task, while varying the number of OpenMP threads.

- #### Strong Scaling (MPI + OpenMP)

    A strong scaling analysis was performed using a fixed global problem size, distributed across an increasing number of nodes.

- #### Weak Scaling (MPI + OpenMP)

    A weak scaling analysis was performed by keeping the problem size per node constant, meaning the global problem size grew proportionally with the number of nodes.