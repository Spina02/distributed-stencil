# Introduction

In this assignment, you are asked to parallelize a stencil computation, beginning either from a complete serial implementation or from a parallel template that provides the basic structure but lacks many components.

A classic example of such a stencil computation is the **Heat Equation**. The heat equation models how heat diffuses through a material over time and is given by:

$$
\frac{\partial u(\vec{x}, t)}{\partial t} = \alpha \nabla^2 u(\vec{x}, t)
$$

where $u(\vec{x}, t)$ is the energy (or temperature) at position $\vec{x}$ and time $t$, and $\alpha$ is the heat diffusivity constant.

To numerically solve this equation, we discretize the spatial domain into a grid and use a finite-difference method. The most common approach is the five-point stencil, where the value at each grid point is ***updated based on its current value and the values of its four immediate neighbors*** (up, down, left, right). This method is both simple and efficient, making it well-suited for parallelization.

Your task is to implement this stencil computation in parallel, taking into account the need to *exchange boundary data between neighboring subdomains* when the grid is partitioned across *multiple processes*. This is a fundamental pattern in scientific computing and is widely used in simulations of physical phenomena.

## Project structure

```c
stencil_template_parallel.c  // Template for parallel code
stencil_template_serial.c    // Template for serial code
```

the template for the parallel code already contains the initialization, some memory allocation, and the skeleton of the main loop. However, it already proposes data structures.

You can start from fresh from the serial version, which is a complete very simple version.

There is also a template batch script for job submission, it is named `go`.

## Problem formulation

The new value of the green point depends on its old values and on the vaues of four surrounding points, i.e. it is “local”.
However, if you are parallelizing the algorithm, the points at the border of the computational domain of task Tk belong to another task Tg, unless they are at the border and the boundary conditions are not periodic.

The new value of the green point depends on its old values and on the vaues of four surrounding points, i.e. it is “local”. However, if you are parallelizing the algorithm, the points at the border of the computational domain of task Tk belong to another task Tg, unless they are at the border and the boundary conditions are not periodic.

Both tasks Tk and Tg need to know the values of the points at the border of their neighbours (for a 5-points stencil only a layer of 1 point is sufficient)

Every task then needs an additional layer that contains the data from its **north**, **east**, **south**, and **west** neighbours. At every timestep, the updated values must be communicated *to/from* the involved tasks.  

Every task may either:
- *embed its local grid patch into a larger one*
-  *have separated buffers*.

When the local grid is embedded into a larger one, checking for the boundaries in the update loop is easier.  
But, if you do not build a derived vector data type, you need to specifically treat *non-contiguous buffers* (in C, the **east** and **west** buffers) when communicating.

Having separated buffers makes the communication easier, but you still have to move data *to/from* the buffers (which, in any case, is still true for **east** and **west**, if you do not use derived data types).

When the boundary conditions are periodic, your plate behaves like it is infinite and with no borders.

<!-- The note above is very important -->

**Note:** We have not seen MPI derived data types during the lectures.  
Hence, you’ll be forced to move data explicitly *from/to* the **east** and **west** buffers.


### Implementation

All in all, the parallel code will perform the following steps:

1. **initialization**

    - read arguments from the command line: `(x,y)`-size of the plate, the number of iterations `N_iter`, the number of heat sources `N_heat`, whether the boundary conditions are periodic.

    - determine the MPI tasks grid: how you decompose the computational domain into a grid of $N_x×N_y$ tasks

    - determine the neighbours of every MPI tasks (depends on the decomposition and on the boundary periodicity)

2. **update loop**

    for N_iter:

    $\quad$ **a)** inject energy from the sources.

    $\quad$ **b)** update the local grid.

    $\quad$ **c)** exchange borders with the neighbours, if any.

    Hint for the last two steps: "is there any possible, at least partial, overlap between these two steps?"

Everything except the a) and c) steps is alreasy present in the templates

## Requirements

**(A)** Build a parallel implementation of the 5-points stencil, based on the provided templates.

- Use the given code structure as a starting point.

**(B)** Parallelization must be with:

- MPI (for distributed memory parallelism)
- OpenMP (for the update loop of the local grid)

**(C)** Instrument your code to measure:

- The time spent on computation
- The time spent on communication

**(D)** Perform a scalability study

### About scalability

- **openmp**

    use 1 MPI task and scale the nr. of threads from 1 up to filling the entire code (there are 56 cores/socket, 2 sockets; it is not necessary to have 112 different runs, a series like 1,2,4,8,16,32,56,84,112 is ok).

- **MPI**

    from the openmp scaling find the best threads-per-task ratios and elect your preferred choice; for instance, let’s say that you choose 8 MPI tasks and 14 openmp threads per node. That given, perform strong and weak scaling using nodes as reference, ranging from 1 up to 16 or 32 nodes (1,2,4,8,16,32)

Each student can use about 50 core-hours maximum of Leonardo Supercomputer

Using the command:
```bash
saldo -b --dcgp
```
 you can inspect how much time has been consumed

As of June 25th, the account is currently on the "booster" partition (GPU-based). However, for this project, we require the "data-centric" (dcgp) partition, which is preferable as it provides a larger amount of core-hours.

- When you run `saldo -b --dcgp`, you should currently see no budget available.
- Running `saldo -b` (without the `--dcgp` flag) will show your current budget on the booster partition.

A request has been made to move the budget to the dcgp partition, and CINECA is processing it. You will know the transfer is complete when your budget appears with `saldo -b --dcgp`.

In the meantime, there are two batch files provided:
- `go_booster`
- `go_dcgp`

Until the budget is moved to the dcgp partition, please use `go_booster` for your jobs (the difference between the two files is minimal).

I remind you that **“scalability”** is the capability of your code to efficiently use different amount of computational resources when

1) the problem’s size is constant, which is called **strong scaling**
2) the problem’s size per resource is constant, which is called **weak scaling**

“***Resources***” may be MPI tasks and/or threads (then basically how many cores you use) or nodes, if you always use entire nodes (which is oftne the case)

#### Strong Scaling

Considering a given problem, we have seen in the lectures that we expect the
time-to-solution to scale as $∝ 1/n$ if $n$ is the amount of computational
resources, for a perfect scaling.

The difference between the real scaling and the perfect scaling is due to either
the intrinsic serial fraction $f$ of the problem, or to parallel overhead $k$.

If your problem has a serial fraction, then your speed-up has an hard limit:
$$
Sp = \dfrac{t(1, s)}{t(n, s)} ⟶ \dfrac1{1-f}
$$
where $t(n,s)$ is the run-time for a size $s$, using $n$ resources (which I remind you is called “*the Amdahls’s law*").

Then, when you test a code for its strong scaling capability, as first analyze your
algorithm and individuate where you loose parallelism by its very nature (and, at
the discussion, try to imagine a different algorithm or a way to express more
parallism).

Then, check whether your scaling results match with your analysis, and how the
run-time scales once you purge the inherently serial part.
That is why requirement c) is about instrumenting your code so that you can
check the timing of different sections and of the communications.

At the end, you should be able to characterize your code, anche discuss
whether any loss of parallelism is expected because it is due to the algorithm, or
it is due to your implementation.

Given that you expect a behaviour, the best way to present you results is not
by plotting absolute times.

In fact, human eye is not good at judging non-linear behaviours, and in case you
need to plot real timings, always superimpose the ideal scaling to your data.
The best way to report about scalability is to plot the Speedup, as defined
before and in the lectures and that it is expected to be linear, or the efficiency,
which is defined as:

$$
\text{Eff}(n,s) = \text{Sp}(n,s) / n
$$

which is expected to be $\sim 1$ for an ideal behaviour

**Plot**

At the end, you end up with a plot that is similar to the one described.

on the x-axis there are the computational resources (in this case the nodes, not the
MPI tasks or the threads), and on the y-axis there is speed up.

Note that the ideal scaling is reported to guide the eye. 

> There will always be some lost of efficiency for a large scaling factor, because for a large-enough amount of resources the parallel overhead will be dominant. There is not a formal definition of this threshold, since that is very much dependent on the problem and the algorithm. However I would say for good code that should happen beyond 16× .

> gap due to inherent non-parallelism, communication, overhead due to inefficiency or to parallelization itself?


#### Weak Scaling

Weak Scaling measure how the time-to-solution varies when you grow the amount of resources $n$ while keeping the workload per resource constant. For instance, for the case under exam, the problem’s size scales as $N_{grid}× N_{grid}$ and then the workload per resource is $N_{grid}^2/n$.

Then, when you enlarge $n$, you have to keep this ratio constant to maintain the workload.

Reporting the weak scaling correctly is easier than for the strong scaling, since the expectation is that the run-time reamains constant. Then, also the absolute timing is “correct”. However, still, reporting the **efficiency** in this case give a better glance of how much the code looses (10%, 20%, 50% ?)
