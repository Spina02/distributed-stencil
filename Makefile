# Compiler settings
CC = gcc
MPICC = mpicc

# Directories
SRCDIR = src
INCDIR = include
BUILDDIR = build
BASELINE_SRCDIR = baseline/src
BASELINE_INCDIR = baseline/include

# Include paths
INCLUDES = -I$(INCDIR)
BASELINE_INCLUDES = -I$(BASELINE_INCDIR)

# Base flags
BASE_CFLAGS = -std=c99 -Wall -Wextra -Wpedantic -Wshadow -Wuninitialized -W
# Disable specific pedantic warnings that are problematic for this code
PEDANTIC_OVERRIDES = -Wno-pointer-sign
MATH_LIBS = -lm

# Optimization flags (your current setup)
OPT_CFLAGS = -Ofast -march=native -flto -fopenmp

# Debug flags
DEBUG_CFLAGS = -g -O0 -DDEBUG

# Release flags (your current setup)
RELEASE_CFLAGS = $(OPT_CFLAGS) -g

# Verbose flags (compile-time verbosity control)
VERBOSE_1_CFLAGS = $(OPT_CFLAGS) -g -DVERBOSE_LEVEL=1
VERBOSE_2_CFLAGS = $(OPT_CFLAGS) -g -DVERBOSE_LEVEL=2
VERBOSE_DEBUG_CFLAGS = -g -O0 -DDEBUG -DVERBOSE_LEVEL=2

# Create build directory
$(BUILDDIR):
	mkdir -p $(BUILDDIR)

# Default target
all: serial stencil_parallel

# Release builds (your current configuration)
serial: $(BUILDDIR)/serial

$(BUILDDIR)/serial: $(SRCDIR)/stencil_serial.c $(INCDIR)/stencil_serial.h | $(BUILDDIR)
	$(CC) $(BASE_CFLAGS) $(PEDANTIC_OVERRIDES) $(RELEASE_CFLAGS) $(INCLUDES) -o $@ $< $(MATH_LIBS)

stencil_parallel: $(BUILDDIR)/stencil_parallel

$(BUILDDIR)/stencil_parallel: $(SRCDIR)/stencil_parallel.c $(INCDIR)/stencil_parallel.h | $(BUILDDIR)
	$(MPICC) $(BASE_CFLAGS) $(PEDANTIC_OVERRIDES) $(RELEASE_CFLAGS) $(INCLUDES) -o $@ $< $(MATH_LIBS)

# Debug builds
debug: serial_debug stencil_parallel_debug

serial_debug: $(BUILDDIR)/serial_debug

$(BUILDDIR)/serial_debug: $(SRCDIR)/stencil_serial.c $(INCDIR)/stencil_serial.h | $(BUILDDIR)
	$(CC) $(BASE_CFLAGS) $(PEDANTIC_OVERRIDES) $(DEBUG_CFLAGS) $(INCLUDES) -o $@ $< $(MATH_LIBS)

stencil_parallel_debug: $(BUILDDIR)/stencil_parallel_debug

$(BUILDDIR)/stencil_parallel_debug: $(SRCDIR)/stencil_parallel.c $(INCDIR)/stencil_parallel.h | $(BUILDDIR)
	$(MPICC) $(BASE_CFLAGS) $(PEDANTIC_OVERRIDES) $(DEBUG_CFLAGS) $(INCLUDES) -o $@ $< $(MATH_LIBS)

# Verbose builds (compile-time verbosity)
verbose1: serial_verbose1 stencil_parallel_verbose1

serial_verbose1: $(BUILDDIR)/serial_verbose1

$(BUILDDIR)/serial_verbose1: $(SRCDIR)/stencil_serial.c $(INCDIR)/stencil_serial.h | $(BUILDDIR)
	$(CC) $(BASE_CFLAGS) $(PEDANTIC_OVERRIDES) $(VERBOSE_1_CFLAGS) $(INCLUDES) -o $@ $< $(MATH_LIBS)

stencil_parallel_verbose1: $(BUILDDIR)/stencil_parallel_verbose1

$(BUILDDIR)/stencil_parallel_verbose1: $(SRCDIR)/stencil_parallel.c $(INCDIR)/stencil_parallel.h | $(BUILDDIR)
	$(MPICC) $(BASE_CFLAGS) $(PEDANTIC_OVERRIDES) $(VERBOSE_1_CFLAGS) $(INCLUDES) -o $@ $< $(MATH_LIBS)

verbose2: serial_verbose2 stencil_parallel_verbose2

serial_verbose2: $(BUILDDIR)/serial_verbose2

$(BUILDDIR)/serial_verbose2: $(SRCDIR)/stencil_serial.c $(INCDIR)/stencil_serial.h | $(BUILDDIR)
	$(CC) $(BASE_CFLAGS) $(PEDANTIC_OVERRIDES) $(VERBOSE_2_CFLAGS) $(INCLUDES) -o $@ $< $(MATH_LIBS)

stencil_parallel_verbose2: $(BUILDDIR)/stencil_parallel_verbose2

$(BUILDDIR)/stencil_parallel_verbose2: $(SRCDIR)/stencil_parallel.c $(INCDIR)/stencil_parallel.h | $(BUILDDIR)
	$(MPICC) $(BASE_CFLAGS) $(PEDANTIC_OVERRIDES) $(VERBOSE_2_CFLAGS) $(INCLUDES) -o $@ $< $(MATH_LIBS)

# Verbose debug builds (maximum verbosity with debug symbols)
verbose_debug: serial_verbose_debug stencil_parallel_verbose_debug

serial_verbose_debug: $(BUILDDIR)/serial_verbose_debug

$(BUILDDIR)/serial_verbose_debug: $(SRCDIR)/stencil_serial.c $(INCDIR)/stencil_serial.h | $(BUILDDIR)
	$(CC) $(BASE_CFLAGS) $(PEDANTIC_OVERRIDES) $(VERBOSE_DEBUG_CFLAGS) $(INCLUDES) -o $@ $< $(MATH_LIBS)

stencil_parallel_verbose_debug: $(BUILDDIR)/stencil_parallel_verbose_debug

$(BUILDDIR)/stencil_parallel_verbose_debug: $(SRCDIR)/stencil_parallel.c $(INCDIR)/stencil_parallel.h | $(BUILDDIR)
	$(MPICC) $(BASE_CFLAGS) $(PEDANTIC_OVERRIDES) $(VERBOSE_DEBUG_CFLAGS) $(INCLUDES) -o $@ $< $(MATH_LIBS)

# Baseline builds (for comparison)
baseline: baseline_serial baseline_parallel

baseline_serial: $(BUILDDIR)/baseline_serial

$(BUILDDIR)/baseline_serial: $(BASELINE_SRCDIR)/stencil_template_serial.c $(BASELINE_INCDIR)/stencil_template_serial.h | $(BUILDDIR)
	$(CC) $(BASE_CFLAGS) $(PEDANTIC_OVERRIDES) $(RELEASE_CFLAGS) $(BASELINE_INCLUDES) -o $@ $< $(MATH_LIBS)

baseline_parallel: $(BUILDDIR)/baseline_parallel

$(BUILDDIR)/baseline_parallel: $(BASELINE_SRCDIR)/stencil_template_parallel.c $(BASELINE_INCDIR)/stencil_template_parallel.h | $(BUILDDIR)
	$(MPICC) $(BASE_CFLAGS) $(PEDANTIC_OVERRIDES) $(RELEASE_CFLAGS) $(BASELINE_INCLUDES) -o $@ $< $(MATH_LIBS)

# Clean targets
clean:
	rm -rf $(BUILDDIR)

clean_all: clean
	rm -f *.o *.out core.*

# Performance testing
test_serial: $(BUILDDIR)/serial
	$(BUILDDIR)/serial -v 1

test_parallel: $(BUILDDIR)/stencil_parallel
	mpirun -np 4 $(BUILDDIR)/stencil_parallel -v 1

# Show current compiler flags
show_flags:
	@echo "Base flags: $(BASE_CFLAGS)"
	@echo "Release flags: $(RELEASE_CFLAGS)"
	@echo "Debug flags: $(DEBUG_CFLAGS)"
	@echo "Verbose1 flags: $(VERBOSE_1_CFLAGS)"
	@echo "Verbose2 flags: $(VERBOSE_2_CFLAGS)"
	@echo "Verbose debug flags: $(VERBOSE_DEBUG_CFLAGS)"
	@echo "Includes: $(INCLUDES)"
	@echo "Math libs: $(MATH_LIBS)"
	@echo "Build dir: $(BUILDDIR)"

# Help target
help:
	@echo "Available targets:"
	@echo "  all              - Build release versions (default)"
	@echo "  serial           - Build serial version"
	@echo "  stencil_parallel - Build parallel version"
	@echo "  debug            - Build debug versions"
	@echo "  verbose1         - Build with VERBOSE_LEVEL=1 (basic output)"
	@echo "  verbose2         - Build with VERBOSE_LEVEL=2 (detailed output)"
	@echo "  verbose_debug    - Build verbose with debug symbols"
	@echo "  baseline         - Build baseline versions"
	@echo "  clean            - Remove build directory"
	@echo "  clean_all        - Remove all generated files"
	@echo "  test_serial      - Run serial test"
	@echo "  test_parallel    - Run parallel test with 4 processes"
	@echo "  show_flags       - Show current compiler flags"
	@echo "  help             - Show this help message"
	@echo ""
	@echo "All executables are built in the $(BUILDDIR)/ directory"
	@echo ""
	@echo "Verbose levels:"
	@echo "  VERBOSE_LEVEL=1  - Basic verbose output (neighbor info, timing)"
	@echo "  VERBOSE_LEVEL=2  - Detailed output (includes grid visualization)"

.PHONY: all debug baseline verbose1 verbose2 verbose_debug clean clean_all test_serial test_parallel show_flags help serial stencil_parallel serial_debug stencil_parallel_debug baseline_serial baseline_parallel serial_verbose1 stencil_parallel_verbose1 serial_verbose2 stencil_parallel_verbose2 serial_verbose_debug stencil_parallel_verbose_debug 