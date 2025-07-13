# Compiler settings
CC = gcc
MPICC = mpicc

# Directories
SRCDIR = src
INCDIR = include
BASELINE_SRCDIR = baseline/src
BASELINE_INCDIR = baseline/include

# Include paths
INCLUDES = -I$(INCDIR)
BASELINE_INCLUDES = -I$(BASELINE_INCDIR)

# Base flags
BASE_CFLAGS = -std=c99 -Wall -Wextra -Wpedantic -Wshadow -Wuninitialized -W
MATH_LIBS = -lm

# Optimization flags (your current setup)
OPT_CFLAGS = -Ofast -march=native -flto -fopenmp

# Debug flags
DEBUG_CFLAGS = -g -O0 -DDEBUG

# Release flags (your current setup)
RELEASE_CFLAGS = $(OPT_CFLAGS) -g

# Default target
all: serial stencil_parallel

# Release builds (your current configuration)
serial: $(SRCDIR)/stencil_serial.c $(INCDIR)/stencil_serial.h
	$(CC) $(BASE_CFLAGS) $(RELEASE_CFLAGS) $(INCLUDES) -o $@ $< $(MATH_LIBS)

stencil_parallel: $(SRCDIR)/stencil_parallel.c $(INCDIR)/stencil_parallel.h
	$(MPICC) $(BASE_CFLAGS) $(RELEASE_CFLAGS) $(INCLUDES) -o $@ $< $(MATH_LIBS)

# Debug builds
debug: serial_debug stencil_parallel_debug

serial_debug: $(SRCDIR)/stencil_serial.c $(INCDIR)/stencil_serial.h
	$(CC) $(BASE_CFLAGS) $(DEBUG_CFLAGS) $(INCLUDES) -o $@ $< $(MATH_LIBS)

stencil_parallel_debug: $(SRCDIR)/stencil_parallel.c $(INCDIR)/stencil_parallel.h
	$(MPICC) $(BASE_CFLAGS) $(DEBUG_CFLAGS) $(INCLUDES) -o $@ $< $(MATH_LIBS)

# Baseline builds (for comparison)
baseline: baseline_serial baseline_parallel

baseline_serial: $(BASELINE_SRCDIR)/stencil_template_serial.c $(BASELINE_INCDIR)/stencil_template_serial.h
	$(CC) $(BASE_CFLAGS) $(RELEASE_CFLAGS) $(BASELINE_INCLUDES) -o $@ $< $(MATH_LIBS)

baseline_parallel: $(BASELINE_SRCDIR)/stencil_template_parallel.c $(BASELINE_INCDIR)/stencil_template_parallel.h
	$(MPICC) $(BASE_CFLAGS) $(RELEASE_CFLAGS) $(BASELINE_INCLUDES) -o $@ $< $(MATH_LIBS)

# Clean targets
clean:
	rm -f serial stencil_parallel serial_debug stencil_parallel_debug
	rm -f baseline_serial baseline_parallel

clean_all: clean
	rm -f *.o *.out core.*

# Performance testing
test_serial: serial
	./serial -v 1

test_parallel: stencil_parallel
	mpirun -np 4 ./stencil_parallel -v 1

# Show current compiler flags
show_flags:
	@echo "Base flags: $(BASE_CFLAGS)"
	@echo "Release flags: $(RELEASE_CFLAGS)"
	@echo "Debug flags: $(DEBUG_CFLAGS)"
	@echo "Includes: $(INCLUDES)"
	@echo "Math libs: $(MATH_LIBS)"

# Help target
help:
	@echo "Available targets:"
	@echo "  all            - Build release versions (default)"
	@echo "  serial         - Build serial version"
	@echo "  stencil_parallel - Build parallel version"
	@echo "  debug          - Build debug versions"
	@echo "  baseline       - Build baseline versions"
	@echo "  clean          - Remove built executables"
	@echo "  clean_all      - Remove all generated files"
	@echo "  test_serial    - Run serial test"
	@echo "  test_parallel  - Run parallel test with 4 processes"
	@echo "  show_flags     - Show current compiler flags"
	@echo "  help           - Show this help message"

.PHONY: all debug baseline clean clean_all test_serial test_parallel show_flags help 