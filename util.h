/**
 * Several utility functions for displaying results.
 */

#pragma once

#ifndef _XOPEN_SOURCE
#if __STDC_VERSION__ >= 199901L
#define _XOPEN_SOURCE 600
#else
#define _XOPEN_SOURCE 500
#endif
#endif

#include <stdlib.h>
#include <time.h>
#include <sys/time.h>

#ifdef __cplusplus
extern "C" {
#endif


/**
 * Prints a positive number with the given number of sigfigs and a unit. The
 * value is scaled to the correct unit (which are mult apart - 1000 for SI and
 * 1024 for digit prefixes).
 */
void print_with_unit(double val, int sigfigs, int mult,
                     const char** units, size_t n_units);

/**
 * Prints a number of bytes after converting to a nicer unit.
 */
void print_bytes(size_t n);

/**
 * Print the time (in seconds) with the right units and 3 significant digits.
 */
void print_time(double seconds);

/**
 * Get the difference between two times.
 */
double get_time_diff(struct timespec* start, struct timespec* end);

#ifdef __cplusplus
}
#endif
