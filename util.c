/**
 * Several utility functions for displaying results.
 */

#include "util.h"

#include <stdlib.h>
#include <stdio.h>
#include <time.h>

/**
 * Prints a positive number with the given number of sigfigs and a unit. The
 * value is scaled to the correct unit (which are mult apart - 1000 for SI and
 * 1024 for digit prefixes).
 */
void print_with_unit(double val, int sigfigs, int mult,
                     const char** units, size_t n_units) {
    size_t i_unit = 0;
    while (i_unit < n_units && val >= mult) { val /= mult; i_unit++; }
    if (i_unit == 0) { sigfigs = 0; }
    else if (val < 10) { sigfigs -= 1; }
    else if (val < 100) { sigfigs -= 2; }
    else { sigfigs -= 3; }
    printf("%.*f %s", sigfigs, val, units[i_unit]);
}

/**
 * Prints a number of bytes after converting to a nicer unit.
 */
void print_bytes(size_t n) {
    static const char* units[4] = {"bytes", "KiB", "MiB", "GiB"};
    print_with_unit(n, 3, 1024, units, 4);
}

/**
 * Print the time (in seconds) with the right units and 3 significant digits.
 */
void print_time(double seconds) {
    static const char* units[4] = {"ns", "us", "ms", "s"};
    print_with_unit(seconds * 1000000000.0, 3, 1000, units, 4);
}

/**
 * Get the difference between two times.
 */
double get_time_diff(struct timespec* start, struct timespec* end) {
    double diff = end->tv_sec - start->tv_sec;
    diff += (end->tv_nsec - start->tv_nsec) / 1000000000.0;
    return diff;
}
