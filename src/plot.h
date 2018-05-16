#ifndef __PLOT__H
#define __PLOT__H

#include <stdlib.h>
#include <stdio.h>
#include <sys/stat.h>

#define MAX_BYTES  16

void generate_combined_plot(char *plot_name, double *x, double *y, int npoints, 
				   int byte, int nbytes, int is_compressible );

void generate_plot( char *plot_name, double *x, double *y, int npoints );

#endif
