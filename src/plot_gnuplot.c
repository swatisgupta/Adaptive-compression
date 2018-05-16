#include "plot.h"

int compressible[MAX_BYTES];

void write_to_file(double *x, double *y, int npoints, int byte) {
  double mean = 0;
  int i = 0; 
  char filename[30];
  sprintf(filename, "datafile%d.dat", byte);
  FILE *datafile = fopen(filename, "w");
  for (i = 0; i < npoints; i++) {
  	fprintf(datafile, " %g %g\n", x[i], y[i]);
	mean += y[i]; 
  }
  mean /= npoints;

  fclose(datafile);

} 

FILE * set_plot(char * plot_name, int npoints) {

  struct stat status = {0};
  if (stat("plots", &status) == -1) {
    mkdir("plots", 0700);
  }
  FILE *gnuplot = popen("gnuplot", "w");
  fprintf(gnuplot, "set term png\n");
  fprintf(gnuplot, "set output \"plots/%s.png\"\n", plot_name);
  fprintf(gnuplot, "set datafile separator \" \"\n");
  fprintf(gnuplot, "set xrange[0:%d] \n", npoints);

  return gnuplot;
}

FILE * set_multi_plot(char * plot_name, int npoints) {

  struct stat status = {0};
  if (stat("plots", &status) == -1) {
    mkdir("plots", 0700);
  }
  FILE *gnuplot = popen("gnuplot", "w");

  fprintf(gnuplot, "set term png\n");
  fprintf(gnuplot, "set output \"plots/%s.png\"\n", plot_name);
  fprintf(gnuplot, "set datafile separator \" \"\n");
  fprintf(gnuplot, "set xrange[0:%d] \n", npoints);
  fprintf(gnuplot, "set xtics rotate by 45 \n");
  fprintf(gnuplot, "set multiplot layout 2,4 title \"Frequency distribution of different byte columns \" \n");
  fprintf(gnuplot, "set lmargin 7\n");
  fprintf(gnuplot, "set rmargin 1\n");
  fprintf(gnuplot, "set tmargin 1\n");
  fprintf(gnuplot, "set bmargin 1\n");
  return gnuplot;
}


void plot_function( FILE* plot, char* function_name, char* function, int type, int compress,  double mean) {

   static double origin = 0;
  //fprintf(plot, "set title \"%s\"\n", function_name);
  fprintf(plot, "set size 0.24, 0.24 \n");
  //fprintf(plot, "set origin %0.1g, %0.2g \n", origin);
  //origin += 0.2;
  switch(type) {
    case 1: /** Points **/
        fprintf(plot, "set style line 1 lc rgb \"green\" lt 1 lw 1  pt 1 ps 0.2\n");
  	fprintf(plot, "set style line 2 lc rgb \"red\" lt 1 lw 1  pt 2 ps 0.2 \n");
#if 0
        fprintf(plot, "f(X) = %g\n", mean); 
	fprintf(plot, "plot \'%s\' with linespoints ls 1, f(X) lw 3 lc rgb \"blue\" \n", function);
#else
	if (compress) {
		//fprintf(plot, "plot \'%s\' with linespoints  ls 1 title \'compressible\' \n", function);
		fprintf(plot, "plot \'%s\' with linespoints  ls 1 notitle\n", function);
	} else {
		//fprintf(plot, "plot \'%s\' with linespoints  ls 2 title \'non-compressible\' \n", function);
		fprintf(plot, "plot \'%s\' with linespoints  ls 2 notitle\n", function);
	} 
	
#endif
	break;
    case 2: /** Horizontal line **/
        fprintf(plot, "f(X) = %g\n", mean); 
	fprintf(plot, "plot f(X) lw 3 lc rgb \"blue\"  \n");
	break;
    case 3: /** Vertical line **/
	fprintf(plot, "set parametric\n");
        fprintf(plot, "plot %g,t lt 1 lw 3 lc rgb \"blue\" \n", mean);
        
   }
}

void close_plot(FILE * plot) {
  fprintf(plot, "replot \n");
  fflush(plot); 
  pclose(plot);
 
}

void generate_combined_plot(char *plot_name, double *x, double *y, int npoints, int byte, int nbytes, int is_compressible )  {

    char filename[50]; 
    char function_name[20];
    int i;
    FILE * plot = NULL;
    compressible[byte] = is_compressible;
    write_to_file(x, y, npoints, byte);
	
    if( byte == nbytes - 1 ) {
	plot = set_multi_plot(plot_name, npoints);
	for (i = 0; i < nbytes; i++) {
		sprintf(filename, "datafile%d.dat", i); 
		sprintf(function_name, "Byte %d", i); 
		plot_function(plot, function_name, filename, 1, compressible[i], 0.0);
	}	
        close_plot(plot);	
	for (i = 0; i < nbytes; i++) {
		sprintf(filename, "datafile%d.dat", i); 
		sprintf(function_name, "Byte %d", i); 
                remove(filename);
	}	
    }
}	

void generate_plot( char *plot_name, double *x, double *y, int npoints ) {

  FILE *datafile = fopen("datafile.dat", "w");

  struct stat status = {0};
  if (stat("plots", &status) == -1) {
    mkdir("plots", 0700);
  }
  
  int i;
  double mean = 0;
  for (i = 0; i < npoints; i++) {
  	fprintf(datafile, " %g %g\n", x[i], y[i]);
	mean += y[i]; 
  }
  mean /= npoints;

  fclose(datafile);

  FILE *gnuplot = popen("gnuplot", "w");
  fprintf(gnuplot, "set term png\n");
  fprintf(gnuplot, "set output \"plots/%s.png\"\n", plot_name);
  fprintf(gnuplot, "set datafile separator \" \"\n");
  fprintf(gnuplot, "set xrange[0:%d] \n", npoints);
  fprintf(gnuplot, "set yrange[0:255] \n");
  fprintf(gnuplot, "set output \"plots/%s.png\"\n", plot_name);
  
  fprintf(gnuplot, "plot \'%s\' with linespoints\n", "datafile.dat");
  fprintf(gnuplot, "f(x) = %g \n", mean);
  fprintf(gnuplot, "plot f(x) lw 3 lc rgb \"blue\"  \n");
  fprintf(gnuplot, "replot\n");
  //fprintf(gnuplot, "set parametric\n");
  //fprintf(gnuplot, "plot %g,t lt 1 lw 3 lc rgb \"red\" \n", mean);
  fflush(gnuplot); 
  pclose(gnuplot);
  remove("datafile.dat");
}

