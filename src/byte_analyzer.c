#include <limits.h>
#include <math.h>
#include <time.h>
#include "util.h"
#include "byte_analyzer.h"
#include "preconditioners.h"
#include "compressor.h"

void entropy_analysis( acomp_size_t *byte_col_freq, acomp_size_t num_elements, 
		       ubyte_t *byte_flags, int *num_compressible, int byte_pos);
void skewness_analysis( acomp_size_t *byte_col_freq, acomp_size_t num_elements, 
		       ubyte_t *byte_flags, int *num_compressible, int byte_pos);
void kurtosis_analysis( acomp_size_t *byte_col_freq, acomp_size_t num_elements, 
		       ubyte_t *byte_flags, int *num_compressible, int byte_pos);
void naive_analysis( acomp_size_t *byte_col_freq, acomp_size_t num_elements, 
		       ubyte_t *byte_flags, int *num_compressible, int byte_pos);

void isobar_analysis( acomp_size_t *byte_col_freq, acomp_size_t num_elements, 
		       ubyte_t *byte_flags, int *num_compressible, int byte_pos);

#ifdef PLOT
 FILE * stat_file = NULL;
#endif

void check_compressibility (byte_t* data, int *byte_col_freq, acomp_size_t num_elements, 
			    int ele_size, ubyte_t *byte_flags, int *num_compressible,
								    int byte_pos ) {
#ifdef PLOT
    time_t current_time;
    int date, hr, min, sec, year;
    char month[30], day[30], plotname[40]; 
    current_time = time(NULL);
    char * str_time = ctime(&current_time);
    double x_val[UCHAR_MAX], y_val[UCHAR_MAX];
    acomp_size_t byte_val;
#endif

    acomp_size_t c_size = num_elements;
    byte_t *prep_data = (byte_t *)malloc(num_elements * ele_size),
            *uncomp_buffer = (byte_t*)malloc(num_elements * ele_size),
    	    *comp_buffer = (byte_t*)malloc(num_elements * ele_size); 

    ubyte_t b_flag = 0;


#ifdef PLOT         
    for ( byte_val = 0; byte_val < UCHAR_MAX; byte_val++ ) {
      x_val[byte_val] = byte_val; 
      y_val[byte_val] = byte_col_freq[byte_val];
    }    

    sscanf(str_time, "%s %s %d %d:%d:%d %d", day, month, &date, &hr, &min,
							     &sec, &year);  
    sprintf(plotname, "%s_%s_%d_%d_time_%d:%d", c_var, month, date, year,
								 hr, min);
#endif

    b_flag = 1 << byte_pos;

    preprocess_data(b_flag, data, prep_data, uncomp_buffer, ele_size, 
							num_elements, 1);


    c_size = compress_prep(prep_data, comp_buffer, num_elements, 0);
    DB_PRINTF("Checking %d :csize %lu and actual size %lu \n", 
					      byte_pos, c_size, num_elements); 	
    if ( c_size < num_elements ) {  
     //       c_size = compress_prep(prep_data, comp_buffer, num_elements, 2);
            if ( c_size < num_elements )  { 
#ifdef PLOT
                generate_combined_plot(plotname, x_val, y_val, UCHAR_MAX, byte_pos, 
								ele_size, 1);  
#endif
                 *byte_flags |=  1 << byte_pos;
            } else {
#ifdef PLOT
               generate_combined_plot(plotname, x_val, y_val, UCHAR_MAX, byte_pos, 
								ele_size, 0);    
#endif
                 (*num_compressible) --;
            }     
    } else {
#ifdef PLOT
            generate_combined_plot(plotname, x_val, y_val, UCHAR_MAX, byte_pos, 
								ele_size, 0);    
#endif
            (*num_compressible) --;
    }   

    free(prep_data); 
    free(comp_buffer); 
    free(uncomp_buffer); 
}

ubyte_t get_compressibility_flags(int dtype, int ele_size, acomp_size_t data_size, 
				   byte_t* input_data, int *num_compressible) 
{
  ubyte_t byte_flags = 0, *data = (ubyte_t*) input_data;
  char * msr;
  int measure = 1;
  acomp_size_t byte_col_freq[UCHAR_MAX], byte_pos, byte_val, 
						num_elements = data_size/ele_size;
  acomp_size_t tmp;

#ifdef PLOT
  ubyte_t byte_f = 0;
  byte_t varname[30];
  int num_comp = ele_size;
  FILE *var_data_file = NULL; 
#endif   

  if( (msr = getenv("MEASURE")) != NULL ) 
	measure = atoi(msr);	

  DB_PRINTF("Inside get_byte_flags_dist: element size [%d], num_elements [%lu], \
			data_size[%lu] Is real? %s\n", ele_size, num_elements, 
				data_size, (dtype == IS_REAL)? "[Yes]": "[No]" );

  *num_compressible = ele_size;

/**
  if ( dtype != IS_REAL ) {
    return (UCHAR_MAX - 1);	
  } 
**/
	
#ifdef PLOT
   stat_file = fopen("statistics.txt", "a+");
#if 0 

  sprintf(varname, "%s.txt", c_var);	

  var_data_file = fopen(c_var, "w");

  for( byte_val = 0; byte_val< num_elements; byte_val++) {
        fprintf(var_data_file, "%d, ", byte_val);
        for( byte_pos = 0; byte_pos < ele_size; byte_pos++) 
	   fprintf(var_data_file, "%lu, ", data[byte_val* ele_size + byte_pos]); 
        fprintf(var_data_file, "\n");
  }
  fclose(var_data_file);
#endif

#endif 

  for( byte_pos = 0; byte_pos < ele_size; byte_pos++) {

    for( byte_val = 0; byte_val < UCHAR_MAX; byte_val++) {
	byte_col_freq[byte_val] = 0;
    }

    for( byte_val = 0; byte_val < num_elements; byte_val++) {
      DB_PRINTF("byte_val %lu, ele_size %d, byte_pos %lu \n", 
				byte_val, ele_size, byte_pos, num_elements);
      if ( byte_val * ele_size  + byte_pos >= data_size) {  
		printf("Error in index, byte_val %lu, ele_size %d, byte_pos %lu \n", 
						      byte_val, ele_size, byte_pos);
		exit(0);
      }
     
      tmp = data[byte_val * ele_size  + byte_pos];
      byte_col_freq[tmp]++;
    }

  
    switch (measure) {
       case 1: /* Entropy */
           entropy_analysis(byte_col_freq, num_elements, &byte_flags, 
							num_compressible, byte_pos );
	   fflush(stdout);
	   break;
      case 2: 
	   skewness_analysis(byte_col_freq, num_elements, &byte_flags, 
							num_compressible, byte_pos );
	   break;
      case 3: 
	   kurtosis_analysis(byte_col_freq, num_elements, &byte_flags, 
							num_compressible, byte_pos );
	   break;
      case 4: 
	   isobar_analysis(byte_col_freq, num_elements, &byte_flags, 
							num_compressible, byte_pos );
	   break;
      case 5: 
	   naive_analysis(byte_col_freq, num_elements, &byte_flags, 
							num_compressible, byte_pos );
	   break;
      default:
           byte_flags |=  1 << byte_pos;   
	   break;
    }

#ifdef PLOT 
   #if 1
	check_compressibility ( input_data, byte_col_freq, num_elements, 
			    ele_size, &byte_f, &num_comp, byte_pos );
   #endif 
#endif 
  }

#if 0
  if (measure == 6) 
 	byte_flags = byte_f;
#endif 

/***
  if (*num_compressible == 0) {
	byte_flags = 2; 
  }	  
***/

  STAT_PRINTF("Byte Flags %u, No of bytes in data: %d, number of compressible bytes: %d\n", 
						 byte_flags, ele_size , *num_compressible);

#ifdef PLOT
  fclose(stat_file);
  stat_file = NULL;
#endif

  return byte_flags;

}

void entropy_analysis( acomp_size_t *byte_col_freq, acomp_size_t num_elements, 
		       ubyte_t *byte_flags, int *num_compressible, int byte_pos) {

    double total_entropy = 0.0f, val = 0.0f,
    	   entropy[UCHAR_MAX], E_thres = CHAR_BIT - 1;
    ubyte_t byte_val;

    for ( byte_val =  0; byte_val < UCHAR_MAX; byte_val++ ) {
	val  =  (double)byte_col_freq[ byte_val]/(double)num_elements;

	if( byte_col_freq[byte_val] == 0 || FLOAT_EQ(val, 0) )
		 entropy[byte_val]  =  0.0f;
	else 
	  entropy[byte_val]  =  -1 * val * log2(val);
        
	 DB_PRINTF("Entropy of byte column %d: entropy state %d (%d), \
		(-1* val* log2(val)  =  -1 * %g * %g)  =  %g \n", 
		byte_pos, byte_val, byte_col_freq[byte_val],  val, 
				     log2(val), entropy[byte_val]);

	total_entropy +=  entropy[byte_val];                                                          
   }
   
   STAT_PRINTF("Entropy of the byte column %d is %g \n", byte_pos, total_entropy);					

   #ifdef PLOT
	fprintf(stat_file, "%g ", total_entropy);
   #endif

   #if 1 
   if(total_entropy > E_thres )   {   
            // required more than 7 bits to encode a byte_tacter
	STAT_PRINTF("Entropy of the byte column %d is %g high\n", byte_pos,
							    total_entropy);					
          (*num_compressible) --;
   } else {
	 *byte_flags |=  1 << byte_pos;
   }
   
   DB_PRINTF("Entropy of the byte column %d is %g\n", byte_pos, total_entropy);					
   #endif
 
}


void skewness_analysis( acomp_size_t *byte_col_freq, acomp_size_t num_elements, 
		       ubyte_t *byte_flags, int *num_compressible, int byte_pos) {

   double skewness = 0.0f, std_devitaion = 0.0f, S_thres = 5.0f,
      	     _mean = 0;
   acomp_size_t high = 0, low = num_elements;
   ubyte_t byte_val;
   acomp_size_t freq_zero = 0;	
  
   for ( byte_val = 0; byte_val < UCHAR_MAX; byte_val ++ ) {
	_mean += byte_col_freq[byte_val];
	if ( byte_col_freq[byte_val] == 0 ) {
		freq_zero ++;
	} 
	if(byte_col_freq[byte_val] > high )
		high = byte_col_freq[byte_val];
         
	if(byte_col_freq[byte_val] < low )
		low = byte_col_freq[byte_val];
    }	 
    _mean /=  UCHAR_MAX ;
    
    for( byte_val = 0; byte_val < UCHAR_MAX; byte_val++) { 	
      	skewness += pow(byte_col_freq[byte_val] - _mean, 3)/(UCHAR_MAX );
      	std_devitaion += pow(fabs(byte_col_freq[byte_val] - _mean), 2);
    }

    std_devitaion = sqrt(std_devitaion/UCHAR_MAX); 

    if (std_devitaion) 
    	skewness /=  pow(std_devitaion, 3); 	
    
    skewness *= sqrt((double)(256 * 255/254)); 	

    STAT_PRINTF(" Skewness of byte column %d is %g mean %g std %g\n", byte_pos,
						 skewness, _mean, std_devitaion );

   #ifdef PLOT
    fprintf(stat_file, "%g ", skewness);
   #endif

   #if 1
    if( fabs(skewness) < S_thres && _mean / std_devitaion > 10 )   {    
         STAT_PRINTF("Skewness of the byte column %d is %g too low\n", byte_pos, skewness);					
        (*num_compressible) --;

    } else {
        *byte_flags |=  1<< byte_pos;
    }

    DB_PRINTF("Skewness of the byte column %d is %g\n", byte_pos, skewness);	

  #endif

}


void kurtosis_analysis( acomp_size_t *byte_col_freq, acomp_size_t num_elements, 
		       ubyte_t *byte_flags, int *num_compressible, int byte_pos) {

  double kurtosis = 0.0f, std_devitaion = 0.0f, K_thres = 2.0f,
  						      _mean = 0.0f;
  ubyte_t byte_val;

  for ( byte_val = 0; byte_val < UCHAR_MAX; byte_val ++ ) {
	_mean +=  byte_col_freq[byte_val];
  }
  _mean /= UCHAR_MAX;
 
  for( byte_val  =  0; byte_val < UCHAR_MAX; byte_val ++) { 	
      	kurtosis +=  pow(byte_col_freq[byte_val] - _mean, 4)/UCHAR_MAX;
      	std_devitaion +=  pow(fabs(byte_col_freq[byte_val] - _mean), 2);
  }

  std_devitaion  =  sqrt(std_devitaion / UCHAR_MAX); 
  if (std_devitaion) 
  	kurtosis /=  pow(std_devitaion, 4); 	
  kurtosis -= 3; 	

  STAT_PRINTF(" Kurtosis of byte column %d is %g\n", byte_pos, kurtosis );
  #ifdef PLOT
	fprintf(stat_file, "%g\n", kurtosis );
  #endif
  #if 1	
      if( fabs(kurtosis) < K_thres )   {  
          STAT_PRINTF("kurtosis of the byte column %d is %g too low\n", 
						     byte_pos, kurtosis);					
          (*num_compressible) --;
      } else {
	 *byte_flags |=  1 << byte_pos;
      }
  #endif		   

}

void isobar_analysis( acomp_size_t *byte_col_freq, acomp_size_t num_elements, 
		       ubyte_t *byte_flags, int *num_compressible, int byte_pos) {
  double IS_theta = 1.38, I_thres = IS_theta * (num_elements/256); 
  ubyte_t byte_val;
  int flag = 0;

  for ( byte_val = 0; byte_val < UCHAR_MAX; byte_val ++ ) {
	if(byte_col_freq[byte_val] > I_thres ) {
		flag = 1; break; 
	}
  }
  #ifdef PLOT
	fprintf(stat_file, "%d\n", flag );
  #endif

  if (!flag) {
	STAT_PRINTF("Isobar method: The byte column %d is not compressible", byte_pos);
              (*num_compressible) --;
  } else {
	*byte_flags |=  1 << byte_pos;
  }  

}

void naive_analysis( acomp_size_t *byte_col_freq, acomp_size_t num_elements, 
		       ubyte_t *byte_flags, int *num_compressible, int byte_pos) {

   acomp_size_t sum_freq = 0, freq, max_freq[3];
   int k;
   ubyte_t byte_val;	
   double _mean = 0.0f;

   for (k = 0; k < 3; k++) {
        max_freq[k] = 0;
   }
  
  for ( byte_val = 0; byte_val < UCHAR_MAX; byte_val ++ ) {
	_mean +=  byte_col_freq[byte_val];
  }
  _mean /= UCHAR_MAX;

   for( byte_val  =  0; byte_val < UCHAR_MAX; byte_val ++) { 	
        freq = byte_col_freq[byte_val];
        for ( k = 0; k < 3; k++ ) {  
      	   if(max_freq[k] < freq) { 
		max_freq[k] = freq;
		freq = max_freq[k];  
    	   }
        }
   } 

  for (k = 0; k < 3; k++) {
        sum_freq += max_freq[k];
  }

  STAT_PRINTF(" Sum of top 3 frequencies of the Byte %d is %lu, \
                mean frequency is %g (creteria is >= 0.60 * %lu  = %g)\n\n", 
                byte_pos, sum_freq, _mean, num_elements, 0.60 * num_elements);
  #if 1
  if ( sum_freq < 0.60 * num_elements )  {    
          STAT_PRINTF("Sum of top 3 frequencies of the byte column %d is \
                                        %lu too low \n", byte_pos, sum_freq);
	(*num_compressible) --;
  } else {
	*byte_flags |=  1 << byte_pos;
  }
 #endif	

}

