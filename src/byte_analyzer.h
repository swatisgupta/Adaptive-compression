#ifndef __BYTE_ANALYZER_H__
#define __BYTE_ANALYZER_H__

#include "util.h"

#define FLOAT_EQ(a,b) fabs((a) - (b)) < (DBL_EPSILON * fabs((a) + (b) ))
 
ubyte_t get_compressibility_flags(int dtype, int ele_size, acomp_size_t data_size, 
			    byte_t* data, int *num_compressible); 

void check_compressibility (byte_t* data, int *byte_col_freq, acomp_size_t num_elements, 
			    int ele_size, ubyte_t *byte_flags, int *num_compressible,
								    int byte_pos );

#endif
