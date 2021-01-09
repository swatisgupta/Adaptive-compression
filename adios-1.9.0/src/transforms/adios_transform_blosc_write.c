
#include <stdint.h>
#include <assert.h>
#include <limits.h>
#include <sys/time.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

#include "core/adios_logger.h"
#include "core/transforms/adios_transforms_common.h"
#include "core/transforms/adios_transforms_write.h"
#include "core/transforms/adios_transforms_hooks_write.h"
#include "core/transforms/adios_transforms_util.h"
#include "core/adios_internals.h" 

#ifdef BLOSC

#include "blosc.h"

int get_datatype_size(enum ADIOS_DATATYPES type);

/** 
int get_datatype_size(enum ADIOS_DATATYPES type) {
        switch(type) {
         case adios_unknown : return -1;             
         case adios_byte: return 1;                  
         case adios_short: return 2;                 
         case adios_integer: return 4;               
         case adios_long: return 8;                    
         case adios_unsigned_byte: return 1;         
         case adios_unsigned_short: return 2;      
         case adios_unsigned_integer: return 4;    
         case adios_unsigned_long: return 8;       
  
         case adios_real: return 4;                 
         case adios_double: return 8;               
         case adios_long_double: return 16;         
  
         case adios_complex: return 8;             
         case adios_double_complex: return 16;      
  
        }    
        return -1;
}
***/

uint16_t adios_transform_blosc_get_metadata_size(struct adios_transform_spec *transform_spec)
{
    return 0; 
}


void adios_transform_blosc_transformed_size_growth(
		const struct adios_var_struct *var, const struct adios_transform_spec *transform_spec,
		unsigned int *constant_factor, double *linear_factor, double *capped_linear_factor, unsigned int *capped_linear_cap)
{
	// Do nothing (defaults to "no transform effect on data size")
}

int adios_transform_blosc_apply(struct adios_file_struct *fd,
                                struct adios_var_struct *var,
                                unsigned int *transformed_len,
                                int use_shared_buffer,
                                int *wrote_to_shared_buffer)
{
    // Assume this function is only called for BLOSC transform type
    assert(var->transform_type == adios_transform_blosc);

    // Get the input data and data length 
    unsigned int input_size = adios_transform_get_pre_transform_var_size(var);
    unsigned int actual_output_size = 0; //output_size;
    int t_size;
    int clevel = 5, shuffle = 1;
    
    enum ADIOS_DATATYPES dtype =var->type;
    int nthreads = 4, pnthreads, rcode;
   
    char *output_buff = NULL;
    char *input_buffer = (char*)var->data;
    
    char * compressor = "blosclz";
    char * cl = NULL, *sh = NULL, *nth = NULL;
    FILE *fp = fopen("stat2.txt", "a"); 
    struct timeval start_time, end_time;
    double total_time = 0;


    if( (cl  = getenv("CLEVEL")) != NULL) {
	     clevel = atof(cl);
    }

    if( (sh  = getenv("SHUFFLE")) != NULL) {
	      shuffle = atof(sh);
    }
    
    if( (nth = getenv("NTHREADS")) != NULL) {
	      nthreads = atof(nth);
    } 
    
    if( getenv("COMRPESSOR") != NULL) {
	      compressor = getenv("COMRPESSOR");
    }
 
   // decide the output buffer
    if (use_shared_buffer)    // If shared buffer is permitted, serialize to there
    {
        *wrote_to_shared_buffer = 1;
        if (!shared_buffer_reserve(fd, input_size))
        {
            log_error("Out of memory allocating %llu bytes for %s for zlib transform\n", input_size, var->name);
            return 0;
        }

        // Write directly to the shared buffer
        output_buff = fd->buffer + fd->offset;
    }
    else    // Else, fall back to var->adata memory allocation
    {
        *wrote_to_shared_buffer = 0;
        output_buff = malloc(input_size);
        if (!output_buff)
        {
            log_error("Out of memory allocating %llu bytes for %s for zlib transform\n", input_size, var->name);
            return 0;
        }
    }

    //get the spatial size of the variable to be written...
    t_size = get_datatype_size(var->pre_transform_type); 


    blosc_init();

    pnthreads = blosc_set_nthreads(nthreads);

    rcode = blosc_set_compressor(compressor);

    if (rcode < 0) {
      printf("Error setting %s compressor.  It really exists?",
	     compressor);
      blosc_set_compressor("blosclz");
    }     

    gettimeofday(&start_time, NULL);    
	
    actual_output_size = blosc_compress(clevel, shuffle, t_size, input_size, input_buffer, output_buff, input_size );

    gettimeofday(&end_time, NULL);    

    blosc_destroy();

    total_time  =  (double)(end_time.tv_usec - start_time.tv_usec) / 1000000 
      + (double)(end_time.tv_sec - start_time.tv_sec);

    fprintf(fp, "Total time %g compressed size %d clevel %d\n", total_time, actual_output_size, clevel);
    fclose(fp);

    if( actual_output_size <= 0 ) {
	actual_output_size = input_size;
         	
    }
	
    if (use_shared_buffer)
    {
        shared_buffer_mark_written(fd, actual_output_size);
    }
    else
    {
        var->adata = output_buff;
        var->data_size = actual_output_size;
        var->free_data = adios_flag_yes;
    }

    *transformed_len = actual_output_size; 
    	
    return 1;
}

#else

DECLARE_TRANSFORM_WRITE_METHOD_UNIMPL(blosc)

#endif
