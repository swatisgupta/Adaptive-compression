
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

#ifdef ACOMP

#include "acomp.h"

int get_datatype_size(enum ADIOS_DATATYPES type) {
        switch(type) {
         case adios_unknown : return -1;             /* (size) */
         case adios_byte: return 1;                  /* (1) */
         case adios_short: return 2;                 /* (2) */
         case adios_integer: return 4;               /* (4) */
         case adios_long: return 8;                  /* (8) */
  
         case adios_unsigned_byte: return 1;         /* (1) */
         case adios_unsigned_short: return 2;      /* (2) */
         case adios_unsigned_integer: return 4;    /* (4) */
         case adios_unsigned_long: return 8;       /* (8) */
  
         case adios_real: return 4;                  /* (4) */
         case adios_double: return 8;               /* (8) */
         case adios_long_double: return 16;          /* (16) */
  
         case adios_complex: return 8;             /* (8) */
         case adios_double_complex: return 16;      /* (16) */
  
        }    
        return -1;
}

uint16_t adios_transform_acomp_get_metadata_size(struct adios_transform_spec *transform_spec)
{
    return acomp_get_max_metadata_size(); 
}


void adios_transform_acomp_transformed_size_growth(
		const struct adios_var_struct *var, const struct adios_transform_spec *transform_spec,
		unsigned int *constant_factor, double *linear_factor, double *capped_linear_factor, unsigned int *capped_linear_cap)
{
	// Do nothing (defaults to "no transform effect on data size")
}

int adios_transform_acomp_apply(struct adios_file_struct *fd,
                                struct adios_var_struct *var,
                                unsigned int *transformed_len,
                                int use_shared_buffer,
                                int *wrote_to_shared_buffer)
{
    // Assume this function is only called for ACOMP transform type
    assert(var->transform_type == adios_transform_acomp);

    // Get the input data and data length 
    unsigned int  spatial_size = 1;
    unsigned int input_size = adios_transform_get_pre_transform_var_size(var);
    unsigned int spatial_dim_len[5] ,  actual_output_size = 0; //output_size;
    struct adios_dimension_struct * var_dim;
    int ndims = 0, t_size;
    double wt_cs=DEFAULT_W_CS, wt_cr=DEFAULT_W_CR, delt=DEFAULT_DELTA;
    
    enum ADIOS_DATATYPES dtype =var->type;
    int is_real = 0;  
  
    char* output_buff = NULL;
    char *input_buffer = (char*)var->data;
    unsigned char * metadata;
    
    char * varname;
    char * wr = NULL, *ws = NULL;
    

/**    // parse the compressiong parameter
    if (var->transform_spec->param_count > 0 ) {
        threshold = atof(var->transform_spec->params[0].key);
    } 
	
    if (var->transform_spec->param_count > 1 ) {
        wt_cr = atof(var->transform_spec->params[1].key);
        wt_cs = atof(var->transform_spec->params[2].key);
    }
**/
    if( (wr  = getenv("WEIGHT_CR")) != NULL) {
	     wt_cr = atof(wr);
    }

    if( (ws  = getenv("WEIGHT_CS")) != NULL) {
	      wt_cs = atof(ws);
    }
  
    //printf(" wr %g, ws %g \n", wt_cr, wt_cs); 
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
        output_buff = malloc(input_size * 20);
        if (!output_buff)
        {
            log_error("Out of memory allocating %llu bytes for %s for zlib transform\n", input_size, var->name);
            return 0;
        }
    }

    //get the spatial size of the variable to be written...
    dtype = var->pre_transform_type; 
    t_size = get_datatype_size(var->pre_transform_type); 
    varname = var->name;
    if(dtype == adios_real || dtype == adios_double) {
         is_real = 1;
    }

    acomp_conf(wt_cs, wt_cr, delt, fd->group->process_id);

    var_dim = var->pre_transform_dimensions;
    while(var_dim != NULL) {
         if(var_dim->dimension.is_time_index == adios_flag_no) {
            spatial_dim_len[ndims] = adios_get_dim_value(&(var_dim->dimension)); 
	        spatial_size *= spatial_dim_len[ndims];
	        ndims++;	
	     }
         var_dim = var_dim->next;
    }
    
    metadata = var->transform_metadata;
    //printf("Is real %d, size %d spatial_size %u input_size %u \n", is_real, spatial_size, input_size);
    actual_output_size = acomp_compress( varname, ndims, spatial_dim_len, t_size, is_real, 
	                                    input_size, spatial_size, output_buff, 
	                                    input_buffer,  metadata);

    //printf("Actual size %u, Input size %u\n", actual_output_size, input_size);       
    // Wrap up, depending on buffer mode
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

DECLARE_TRANSFORM_WRITE_METHOD_UNIMPL(acomp)

#endif
