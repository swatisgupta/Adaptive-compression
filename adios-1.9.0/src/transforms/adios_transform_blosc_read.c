
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <assert.h>

#include "core/util.h"
#include "core/adios_logger.h"
#include "public/adios_types.h"
#include "core/transforms/adios_transforms_hooks_read.h"
#include "core/transforms/adios_transforms_reqgroup.h"
#include "core/adios_internals.h" // adios_get_type_ize()

#ifdef BLOSC

#include "blosc.h"

int adios_transform_blosc_is_implemented (void) {return 1;}



int adios_transform_blosc_generate_read_subrequests(adios_transform_read_request *reqgroup,
                                                    adios_transform_pg_read_request *pg_reqgroup)
{
    void *buf = malloc(pg_reqgroup->raw_var_length);
    assert(buf);
    adios_transform_raw_read_request *subreq = adios_transform_raw_read_request_new_whole_pg(pg_reqgroup, buf);
    adios_transform_raw_read_request_append(pg_reqgroup, subreq);
    return 0;
}

// Do nothing for individual subrequest
adios_datablock * adios_transform_blosc_subrequest_completed(adios_transform_read_request *reqgroup,
                                                            adios_transform_pg_read_request *pg_reqgroup,
                                                            adios_transform_raw_read_request *completed_subreq)
{
    return NULL;
}



adios_datablock * adios_transform_blosc_pg_reqgroup_completed(adios_transform_read_request *reqgroup,
                                                             adios_transform_pg_read_request *completed_pg_reqgroup)
{
    unsigned int compressed_size = (unsigned int)completed_pg_reqgroup->raw_var_length;
 
    unsigned char* compressed_data = (unsigned char*) completed_pg_reqgroup->subreqs->data;
    
    unsigned int t_size = adios_get_type_size(reqgroup->transinfo->orig_type, "");

    unsigned int uncompressed_size; 
 
    int d = 0;
    
    unsigned char * output_data = NULL, *nth = NULL;

    int nthreads = 4;
   
    if ( ( nth = getenv("NTHREADS")) != NULL ) {
  		nthreads = atoi(nth); 
    }	  	

    blosc_set_nthreads(nthreads);    

    uncompressed_size = t_size;

    for(d = 0; d < reqgroup->transinfo->orig_ndim; d++) {
        uncompressed_size *= 
            (unsigned int)(completed_pg_reqgroup->orig_varblock->count[d]);

    }
   
    output_data = malloc(uncompressed_size); 
    
    if ( blosc_decompress( compressed_data, output_data, uncompressed_size ) < 0) {
       output_data = compressed_data;
    } 


    return adios_datablock_new_whole_pg(reqgroup, completed_pg_reqgroup, output_data);
}

adios_datablock * adios_transform_blosc_reqgroup_completed(adios_transform_read_request *completed_reqgroup)
{
    return NULL;
}


#else

DECLARE_TRANSFORM_READ_METHOD_UNIMPL(blosc);

#endif

