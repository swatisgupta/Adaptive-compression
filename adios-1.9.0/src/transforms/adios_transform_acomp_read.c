
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

#ifdef COMP

#include "zlib.h"
#include "acomp.h"

int adios_transform_acomp_is_implemented (void) {return 1;}



int adios_transform_acomp_generate_read_subrequests(adios_transform_read_request *reqgroup,
                                                    adios_transform_pg_read_request *pg_reqgroup)
{
    void *buf = malloc(pg_reqgroup->raw_var_length);
    assert(buf);
    adios_transform_raw_read_request *subreq = adios_transform_raw_read_request_new_whole_pg(pg_reqgroup, buf);
    adios_transform_raw_read_request_append(pg_reqgroup, subreq);
    return 0;
}

// Do nothing for individual subrequest
adios_datablock * adios_transform_acomp_subrequest_completed(adios_transform_read_request *reqgroup,
                                                            adios_transform_pg_read_request *pg_reqgroup,
                                                            adios_transform_raw_read_request *completed_subreq)
{
    return NULL;
}



adios_datablock * adios_transform_acomp_pg_reqgroup_completed(adios_transform_read_request *reqgroup,
                                                             adios_transform_pg_read_request *completed_pg_reqgroup)
{
    unsigned int compressed_size = (unsigned int)completed_pg_reqgroup->raw_var_length;
 
    unsigned char* metadata = (unsigned char*)completed_pg_reqgroup->transform_metadata;
    unsigned char* compressed_data = (unsigned char*) completed_pg_reqgroup->subreqs->data;
    
    unsigned int uncompressed_size_meta;

    unsigned int t_size = adios_get_type_size(reqgroup->transinfo->orig_type, "");

    unsigned int uncompressed_size, uncompressed_subreq_size = t_size; 
 
    int ndims = 0;
   
    
    unsigned char * output_data = NULL;
    enum ADIOS_DATATYPES var_type = adios_real;
    
    int d = 0, safe_size = 1024;
    unsigned int dimLen[5], req_start[5], req_count[5]; 	 

    uncompressed_size = t_size;
    for(d = 0; d < reqgroup->transinfo->orig_ndim; d++)
    {   dimLen[d] = completed_pg_reqgroup->orig_varblock->count[d];
	      req_start[d] = completed_pg_reqgroup->pg_bounds_sel->u.bb.start[d];
	      req_count[d] = completed_pg_reqgroup->pg_bounds_sel->u.bb.count[d];
        uncompressed_size *= 
            (unsigned int)(completed_pg_reqgroup->orig_varblock->count[d]);

    }

    acomp_get_actual_size( metadata, &uncompressed_size_meta);
    
    if(uncompressed_size_meta != uncompressed_size)
    {
        printf("WARNING: possible wrong data size or corrupted metadata\n");
    }
    
    ndims = reqgroup->transinfo->orig_ndim;
    
    if ( acomp_decompress( t_size, ndims, metadata, compressed_data,
                        output_data, 
                    dimLen, req_start, req_count ) == 0) {
       output_data = compressed_data;
    } 


    return adios_datablock_new_whole_pg(reqgroup, completed_pg_reqgroup, output_data);
}

adios_datablock * adios_transform_acomp_reqgroup_completed(adios_transform_read_request *completed_reqgroup)
{
    return NULL;
}


#else

DECLARE_TRANSFORM_READ_METHOD_UNIMPL(acomp);

#endif

