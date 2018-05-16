#ifndef __ACOMP_CHUNKING_H__
#define __ACOMP_CHUNKING_H__

#include <assert.h>
#include "util.h"

#define MAX_CHUNKS 6
//#define DEFAULT_CHUNK_SIZE 2097152
#define DEFAULT_CHUNK_SIZE 3145728

struct chunk_info_t 
{
  int  chunk_id;
  byte_t     *data;
  acomp_size_t  chunk_size;
  acomp_size_t  chunk_offset; 
  acomp_size_t  compSize; 

};


int get_chunks_info_1D (struct chunk_info_t **dchunks, acomp_size_t* dim_size, 
					  int ndims, int ele_size);

int get_overlap_chunks_1D( struct chunk_info_t ** chunks, int ndims, 
			   acomp_size_t* dimLen, int num_chunks, 
			   acomp_size_t* subreq_start, 
			   acomp_size_t* subreq_count, int ele_size);
 
#endif
