#include "chunking.h" 

int get_chunks_info_1D (struct chunk_info_t **dchunks, acomp_size_t* dim_size, 
					  int ndims, int ele_size) {

  int i, dim = ndims -1, num_chunks;
  acomp_size_t total_size = 1, chunk_size = DEFAULT_CHUNK_SIZE / ele_size;

  struct chunk_info_t * chunks  =  NULL;	

  for ( ; dim >= 0; dim-- ) { 
    //printf("Dimlen %d  = %u \n", dim, dim_size[dim] );
    total_size *=  dim_size[dim];
  }
	
  num_chunks  =  total_size / chunk_size;
       
  if(total_size % chunk_size) 
    num_chunks ++;
 
  if( num_chunks > MAX_CHUNKS) { 
    chunk_size  =  total_size / MAX_CHUNKS;
    num_chunks  =   total_size / chunk_size;
  }

  if( num_chunks  ==  0) 
    num_chunks  =  1;
	
  *dchunks  =  (struct chunk_info_t *) malloc(sizeof(struct chunk_info_t) * num_chunks);

  assert(*dchunks !=  NULL );
	
  chunks = *dchunks;

  chunk_size *= ele_size;
	
  DB_PRINTF("chunk_size %lu, total_size %lu, num_chunks %d ndims %d\n", chunk_size, 
						total_size * ele_size, num_chunks, ndims);
   
  for (i  =  0; i < num_chunks ; i++) {
    chunks[i].chunk_offset  =  i * chunk_size;
    chunks[i].chunk_id  =  i;
    if( i  ==  num_chunks - 1) {
      chunks[i].chunk_size  =  (total_size * ele_size) - i * chunk_size;
      chunks[i].compSize  =  (total_size * ele_size) - i * chunk_size;
    } else  {
      chunks[i].chunk_size  =  chunk_size;
      chunks[i].compSize  =  chunk_size;
    }
    chunks[i].data  =  NULL;
    DB_PRINTF("chunk %d  => offset %lu, size %lu\n", chunks[i].chunk_id, 
					   chunks[i].chunk_offset, chunks[i].chunk_size);
  }

  return num_chunks;
}

int get_overlap_chunks_1D( struct chunk_info_t **dchunks, int ndims, acomp_size_t* dimLen, int num_chunks, 
			   acomp_size_t* subreq_start, acomp_size_t* subreq_count, int ele_size) {

  int i =  ndims - 1;
  acomp_size_t treq_size  =  ele_size, total_size  =  1, offset  = 0;
  int chunk_start, chunk_end;
  acomp_size_t chunk_size = DEFAULT_CHUNK_SIZE / ele_size; 
  struct chunk_info_t * chunks  =  NULL;

  for(; i>= 0; i--) { 
    offset *=  (subreq_start[i]* ele_size);
    treq_size *=  subreq_count[i];   
    total_size *=  dimLen[i];   
  }
	
  if( num_chunks == MAX_CHUNKS ) {
      chunk_size = total_size / MAX_CHUNKS;	
  }	

  chunk_size *= ele_size;

  chunk_start  =  offset / chunk_size;

  chunk_end  =  (offset + treq_size) / chunk_size; 
	
  if( chunk_end >=  num_chunks) 
    chunk_end  =  num_chunks -1;

  *dchunks  =  malloc(sizeof(struct chunk_info_t)*(chunk_end - chunk_start + 1));

  assert(*dchunks !=  NULL);
   
  chunks  =  *dchunks;
  DB_PRINTF("Noverlap %d, num_chunks %d,  chunk_start %d, chunk_end %d, \
		treq_size %lu, total_size %lu\n", (chunk_end - chunk_start + 1), 
		num_chunks, chunk_start, chunk_end, treq_size, total_size * ele_size);
        

  for( i  =  0 ; chunk_start <=  chunk_end; chunk_start++) {
    chunks[i].chunk_id  =  chunk_start;
    chunks[i].chunk_offset  =  chunk_start * chunk_size;
    if(chunk_start  ==  chunk_end) { 
      chunks[i].chunk_size  =  (total_size * ele_size) - chunk_start * chunk_size;
      chunks[i].compSize  =  chunks[i].chunk_size;
    } else {
      chunks[i].chunk_size  =  chunk_size;
      chunks[i].compSize  =  chunk_size;
    }
    chunks[i].data  =  NULL; 
    DB_PRINTF("chunk %d  => chunk_offset %lu, chunk_size %lu\n", chunks[i].chunk_id, 
					chunks[i].chunk_offset, chunks[i].chunk_size);
    i++;
  }
	
  return i; 
}

#if 0
acomp_size_t restore_chunk_2D(byte_t* input_data, byte_t* output_data, acomp_size_t offset, 
		      acomp_size_t count, int ele_size) {
    
  byte_t* out_data  =  output_data;
  byte_t* in_data  =  input_data;
  acomp_size_t datalen  =  0; 

  memcpy(out_data + offset * ele_size, in_data, count * ele_size);
  datalen  =  count * ele_size;
  return datalen;

}


acomp_size_t extract_chunk_2D(byte_t* input_data, byte_t* output_data, acomp_size_t offset, 
			  acomp_size_t count, int ndims, int ele_size) {
 
  byte_t* out_data  =  output_data;
  byte_t* in_data  =  input_data;
  acomp_size_t datalen  =  0; 
    
  memcpy(out_data, in_data + offset * ele_size, count * ele_size);
  datalen  =  count * ele_size;
  return datalen;
}
#endif


