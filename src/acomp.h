#ifndef ACOMP_H
#define ACOMP_H

#include "acomp_types.h"

#define DEFAULT_W_CS 0.0f
#define DEFAULT_W_CR 1.0f
#define DEFAULT_DELTA 100.00f
#define IO_BANDWIDTH 2.7E-09

int acomp_decompress(int ele_size, int ndims, unsigned char* metadata, char* compressed_data, 
	                char* output_data, 
	                acomp_size_t* dimLen, acomp_size_t* subset_start, acomp_size_t * subset_count );

int acomp_get_max_metadata_size();

acomp_size_t acomp_compress( char * varname, int ndims, acomp_size_t* dim_len, int ele_size, int is_real, 
	                acomp_size_t input_size_bytes, acomp_size_t n_elements,  
	                char* output_buff,  char* input_buff, 
	                unsigned char * metadata);


void acomp_get_actual_size( unsigned char* metadata, acomp_size_t *uncompressed_size_meta );

void acomp_conf(double w_cs, double w_cr, double delta);

int acomp_get_nchunks(unsigned char* metadata, acomp_size_t** chunk_meta);

unsigned char acomp_get_byteflags(unsigned char* metadata);

int acomp_is_compressed(unsigned char* metadata, int* compression_method);

#endif
